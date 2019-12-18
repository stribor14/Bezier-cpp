#include "BezierCpp/bezier.h"
#include "BezierCpp/legendre_gauss.h"

#include <exception>
#include <numeric>

#include <iostream>

#include <unsupported/Eigen/MatrixFunctions>

inline double factorial(uint k) { return std::tgamma(k + 1); }
inline double binomial(uint n, uint k) { return factorial(n) / (factorial(k) * factorial(n - k)); }

namespace Bezier
{

Curve::CoeffsMap Curve::bernstein_coeffs_ = CoeffsMap();
Curve::CoeffsMap Curve::splitting_coeffs_left_ = CoeffsMap();
Curve::CoeffsMap Curve::splitting_coeffs_right_ = CoeffsMap();
Curve::CoeffsMap Curve::elevate_order_coeffs_ = CoeffsMap();
Curve::CoeffsMap Curve::lower_order_coeffs_ = CoeffsMap();

void Curve::resetCache()
{
  cached_derivative_.reset();
  cached_roots_.reset();
  cached_bounding_box_tight_.reset();
  cached_bounding_box_relaxed_.reset();
  cached_polyline_.reset();
}

Curve::Coeffs Curve::bernsteinCoeffs() const
{
  if (bernstein_coeffs_.find(N_) == bernstein_coeffs_.end())
  {
    bernstein_coeffs_.insert(std::make_pair(N_, Coeffs::Zero(N_, N_)));
    bernstein_coeffs_.at(N_).diagonal(-1) = -Eigen::ArrayXd::LinSpaced(N_ - 1, 1, N_ - 1);
    bernstein_coeffs_.at(N_) = bernstein_coeffs_.at(N_).exp();
    for (uint k = 0; k < N_; k++)
      bernstein_coeffs_.at(N_).row(k) *= binomial(N_ - 1, k);
  }
  return bernstein_coeffs_.at(N_);
}

Curve::Coeffs Curve::splittingCoeffsLeft(double z) const
{
  if (z == 0.5)
  {
    if (splitting_coeffs_left_.find(N_) == splitting_coeffs_left_.end())
    {
      splitting_coeffs_left_.insert(std::make_pair(N_, Coeffs::Zero(N_, N_)));
      splitting_coeffs_left_.at(N_).diagonal() = Eigen::pow(0.5, Eigen::ArrayXd::LinSpaced(N_, 0, N_ - 1));
      splitting_coeffs_left_.at(N_) = bernsteinCoeffs().inverse() * splitting_coeffs_left_.at(N_) * bernsteinCoeffs();
    }
    return splitting_coeffs_left_.at(N_);
  }
  else
  {
    Curve::Coeffs coeffs(Coeffs::Zero(N_, N_));
    coeffs.diagonal() = Eigen::pow(z, Eigen::ArrayXd::LinSpaced(N_, 0, N_ - 1));
    coeffs = bernsteinCoeffs().inverse() * coeffs * bernsteinCoeffs();
    return coeffs;
  }
}

Curve::Coeffs Curve::splittingCoeffsRight(double z) const
{
  if (z == 0.5)
  {
    if (splitting_coeffs_right_.find(N_) == splitting_coeffs_right_.end())
    {
      splitting_coeffs_right_.insert(std::make_pair(N_, Coeffs::Zero(N_, N_)));
      Curve::Coeffs temp_splitting_coeffs_left = splittingCoeffsLeft();
      for (uint k = 0; k < N_; k++)
        splitting_coeffs_right_.at(N_).block(k, k, 1, N_ - k) =
            temp_splitting_coeffs_left.block(N_ - 1 - k, 0, 1, N_ - k);
    }
    return splitting_coeffs_right_.at(N_);
  }
  else
  {
    Curve::Coeffs coeffs(Coeffs::Zero(N_, N_));
    Curve::Coeffs temp_splitting_coeffs_left = splittingCoeffsLeft(z);
    for (uint k = 0; k < N_; k++)
      coeffs.block(k, k, 1, N_ - k) = temp_splitting_coeffs_left.block(N_ - 1 - k, 0, 1, N_ - k);
    return coeffs;
  }
}

Curve::Coeffs Curve::elevateOrderCoeffs(uint n) const
{
  if (elevate_order_coeffs_.find(n) == elevate_order_coeffs_.end())
  {
    elevate_order_coeffs_.insert(std::make_pair(n, Coeffs::Zero(n + 1, n)));
    elevate_order_coeffs_.at(n).diagonal() = 1 - Eigen::ArrayXd::LinSpaced(n, 0, n - 1) / n;
    elevate_order_coeffs_.at(n).diagonal(-1) = Eigen::ArrayXd::LinSpaced(n, 1, n) / n;
  }
  return elevate_order_coeffs_.at(n);
}

Curve::Coeffs Curve::lowerOrderCoeffs(uint n) const
{
  if (lower_order_coeffs_.find(n) == lower_order_coeffs_.end())
  {
    lower_order_coeffs_.insert(std::make_pair(n, Coeffs::Zero(n - 1, n)));
    lower_order_coeffs_.at(n).noalias() =
        (elevateOrderCoeffs(n - 1).transpose() * elevateOrderCoeffs(n - 1)).inverse() *
        elevateOrderCoeffs(n - 1).transpose();
  }
  return lower_order_coeffs_.at(n);
}

Curve::Curve(const Eigen::MatrixX2d& points)
{
  N_ = static_cast<uint>(points.rows());
  control_points_ = points;
}

Curve::Curve(const PointVector& points)
{
  N_ = static_cast<uint>(points.size());
  control_points_.resize(N_, 2);
  for (uint k = 0; k < N_; k++)
    control_points_.row(k) = points.at(k);
}

Curve::Curve(const Curve& curve) : Curve(curve.getControlPoints()) {}

uint Curve::getOrder() { return N_ - 1; }

PointVector Curve::getControlPoints() const
{
  PointVector points(N_);
  for (uint k = 0; k < N_; k++)
    points.at(k) = control_points_.row(k);
  return points;
}

std::pair<Point, Point> Curve::getEndPoints() const
{
  return std::make_pair(control_points_.row(0), control_points_.row(N_ - 1));
}

PointVector Curve::getPolyline(double smoothness, double precision) const
{
  if (!cached_polyline_ || cached_polyline_params_ != std::make_tuple(smoothness, precision))
  {
    PointVector* polyline = new PointVector;
    std::vector<Eigen::MatrixX2d> subcurves;
    subcurves.push_back(control_points_);
    polyline->push_back(control_points_.row(0));
    while (!subcurves.empty())
    {
      auto cp = subcurves.back();
      subcurves.pop_back();

      double string_length = (cp.row(0) - cp.row(N_ - 1)).norm();
      double hull_length = 0.0;
      for (uint k = 1; k < N_; k++)
        hull_length += (cp.row(k) - cp.row(k - 1)).norm();

      if (hull_length <= smoothness * string_length || string_length <= precision)
      {
        polyline->push_back(cp.row(N_ - 1));
      }
      else
      {
        subcurves.push_back(splittingCoeffsRight(0.5) * cp);
        subcurves.push_back(splittingCoeffsLeft(0.5) * cp);
      }
    }
    (const_cast<Curve*>(this))->cached_polyline_params_ = {smoothness, precision};
    (const_cast<Curve*>(this))->cached_polyline_.reset(polyline);
  }
  return *cached_polyline_;
}

double Curve::getLength() const { return getLength(0.0, 1.0); }

double Curve::getLength(double t) const { return getLength(0.0, t); }

double Curve::getLength(double t1, double t2) const
{
  double sum = 0;

  for (uint k = 0; k < LegendreGauss::N; k++)
    sum += LegendreGauss::weights.at(k) *
           getDerivativeAt(LegendreGauss::abcissae.at(k) * (t2 - t1) / 2 + (t1 + t2) / 2).norm();

  return sum * (t2 - t1) / 2;
}

double Curve::iterateByLength(double t, double s, double epsilon, std::size_t max_iter) const
{
  const double s_t = getLength(t);
  //  if (s_t + s < 0 || s_t + s > getLength())
  //    throw std::out_of_range{"Resulting parameter t not in [0, 1] range."};
  if (s_t + s < 0)
    return 0;
  if (s_t + s > getLength())
    return 1;

  std::size_t current_iter = 0;
  while (current_iter < max_iter)
  {
    // Newton-Raphson
    double f = (getLength(t) - s_t - s);
    double f_d = getDerivativeAt(t).norm();
    t -= f / f_d;

    // if there is no change to t
    if (std::fabs(f) < epsilon)
      break;

    current_iter++;
  }

  return t;
}

void Curve::reverse()
{
  control_points_ = control_points_.colwise().reverse().eval();
  resetCache();
}

void Curve::manipulateControlPoint(uint idx, const Point& point)
{
  control_points_.row(idx) = point;
  resetCache();
}

void Curve::manipulateCurvature(double t, const Point& point)
{
  if (N_ < 3 || N_ > 4)
    throw std::logic_error{"Only quadratic and cubic curves can be manipulated"};

  double r =
      std::fabs((std::pow(t, N_ - 1) + std::pow(1 - t, N_ - 1) - 1) / (std::pow(t, N_ - 1) + std::pow(1 - t, N_ - 1)));
  double u = std::pow(1 - t, N_ - 1) / (std::pow(t, N_ - 1) + std::pow(1 - t, N_ - 1));
  Point C = u * control_points_.row(0) + (1 - u) * control_points_.row(N_ - 1);
  Point B = point;
  Point A = B - (C - B) / r;

  switch (N_)
  {
  case 3:
    control_points_.row(1) = A;
    break;
  case 4:
    Point e1 = control_points_.row(0) * std::pow(1 - t, 2) + control_points_.row(1) * 2 * t * (1 - t) +
               control_points_.row(2) * std::pow(t, 2);
    Point e2 = control_points_.row(1) * std::pow(1 - t, 2) + control_points_.row(2) * 2 * t * (1 - t) +
               control_points_.row(3) * std::pow(t, 2);
    e1 = B + e1 - valueAt(t);
    e2 = B + e2 - valueAt(t);
    Point v1 = A - (A - e1) / (1 - t);
    Point v2 = A + (e2 - A) / t;
    control_points_.row(1).noalias() = control_points_.row(0) + (v1.transpose() - control_points_.row(0)) / t;
    control_points_.row(2).noalias() = control_points_.row(3) - (control_points_.row(3) - v2.transpose()) / (1 - t);
  }
  resetCache();
}

void Curve::elevateOrder()
{
  Eigen::MatrixXd new_points = elevateOrderCoeffs(N_) * control_points_;
  control_points_.resize(++N_, 2);
  control_points_ = new_points;
  resetCache();
}

void Curve::lowerOrder()
{
  if (N_ == 2)
    throw std::logic_error{"Cannot further reduce the order of curve."};
  Eigen::MatrixXd new_points = lowerOrderCoeffs(N_) * control_points_;
  control_points_.resize(--N_, 2);
  control_points_ = new_points;
  resetCache();
}

Point Curve::valueAt(double t) const
{
  if (N_ == 0)
    return Bezier::Point(0, 0);
  Eigen::VectorXd power_basis = Eigen::pow(t, Eigen::ArrayXd::LinSpaced(N_, 0, N_ - 1));
  return (power_basis.transpose() * bernsteinCoeffs() * control_points_).transpose();
}

double Curve::curvatureAt(double t) const
{
  Point d1 = getDerivativeAt(t);
  Point d2 = getDerivativeAt(2, t);

  return (d1.x() * d2.y() - d1.y() * d2.x()) / std::pow(d1.norm(), 3);
}

double Curve::curvatureDerivativeAt(double t) const
{
  auto d1 = getDerivativeAt(t);
  auto d2 = getDerivativeAt(2, t);
  auto d3 = getDerivativeAt(3, t);

  return (d1.x() * d3.y() - d1.y() * d3.x()) / std::pow(d1.norm(), 3) -
         3 * d1.dot(d2) * (d1.x() * d2.y() - d1.y() * d2.x()) / std::pow(d1.norm(), 5);
}

Vec2 Curve::tangentAt(double t, bool normalize) const
{
  Point p(getDerivativeAt(t));
  if (normalize && p.norm() > 0)
    p.normalize();
  return p;
}

Vec2 Curve::normalAt(double t, bool normalize) const
{
  Point tangent = tangentAt(t, normalize);
  return Vec2(-tangent.y(), tangent.x());
}

ConstCurvePtr Curve::getDerivative() const
{
  if (!cached_derivative_)
  {
    Eigen::MatrixX2d new_points = (N_ - 1) * (control_points_.bottomRows(N_ - 1) - control_points_.topRows(N_ - 1));
    (const_cast<Curve*>(this))->cached_derivative_ = std::make_shared<const Curve>(new_points);
  }
  return cached_derivative_;
}

Point Curve::getDerivativeAt(double t) const { return getDerivative()->valueAt(t); }

ConstCurvePtr Curve::getDerivative(uint n) const
{
  if (n == 0)
    throw std::invalid_argument{"Parameter 'n' cannot be zero."};
  ConstCurvePtr nth_derivative = getDerivative();
  for (uint k = 1; k < n; k++)
    nth_derivative = nth_derivative->getDerivative();
  return nth_derivative;
}

Point Curve::getDerivativeAt(uint n, double t) const { return getDerivative(n)->valueAt(t); }

PointVector Curve::getRoots(double step, double epsilon, std::size_t max_iter) const
{
  if (!cached_roots_ || cached_roots_params_ != std::make_tuple(step, epsilon, max_iter))
  {
    (const_cast<Curve*>(this))->cached_roots_params_ = {step, epsilon, max_iter};
    (const_cast<Curve*>(this))->cached_roots_ = std::make_shared<PointVector>();
    std::vector<double> added_t;

    // check both axes
    for (uint k = 0; k < 2; k++)
    {
      double t = 0;
      while (t <= 1.0)
      {
        double t_newton = t;
        std::size_t current_iter = 0;

        // it has to converge in max_iter steps
        while (current_iter < max_iter)
        {
          // Newton-Raphson
          double f = getDerivativeAt(t_newton)[k];
          double f_d = getDerivativeAt(2, t_newton)[k];
          t_newton -= f / f_d;
          // if there is no change to t_current
          if (std::fabs(f) < epsilon)
          {
            // check if between [0, 1]
            if (t_newton >= 0.0 && t_newton <= 1.0)
            {
              // check if same value wasn't found before
              if (added_t.end() == std::find_if(added_t.begin(), added_t.end(), [t_newton, epsilon](const double& val) {
                    return std::fabs(val - t_newton) < epsilon;
                  }))
              {
                // add new value and point
                added_t.push_back(t_newton);
                (const_cast<Curve*>(this))->cached_roots_->push_back(valueAt(t_newton));
              }
            }

            // this t_newton converged
            break;
          }

          current_iter++;
        }

        t += step;
      }
    }
  }
  return *cached_roots_;
}

BBox Curve::getBBox(bool use_roots) const
{
  if (!(use_roots ? cached_bounding_box_tight_ : cached_bounding_box_relaxed_))
  {
    PointVector extremes;
    if (use_roots)
    {
      extremes = getRoots();
      extremes.push_back(control_points_.row(0));
      extremes.push_back(control_points_.row(N_ - 1));
    }
    else
    {
      for (uint k = 0; k < control_points_.rows(); k++)
        extremes.push_back(control_points_.row(k));
    }

    // find mininum and maximum along each axis
    auto x_extremes = std::minmax_element(extremes.begin(), extremes.end(),
                                          [](const Point& lhs, const Point& rhs) { return lhs.x() < rhs.x(); });
    auto y_extremes = std::minmax_element(extremes.begin(), extremes.end(),
                                          [](const Point& lhs, const Point& rhs) { return lhs.y() < rhs.y(); });
    (use_roots ? (const_cast<Curve*>(this))->cached_bounding_box_tight_
               : (const_cast<Curve*>(this))->cached_bounding_box_relaxed_) =
        std::make_shared<BBox>(Point(x_extremes.first->x(), y_extremes.first->y()),
                               Point(x_extremes.second->x(), y_extremes.second->y()));
  }
  return *(use_roots ? cached_bounding_box_tight_ : cached_bounding_box_relaxed_);
}

std::pair<Curve, Curve> Curve::splitCurve(double z) const
{
  return std::make_pair(Curve(splittingCoeffsLeft(z) * control_points_),
                        Curve(splittingCoeffsRight(z) * control_points_));
}

PointVector Curve::getPointsOfIntersection(const Curve& curve, bool stop_at_first, double epsilon) const
{
  PointVector points_of_intersection;
  std::vector<std::pair<ConstCurvePtr, ConstCurvePtr>> subcurve_pairs;

  if (this != &curve)
  {
    // we don't know if shared_ptr of "this" and "curve" exists
    // if not, make temporary copies and use them
    subcurve_pairs.push_back(std::make_pair(std::shared_ptr<Curve>(), std::shared_ptr<Curve>()));
    try // for "this"
    {
      subcurve_pairs.front().first = this->shared_from_this();
    }
    catch (std::bad_weak_ptr const&)
    {
      subcurve_pairs.front().first = std::make_shared<const Curve>(*this);
    }
    try // for "curve"
    {
      subcurve_pairs.front().second = curve.shared_from_this();
    }
    catch (std::bad_weak_ptr const&)
    {
      subcurve_pairs.front().second = std::make_shared<const Curve>(curve);
    }
  }
  else
  {
    // self intersections

    // get all inflection points (roots)
    auto roots = getRoots();
    std::map<double, Point> t_point_pair;
    for (const auto& root : roots)
      t_point_pair.insert(std::make_pair(projectPoint(root), root));

    // divide curve into subcurves at inflection points
    std::vector<ConstCurvePtr> subcurves;
    for (auto& root_pair : t_point_pair)
    {
      if (subcurves.empty())
      {
        subcurves.push_back(
            std::make_shared<Curve>(splittingCoeffsLeft(root_pair.first - epsilon / 2) * control_points_));
        subcurves.push_back(
            std::make_shared<Curve>(splittingCoeffsRight(root_pair.first + epsilon / 2) * control_points_));
      }
      else
      {
        double new_t = subcurves.back()->projectPoint(root_pair.second);
        auto new_cp = subcurves.back()->control_points_;
        subcurves.pop_back();
        subcurves.push_back(std::make_shared<Curve>(splittingCoeffsLeft(new_t - epsilon / 2) * new_cp));
        subcurves.push_back(std::make_shared<Curve>(splittingCoeffsRight(new_t + epsilon / 2) * new_cp));
      }
    }

    // create all pairs of subcurves
    for (uint k = 0; k < subcurves.size(); k++)
      for (uint i = k + 1; i < subcurves.size(); i++)
        subcurve_pairs.push_back(std::make_pair(subcurves.at(k), subcurves.at(i)));
  }

  while (!subcurve_pairs.empty())
  {
    ConstCurvePtr part_a = std::get<0>(subcurve_pairs.back());
    ConstCurvePtr part_b = std::get<1>(subcurve_pairs.back());
    subcurve_pairs.pop_back();

    BBox bbox1 = part_a->getBBox(false); // very slow with tight BBox (roots)
    BBox bbox2 = part_b->getBBox(false); // very slow with tight BBox (roots)
    if (!bbox1.intersects(bbox2))
    {
      // no intersection
      continue;
    }

    if (bbox1.diagonal().norm() < epsilon && bbox2.diagonal().norm() < epsilon)
    {
      // segments converged, check if not already found and add new
      Point new_point = part_a->valueAt(0.5);
      if (points_of_intersection.end() ==
          std::find_if(points_of_intersection.begin(), points_of_intersection.end(),
                       [new_point, epsilon](const Point& point) { return (point - new_point).norm() < epsilon; }))
      {
        points_of_intersection.push_back(new_point);

        // if only first point is needed, stop
        if (stop_at_first)
          return points_of_intersection;
      }
      continue;
    }

    // intersection exists, but segments are still too large
    // divide both segments in half and new pairs
    // LIFO : we want to first discover closest intersection (smallest t on this curve)
    // so it is important which pair of subcurves is inserted first
    std::vector<ConstCurvePtr> subcurves_a;
    std::vector<ConstCurvePtr> subcurves_b;

    if (bbox1.diagonal().norm() < epsilon)
    {
      // if small enough, do not divide it further
      subcurves_a.push_back(part_a);
    }
    else
    {
      // divide into two subcurves
      // first insert 2nd subcurve t = [0.5 to 1]
      subcurves_a.push_back(std::make_shared<Curve>(part_a->splittingCoeffsRight() * part_a->control_points_));
      subcurves_a.push_back(std::make_shared<Curve>(part_a->splittingCoeffsLeft() * part_a->control_points_));
    }

    if (bbox2.diagonal().norm() < epsilon)
    {
      // if small enough, do not divide it further
      // first insert 2nd subcurve t = [0.5 to 1]
      subcurves_b.push_back(part_b);
    }
    else
    {
      // divide into two subcurves
      subcurves_b.push_back(std::make_shared<Curve>(part_b->splittingCoeffsRight() * part_b->control_points_));
      subcurves_b.push_back(std::make_shared<Curve>(part_b->splittingCoeffsLeft() * part_b->control_points_));
    }

    // insert all combinations for next iteration
    // last pair is one where both subcurves have smalles t ranges
    for (auto&& subcurve_b : subcurves_b)
      for (auto&& subcurve_a : subcurves_a)
        subcurve_pairs.push_back(std::make_pair(subcurve_a, subcurve_b));
  }
  return points_of_intersection;
}

double Curve::projectPoint(const Point& point, double step, double epsilon, std::size_t max_iter) const
{
  double t = 0;
  double t_dist = (valueAt(t) - point).norm();

  // Coarse search
  for (double k = step; k < 1 + step; k += step)
  {
    double new_dist = (valueAt(k) - point).norm();
    if (new_dist < t_dist)
    {
      t_dist = new_dist;
      t = k;
    }
  }

  // Fine search - Newton-Raphson
  // function to minimize is a dot product between projection vector and tangent
  // - projection vector is a vector between point we are projecting and our current guess
  double t_old = t;
  std::size_t current_iter = 0;
  while (current_iter < max_iter)
  {
    Point P = valueAt(t);
    Point d1 = getDerivativeAt(t);
    Point d2 = getDerivativeAt(2, t);
    double f = (point - P).dot(d1);
    double f_d = (point - P).dot(d2) - d1.dot(d1);
    t -= f / f_d;
    if (t < 0 || t > 1)
    {
      return t_old;
    }
    if (std::fabs(f) < epsilon)
      return t;

    current_iter++;
  }
}

void Curve::applyContinuity(const Curve& source_curve, std::vector<double>& beta_coeffs)
{
  ulong c_order = beta_coeffs.size();

  Eigen::MatrixXd pascal_alterating_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  pascal_alterating_matrix.diagonal(-1) = -Eigen::ArrayXd::LinSpaced(c_order, 1, c_order);
  pascal_alterating_matrix = pascal_alterating_matrix.exp();

  Eigen::MatrixXd bell_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  bell_matrix(0, c_order) = 1;

  for (uint i = 0; i < c_order; i++)
  {
    bell_matrix.block(1, c_order - i - 1, i + 1, 1) =
        bell_matrix.block(0, c_order - i, i + 1, i + 1) *
        pascal_alterating_matrix.block(i, 0, 1, i + 1)
            .cwiseAbs()
            .transpose()
            .cwiseProduct(Eigen::Map<Eigen::MatrixXd>(beta_coeffs.data(), i + 1, 1));
  }

  Eigen::MatrixXd factorial_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  for (uint i = 0; i < c_order + 1; i++)
  {
    factorial_matrix(i, i) = factorial(N_ - 1) / factorial(N_ - 1 - i);
  }

  Eigen::MatrixXd derivatives(Eigen::MatrixXd::Zero(2, c_order + 1));
  derivatives.col(0) = source_curve.valueAt(1);
  for (uint i = 1; i < c_order + 1; i++)
    derivatives.col(i) = source_curve.getDerivativeAt(i, 1);

  Eigen::MatrixXd derivatives_wanted = (derivatives * bell_matrix).rowwise().reverse().transpose();

  Eigen::MatrixXd control_points = (factorial_matrix * pascal_alterating_matrix).inverse() * derivatives_wanted;

  for (uint i = 0; i < c_order + 1; i++)
  {
    manipulateControlPoint(i, control_points.row(i));
  }
}

} // namespace Bezier