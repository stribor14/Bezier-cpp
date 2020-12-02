#include "Bezier/bezier.h"
#include "Bezier/legendre_gauss.h"

#include <numeric>

#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/Polynomials>

using namespace Bezier;

inline double factorial(uint k) { return std::tgamma(k + 1); }
inline double binomial(uint n, uint k) { return factorial(n) / (factorial(k) * factorial(n - k)); }
inline Eigen::VectorXd trimZeroes(const Eigen::VectorXd& vec)
{
  auto idx = vec.size();
  while (vec(idx - 1) == 0.0)
    --idx;
  return vec.head(idx);
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
    control_points_.row(k) = points[k];
}

Curve::Curve(const Curve& curve) : Curve(curve.controlPoints()) {}

uint Curve::order() { return N_ - 1; }

PointVector Curve::controlPoints() const
{
  PointVector points(N_);
  for (uint k = 0; k < N_; k++)
    points[k] = control_points_.row(k);
  return points;
}

std::pair<Point, Point> Curve::endPoints() const { return {control_points_.row(0), control_points_.row(N_ - 1)}; }

PointVector Curve::polyline(double smoothness, double precision) const
{
  if (!cached_polyline_ || cached_polyline_params_ != std::make_pair(smoothness, precision))
  {
    auto* polyline = new PointVector;
    std::vector<Eigen::MatrixX2d> subcurves;
    subcurves.emplace_back(control_points_);
    polyline->emplace_back(control_points_.row(0));
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
        polyline->emplace_back(cp.row(N_ - 1));
      }
      else
      {
        subcurves.emplace_back(splittingCoeffsRight(N_) * cp);
        subcurves.emplace_back(splittingCoeffsLeft(N_) * cp);
      }
    }

    (const_cast<Curve*>(this))->cached_polyline_params_ = {smoothness, precision};
    (const_cast<Curve*>(this))->cached_polyline_.reset(polyline);
  }
  return *cached_polyline_;
}

double Curve::length() const { return length(0.0, 1.0); }

double Curve::length(Parameter t) const { return length(0.0, t); }

double Curve::length(Parameter t1, Parameter t2) const
{
  double sum = 0;

  for (uint k = 0; k < LegendreGauss::N; k++)
    sum += LegendreGauss::weights[k] * derivativeAt(LegendreGauss::abcissae[k] * (t2 - t1) / 2 + (t1 + t2) / 2).norm();

  return sum * (t2 - t1) / 2;
}

double Curve::iterateByLength(Parameter t, double s, double epsilon) const
{
  const double s_t = length(t);

  if (s_t + s < 0)
    return 0;
  if (s_t + s > length())
    return 1;

  double f = (length(t) - s_t - s);
  while (std::fabs(f) > epsilon)
  {
    // Halley
    f = (length(t) - s_t - s);
    double f_d = derivativeAt(t).norm();
    double f_d2 = derivativeAt(2, t).norm();

    t -= (2 * f * f_d) / (2 * f_d * f_d - f * f_d2);
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

void Curve::manipulateCurvature(Parameter t, const Point& point)
{
  if (N_ < 3 || N_ > 4)
    throw std::logic_error{"Only quadratic and cubic curves can be manipulated"};

  double r =
      std::fabs((std::pow(t, N_ - 1) + std::pow(1 - t, N_ - 1) - 1) / (std::pow(t, N_ - 1) + std::pow(1 - t, N_ - 1)));
  double u = std::pow(1 - t, N_ - 1) / (std::pow(t, N_ - 1) + std::pow(1 - t, N_ - 1));
  Point C = u * control_points_.row(0) + (1 - u) * control_points_.row(N_ - 1);
  const Point& B = point;
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
  control_points_ = elevateOrderCoeffs(N_++) * control_points_;
  resetCache();
}

void Curve::lowerOrder()
{
  if (N_ == 2)
    throw std::logic_error{"Cannot further reduce the order of curve."};
  control_points_ = lowerOrderCoeffs(N_--) * control_points_;
  resetCache();
}

Point Curve::valueAt(Parameter t) const
{
  if (N_ == 0)
    return {0, 0};
  Eigen::VectorXd power_basis = Eigen::pow(t, Eigen::ArrayXd::LinSpaced(N_, 0, N_ - 1));
  return (power_basis.transpose() * bernsteinCoeffs(N_) * control_points_).transpose();
}

PointVector Curve::valueAt(ParameterVector t_vector) const
{
  PointVector points;
  points.reserve(t_vector.size());

  auto t_matrix = Eigen::Map<Eigen::VectorXd>(t_vector.data(), static_cast<int>(t_vector.size())).replicate(1, N_);
  auto p_matrix = Eigen::ArrayXd::LinSpaced(N_, 0, N_ - 1).transpose().replicate(static_cast<int>(t_vector.size()), 1);
  Eigen::MatrixXd power_basis = Eigen::pow(t_matrix.array(), p_matrix.array());
  Eigen::MatrixXd points_eigen = (power_basis * bernsteinCoeffs(N_) * control_points_);

  for (uint k = 0; k < points_eigen.rows(); k++)
    points.emplace_back(points_eigen.row(k));

  return points;
}

double Curve::curvatureAt(Parameter t) const
{
  Point d1 = derivativeAt(t);
  Point d2 = derivativeAt(2, t);

  return (d1.x() * d2.y() - d1.y() * d2.x()) / std::pow(d1.norm(), 3);
}

double Curve::curvatureDerivativeAt(Parameter t) const
{
  auto d1 = derivativeAt(t);
  auto d2 = derivativeAt(2, t);
  auto d3 = derivativeAt(3, t);

  return (d1.x() * d3.y() - d1.y() * d3.x()) / std::pow(d1.norm(), 3) -
         3 * d1.dot(d2) * (d1.x() * d2.y() - d1.y() * d2.x()) / std::pow(d1.norm(), 5);
}

Vector Curve::tangentAt(Parameter t, bool normalize) const
{
  Point p(derivativeAt(t));
  if (normalize && p.norm() > 0)
    p.normalize();
  return p;
}

Vector Curve::normalAt(Parameter t, bool normalize) const
{
  Point tangent = tangentAt(t, normalize);
  return {-tangent.y(), tangent.x()};
}

std::shared_ptr<const Curve> Curve::derivative() const
{
  if (!cached_derivative_)
  {
    (const_cast<Curve*>(this))->cached_derivative_ =
        N_ == 1 ? std::make_shared<const Curve>(PointVector{Point(0, 0)})
                : std::make_shared<const Curve>(
                      ((N_ - 1) * (control_points_.bottomRows(N_ - 1) - control_points_.topRows(N_ - 1))).eval());
  }
  return cached_derivative_;
}

std::shared_ptr<const Curve> Curve::derivative(uint n) const
{
  if (n == 0)
    throw std::invalid_argument{"Parameter 'n' cannot be zero."};
  std::shared_ptr<const Curve> nth_derivative = derivative();
  for (uint k = 1; k < n; k++)
    nth_derivative = nth_derivative->derivative();
  return nth_derivative;
}

Vector Curve::derivativeAt(Parameter t) const { return derivative()->valueAt(t); }

Vector Curve::derivativeAt(uint n, Parameter t) const { return derivative(n)->valueAt(t); }

ParameterVector Curve::roots() const
{
  if (!cached_roots_)
  {
    ParameterVector* roots = new ParameterVector();
    if (N_ > 1)
    {
      std::vector<double> roots_X, roots_Y;
      Eigen::MatrixXd bezier_polynomial = bernsteinCoeffs(N_) * control_points_;
      Eigen::PolynomialSolver<double, Eigen::Dynamic> poly_solver;
      poly_solver.compute(trimZeroes(bezier_polynomial.col(0)));
      poly_solver.realRoots(roots_X);
      poly_solver.compute(trimZeroes(bezier_polynomial.col(1)));
      poly_solver.realRoots(roots_Y);
      roots->reserve(roots_X.size() + roots_Y.size());
      std::copy_if(std::make_move_iterator(roots_X.begin()), std::make_move_iterator(roots_X.end()),
                   std::back_inserter(*roots), [](double t) { return t >= 0 && t <= 1; });
      std::copy_if(std::make_move_iterator(roots_Y.begin()), std::make_move_iterator(roots_Y.end()),
                   std::back_inserter(*roots), [](double t) { return t >= 0 && t <= 1; });
    }
    const_cast<Curve*>(this)->cached_roots_.reset(roots);
  }
  return *cached_roots_;
}

ParameterVector Curve::extrema() const { return derivative()->roots(); }

BoundingBox Curve::boundingBox() const
{
  if (!cached_bounding_box_)
  {
    PointVector extremes = valueAt(extrema());

    extremes.emplace_back(control_points_.row(0));
    extremes.emplace_back(control_points_.row(N_ - 1));

    // find mininum and maximum along each axis
    auto x_extremes = std::minmax_element(extremes.begin(), extremes.end(),
                                          [](const Point& lhs, const Point& rhs) { return lhs.x() < rhs.x(); });
    auto y_extremes = std::minmax_element(extremes.begin(), extremes.end(),
                                          [](const Point& lhs, const Point& rhs) { return lhs.y() < rhs.y(); });
    const_cast<Curve*>(this)->cached_bounding_box_.reset(new BoundingBox(
        Point(x_extremes.first->x(), y_extremes.first->y()), Point(x_extremes.second->x(), y_extremes.second->y())));
  }
  return *cached_bounding_box_;
}

std::pair<Curve, Curve> Curve::splitCurve(double z) const
{
  return {Curve(splittingCoeffsLeft(N_, z) * control_points_), Curve(splittingCoeffsRight(N_, z) * control_points_)};
}

PointVector Curve::intersection(const Curve& curve, bool stop_at_first, double epsilon) const
{
  PointVector points_of_intersection;

  std::vector<std::pair<Eigen::MatrixX2d, Eigen::MatrixX2d>> subcurve_pairs;

  if (this != &curve)
  {
    subcurve_pairs.emplace_back(control_points_, curve.control_points_);
  }
  else
  {
    // self intersections

    // get all extrema of the curve
    std::map<double, Point> t_point_pair;
    for (const auto& root : valueAt(extrema()))
      t_point_pair.insert({projectPoint(root), root});

    // divide curve into subcurves at roots
    std::vector<Eigen::MatrixX2d> subcurves;
    for (const auto& root_pair : t_point_pair)
    {
      if (subcurves.empty())
      {
        subcurves.emplace_back(splittingCoeffsLeft(N_, root_pair.first - epsilon / 2) * control_points_);
        subcurves.emplace_back(splittingCoeffsRight(N_, root_pair.first + epsilon / 2) * control_points_);
      }
      else
      {
        Curve temp_curve(subcurves.back());
        double new_t = temp_curve.projectPoint(root_pair.second);
        auto new_cp = subcurves.back();
        subcurves.pop_back();
        subcurves.emplace_back(splittingCoeffsLeft(N_, new_t - epsilon / 2) * new_cp);
        subcurves.emplace_back(splittingCoeffsRight(N_, new_t + epsilon / 2) * new_cp);
      }
    }

    // create all pairs of subcurves
    for (uint k = 0; k < subcurves.size(); k++)
      for (uint i = k + 1; i < subcurves.size(); i++)
        subcurve_pairs.emplace_back(subcurves[k], subcurves[i]);
  }

  auto bbox = [](Eigen::MatrixX2d cp) {
    return BoundingBox(Point(cp.col(0).minCoeff(), cp.col(1).minCoeff()),
                       Point(cp.col(0).maxCoeff(), cp.col(1).maxCoeff()));
  };

  while (!subcurve_pairs.empty())
  {
    Eigen::MatrixX2d part_a = std::get<0>(subcurve_pairs.back());
    Eigen::MatrixX2d part_b = std::get<1>(subcurve_pairs.back());
    subcurve_pairs.pop_back();

    BoundingBox bbox1 = bbox(part_a);
    BoundingBox bbox2 = bbox(part_b);
    if (!bbox1.intersects(bbox2))
    {
      // no intersection
      continue;
    }

    if (bbox1.diagonal().norm() < epsilon && bbox2.diagonal().norm() < epsilon)
    {
      // segments converged, check if not already found and add new
      Point new_point = bbox1.center();
      if (points_of_intersection.end() ==
          std::find_if(points_of_intersection.begin(), points_of_intersection.end(),
                       [new_point, epsilon](const Point& point) { return (point - new_point).norm() < epsilon; }))
      {
        points_of_intersection.emplace_back(new_point);

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
    std::vector<Eigen::MatrixX2d> subcurves_a;
    std::vector<Eigen::MatrixX2d> subcurves_b;

    if (bbox1.diagonal().norm() < epsilon)
    {
      // if small enough, do not divide it further
      subcurves_a.emplace_back(part_a);
    }
    else
    {
      // divide into two subcurves
      // first insert 2nd subcurve t = [0.5 to 1]
      subcurves_a.emplace_back(splittingCoeffsRight(N_) * part_a);
      subcurves_a.emplace_back(splittingCoeffsLeft(N_) * part_a);
    }

    if (bbox2.diagonal().norm() < epsilon)
    {
      // if small enough, do not divide it further
      // first insert 2nd subcurve t = [0.5 to 1]
      subcurves_b.emplace_back(part_b);
    }
    else
    {
      // divide into two subcurves
      subcurves_b.emplace_back(curve.splittingCoeffsRight(N_) * part_b);
      subcurves_b.emplace_back(curve.splittingCoeffsLeft(N_) * part_b);
    }

    // insert all combinations for next iteration
    // last pair is one where both subcurves have smalles t ranges
    for (auto& subcurve_b : subcurves_b)
      for (auto& subcurve_a : subcurves_a)
        subcurve_pairs.emplace_back(subcurve_a, subcurve_b);
  }

  return points_of_intersection;
}

Parameter Curve::projectPoint(const Point& point) const
{
  if (!cached_projection_polynomial_part_)
  {
    Eigen::MatrixXd curve_polynomial = (bernsteinCoeffs(N_) * control_points_);
    Eigen::MatrixXd derivate_polynomial = (bernsteinCoeffs(N_ - 1) * derivative()->control_points_);

    Eigen::VectorXd polynomial_part = Eigen::VectorXd::Zero(curve_polynomial.rows() + derivate_polynomial.rows() - 1);
    for (uint k = 0; k < curve_polynomial.rows(); k++)
      polynomial_part.middleRows(k, derivate_polynomial.rows()) +=
          derivate_polynomial * curve_polynomial.row(k).transpose();

    const_cast<Curve*>(this)->cached_projection_polynomial_part_.reset(new Eigen::VectorXd(std::move(polynomial_part)));
    const_cast<Curve*>(this)->cached_projection_polynomial_derivative_ = std::move(derivate_polynomial);
  }

  Eigen::VectorXd polynomial = *cached_projection_polynomial_part_;
  polynomial.topRows(cached_projection_polynomial_derivative_.rows()) -=
      cached_projection_polynomial_derivative_ * point;

  std::vector<double> candidates;
  Eigen::PolynomialSolver<double, Eigen::Dynamic> poly_solver(trimZeroes(polynomial));
  poly_solver.realRoots(candidates);

  double projection = (point - valueAt(0.0)).norm() < (point - valueAt(1.0)).norm() ? 0.0 : 1.0;
  double min = (point - valueAt(projection)).norm();

  for (auto candidate : candidates)
  {
    if (candidate < 0 || candidate > 1)
      continue;

    double dist = (point - valueAt(candidate)).norm();
    if (dist < min)
    {
      projection = candidate;
      min = dist;
    }
  }
  return projection;
}

ParameterVector Curve::projectPoint(PointVector point_vector) const
{
  ParameterVector t_vector;
  t_vector.reserve(point_vector.size());
  for (auto point : point_vector)
    t_vector.emplace_back(projectPoint(point));
  return t_vector;
}

double Curve::distance(const Point& point) const { return (point - valueAt(projectPoint(point))).norm(); }

std::vector<double> Curve::distance(PointVector point_vector) const
{
  std::vector<double> dist_vector;
  dist_vector.reserve(point_vector.size());
  for (auto t : point_vector)
    dist_vector.emplace_back(distance(t));
  return dist_vector;
}

void Curve::applyContinuity(const Curve& source_curve, const std::vector<double>& beta_coeffs)
{
  uint c_order = static_cast<uint>(beta_coeffs.size());

  Eigen::MatrixXd pascal_alterating_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  pascal_alterating_matrix.diagonal(-1) = -Eigen::ArrayXd::LinSpaced(c_order, 1, c_order);
  pascal_alterating_matrix = pascal_alterating_matrix.exp();

  Eigen::MatrixXd bell_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  bell_matrix(0, c_order) = 1;
  for (uint i = 0; i < c_order; i++)
    bell_matrix.block(1, c_order - i - 1, i + 1, 1) =
        bell_matrix.block(0, c_order - i, i + 1, i + 1) *
        pascal_alterating_matrix.block(i, 0, 1, i + 1)
            .cwiseAbs()
            .transpose()
            .cwiseProduct(Eigen::Map<const Eigen::MatrixXd>(beta_coeffs.data(), i + 1, 1));

  Eigen::MatrixXd factorial_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  for (uint i = 0; i < c_order + 1; i++)
    factorial_matrix(i, i) = factorial(N_ - 1) / factorial(N_ - 1 - i);

  Eigen::MatrixXd derivatives;
  derivatives.resize(2, c_order + 1);
  derivatives.col(0) = source_curve.control_points_.bottomRows(1);
  for (uint i = 1; i < c_order + 1; i++)
    derivatives.col(i) = source_curve.derivative(i)->control_points_.bottomRows(1);

  Eigen::MatrixXd derivatives_wanted = (derivatives * bell_matrix).rowwise().reverse().transpose();

  control_points_.topRows(c_order + 1).noalias() =
      (factorial_matrix * pascal_alterating_matrix).inverse() * derivatives_wanted;
  resetCache();
}

void Curve::resetCache()
{
  cached_derivative_.reset();
  cached_roots_.reset();
  cached_bounding_box_.reset();
  cached_polyline_.reset();
  cached_projection_polynomial_part_.reset();
}

Curve::CoeffsMap Curve::bernstein_coeffs_ = CoeffsMap();
Curve::CoeffsMap Curve::splitting_coeffs_left_ = CoeffsMap();
Curve::CoeffsMap Curve::splitting_coeffs_right_ = CoeffsMap();
Curve::CoeffsMap Curve::elevate_order_coeffs_ = CoeffsMap();
Curve::CoeffsMap Curve::lower_order_coeffs_ = CoeffsMap();

Curve::Coeffs Curve::bernsteinCoeffs(uint n)
{
  if (bernstein_coeffs_.find(n) == bernstein_coeffs_.end())
  {
    bernstein_coeffs_.insert({n, Coeffs::Zero(n, n)});
    bernstein_coeffs_[n].diagonal(-1) = -Eigen::ArrayXd::LinSpaced(n - 1, 1, n - 1);
    bernstein_coeffs_[n] = bernstein_coeffs_[n].exp();
    for (uint k = 0; k < n; k++)
      bernstein_coeffs_[n].row(k) *= binomial(n - 1, k);
  }
  return bernstein_coeffs_[n];
}

Curve::Coeffs Curve::splittingCoeffsLeft(uint n, Parameter z)
{
  if (z == 0.5)
  {
    if (splitting_coeffs_left_.find(n) == splitting_coeffs_left_.end())
    {
      splitting_coeffs_left_.insert({n, Coeffs::Zero(n, n)});
      splitting_coeffs_left_[n].diagonal() = Eigen::pow(0.5, Eigen::ArrayXd::LinSpaced(n, 0, n - 1));
      splitting_coeffs_left_[n] = bernsteinCoeffs(n).inverse() * splitting_coeffs_left_[n] * bernsteinCoeffs(n);
    }
    return splitting_coeffs_left_[n];
  }
  else
  {
    Curve::Coeffs coeffs(Coeffs::Zero(n, n));
    coeffs.diagonal() = Eigen::pow(z, Eigen::ArrayXd::LinSpaced(n, 0, n - 1));
    coeffs = bernsteinCoeffs(n).inverse() * coeffs * bernsteinCoeffs(n);
    return coeffs;
  }
}

Curve::Coeffs Curve::splittingCoeffsRight(uint n, Parameter z)
{
  if (z == 0.5)
  {
    if (splitting_coeffs_right_.find(n) == splitting_coeffs_right_.end())
    {
      splitting_coeffs_right_.insert({n, Coeffs::Zero(n, n)});
      Curve::Coeffs temp_splitting_coeffs_left = splittingCoeffsLeft(n);
      for (uint k = 0; k < n; k++)
        splitting_coeffs_right_[n].block(k, k, 1, n - k) = temp_splitting_coeffs_left.block(n - 1 - k, 0, 1, n - k);
    }
    return splitting_coeffs_right_[n];
  }
  else
  {
    Curve::Coeffs coeffs(Coeffs::Zero(n, n));
    Curve::Coeffs temp_splitting_coeffs_left = splittingCoeffsLeft(n, z);
    for (uint k = 0; k < n; k++)
      coeffs.block(k, k, 1, n - k) = temp_splitting_coeffs_left.block(n - 1 - k, 0, 1, n - k);
    return coeffs;
  }
}

Curve::Coeffs Curve::elevateOrderCoeffs(uint n)
{
  if (elevate_order_coeffs_.find(n) == elevate_order_coeffs_.end())
  {
    elevate_order_coeffs_.insert({n, Coeffs::Zero(n + 1, n)});
    elevate_order_coeffs_[n].diagonal() = 1 - Eigen::ArrayXd::LinSpaced(n, 0, n - 1) / n;
    elevate_order_coeffs_[n].diagonal(-1) = Eigen::ArrayXd::LinSpaced(n, 1, n) / n;
  }
  return elevate_order_coeffs_[n];
}

Curve::Coeffs Curve::lowerOrderCoeffs(uint n)
{
  if (lower_order_coeffs_.find(n) == lower_order_coeffs_.end())
  {
    lower_order_coeffs_.insert({n, Coeffs::Zero(n - 1, n)});
    lower_order_coeffs_[n].noalias() = (elevateOrderCoeffs(n - 1).transpose() * elevateOrderCoeffs(n - 1)).inverse() *
                                       elevateOrderCoeffs(n - 1).transpose();
  }
  return lower_order_coeffs_[n];
}
