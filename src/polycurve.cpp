#include "Bezier/polycurve.h"
#include "Bezier/bezier.h"

#include <numeric>
#include <utility>

inline double binomial(unsigned n, unsigned k) { return tgamma(n + 1) / (tgamma(k + 1) * tgamma(n - k + 1)); }

using namespace Bezier;

PolyCurve::PolyCurve(std::deque<std::shared_ptr<Curve>> curve_list) : curves_(std::move(curve_list)) {}

PolyCurve::PolyCurve(std::shared_ptr<Curve>& curve) { curves_.push_back(curve); }

PolyCurve::PolyCurve(std::vector<std::shared_ptr<Curve>>& curve_list)
{
  for (auto& curve_ptr : curve_list)
    insertBack(curve_ptr);
}

PolyCurve::PolyCurve(const PolyCurve& poly_curve) : PolyCurve(poly_curve.curves_) {}

void PolyCurve::insertAt(unsigned idx, std::shared_ptr<Curve>& curve)
{
  Point s_1, s_2, e_1, e_2;
  std::tie(s_1, e_1) = curve->endPoints();
  if (idx > 0) // check with curve before
  {
    std::tie(s_2, e_2) = curves_[idx - 1]->endPoints();

    double s_e = (s_1 - e_2).norm();
    double e_e = (e_1 - e_2).norm();

    if (e_e < s_e) // we need to reverse the curve
    {
      curve->reverse();
      std::tie(s_1, e_1) = curve->endPoints();
    }

    curve->manipulateControlPoint(0, (s_1 + e_2) / 2);
    curves_[idx - 1]->manipulateControlPoint(curves_[idx - 1]->order(), (s_1 + e_2) / 2);
  }
  if (idx + 1 < size()) // check with curve after
  {
    std::tie(s_2, e_2) = curves_[idx]->endPoints();

    double e_s = (e_1 - s_2).norm();
    double s_s = (s_1 - s_2).norm();

    if (s_s < e_s) // we ned to reverse the curve
    {
      curve->reverse();
      std::tie(s_1, e_1) = curve->endPoints();
    }

    curve->manipulateControlPoint(curve->order(), (e_1 + s_2) / 2);
    curves_[idx + 1]->manipulateControlPoint(0, (e_1 + s_2) / 2);
  }

  curves_.insert(curves_.begin() + idx, curve);
}

void PolyCurve::insertFront(std::shared_ptr<Curve>& curve) { insertAt(0, curve); }

void PolyCurve::insertBack(std::shared_ptr<Curve>& curve) { insertAt(size(), curve); }

void PolyCurve::removeAt(unsigned idx)
{
  if (idx == 0)
    removeFirst();
  else if (idx == size() - 1)
    removeBack();
  else
  {
    Point s_1, s_2, e_1, e_2;
    std::tie(s_1, e_1) = curves_[idx - 1]->endPoints();
    std::tie(s_2, e_2) = curves_[idx + 1]->endPoints();
    curves_[idx - 1]->manipulateControlPoint(curves_[idx - 1]->order(), (e_1 + s_2) / 2);
    curves_[idx + 1]->manipulateControlPoint(0, (e_1 + s_2) / 2);
    curves_.erase(curves_.begin() + idx);
  }
}

void PolyCurve::removeFirst() { curves_.pop_front(); }

void PolyCurve::removeBack() { curves_.pop_back(); }

PolyCurve PolyCurve::subPolyCurve(unsigned idx_l, unsigned idx_r) const
{
  return PolyCurve(std::deque<std::shared_ptr<Curve>>(curves_.begin() + idx_l, curves_.begin() + idx_r));
}

unsigned PolyCurve::size() const { return static_cast<unsigned>(curves_.size()); }

unsigned PolyCurve::curveIdx(double t) const
{
  unsigned idx = static_cast<unsigned>(t);
  return idx - (idx == size());
}

std::shared_ptr<Curve> PolyCurve::curvePtr(unsigned idx) const { return curves_[idx]; }

std::vector<std::shared_ptr<Curve>> PolyCurve::curveList() const
{
  return std::vector<std::shared_ptr<Curve>>(curves_.begin(), curves_.end());
}

PointVector PolyCurve::polyline(double smoothness, double precision) const
{
  PointVector polyline;
  for (unsigned k = 0; k < size(); k++)
  {
    auto new_poly = curvePtr(k)->polyline(smoothness, precision);
    polyline.reserve(polyline.size() + new_poly.size() - (k ? 1 : 0));
    polyline.insert(polyline.end(), new_poly.begin() + (k ? 1 : 0), new_poly.end());
  }
  return polyline;
}

double PolyCurve::length() const { return length(0, size()); }

double PolyCurve::length(double t) const { return length(0, t); }

double PolyCurve::length(double t1, double t2) const
{
  unsigned idx1 = curveIdx(t1);
  unsigned idx2 = curveIdx(t2);

  if (idx1 == idx2)
    return curves_[idx1]->length(t1 - idx1, t2 - idx2);
  else if (idx1 + 1 == idx2)
    return curves_[idx1]->length(t1 - idx1, 1.0) + curves_[idx2]->length(0.0, t2 - idx2);
  else
    return std::accumulate(begin(curves_) + idx1 + 1, begin(curves_) + idx2,
                           curves_[idx1]->length(t1 - idx1, 1.0) + curves_[idx2]->length(0.0, t2 - idx2),
                           [](double sum, std::shared_ptr<Curve> curve) { return sum + curve->length(); });
}

double PolyCurve::iterateByLength(double t, double s, double epsilon, std::size_t max_iter) const
{
  double s_t = length(t);
  //  if (s_t + s < 0 || s_t + s > length())
  //    throw std::out_of_range{"Resulting parameter t not in [0, n] range."};
  if (s_t + s < 0)
    return 0;
  if (s_t + s > length())
    return size();

  unsigned idx = curveIdx(t);
  t -= idx;
  s_t -= length(idx);

  bool first_iteration = true;
  while (curvePtr(idx)->length() - s_t < s)
  {
    s -= curvePtr(idx++)->length() - s_t;
    if (first_iteration)
    {
      first_iteration = false;
      s_t = 0;
      t = 0;
    }
  }

  return idx + curvePtr(idx)->iterateByLength(t, s, epsilon, max_iter);
}

std::pair<Point, Point> PolyCurve::endPoints() const
{
  return std::make_pair(curves_.front()->endPoints().first, curves_.back()->endPoints().second);
}

PointVector PolyCurve::controlPoints() const
{
  PointVector cp;
  for (auto& curve_ptr : curves_)
  {
    auto cp_c = curve_ptr->controlPoints();
    cp.reserve(cp.size() + cp_c.size());
    cp.insert(cp.end(), cp_c.begin(), cp_c.end());
  }
  return cp;
}

void PolyCurve::manipulateControlPoint(unsigned idx, const Point& point)
{
  for (auto& curve_ptr : curves_)
    if (idx <= curve_ptr->order())
    {
      curve_ptr->manipulateControlPoint(idx, point);
      break;
    }
    else
      --idx -= curve_ptr->order();
}

Point PolyCurve::valueAt(double t) const
{
  unsigned idx = curveIdx(t);
  return curvePtr(idx)->valueAt(t - idx);
}

double PolyCurve::curvatureAt(double t) const
{
  unsigned idx = curveIdx(t);
  return curvePtr(idx)->curvatureAt(t - idx);
}

double PolyCurve::curvatureDerivativeAt(double t) const
{
  unsigned idx = curveIdx(t);
  return curvePtr(idx)->curvatureDerivativeAt(t - idx);
}

Vector PolyCurve::tangentAt(double t, bool normalize) const
{
  unsigned idx = curveIdx(t);
  return curvePtr(idx)->tangentAt(t - idx, normalize);
}

Vector PolyCurve::normalAt(double t, bool normalize) const
{
  unsigned idx = curveIdx(t);
  return curvePtr(idx)->normalAt(t - idx, normalize);
}

Point PolyCurve::derivativeAt(double t) const
{
  unsigned idx = curveIdx(t);
  return curvePtr(idx)->derivativeAt(t - idx);
}

Point PolyCurve::derivativeAt(unsigned n, double t) const
{
  unsigned idx = curveIdx(t);
  return curvePtr(idx)->derivativeAt(n, t - idx);
}

BoundingBox PolyCurve::boundingBox(bool use_roots) const
{
  BoundingBox bbox;
  for (auto& curve_ptr : curves_)
    bbox.extend(curve_ptr->boundingBox(use_roots));
  return bbox;
}

// namespace is a workaround for a bug on old gcc versions:
// https://stackoverflow.com/questions/25311512/specialization-of-template-in-different-namespace
namespace Bezier
{

template <>
PointVector PolyCurve::pointsOfIntersection<Curve>(const Curve& curve, bool stop_at_first, double epsilon) const
{
  PointVector points;
  for (auto& curve_ptr : curves_)
  {
    auto new_points = curve_ptr->pointsOfIntersection(curve, stop_at_first, epsilon);
    points.reserve(points.size() + new_points.size());
    points.insert(points.end(), new_points.begin(), new_points.end());
    if (!points.empty() && stop_at_first)
      break;
  }
  return points;
}

template <>
PointVector PolyCurve::pointsOfIntersection<PolyCurve>(const PolyCurve& poly_curve, bool stop_at_first,
                                                       double epsilon) const
{
  PointVector points;
  for (auto& curve_ptr : curves_)
  {
    auto new_points = poly_curve.pointsOfIntersection(*curve_ptr, stop_at_first, epsilon);
    points.reserve(points.size() + new_points.size());
    points.insert(points.end(), new_points.begin(), new_points.end());
    if (!points.empty() && stop_at_first)
      break;
  }
  return points;
}

} // namespace Bezier

double PolyCurve::projectPoint(const Point& point, double step, double epsilon) const
{
  double min_t = curves_.front()->projectPoint(point, step, epsilon);
  double min_dist = (point - curves_.front()->valueAt(min_t)).norm();

  for (unsigned k = 1; k < size(); k++)
  {
    double t = curves_[k]->projectPoint(point, step, epsilon);
    double dist = (point - curves_[k]->valueAt(t)).norm();
    if (dist < min_dist)
    {
      min_dist = dist;
      min_t = k + t;
    }
  }
  return min_t;
}
