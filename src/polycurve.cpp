#include "Bezier/polycurve.h"
#include "Bezier/bezier.h"

#include <numeric>
#include <utility>

inline double binomial(uint n, uint k) { return tgamma(n + 1) / (tgamma(k + 1) * tgamma(n - k + 1)); }

using namespace Bezier;

PolyCurve::PolyCurve(std::deque<std::shared_ptr<Curve>> curve_list) : curves_(std::move(curve_list)) {}

PolyCurve::PolyCurve(std::shared_ptr<Curve>& curve) { curves_.push_back(curve); }

PolyCurve::PolyCurve(std::vector<std::shared_ptr<Curve>>& curve_list)
{
  for (auto& curve_ptr : curve_list)
    insertBack(curve_ptr);
}

PolyCurve::PolyCurve(const PolyCurve& poly_curve) : PolyCurve(poly_curve.curves_) {}

void PolyCurve::insertAt(uint idx, std::shared_ptr<Curve>& curve)
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

void PolyCurve::removeAt(uint idx)
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

PolyCurve PolyCurve::subPolyCurve(uint idx_l, uint idx_r) const
{
  return PolyCurve(std::deque<std::shared_ptr<Curve>>(curves_.begin() + idx_l, curves_.begin() + idx_r));
}

uint PolyCurve::size() const { return static_cast<uint>(curves_.size()); }

uint PolyCurve::curveIdx(Parameter t) const
{
  uint idx = static_cast<uint>(t);
  return idx - (idx == size());
}

std::shared_ptr<Curve> PolyCurve::curvePtr(uint idx) const { return curves_[idx]; }

std::vector<std::shared_ptr<Curve>> PolyCurve::curveList() const
{
  return std::vector<std::shared_ptr<Curve>>(curves_.begin(), curves_.end());
}

PointVector PolyCurve::polyline(double smoothness, double precision) const
{
  PointVector polyline;
  for (uint k = 0; k < size(); k++)
  {
    auto new_poly = curvePtr(k)->polyline(smoothness, precision);
    polyline.reserve(polyline.size() + new_poly.size() - (k ? 1 : 0));
    polyline.insert(polyline.end(), new_poly.begin() + (k ? 1 : 0), new_poly.end());
  }
  return polyline;
}

double PolyCurve::length() const { return length(0, size()); }

double PolyCurve::length(Parameter t) const { return length(0, t); }

double PolyCurve::length(Parameter t1, Parameter t2) const
{
  uint idx1 = curveIdx(t1);
  uint idx2 = curveIdx(t2);

  if (idx1 == idx2)
    return curves_[idx1]->length(t1 - idx1, t2 - idx2);
  else if (idx1 + 1 == idx2)
    return curves_[idx1]->length(t1 - idx1, 1.0) + curves_[idx2]->length(0.0, t2 - idx2);
  else
    return std::accumulate(begin(curves_) + idx1 + 1, begin(curves_) + idx2,
                           curves_[idx1]->length(t1 - idx1, 1.0) + curves_[idx2]->length(0.0, t2 - idx2),
                           [](double sum, std::shared_ptr<Curve> curve) { return sum + curve->length(); });
}

double PolyCurve::iterateByLength(Parameter t, double s, double epsilon) const
{
  double s_t = length(t);
  //  if (s_t + s < 0 || s_t + s > length())
  //    throw std::out_of_range{"Resulting parameter t not in [0, n] range."};
  if (s_t + s < 0)
    return 0;
  if (s_t + s > length())
    return size();

  uint idx = curveIdx(t);
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

  return idx + curvePtr(idx)->iterateByLength(t, s, epsilon);
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

void PolyCurve::manipulateControlPoint(uint idx, const Point& point)
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

Point PolyCurve::valueAt(Parameter t) const
{
  uint idx = curveIdx(t);
  return curvePtr(idx)->valueAt(t - idx);
}

PointVector PolyCurve::valueAt(ParameterVector t_vector) const
{
  PointVector points;
  points.reserve(t_vector.size());
  for (auto t : t_vector)
    points.emplace_back(valueAt(t));
  return points;
}

double PolyCurve::curvatureAt(Parameter t) const
{
  uint idx = curveIdx(t);
  return curvePtr(idx)->curvatureAt(t - idx);
}

double PolyCurve::curvatureDerivativeAt(Parameter t) const
{
  uint idx = curveIdx(t);
  return curvePtr(idx)->curvatureDerivativeAt(t - idx);
}

Vector PolyCurve::tangentAt(Parameter t, bool normalize) const
{
  uint idx = curveIdx(t);
  return curvePtr(idx)->tangentAt(t - idx, normalize);
}

Vector PolyCurve::normalAt(Parameter t, bool normalize) const
{
  uint idx = curveIdx(t);
  return curvePtr(idx)->normalAt(t - idx, normalize);
}

Point PolyCurve::derivativeAt(Parameter t) const
{
  uint idx = curveIdx(t);
  return curvePtr(idx)->derivativeAt(t - idx);
}

Point PolyCurve::derivativeAt(uint n, Parameter t) const
{
  uint idx = curveIdx(t);
  return curvePtr(idx)->derivativeAt(n, t - idx);
}

BoundingBox PolyCurve::boundingBox() const
{
  BoundingBox bbox;
  for (auto& curve_ptr : curves_)
    bbox.extend(curve_ptr->boundingBox());
  return bbox;
}

template <> PointVector PolyCurve::intersection<Curve>(const Curve& curve, bool stop_at_first, double epsilon) const
{
  PointVector points;
  for (auto& curve_ptr : curves_)
  {
    auto new_points = curve_ptr->intersection(curve, stop_at_first, epsilon);
    points.reserve(points.size() + new_points.size());
    points.insert(points.end(), new_points.begin(), new_points.end());
    if (!points.empty() && stop_at_first)
      break;
  }
  return points;
}

template <>
PointVector PolyCurve::intersection<PolyCurve>(const PolyCurve& poly_curve, bool stop_at_first, double epsilon) const
{
  PointVector points;
  for (auto& curve_ptr : curves_)
  {
    auto new_points = poly_curve.intersection(*curve_ptr, stop_at_first, epsilon);
    points.reserve(points.size() + new_points.size());
    points.insert(points.end(), new_points.begin(), new_points.end());
    if (!points.empty() && stop_at_first)
      break;
  }
  return points;
}

Parameter PolyCurve::projectPoint(const Point& point) const
{
  double min_t = curves_.front()->projectPoint(point);
  double min_dist = (point - curves_.front()->valueAt(min_t)).norm();

  for (uint k = 1; k < size(); k++)
  {
    Parameter t = curves_[k]->projectPoint(point);
    double dist = (point - curves_[k]->valueAt(t)).norm();
    if (dist < min_dist)
    {
      min_dist = dist;
      min_t = k + t;
    }
  }
  return min_t;
}

ParameterVector PolyCurve::projectPoint(PointVector point_vector) const
{
  ParameterVector t_vector;
  t_vector.reserve(point_vector.size());
  for (auto point : point_vector)
    t_vector.emplace_back(projectPoint(point));
  return t_vector;
}
