/*
 * Copyright 2019 Mirko Kokot
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef BEZIER_H
#define BEZIER_H

#include <map>

#include "declarations.h"

namespace Bezier
{
/*!
 * \brief A Bezier curve class
 *
 * A class for storing and using any-order Bezier curve.
 * It uses private and static caching for storing often accessed data.
 * Private caching is used for data concerning individual curve, while
 * static caching is used for common data (coefficient matrices)
 */
class Curve
{
public:
  /*!
   * \brief Create the Bezier curve
   * \param points Nx2 matrix where each row is one of N control points that define the curve
   */
  Curve(const Eigen::MatrixX2d& points);

  /*!
   * \brief Create the Bezier curve
   * \param points A vector of control points that define the curve
   */
  Curve(const PointVector& points);

  /*!
   * \brief Create the Bezier curve copy
   * \param curve A Bezier curve to copy
   */
  Curve(const Curve& curve);

  /*!
   * \brief Get order of the curve (Nth order curve is described with N+1 points);
   * \return Order of curve
   */
  uint order();

  /*!
   * \brief Get the control points
   * \return A vector of control points
   */
  PointVector controlPoints() const;

  /*!
   * \brief Get first and last control points
   * \return A pair of end points
   */
  std::pair<Point, Point> endPoints() const;

  /*!
   * \brief Get a polyline representation of the curve as a vector of points on curve
   * \param smoothness Smoothness factor > 1 (more resulting points when closer to 1)
   * \param precision Minimal distance between two subsequent points
   * \return A vector of polyline vertices
   */
  PointVector polyline(double smoothness = 1.0001, double precision = 1.0) const;

  /*!
   * \brief Compute exaxt arc length with Legendre-Gauss quadrature
   * \return Arc length
   * \warning Precision depends on value of LEGENDRE_GAUSS_N at compile time
   */
  double length() const;

  /*!
   * \brief Compute exact arc length with Legendre-Gauss quadrature
   * \param t Curve parameter to which length is computed
   * \return Arc length from start to parameter t
   * \warning Precision depends on value of LEGENDRE_GAUSS_N at compile time
   */
  double length(Parameter t) const;

  /*!
   * \brief Compute exact arc length with Legendre-Gauss quadrature
   * \param t1 Curve parameter from which length is computed
   * \param t2 Curve parameter to which length is computed
   * \return Arc length between paramaters t1 and t2
   * \warning Precision depends on value of LEGENDRE_GAUSS_N at compile time
   */
  double length(Parameter t1, Parameter t2) const;

  /*!
   * \brief Compute parameter t which is S distance from given t
   * \param t Curve parameter
   * \param s Distance to iterate
   * \param epsilon Precision of resulting t
   * \return New parameter t
   */
  Parameter iterateByLength(Parameter t, double s, double epsilon = 0.001) const;

  /*!
   * \brief Reverse order of control points
   */
  void reverse();

  /*!
   * \brief Set the new coordinates to a control point
   * \param index Index of chosen control point
   * \param point New control point
   */
  void manipulateControlPoint(uint idx, const Point& point);

  /*!
   * \brief Manipulate the curve so that it passes through wanted point for given 't'
   * \param t Curve parameter
   * \param point Point where curve should pass for a given t
   *
   * \warning Only works for quadratic and cubic curves
   * \warning Resets cached data
   */
  void manipulateCurvature(Parameter t, const Point& point);

  /*!
   * \brief Raise the curve order by 1
   *
   * Curve will always retain its shape
   * \warning Resets cached data
   */
  void elevateOrder();

  /*!
   * \brief Lower the curve order by 1
   *
   * If current shape cannot be described by lower order, it will be best aproximation
   * \warning Resets cached data
   */
  void lowerOrder();

  /*!
   * \brief Get the point on curve for a given t
   * \param t Curve parameter
   * \return Point on a curve for a given t
   */
  Point valueAt(Parameter t) const;

  /*!
   * \brief Get the point vector on curve for given parameters
   * \param t_vector Curve parameters
   * \return Vector of points on a curve for given parameters
   */
  PointVector valueAt(ParameterVector t_vector) const;

  /*!
   * \brief Get curvature of the curve for a given t
   * \param t Curve parameter
   * \return Curvature of a curve for a given t
   */
  double curvatureAt(Parameter t) const;

  /*!
   * \brief Get curvature derivative of the curve for a given t
   * \param t Curve parameter
   * \return Curvature derivative of a curve for a given t
   */
  double curvatureDerivativeAt(Parameter t) const;

  /*!
   * \brief Get the tangent of the curve for a given t
   * \param t Curve parameter
   * \param normalize If the resulting tangent should be normalized
   * \return Tangent of a curve for a given t
   */
  Vector tangentAt(Parameter t, bool normalize = true) const;

  /*!
   * \brief Get the normal of the curve for a given t
   * \param t Curve parameter
   * \param normalize If the resulting normal should be normalized
   * \return Normal of a curve for given t
   */
  Vector normalAt(Parameter t, bool normalize = true) const;

  /*!
   * \brief Get the derivative of a curve
   * \return Derivative curve
   */
  std::shared_ptr<const Curve> derivative() const;

  /*!
   * \brief Get the nth derivative of a curve
   * \param n Desired number of derivative
   * \return Derivative curve
   * \warning Parameter n cannot be zero
   */
  std::shared_ptr<const Curve> derivative(uint n) const;

  /*!
   * \brief Get value of a derivative for a given t
   * \param t Curve parameter
   * \return Curve derivative at t
   */
  Vector derivativeAt(Parameter t) const;

  /*!
   * \brief Get value of an nth derivative for a given t
   * \param n Desired number of derivative
   * \param t Curve parameter
   * \return nth curve derivative at t
   */
  Vector derivativeAt(uint n, Parameter t) const;

  /*!
   * \brief Get roots of the curve on both axes
   * \return A vector of parameters where curve passes through axes
   */
  ParameterVector roots() const;

  /*!
   * \brief Get all extrema of the curve
   * \return A vector of parameters where maxima are
   */
  ParameterVector extrema() const;

  /*!
   * \brief Get the bounding box of curve
   * \return Bounding box (if use_roots is false, returns the bounding box of control points)
   */
  BoundingBox boundingBox() const;

  /*!
   * \brief Split the curve into two subcurves
   * \param z Parameter t at which to split the curve
   * \return Pair of two subcurves
   */
  std::pair<Curve, Curve> splitCurve(double z = 0.5) const;

  /*!
   * \brief Get the points of intersection with another curve
   * \param curve Curve to intersect with
   * \param stop_at_first If first point of intersection is enough
   * \param epsilon Precision of resulting intersection
   * \return A vector af points of intersection between curves
   */
  PointVector intersection(const Curve& curve, bool stop_at_first = false, double epsilon = 0.001) const;

  /*!
   * \brief Get the parameter t where curve is closest to given point
   * \param point Point to project on curve
   * \return Parameter t
   */
  Parameter projectPoint(const Point& point) const;

  /*!
   * \brief Get the parameter t vector where curve is closest to given points
   * \param point_vector Points to project on curve
   * \return Vector of parameters t
   */
  ParameterVector projectPoint(PointVector point_vector) const;

  /*!
   * \brief Get distance of the point to the curve
   * \param point Point to project on curve
   * \return Distance to the curve
   */
  double distance(const Point& point) const;

  /*!
   * \brief Get the distance vector of points to the curve
   * \param point_vector Points to project on curve
   * \return Vector of distances
   */
  std::vector<double> distance(PointVector point_vector) const;

  /*!
   * \brief applyContinuity Apply geometric continuity based on the another curve.
   * \param locked_curve Curve on which calculation are based.
   * \param beta_coeffs Beta-constraints used to calculate continuity. Size defines continuity order.
   */
  void applyContinuity(const Curve& source_curve, const std::vector<double> &beta_coeffs);

private:
  /*!
   * \brief Coefficients for matrix operations
   */
  using Coeffs = Eigen::MatrixXd;
  /*!
   * \brief Map of different coefficient matrices, depending on the order of the curve
   */
  using CoeffsMap = std::map<uint, Coeffs>;

  /// Number of control points (order + 1)
  uint N_;
  /// N x 2 matrix where each row corresponds to control Point
  Eigen::MatrixX2d control_points_;

  // private caching
  std::shared_ptr<const Curve> cached_derivative_;         /*! If generated, stores derivative for later use */
  std::unique_ptr<ParameterVector> cached_roots_;          /*! If generated, stores roots for later use */
  std::unique_ptr<BoundingBox> cached_bounding_box_;       /*! If generated, stores bounding box for later use */
  std::unique_ptr<PointVector> cached_polyline_;           /*! If generated, stores polyline for later use */
  std::pair<double, double> cached_polyline_params_{0, 0}; /*! Smootheness and precision of cached polyline */
  std::unique_ptr<Eigen::VectorXd>
      cached_projection_polynomial_part_;                   /*! Constant part of point projection polynomial */
  Eigen::MatrixXd cached_projection_polynomial_derivative_; /*! Polynomial representation of the curve derivative */

  /// Reset all privately cached data
  inline void resetCache();

  // static caching
  static CoeffsMap bernstein_coeffs_;       /*! Map of Bernstein coefficients */
  static CoeffsMap splitting_coeffs_left_;  /*! Map of coefficients to get subcurve for t = [0, 0.5] */
  static CoeffsMap splitting_coeffs_right_; /*! Map of coefficients to get subcurve for t = [0.5, 1] */
  static CoeffsMap elevate_order_coeffs_;   /*! Map of coefficients for elevating the order of the curve */
  static CoeffsMap lower_order_coeffs_;     /*! Map of coefficients for lowering the order of the curve */

  /// Static getter function for Bernstein coefficients
  static Coeffs bernsteinCoeffs(uint n);
  /// Static getter function for coefficients to get a subcurve t = [0, z];
  static Coeffs splittingCoeffsLeft(uint n, Parameter z = 0.5);
  /// Static getter function for coefficients to get a subcurve t = [z, 1];
  static Coeffs splittingCoeffsRight(uint n, Parameter z = 0.5);
  /// Static getter function for coefficients to elevate order of curve
  static Coeffs elevateOrderCoeffs(uint n);
  /// Static getter function for coefficients to lower order of curve
  static Coeffs lowerOrderCoeffs(uint n);
};

} // namespace Bezier

#endif // BEZIER_H
