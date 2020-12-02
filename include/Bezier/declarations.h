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

#ifndef DECLARATIONS_H
#define DECLARATIONS_H

#include <Eigen/Dense>
#include <memory>
#include <vector>

/*!
 * Nominal namespace containing class pre-definitions and typedefs
 */
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
class Curve;

/*!
 * \brief A polyline class
 *
 * A class for operations on polyline, mainly used for
 * polyline representation of curves. Underlying structure
 * is a vector of points.
 */
class PolyLine;

/*!
 * \brief A Bezier polycurve class
 *
 * A class for linking multiple Bezier curves with at least
 * C0 continuity. It allows subcurve and continuity manipulation.
 * Both parametric and geometric continuity are supported.
 */
class PolyCurve;

/*!
 * \brief Parameter of the curve function
 */
using Parameter = double;

/*!
 * \brief A vector of Parameters
 */
using ParameterVector = std::vector<Parameter>;

/*!
 * \brief Point in xy plane
 */
using Point = Eigen::Vector2d;

/*!
 * \brief A vector of Points
 */
using PointVector = std::vector<Point>;

/*!
 * \brief A Vector in xy plane
 */
using Vector = Eigen::Vector2d;

/*!
 * \brief Bounding box class
 */
using BoundingBox = Eigen::AlignedBox2d;
} // namespace Bezier
#endif // DECLARATIONS_H
