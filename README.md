# Bezier
[![Build Status](https://api.travis-ci.com/romb-technologies/Bezier.svg?branch=master)](https://api.travis-ci.com/romb-technologies/Bezier)
![v0.2](https://img.shields.io/badge/version-0.2-blue.svg)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/aceb46ce7de1407abd56cfc127dba5f1)](https://www.codacy.com/manual/stribor14/Bezier?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=stribor14/Bezier-cpp&amp;utm_campaign=Badge_Grade)

Fast and lightweight class for using the Bezier curves of any order in C++

*Algorithm implementations are based on [A Primer on Bezier Curves](https://pomax.github.io/bezierinfo/) by Pomax*

## Key Features
  - Any number of control points
  - Fast operations on curves
  - Dynamic manipulation
  - Composite Bezier curves (polycurves)

CMake *find_package()* compatible!
```
find_package(Bezier)
target_link_libraries(target bezier)
```

## Implemented methods
  - Get value, derivative, curvature, tangent and normal for parameter *t*
  - Get t from projection any point onto a curve
  - Get precise length for any part of curve
  - Get a derivative curve (hodograph)
  - Split into two subcurves
  - Find curve roots and bounding box
  - Find points of intersection with another curve
  - Elevate/lower order
  - Apply parametric and geometric continuities
  - etc.
  
## Wish list
    - [ ] Polycurve - oversee continuities between consecutive sub-curves
    - [ ] Polycurve - propagation of sub-curve manipulation depending on continutiy
    - [ ] Bezier shapes
    - [ ] More sophisticated example

## Dependencies
  - c++11
  - Eigen3

## Instalation
### System-wide installation
```
git clone https://github.com/romb-technologies/Bezier
mkdir Bezier/build
cd Bezier/build
cmake ..
make
make install
```
### ROS
- for use within a ROS workspace without the system-wide installation, clone the repo to src folder in you catkin workspace 

## Example program __[OUTDATED]__ 
A small Qt5 based program written as a playground for manipulating Bezier curves.
- press *__H__* for a list of possible actions

### Additional dependencies
 - qt5-default 

## Licence
Apache License Version 2.0

## Credit

The development of this software was in part supported by the European Cohesion Fund, through grant number KK.03.2.2.04.314 "Software modules for advanced autonomous motion of load-transportation vehicles (Soft4AGV)" of the "Innovations of newly established SMEs - phase 2" open call.
