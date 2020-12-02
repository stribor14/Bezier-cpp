
# Bezier
[![Build Status](https://api.travis-ci.com/romb-technologies/Bezier.svg?branch=master)](https://api.travis-ci.com/romb-technologies/Bezier)
![v0.3](https://img.shields.io/badge/version-0.3-blue.svg)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/aceb46ce7de1407abd56cfc127dba5f1)](https://www.codacy.com/manual/stribor14/Bezier?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=stribor14/Bezier-cpp&amp;utm_campaign=Badge_Grade)

Fast and lightweight class for using the Bezier curves of any order in C++

*Some algorithm implementations are based on [A Primer on Bezier Curves](https://pomax.github.io/bezierinfo/) by Pomax*

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
  - Get value, derivative, curvature, tangent and normal for a parameter *t*
  - Point projection onto curve
  - Get curve length, parameter iteration by length
  - Get a derivative curve (hodograph)
  - Split into subcurves
  - Find curve roots, extrema and bounding box
  - Find points of intersection
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

## Example program
A small Qt5 based program written as a playground for manipulating Bezier curves.
- press *__H__* for a list of possible actions
- *outdated* - newer features are not used in the example application

### Additional dependencies
 - qt5-default 

## Licence
Apache License Version 2.0
