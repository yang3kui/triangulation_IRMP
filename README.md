# triangulation_IRMP

This repo contains an implementation of the Iteratively Reweighted MidPoint method (IRMP) for fast multiple view triangulation.  It also contains a comparsion between the IRMP method and doing multiple view triangulation using Ceres. The preprint version of the corresponding paper is now available at [IEEE Xplore](https://ieeexplore.ieee.org/document/8611369).

If you use this code in an academic context, please cite the following paper:

```
@ARTICLE{yangkui_irmp,
author={K. Yang and W. Fang and Y. Zhao and N. Deng},
journal={IEEE Robotics and Automation Letters},
title={Iteratively Reweighted Midpoint Method for Fast Multiple View Triangulation},
year={2019},
volume={4},
number={2},
pages={708-715},
doi={10.1109/LRA.2019.2893022},
ISSN={2377-3766},
month={April},}
```

## Dependencies 
* Ceres 
* Eigen

## Build
```bash
mkdir build
cd build 
cmake ..
make -j4
```

## Run:
```
./triangulation
```
