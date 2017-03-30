## camera-resection

Computes camera position, rotations, and pixel scale parameters, given points' 3d and screen coordinates. Uses a
modified Gauss-Newton non-linear solver.

## Motivation

A useful library for obtaining camera calibration data.

## Github page
http://joel1fx.github.io/camera-resection/

## Installation

Make sure Eigen is installed. If it is installed in a non-standard path, you will need to edit the Makefile.

For more information about Eigen and its license see:

http://eigen.tuxfamily.org/index.php?title=Main_Page

Download the repo and
```
make
```

The code is std C++, so it should compile readily.

## Command Line Reference

Usage:

```
camera_resection
```

## License

MIT license

