# Tutte Parameterization

## Dependency

* [Eigen3.4.0](http://eigen.tuxfamily.org/)

## Usage(Linux)

```SHELL
git clone https://github.com/FlynnHHH/TutteParameterization.git
```

Edit line **9** of CmakeLists.txt to set the values of **EIGEN_PATH**

```SHELL
mkdir build && cd build
cmake ../
make
./TUTTE_PARAMETERIZATION path/to/input_mesh.obj path/to/output_mesh.obj
```
