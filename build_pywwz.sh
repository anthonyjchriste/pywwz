#!/bin/bash


python3 -m numpy.f2py -c fortran/wwz11_f2py.f90 -m pywwz
mkdir -p build
mv pywwz.cpython* build/.