#!/bin/bash
g++ `root-config --cflags` -g -O2 -o HDF5Reader HDF5Reader.cpp `root-config --libs` -I /usr/include/hdf5/serial/ -lhdf5_serial -lhdf5_cpp
