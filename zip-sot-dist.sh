#!/bin/bash

# version of the LPFG-Self-Organizing-Trees code
NAME=sot
VER=1.2.8

# Source files to copy
CPPS="methods.cpp mtg.cpp randist.cpp scatter.cpp parameters.inc"
# Header files to copy
HPPS="user.h quaternion.hpp hlu.hpp env.hpp cyl.hpp default_parameters.inc"
# LPFG files
LPFGS="lsystem.l view.v material.mat bark3-64.rgb leaf.s leaft.rgb" 
# Helping files
HELPS="README.md NEWS"

# Create temporary directory
DISTDIR=$NAME-dist-$VER
mkdir $DISTDIR

# Copy the file into the directory
cp $CPPS $HPPS $LPFGS $HELPS $DISTDIR

# Create an archive
ZIPARCH=$NAME-dist-$VER.zip
zip -r $ZIPARCH $DISTDIR

# Remove the directory
rm -rf $DISTDIR
