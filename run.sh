#!/bin/bash

L=20
J=0.50

mkdir "data"
cd data

#for L in 20 30 40 50
#do
#  for J in 0.40 0.42 0.43 0.44 0.45 0.46 0.50
#  do
    dir="L"$L"-J"$J
    mkdir $dir
    cd $dir
    ../../susIsing -L $L -J $J
    cd ..
#  done
#done

cd ..
