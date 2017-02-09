#!/bin/bash
input=$1
natoms=$(head -n1 $input)
nlines=$(wc -l $input)
echo $natoms $nlines
./com.f $input $natoms $nlines
