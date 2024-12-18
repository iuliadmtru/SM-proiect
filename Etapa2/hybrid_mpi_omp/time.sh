#!/bin/bash

module load mpi/openmpi-x86_64

make -f Makefile.fep

ldd gaussian_blur

threads=$1

echo -e "**$threads threads**\n"

echo "| Size | Time (s) |"
echo "|------|----------|"

for img in ../tests/*; do
    img=$(realpath -s --relative-to=../tests/ $img)
    size=$(basename $img .jpg)

    { time mpirun ./gaussian_blur ../tests/$img out/$img 10 $threads; } 2> time.log
    echo "|" $size "|" $(cat time.log | tail -n +2 | awk -F "[ \t]*" '{print $2}' | head -n 1) "|"
done
echo ""