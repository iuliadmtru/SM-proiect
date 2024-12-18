#!/bin/bash

make

PROGNAME=$(basename $0)
if [ $# -lt 3 ]; then
    echo "Usage: ./get_time [start] [incr] [end]"
    exit 1
fi

start=$1
incr=$2
end=$3

for ((i=start;i<=end;i*=incr)); do
    echo -e "**$i threads**\n"

    echo "| Size | Time (s) |"
    echo "|------|----------|"

    for img in ../tests/*; do
        img=$(realpath -s --relative-to=../tests/ $img)
        size=$(basename $img .jpg)

        { time ./gaussian_blur ../tests/$img out/$img 10 $i; } 2> time.log
        echo "|" $size "|" $(cat time.log | tail -n +2 | awk -F "[ \t]*" '{print $2}' | head -n 1) "|"
    done
    echo ""
done
