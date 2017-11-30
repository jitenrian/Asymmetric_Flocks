#!/bin/bash
rm gnuplot

for((i=0; i<2400;i++))

do

echo "set xrange [-0:40]">>gnuplot
echo "set yrange [-0:40]">>gnuplot
echo "set title 'Flocking'" >>gnuplot
echo "set term png" >>gnuplot
echo "set output '$i.png'" >>gnuplot	
echo "p '$i.dat' u 1:2:3:4 with vectors" >>gnuplot
#echo "set term x11"

done

gnuplot gnuplot
ffmpeg -f image2 -r 25 -i %d.png -vb 20M -vcodec mpeg4 -y movie.mp4

rm *.png
