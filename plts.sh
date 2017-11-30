#!/bin/bash
rm plot

for((i=2; i<10;i++))
do

echo "set title '$i'" >>plot
echo "set term png" >>plot
echo "set output '$i.png'" >>plot	
echo "p 'order_parameter.dat' u 1:$i" >>plot


done

gnuplot plot
