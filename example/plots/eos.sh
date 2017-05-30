#!/bin/bash
gnuplot eos.plot
convert -density 300 eos.eps eos.png 
