#!/usr/bin/env python3                          
  
import sys
import os

def str_split(lines):
list2 = lines.split()
return list2

pathname = os.path.dirname(sys.argv[0])
WORKING_DIR = os.path.abspath(pathname)
FRiP = float(sys.argv[1])
input_bedGraph = sys.argv[2]
output_bedGraph = sys.argv[3]

f_bedGraph = open(WORKING_DIR + '/' + input_bedGraph, 'r')
lines_bedGraph = f_bedGraph.readlines()
f_bedGraph.close()

fo = open(WORKING_DIR + '/' + output_bedGraph, 'w')

for i in range(len(lines_bedGraph)):
if '#' in lines_bedGraph[i]:
fo.write(lines_bedGraph[i])
else:
fo.write(str_split(lines_bedGraph[i])[0] + '\t' +
str_split(lines_bedGraph[i])[1] + '\t' +
str_split(lines_bedGraph[i])[2] + '\t' +
str(round(float(str_split(lines_bedGraph[i])[3])/FRiP,2)) + '\n'
)
fo.close()
