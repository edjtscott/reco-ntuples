#!/bin/bash

if [ $# -ne 3 ]
then 
  echo "ERROR - needs three arguments"
  echo "1. Sample name, e.g. Photon_Pt35"
  echo "2. Conversion type, e.g. Unconverted"
  echo "3. Extension, e.g. Old"
  exit 1
fi

SAMPLE=$1
echo "Sample name is $SAMPLE"
  
TYPE=$2
echo "Type of photon is $TYPE"

EXT=$3
echo "Ext is $EXT"

DOSEL="false"
#DOSEL="true"

ANALYSER="flexibleAnalyser"
echo "Code to be used is ${ANALYSER}.cxx"

#COMBINER="efficientCombineGraphs"
#echo "Multigraph plotter to be used is ${COMBINER}.cxx"

if [ -e $ANALYSER ]
then 
  rm $ANALYSER
fi
echo "compiling..."
g++ -o $ANALYSER ${ANALYSER}.cxx `root-config --cflags --libs`
echo "done"
echo "./$ANALYSER $SAMPLE $TYPE $EXT $DOSEL"
./$ANALYSER $SAMPLE $TYPE $EXT $DOSEL

OUT1="output_"
OUT2=".root"
OUT3=""
OUTFILE=$OUT1$SAMPLE$OUT2
echo "Outfile is $OUT1$SAMPLE$OUT2"

PASS="9"

# check directory exists then copy plots across
if [ ! -d ~/www/HGCclustering/Pass${PASS}/$SAMPLE ]
then
  mkdir ~/www/HGCclustering/Pass${PASS}/$SAMPLE
fi
if [ ! -d ~/www/HGCclustering/Pass${PASS}/$SAMPLE/$TYPE ]
then
  mkdir ~/www/HGCclustering/Pass${PASS}/$SAMPLE/$TYPE
  cp ~/www/ICHEP16/ValidationRecommendedSingle/index.php ~/www/HGCclustering/Pass${PASS}/$SAMPLE/$TYPE
fi
if [ ! -d ~/www/HGCclustering/Pass${PASS}/$SAMPLE/$TYPE/$EXT ]
then
  mkdir ~/www/HGCclustering/Pass${PASS}/$SAMPLE/$TYPE/$EXT
  cp ~/www/ICHEP16/ValidationRecommendedSingle/index.php ~/www/HGCclustering/Pass${PASS}/$SAMPLE/$TYPE/$EXT
fi
cp HGCPlots/* ~/www/HGCclustering/Pass${PASS}/$SAMPLE/$TYPE/$EXT
rm HGCPlots/*
#rm $OUTFILE
