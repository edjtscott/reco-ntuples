#!/bin/bash

#PART="Photon"
PART="Pion"
PT="35"
INFILE="/afs/cern.ch/work/e/escott/public/HGCStudies/Ntuples/FromClemens_${PART}_Pt${PT}.root"
WEB="/afs/cern.ch/user/e/escott/www/HGCclustering/FixClusteringTest/${PART}_Pt${PT}/"
#CONV=1

#echo "python plotHGCal.py -f $INFILE -w $WEB -p $PART -m $PT -c $CONV"
#echo ""
#python plotHGCal.py -f $INFILE -w $WEB -p $PART -m $PT -c $CONV
#python plotHGCal.py -f $INFILE -w $WEB -p $PART -m $PT -c $CONV -b

python plotHGCal.py -f $INFILE -w $WEB -p $PART -m $PT -c 0
python plotHGCal.py -f $INFILE -w $WEB -p $PART -m $PT -c 1
python plotHGCal.py -f $INFILE -w $WEB -p $PART -m $PT -c 2
