#!/bin/bash

#PART="Photon"
PART="Pion"

PT="35"

#NAME=""
#NAME="AcceptedPR"
#NAME="AcceptedPR_EE2FH5BH5"
#NAME="PU200_AcceptedPR"
#NAME="PU200_AcceptedPR_EE2FH5BH5"
#NAME="D17"
#NAME="PU200_D17"
NAME="D17_255_225"
#NAME="PU200_D17_255_225"

#INFILE="/afs/cern.ch/work/e/escott/public/HGCStudies/Ntuples/FromClemens_${PART}_Pt${PT}.root"
INFILE="/afs/cern.ch/work/e/escott/public/HGCStudies/Ntuples/partGun_${PART}_Pt${PT}_${NAME}.root"

#WEB="/afs/cern.ch/user/e/escott/www/HGCclustering/FixClusteringTest/${PART}_Pt${PT}/"
WEB="/afs/cern.ch/user/e/escott/www/HGCclustering/${NAME}/${PART}_Pt${PT}/"

#CONV=1

#echo "python plotHGCal.py -f $INFILE -w $WEB -p $PART -m $PT -c $CONV"
#echo ""
#python plotHGCal.py -f $INFILE -w $WEB -p $PART -m $PT -c $CONV
#python plotHGCal.py -f $INFILE -w $WEB -p $PART -m $PT -c $CONV -b

python plotHGCal.py -f $INFILE -w $WEB -p $PART -m $PT -c 0
#python plotHGCal.py -f $INFILE -w $WEB -p $PART -m $PT -c 1
#python plotHGCal.py -f $INFILE -w $WEB -p $PART -m $PT -c 2
