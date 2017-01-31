#!/bin/bash

#./makePlots.sh Photon_Pt35 Unconverted SensorDependent_MoreConfigurableTest
#./makePlots.sh Photon_Pt35 Converted SensorDependent_MoreConfigurableTest
./makePlots.sh Photon_Pt35 Unconverted SensorDependent_MoreConfigurable_Kappa5
./makePlots.sh Photon_Pt35 Converted SensorDependent_MoreConfigurable_Kappa5

#./makePlots.sh Photon_Pt35 Unconverted SensorDependent
#./makePlots.sh Photon_Pt35 Converted SensorDependent

#./makePlots.sh Photon_Pt35 Unconverted SensorDependent
#./makePlots.sh Photon_Pt10 Unconverted SensorDependent
#./makePlots.sh Photon_Pt5 Unconverted SensorDependent
#./makePlots.sh Photon_Pt35 Converted SensorDependent
#./makePlots.sh Photon_Pt10 Converted SensorDependent
#./makePlots.sh Photon_Pt5 Converted SensorDependent

#./makePlots.sh Pion_Pt35 All SensorDependent
#./makePlots.sh Pion_Pt10 All SensorDependent
#./makePlots.sh Pion_Pt5 All SensorDependent


#EXT="Old"
#EXT="pre15"

#./makePhotonPlots.sh Photon_Pt35 Unconverted $EXT
#./makePhotonPlots.sh Photon_Pt35 Converted $EXT
#./makePhotonPlots.sh Photon_Pt10 Unconverted $EXT
#./makePhotonPlots.sh Photon_Pt10 Converted $EXT
#./makePhotonPlots.sh Photon_Pt5 Unconverted $EXT
#./makePhotonPlots.sh Photon_Pt5 Converted $EXT
#./makePhotonPlots.sh Photon_Pt2 Unconverted $EXT
#./makePhotonPlots.sh Photon_Pt2 Converted $EXT

#PASS="3"
#COMBINER="efficientCombineGraphs"

# check and make combined plots
#if [ ! -d ~/www/HGCclustering/Pass${PASS}/Other ]
#then
#  mkdir ~/www/HGCclustering/Pass${PASS}/Other
#  cp ~/www/ICHEP16/ValidationRecommendedSingle/index.php ~/www/HGCclustering/Pass${PASS}/Other
#fi
#g++ -o $COMBINER ${COMBINER}.cxx `root-config --cflags --libs`
#./$COMBINER
#cp HGCPlots/* ~/www/HGCclustering/Pass${PASS}/Other
#rm HGCPlots/*
