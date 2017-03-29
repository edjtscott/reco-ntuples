#!/bin/bash

#./makePlots.sh Pion_Pt35 All HigherStats_4cm

#./makePlots.sh Pion_Pt25 All RelValNoPU
#./makePlots.sh Pion_Pt25 All RelVal140PU

./makePlots.sh Pion_Pt35 All UpdatedEScale_40mm
./makePlots.sh Pion_Pt35 All UpdatedEScale_50mm
./makePlots.sh Pion_Pt35 All UpdatedEScale_EE2FH4BH4
./makePlots.sh Pion_Pt35 All UpdatedEScale_EE2FH5BH5

#./makePlots.sh Photon_Pt35 Unconverted UpdatedEScale_20mm
#./makePlots.sh Photon_Pt35 Converted UpdatedEScale_20mm
#./makePlots.sh Photon_Pt35 Unconverted KDmultis_20mm
#./makePlots.sh Photon_Pt35 Converted KDmultis_20mm

#./makePlots.sh Photon_Pt35 Unconverted CartesianCorrected_20mm
#./makePlots.sh Photon_Pt35 Converted CartesianCorrected_20mm

#./makePlots.sh Pion_Pt35 All CartesianCorrected_20mm
#./makePlots.sh Pion_Pt35 All CartesianCorrected_30mm
#./makePlots.sh Pion_Pt35 All CartesianCorrected_40mm
#./makePlots.sh Pion_Pt35 All CartesianCorrected_50mm
#./makePlots.sh Pion_Pt35 All CartesianCorrected_60mm

#./makePlots.sh Photon_Pt35 Unconverted SuperClusteringOneEight
#./makePlots.sh Photon_Pt35 Converted SuperClusteringOneEight
#./makePlots.sh Photon_Pt35 Unconverted SuperClusteringTwoTwo
#./makePlots.sh Photon_Pt35 Converted SuperClusteringTwoTwo
#./makePlots.sh Photon_Pt35 Unconverted SuperClusteringTwoSix
#./makePlots.sh Photon_Pt35 Converted SuperClusteringTwoSix
#
#./makePlots.sh Pion_Pt35 All SensorDependent
#./makePlots.sh Pion_Pt35 All SensorDependent_FH2BH2R03
#./makePlots.sh Pion_Pt35 All SensorDependent_FH3BH5R03
#./makePlots.sh Pion_Pt35 All SensorDependent_FH5BH7R03
#./makePlots.sh Pion_Pt35 All SensorDependent_FH3BH5
#./makePlots.sh Pion_Pt35 All SensorDependent_FH5BH7


#./makePlots.sh Photon_Pt35 Unconverted SensorDependent_MoreConfigurableTest
#./makePlots.sh Photon_Pt35 Converted SensorDependent_MoreConfigurableTest
#./makePlots.sh Photon_Pt35 Unconverted SensorDependent_MoreConfigurable_Kappa5
#./makePlots.sh Photon_Pt35 Converted SensorDependent_MoreConfigurable_Kappa5

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
