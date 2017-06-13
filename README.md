# TrackerMaterialAnalysis
Code for running the analysis for material of the CMS tracker

INSTRUCTIONS TO RUN THE CODE

1) Run the file MinBias.cpp. This code must be run twice.
  It depends on two bools: the first one is called "isFirstRun" and must be set to "true" if it is the first time you run it. It produces a .root file with some 1D histograms. The second bool "isMC" must be set to "true" to run on MC or "false" to run on data.

  To run the first time:
    on data:  root -b -q MinBias.cpp+\('true','false'\)
    on MC:    root -b -q MinBias.cpp+\('true','true'\)

  The second run produces another .root file with some 3D histograms.
  To run the second time:
  on data:  root -b -q MinBias.cpp+\('false','false'\)
  on MC:    root -b -q MinBias.cpp+\('false','true'\)

2) Run the file MatPlotter.py to produce a .root file with 1D histograms of the radiation length as function of local eta in the track.
