cd /afs/cern.ch/user/b/bmahakud/work/NewBaconFramework/CMSSW_5_3_13/src/BaconAna/DataFormats/tst
eval `scramv1 runtime -sh`
cp setup.C ~/
cp readBacon.C ~/
cp /afs/cern.ch/user/b/bmahakud/work/NewBaconFramework/CMSSW_5_3_13/src/BaconAna/DataFormats/tst/Input.root ~/
cd /afs/cern.ch/user/b/bmahakud/work/NewBaconFramework/CMSSW_5_3_13/src/BaconAna/DataFormats/tst 
root -l -b launchjob.C\(\"Input.root\",\"3\",\"30\"\)
