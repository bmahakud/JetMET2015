SHELL = /bin/bash
FASTJETLOCATION=/afs/cern.ch/cms/slc5_amd64_gcc462/external/fastjet/3.0.1-cms3
ROOTINC := `root-config --glibs --cflags`
CC=g++
FJCONTRIBLOCATION=/afs/cern.ch/user/b/bmahakud/work/NewBaconFramework/CMSSW_5_3_13/src/fjcontrib-1.009
BINDIR=bin
INCDIR=/afs/cern.ch/user/b/bmahakud/work/NewBaconFramework/CMSSW_5_3_13/src/fjcontrib-1.009/include/fastjet/contrib
CMSSWLIB=/afs/cern.ch/user/b/bmahakud/work/NewBaconFramework/CMSSW_5_3_13
PATHCOND=/afs/cern.ch/user/b/bmahakud/work/NewBaconFramework/CMSSW_5_3_13/src/CondFormats/JetMETObjects/interface

readingBacon: $(FJCONTRIBLOCATION)/lib/libSoftDrop.a $(CMSSWLIB)/lib/slc5_amd64_gcc462/*.so
        @mkdir -p $(BINDIR)
   # Note: $(CXXFLAGS) is after Fastjet flags as Fastjet includes
        $(CXX) -I$(FJCONTRIBLOCATION)/$(INCDIR) $@.cc \
        `$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
        $(ROOTINC)\
        $(CXXFLAGS) -Wno-shadow \
        -o $(BINDIR)/$@.exe \
        -L$(FJCONTRIBLOCATION)/lib -lSoftDrop \
        -L$(CMSSWLIB)/lib/slc5_amd64_gcc462 \
        -L$(FASTJETLOCATION)/lib \
        `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins`
        @ln -fs $(BINDIR)/$@.exe $@.exe
        @rm -f $@.o
