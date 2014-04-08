#include<iostream>
using namespace std;


void launchjob(char *rs=NULL,std::string start1="0",std::string end1="0"){


         char puff[50]; // buff is large enough to hold the entire formatted string

        sprintf(puff, "%s", rs);


gROOT->TROOT::ProccessLine(".x setup.C");
gROOT->TROOT::ProccessLine(".L readBacon.C++");
gSystem->Load("readBacon_C.so");
gROOT->TROOT::ProccessLine(readBacon(puff,start1,end1));
