#define AnalysisAvJetMass2_cxx
#include "AnalysisAvJetMass2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnalysisAvJetMass2::Loop()
{
int noEvents=20000;
int JetpTMin=200;
int JetpTMax=350;
char pTRangeAK8[200];
char pTRangeAK10[200];
char pTRangeAK12[200];

sprintf(pTRangeAK8,"Pile up dependence of AK8 Jet mass (%i GeV < Jet pT < %i GeV);Number of primary vertices ;Average Jet mass in GeV",JetpTMin,JetpTMax);

sprintf(pTRangeAK10,"Pile up dependence of AK10 Jet mass (%i GeV < Jet pT < %i GeV);Number of primary vertices ;Average Jet mass in GeV",JetpTMin,JetpTMax);

sprintf(pTRangeAK12,"Pile up dependence of AK12 Jet mass (%i GeV < Jet pT < %i GeV);Number of primary vertices ;Average Jet mass in GeV",JetpTMin,JetpTMax);



TH1F *JetmassAK8_h[41][8];//firstone is no of diff PV,second one is no of grooming 
TH1F *JetmassAK10_h[41][8];
TH1F *JetmassAK12_h[41][8];

   int nPVi=20;
   int nPVf=40;
  int NoGrooming=8;

  int NoOfhistos=nPVf-nPVi+1;
  static const int anPVi=20;
  static const int bnPVf=40; 
  static const int siz=bnPVf-anPVi+1;
  //cout<<"siz = "<<siz<<endl;
  char histname[100];
  double rangeMin=0;
  double rangeMax=1000;
  int NoOfbins=100;


  for(int g0=0;g0<NoGrooming;g0++){
for (int jj=nPVi;jj<nPVf+1;jj++) {
        sprintf(histname,"AveargejetmassAK8%i_%i",jj,g0);
        JetmassAK8_h[jj][g0] = new TH1F(histname,histname,NoOfbins,rangeMin,rangeMax);
        JetmassAK8_h[jj][g0]->SetLineColor(3);
        JetmassAK8_h[jj][g0]->SetLineWidth(3);
        JetmassAK8_h[jj][g0]->GetXaxis()->SetTitle("Number of primary Vertices");
        JetmassAK8_h[jj][g0]->GetXaxis()->SetTitleSize(0.055);
        JetmassAK8_h[jj][g0]->GetXaxis()->SetLabelSize(0.055);
        JetmassAK8_h[jj][g0]->GetYaxis()->SetTitle("Events");
        JetmassAK8_h[jj][g0]->GetYaxis()->SetTitleSize(0.055);
        JetmassAK8_h[jj][g0]->GetYaxis()->SetLabelSize(0.055);
      }
      }

  for(int g0=0;g0<NoGrooming;g0++){
for (int jj=nPVi;jj<nPVf+1;jj++) {
        sprintf(histname,"AveargejetmassAK10%i_%i",jj,g0);
        JetmassAK10_h[jj][g0] = new TH1F(histname,histname,NoOfbins,rangeMin,rangeMax);
        JetmassAK10_h[jj][g0]->SetLineColor(3);
        JetmassAK10_h[jj][g0]->SetLineWidth(3);
        JetmassAK10_h[jj][g0]->GetXaxis()->SetTitle("Number of primary Vertices");
        JetmassAK10_h[jj][g0]->GetXaxis()->SetTitleSize(0.055);
        JetmassAK10_h[jj][g0]->GetXaxis()->SetLabelSize(0.055);
        JetmassAK10_h[jj][g0]->GetYaxis()->SetTitle("Events");
        JetmassAK10_h[jj][g0]->GetYaxis()->SetTitleSize(0.055);
        JetmassAK10_h[jj][g0]->GetYaxis()->SetLabelSize(0.055);
      }
      }

  for(int g0=0;g0<NoGrooming;g0++){
for (int jj=nPVi;jj<nPVf+1;jj++) {
        sprintf(histname,"AveargejetmassAK12%i_%i",jj,g0);
        JetmassAK12_h[jj][g0] = new TH1F(histname,histname,NoOfbins,rangeMin,rangeMax);
        JetmassAK12_h[jj][g0]->SetLineColor(3);
        JetmassAK12_h[jj][g0]->SetLineWidth(3);
        JetmassAK12_h[jj][g0]->GetXaxis()->SetTitle("Number of primary Vertices");
        JetmassAK12_h[jj][g0]->GetXaxis()->SetTitleSize(0.055);
        JetmassAK12_h[jj][g0]->GetXaxis()->SetLabelSize(0.055);
        JetmassAK12_h[jj][g0]->GetYaxis()->SetTitle("Events");
        JetmassAK12_h[jj][g0]->GetYaxis()->SetTitleSize(0.055);
        JetmassAK12_h[jj][g0]->GetYaxis()->SetLabelSize(0.055);
      }
      }











   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
     ///////////////////////////////////////noOfevents to run
    if(jentry==noEvents){
       break;
     } 
     //////////////////////////////////////noOfevents to run
  int nPV=n_PVertex;
/////////////////////////////////////////////////////////////////////////////AK8
  int sizak8=GroomedJet_AK0_8_pt_sd1_corr->size();
      for(int g0=0;g0<NoGrooming;g0++){//no of grooming loop
      for (int jj=0;jj<sizak8;jj++) {//ak8loop
      
      if(nPV >= nPVi && nPV <= nPVf ){//nPV criteria
      if(jj==0 && g0==0){//sd1
            if(GroomedJet_AK0_8_pt_sd1_corr->at(jj) > JetpTMin && GroomedJet_AK0_8_pt_sd1_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_sd1_corr->at(jj));
              }//pt cut              
           }//sd1
         

       if(jj==0 && g0==1){//sd2
       if(GroomedJet_AK0_8_pt_sd2_corr->at(jj) > JetpTMin && GroomedJet_AK0_8_pt_sd2_corr->at(jj) < JetpTMax){//pt cut
       JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_sd2_corr->at(jj));
         }//pt cut                       
         }//sd2
  
     if(jj==0 && g0==2){//sd3
       if(GroomedJet_AK0_8_pt_sd3_corr->at(jj) > JetpTMin && GroomedJet_AK0_8_pt_sd3_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_sd3_corr->at(jj));
                           }//pt cut
                         }//sd3

     if(jj==0 && g0==3){//sd4
        if(GroomedJet_AK0_8_pt_sd4_corr->at(jj) > JetpTMin && GroomedJet_AK0_8_pt_sd4_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_sd4_corr->at(jj));
           }//pt cut             
            }//sd4

     if(jj==0 && g0==4){//tr
             if(GroomedJet_AK0_8_pt_ft_corr->at(jj) > JetpTMin && GroomedJet_AK0_8_pt_ft_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_ft_corr->at(jj));
                    }//pt cut    
                     }//tr

      if(jj==0 && g0==5){//tr
       if(GroomedJet_AK0_8_pt_tr_corr->at(jj) > JetpTMin && GroomedJet_AK0_8_pt_tr_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_tr_corr->at(jj));
                     }
                         }//tr

      if(jj==0 && g0==6){//pr
    if(GroomedJet_AK0_8_pt_pr_corr->at(jj) > JetpTMin && GroomedJet_AK0_8_pt_pr_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_pr_corr->at(jj));
                         }
                         }//pr

     if(jj==0 && g0==7){//ungroomed
          if(GroomedJet_AK0_8_pt_corr->at(jj) > JetpTMin && GroomedJet_AK0_8_pt_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_corr->at(jj));
                         }
                         }//ungroomed


          
       }//nPV criteria
}//ak8 loop
}//no of grooming loop
///////////////////////////////////////////////////////////////////////////AK8


///////////////////////////////////////////////////////////////////////////AK10

int sizak10=GroomedJet_AK1_pt_sd1_corr->size();
      for(int g0=0;g0<NoGrooming;g0++){//no of grooming loop
      for (int jj=0;jj<sizak10;jj++) {//ak10loop
      
      if(nPV >= nPVi && nPV <= nPVf ){//nPV criteria
      if(jj==0 && g0==0){//sd1
       if(GroomedJet_AK1_pt_sd1_corr->at(jj) > JetpTMin && GroomedJet_AK1_pt_sd1_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_sd1_corr->at(jj));
                        } 
                        }//sd1
         

       if(jj==0 && g0==1){//sd2
       if(GroomedJet_AK1_pt_sd2_corr->at(jj) > JetpTMin && GroomedJet_AK1_pt_sd2_corr->at(jj) < JetpTMax){//pt cut
       JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_sd2_corr->at(jj));
               }                         
             }//sd2
  
     if(jj==0 && g0==2){//sd3
     if(GroomedJet_AK1_pt_sd2_corr->at(jj) > JetpTMin && GroomedJet_AK1_pt_sd2_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_sd3_corr->at(jj));
                             }
                         }//sd3

     if(jj==0 && g0==3){//sd4
       if(GroomedJet_AK1_pt_sd4_corr->at(jj) > JetpTMin && GroomedJet_AK1_pt_sd4_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_sd4_corr->at(jj));
                 }//pt cut        
                 }//sd4

     if(jj==0 && g0==4){//tr
       if(GroomedJet_AK1_pt_ft_corr->at(jj) > JetpTMin && GroomedJet_AK1_pt_ft_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_ft_corr->at(jj));
                         }
                         }//tr

      if(jj==0 && g0==5){//tr
      if(GroomedJet_AK1_pt_tr_corr->at(jj) > JetpTMin && GroomedJet_AK1_pt_tr_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_tr_corr->at(jj));
                           }
                         }//tr

      if(jj==0 && g0==6){//pr
      if(GroomedJet_AK1_pt_pr_corr->at(jj) > JetpTMin && GroomedJet_AK1_pt_pr_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_pr_corr->at(jj));
                         }//pt cut
                         }//pr

     if(jj==0 && g0==7){//ungroomed
      if(GroomedJet_AK1_pt_corr->at(jj) > JetpTMin && GroomedJet_AK1_pt_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_corr->at(jj));
            }                         
         }//ungroomed



       
       }//nPV criteria
}//ak10 loop
}//no of grooming loop

///////////////////////////////////////////////////////////////////////////AK10


///////////////////////////////////////////////////////////////////////////////AK12

int sizak12=GroomedJet_AK1_2_pt_sd1_corr->size();
      for(int g0=0;g0<NoGrooming;g0++){//no of grooming loop
      for (int jj=0;jj<sizak12;jj++) {//ak12loop
      
      if(nPV >= nPVi && nPV <= nPVf ){//nPV criteria
      if(jj==0 && g0==0){//sd1
      if(GroomedJet_AK1_2_pt_sd1_corr->at(jj) > JetpTMin && GroomedJet_AK1_2_pt_sd1_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_sd1_corr->at(jj));
              }//pt cut                          
           }//sd1
         

       if(jj==0 && g0==1){//sd2
       if(GroomedJet_AK1_2_pt_sd2_corr->at(jj) > JetpTMin && GroomedJet_AK1_2_pt_sd2_corr->at(jj) < JetpTMax){//pt cut
       JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_sd2_corr->at(jj));
                           }//pt cut
                         }//sd2
  
     if(jj==0 && g0==2){//sd3
     if(GroomedJet_AK1_2_pt_sd3_corr->at(jj) > JetpTMin && GroomedJet_AK1_2_pt_sd3_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_sd3_corr->at(jj));
                           }//pt cut
                         }//sd3

     if(jj==0 && g0==3){//sd4
       if(GroomedJet_AK1_2_pt_sd4_corr->at(jj) > JetpTMin && GroomedJet_AK1_2_pt_sd4_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_sd4_corr->at(jj));
                            }//pt cut
                         }//sd4

     if(jj==0 && g0==4){//tr
     if(GroomedJet_AK1_2_pt_ft_corr->at(jj) > JetpTMin && GroomedJet_AK1_2_pt_ft_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_ft_corr->at(jj));
                             }// pt cut
                         }//tr

      if(jj==0 && g0==5){//tr
       if(GroomedJet_AK1_2_pt_tr_corr->at(jj) > JetpTMin && GroomedJet_AK1_2_pt_tr_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_tr_corr->at(jj));
               }//pt cut
                         }//tr

      if(jj==0 && g0==6){//pr
      if(GroomedJet_AK1_2_pt_pr_corr->at(jj) > JetpTMin && GroomedJet_AK1_2_pt_pr_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_pr_corr->at(jj));
                            }
                         }//pr

     if(jj==0 && g0==7){//ungroomed
     if(GroomedJet_AK1_2_pt_corr->at(jj) > JetpTMin && GroomedJet_AK1_2_pt_corr->at(jj) < JetpTMax){//pt cut
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_corr->at(jj));
                             }//pt cut
                         }//ungroomed








           
       }//nPV criteria
}//ak12loop
}//no of grooming loop


///////////////////////////////////////////////////////////////////////////////ak12

   }//event loop

TCanvas *c1=new TCanvas("mycanvas1","My Canvas1");
TCanvas *c2=new TCanvas("mycanvas2","My Canvas2");
TCanvas *c3=new TCanvas("mycanvas3","My Canvas3");

gStyle->SetOptFit(111);
TMultiGraph *mg = new TMultiGraph();
TMultiGraph *mg2 = new TMultiGraph();
TMultiGraph *mg3 = new TMultiGraph();


mg->SetTitle(pTRangeAK8);
//mg->GetYaxis()->SetRangeUser(0, 360.0);
mg2->SetTitle(pTRangeAK10);
mg3->SetTitle(pTRangeAK12);



double x_corr_mc_nPV[siz];
double y_corr_mc_ak8_sd1[siz];
double y_corr_mc_ak8_sd1_err[siz];

double y_corr_mc_ak8_sd2[siz];
double y_corr_mc_ak8_sd2_err[siz];

double y_corr_mc_ak8_sd3[siz];
double y_corr_mc_ak8_sd3_err[siz];

double y_corr_mc_ak8_sd4[siz];
double y_corr_mc_ak8_sd4_err[siz];

double y_corr_mc_ak8_ft[siz];
double y_corr_mc_ak8_ft_err[siz];

double y_corr_mc_ak8_tr[siz];
double y_corr_mc_ak8_tr_err[siz];

double y_corr_mc_ak8_pr[siz];
double y_corr_mc_ak8_pr_err[siz];

double y_corr_mc_ak8_ung[siz];
double y_corr_mc_ak8_ung_err[siz];
//////////////////////////////////ak10start
double y_corr_mc_ak10_sd1[siz];
double y_corr_mc_ak10_sd1_err[siz];

double y_corr_mc_ak10_sd2[siz];
double y_corr_mc_ak10_sd2_err[siz];

double y_corr_mc_ak10_sd3[siz];
double y_corr_mc_ak10_sd3_err[siz];

double y_corr_mc_ak10_sd4[siz];
double y_corr_mc_ak10_sd4_err[siz];

double y_corr_mc_ak10_ft[siz];
double y_corr_mc_ak10_ft_err[siz];

double y_corr_mc_ak10_tr[siz];
double y_corr_mc_ak10_tr_err[siz];

double y_corr_mc_ak10_pr[siz];
double y_corr_mc_ak10_pr_err[siz];

double y_corr_mc_ak10_ung[siz];
double y_corr_mc_ak10_ung_err[siz];
////////////////////////////////ak10end

////////////////////////////////ak12
double y_corr_mc_ak12_sd1[siz];
double y_corr_mc_ak12_sd1_err[siz];

double y_corr_mc_ak12_sd2[siz];
double y_corr_mc_ak12_sd2_err[siz];

double y_corr_mc_ak12_sd3[siz];
double y_corr_mc_ak12_sd3_err[siz];

double y_corr_mc_ak12_sd4[siz];
double y_corr_mc_ak12_sd4_err[siz];

double y_corr_mc_ak12_ft[siz];
double y_corr_mc_ak12_ft_err[siz];

double y_corr_mc_ak12_tr[siz];
double y_corr_mc_ak12_tr_err[siz];

double y_corr_mc_ak12_pr[siz];
double y_corr_mc_ak12_pr_err[siz];

double y_corr_mc_ak12_ung[siz];
double y_corr_mc_ak12_ung_err[siz];

////////////////////////////////ak12


/////////////////////////////////////////////ak8,ak10,ak12
for(int i0=0;i0<nPVf-nPVi+1;i0++){
x_corr_mc_nPV[i0]=i0+nPVi;

////ak8
y_corr_mc_ak8_sd1[i0]=JetmassAK8_h[i0+nPVi][0]->GetMean();
y_corr_mc_ak8_sd1_err[i0]=JetmassAK8_h[i0+nPVi][0]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][0]->GetEntries());


y_corr_mc_ak8_sd2[i0]=JetmassAK8_h[i0+nPVi][1]->GetMean();
y_corr_mc_ak8_sd2_err[i0]=JetmassAK8_h[i0+nPVi][1]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][1]->GetEntries());

y_corr_mc_ak8_sd3[i0]=JetmassAK8_h[i0+nPVi][2]->GetMean();
y_corr_mc_ak8_sd3_err[i0]=JetmassAK8_h[i0+nPVi][2]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][2]->GetEntries());

y_corr_mc_ak8_sd4[i0]=JetmassAK8_h[i0+nPVi][3]->GetMean();
y_corr_mc_ak8_sd4_err[i0]=JetmassAK8_h[i0+nPVi][3]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][3]->GetEntries());

y_corr_mc_ak8_ft[i0]=JetmassAK8_h[i0+nPVi][4]->GetMean();
y_corr_mc_ak8_ft_err[i0]=JetmassAK8_h[i0+nPVi][4]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][4]->GetEntries());


y_corr_mc_ak8_tr[i0]=JetmassAK8_h[i0+nPVi][5]->GetMean();
y_corr_mc_ak8_tr_err[i0]=JetmassAK8_h[i0+nPVi][5]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][5]->GetEntries());

y_corr_mc_ak8_pr[i0]=JetmassAK8_h[i0+nPVi][6]->GetMean();
y_corr_mc_ak8_pr_err[i0]=JetmassAK8_h[i0+nPVi][6]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][6]->GetEntries());

y_corr_mc_ak8_ung[i0]=JetmassAK8_h[i0+nPVi][7]->GetMean();
y_corr_mc_ak8_ung_err[i0]=JetmassAK8_h[i0+nPVi][7]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][7]->GetEntries());
//ak8

//ak10
y_corr_mc_ak10_sd1[i0]=JetmassAK10_h[i0+nPVi][0]->GetMean();
y_corr_mc_ak10_sd1_err[i0]=JetmassAK10_h[i0+nPVi][0]->GetRMS()/sqrt(JetmassAK10_h[i0+nPVi][0]->GetEntries());


y_corr_mc_ak10_sd2[i0]=JetmassAK10_h[i0+nPVi][1]->GetMean();
y_corr_mc_ak10_sd2_err[i0]=JetmassAK10_h[i0+nPVi][1]->GetRMS()/sqrt(JetmassAK10_h[i0+nPVi][1]->GetEntries());

y_corr_mc_ak10_sd3[i0]=JetmassAK10_h[i0+nPVi][2]->GetMean();
y_corr_mc_ak10_sd3_err[i0]=JetmassAK10_h[i0+nPVi][2]->GetRMS()/sqrt(JetmassAK10_h[i0+nPVi][2]->GetEntries());

y_corr_mc_ak10_sd4[i0]=JetmassAK10_h[i0+nPVi][3]->GetMean();
y_corr_mc_ak10_sd4_err[i0]=JetmassAK10_h[i0+nPVi][3]->GetRMS()/sqrt(JetmassAK10_h[i0+nPVi][3]->GetEntries());

y_corr_mc_ak10_ft[i0]=JetmassAK10_h[i0+nPVi][4]->GetMean();
y_corr_mc_ak10_ft_err[i0]=JetmassAK10_h[i0+nPVi][4]->GetRMS()/sqrt(JetmassAK10_h[i0+nPVi][4]->GetEntries());


y_corr_mc_ak10_tr[i0]=JetmassAK10_h[i0+nPVi][5]->GetMean();
y_corr_mc_ak10_tr_err[i0]=JetmassAK10_h[i0+nPVi][5]->GetRMS()/sqrt(JetmassAK10_h[i0+nPVi][5]->GetEntries());

y_corr_mc_ak10_pr[i0]=JetmassAK10_h[i0+nPVi][6]->GetMean();
y_corr_mc_ak10_pr_err[i0]=JetmassAK10_h[i0+nPVi][6]->GetRMS()/sqrt(JetmassAK10_h[i0+nPVi][6]->GetEntries());

y_corr_mc_ak10_ung[i0]=JetmassAK10_h[i0+nPVi][7]->GetMean();
y_corr_mc_ak10_ung_err[i0]=JetmassAK10_h[i0+nPVi][7]->GetRMS()/sqrt(JetmassAK10_h[i0+nPVi][7]->GetEntries());


//ak10


//ak12
y_corr_mc_ak12_sd1[i0]=JetmassAK12_h[i0+nPVi][0]->GetMean();
y_corr_mc_ak12_sd1_err[i0]=JetmassAK12_h[i0+nPVi][0]->GetRMS()/sqrt(JetmassAK12_h[i0+nPVi][0]->GetEntries());


y_corr_mc_ak12_sd2[i0]=JetmassAK12_h[i0+nPVi][1]->GetMean();
y_corr_mc_ak12_sd2_err[i0]=JetmassAK12_h[i0+nPVi][1]->GetRMS()/sqrt(JetmassAK12_h[i0+nPVi][1]->GetEntries());

y_corr_mc_ak12_sd3[i0]=JetmassAK12_h[i0+nPVi][2]->GetMean();
y_corr_mc_ak12_sd3_err[i0]=JetmassAK12_h[i0+nPVi][2]->GetRMS()/sqrt(JetmassAK12_h[i0+nPVi][2]->GetEntries());

y_corr_mc_ak12_sd4[i0]=JetmassAK12_h[i0+nPVi][3]->GetMean();
y_corr_mc_ak12_sd4_err[i0]=JetmassAK12_h[i0+nPVi][3]->GetRMS()/sqrt(JetmassAK12_h[i0+nPVi][3]->GetEntries());

y_corr_mc_ak12_ft[i0]=JetmassAK12_h[i0+nPVi][4]->GetMean();
y_corr_mc_ak12_ft_err[i0]=JetmassAK12_h[i0+nPVi][4]->GetRMS()/sqrt(JetmassAK12_h[i0+nPVi][4]->GetEntries());


y_corr_mc_ak12_tr[i0]=JetmassAK12_h[i0+nPVi][5]->GetMean();
y_corr_mc_ak12_tr_err[i0]=JetmassAK12_h[i0+nPVi][5]->GetRMS()/sqrt(JetmassAK12_h[i0+nPVi][5]->GetEntries());

y_corr_mc_ak12_pr[i0]=JetmassAK12_h[i0+nPVi][6]->GetMean();
y_corr_mc_ak12_pr_err[i0]=JetmassAK12_h[i0+nPVi][6]->GetRMS()/sqrt(JetmassAK12_h[i0+nPVi][6]->GetEntries());

y_corr_mc_ak12_ung[i0]=JetmassAK12_h[i0+nPVi][7]->GetMean();
y_corr_mc_ak12_ung_err[i0]=JetmassAK12_h[i0+nPVi][7]->GetRMS()/sqrt(JetmassAK12_h[i0+nPVi][7]->GetEntries());

//ak12


}

//////////////////////////////////////////////////////////////////ak8,ak10,ak12


///////////////////////////////////ak8 start

TGraphErrors *graph0  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak8_pr,0,y_corr_mc_ak8_pr_err);
        graph0->SetLineColor(2);
        graph0->SetLineWidth(3);
        graph0->GetXaxis()->SetTitle("Number of primary Vertices");
        graph0->GetXaxis()->SetTitleSize(0.055);
        graph0->GetXaxis()->SetLabelSize(0.055);
        graph0->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph0->GetYaxis()->SetTitleSize(0.055);
        graph0->GetYaxis()->SetLabelSize(0.055);
        graph0->SetMarkerStyle(20);
        graph0->GetYaxis()->SetLimits(0, 360.0);





TGraphErrors *graph1  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak8_sd1,0,y_corr_mc_ak8_sd1_err);
        graph1->SetLineColor(3);
        graph1->SetLineWidth(3);
        graph1->GetXaxis()->SetTitle("Number of primary Vertices");
        graph1->GetXaxis()->SetTitleSize(0.055);
        graph1->GetXaxis()->SetLabelSize(0.055);
        graph1->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph1->GetYaxis()->SetTitleSize(0.055);
        graph1->GetYaxis()->SetLabelSize(0.055);
        graph1->SetMarkerStyle(21);
        graph1->GetYaxis()->SetLimits(0, 360.0);
 
TGraphErrors *graph2  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak8_sd2,0,y_corr_mc_ak8_sd2_err);
        graph2->SetLineColor(4);
        graph2->SetLineWidth(3);
        graph2->GetXaxis()->SetTitle("Number of primary Vertices");
        graph2->GetXaxis()->SetTitleSize(0.055);
        graph2->GetXaxis()->SetLabelSize(0.055);
        graph2->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph2->GetYaxis()->SetTitleSize(0.055);
        graph2->GetYaxis()->SetLabelSize(0.055);
        graph2->SetMarkerStyle(22);
        graph2->GetYaxis()->SetLimits(0, 360.0);
        
TGraphErrors *graph3  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak8_sd3,0,y_corr_mc_ak8_sd3_err);
        graph3->SetLineColor(5);
        graph3->SetLineWidth(3);
        graph3->GetXaxis()->SetTitle("Number of primary Vertices");
        graph3->GetXaxis()->SetTitleSize(0.055);
        graph3->GetXaxis()->SetLabelSize(0.055);
        graph3->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph3->GetYaxis()->SetTitleSize(0.055);
        graph3->GetYaxis()->SetLabelSize(0.055);
        graph3->SetMarkerStyle(23);
        graph3->GetYaxis()->SetLimits(0, 360.0); 
TGraphErrors *graph4  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak8_sd4,0,y_corr_mc_ak8_sd4_err);
        graph4->SetLineColor(6);
        graph4->SetLineWidth(3);
        graph4->GetXaxis()->SetTitle("Number of primary Vertices");
        graph4->GetXaxis()->SetTitleSize(0.055);
        graph4->GetXaxis()->SetLabelSize(0.055);
        graph4->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph4->GetYaxis()->SetTitleSize(0.055);
        graph4->GetYaxis()->SetLabelSize(0.055);
        graph4->SetMarkerStyle(24); 
        graph4->GetYaxis()->SetLimits(0, 360.0);        

TGraphErrors *graph5  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak8_ung,0,y_corr_mc_ak8_ung_err);
        graph5->SetLineColor(7);
        graph5->SetLineWidth(3);
        graph5->GetXaxis()->SetTitle("Number of primary Vertices");
        graph5->GetXaxis()->SetTitleSize(0.055);
        graph5->GetXaxis()->SetLabelSize(0.055);
        graph5->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph5->GetYaxis()->SetTitleSize(0.055);
        graph5->GetYaxis()->SetLabelSize(0.055);
        graph5->SetMarkerStyle(25);
        graph5->SetMarkerSize(1);
        graph5->GetYaxis()->SetLimits(0, 360.0);

////////////////////////////////ak8 start



////////////////////////////ak10

TGraphErrors *graph0t  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak10_pr,0,y_corr_mc_ak10_pr_err);
        graph0t->SetLineColor(2);
        graph0t->SetLineWidth(3);
        graph0t->GetXaxis()->SetTitle("Number of primary Vertices");
        graph0t->GetXaxis()->SetTitleSize(0.055);
        graph0t->GetXaxis()->SetLabelSize(0.055);
        graph0t->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph0t->GetYaxis()->SetTitleSize(0.055);
        graph0t->GetYaxis()->SetLabelSize(0.055);
        graph0t->SetMarkerStyle(20);


TGraphErrors *graph1t  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak10_sd1,0,y_corr_mc_ak10_sd1_err);
        graph1t->SetLineColor(3);
        graph1t->SetLineWidth(3);
        graph1t->GetXaxis()->SetTitle("Number of primary Vertices");
        graph1t->GetXaxis()->SetTitleSize(0.055);
        graph1t->GetXaxis()->SetLabelSize(0.055);
        graph1t->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph1t->GetYaxis()->SetTitleSize(0.055);
        graph1t->GetYaxis()->SetLabelSize(0.055);
        graph1t->SetMarkerStyle(21);
TGraphErrors *graph2t  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak10_sd2,0,y_corr_mc_ak10_sd2_err);
        graph2t->SetLineColor(4);
        graph2t->SetLineWidth(3);
        graph2t->GetXaxis()->SetTitle("Number of primary Vertices");
        graph2t->GetXaxis()->SetTitleSize(0.055);
        graph2t->GetXaxis()->SetLabelSize(0.055);
        graph2t->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph2t->GetYaxis()->SetTitleSize(0.055);
        graph2t->GetYaxis()->SetLabelSize(0.055);
        graph2t->SetMarkerStyle(22);
TGraphErrors *graph3t  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak10_sd3,0,y_corr_mc_ak10_sd3_err);
        graph3t->SetLineColor(5);
        graph3t->SetLineWidth(3);
        graph3t->GetXaxis()->SetTitle("Number of primary Vertices");
        graph3t->GetXaxis()->SetTitleSize(0.055);
        graph3t->GetXaxis()->SetLabelSize(0.055);
        graph3t->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph3t->GetYaxis()->SetTitleSize(0.055);
        graph3t->GetYaxis()->SetLabelSize(0.055);
        graph3t->SetMarkerStyle(23);
TGraphErrors *graph4t  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak10_sd4,0,y_corr_mc_ak10_sd4_err);
        graph4t->SetLineColor(6);
        graph4t->SetLineWidth(3);
        graph4t->GetXaxis()->SetTitle("Number of primary Vertices");
        graph4t->GetXaxis()->SetTitleSize(0.055);
        graph4t->GetXaxis()->SetLabelSize(0.055);
        graph4t->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph4t->GetYaxis()->SetTitleSize(0.055);
        graph4t->GetYaxis()->SetLabelSize(0.055);
        graph4t->SetMarkerStyle(24); 

TGraphErrors *graph5t  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak10_ung,0,y_corr_mc_ak10_ung_err);
        graph5t->SetLineColor(7);
        graph5t->SetLineWidth(3);
        graph5t->GetXaxis()->SetTitle("Number of primary Vertices");
        graph5t->GetXaxis()->SetTitleSize(0.055);
        graph5t->GetXaxis()->SetLabelSize(0.055);
        graph5t->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph5t->GetYaxis()->SetTitleSize(0.055);
        graph5t->GetYaxis()->SetLabelSize(0.055);
        graph5t->SetMarkerStyle(25);
        graph5t->SetMarkerSize(1);

///////////////////////////ak10



////////////////////////////////////ak12

TGraphErrors *graph0te  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak12_pr,0,y_corr_mc_ak12_pr_err);
        graph0te->SetLineColor(2);
        graph0te->SetLineWidth(3);
        graph0te->GetXaxis()->SetTitle("Number of primary Vertices");
        graph0te->GetXaxis()->SetTitleSize(0.055);
        graph0te->GetXaxis()->SetLabelSize(0.055);
        graph0te->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph0te->GetYaxis()->SetTitleSize(0.055);
        graph0te->GetYaxis()->SetLabelSize(0.055);
        graph0te->SetMarkerStyle(20);


TGraphErrors *graph1te  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak12_sd1,0,y_corr_mc_ak12_sd1_err);
        graph1te->SetLineColor(3);
        graph1te->SetLineWidth(3);
        graph1te->GetXaxis()->SetTitle("Number of primary Vertices");
        graph1te->GetXaxis()->SetTitleSize(0.055);
        graph1te->GetXaxis()->SetLabelSize(0.055);
        graph1te->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph1te->GetYaxis()->SetTitleSize(0.055);
        graph1te->GetYaxis()->SetLabelSize(0.055);
        graph1te->SetMarkerStyle(21);
TGraphErrors *graph2te  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak12_sd2,0,y_corr_mc_ak12_sd2_err);
        graph2te->SetLineColor(4);
        graph2te->SetLineWidth(3);
        graph2te->GetXaxis()->SetTitle("Number of primary Vertices");
        graph2te->GetXaxis()->SetTitleSize(0.055);
        graph2te->GetXaxis()->SetLabelSize(0.055);
        graph2te->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph2te->GetYaxis()->SetTitleSize(0.055);
        graph2te->GetYaxis()->SetLabelSize(0.055);
        graph2te->SetMarkerStyle(22);
TGraphErrors *graph3te  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak12_sd3,0,y_corr_mc_ak12_sd3_err);
        graph3te->SetLineColor(5);
        graph3te->SetLineWidth(3);
        graph3te->GetXaxis()->SetTitle("Number of primary Vertices");
        graph3te->GetXaxis()->SetTitleSize(0.055);
        graph3te->GetXaxis()->SetLabelSize(0.055);
        graph3te->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph3te->GetYaxis()->SetTitleSize(0.055);
        graph3te->GetYaxis()->SetLabelSize(0.055);
        graph3te->SetMarkerStyle(23);
TGraphErrors *graph4te  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak12_sd4,0,y_corr_mc_ak12_sd4_err);
        graph4te->SetLineColor(6);
        graph4te->SetLineWidth(3);
        graph4te->GetXaxis()->SetTitle("Number of primary Vertices");
        graph4te->GetXaxis()->SetTitleSize(0.055);
        graph4te->GetXaxis()->SetLabelSize(0.055);
        graph4te->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph4te->GetYaxis()->SetTitleSize(0.055);
        graph4te->GetYaxis()->SetLabelSize(0.055);
        graph4te->SetMarkerStyle(24); 

TGraphErrors *graph5te  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak12_ung,0,y_corr_mc_ak12_ung_err);
        graph5te->SetLineColor(7);
        graph5te->SetLineWidth(3);
        graph5te->GetXaxis()->SetTitle("Number of primary Vertices");
        graph5te->GetXaxis()->SetTitleSize(0.055);
        graph5te->GetXaxis()->SetLabelSize(0.055);
        graph5te->GetYaxis()->SetTitle("Average Jet mass in GeV");
        graph5te->GetYaxis()->SetTitleSize(0.055);
        graph5te->GetYaxis()->SetLabelSize(0.055);
        graph5te->SetMarkerStyle(25);
        graph5te->SetMarkerSize(1);


//////////////////////////////////ak12

//graphs.back()->GetXaxis()->SetRangeUser(0, 360.0);

//////////////////////////////////ak8
mg->Add(graph0,"lp");//lp
mg->Add(graph1,"lp");//lp
mg->Add(graph2,"lp");//lp
mg->Add(graph3,"lp");//lp
mg->Add(graph4,"lp");//lp
mg->Add(graph5,"lp");//lp

c1->cd();
mg->Draw("AP");
leg = new TLegend(0.6,0.7,0.89,0.89);
leg->AddEntry(graph0,"AK8 prunned","lp");
leg->AddEntry(graph1,"AK8 softdrop #beta = 0","lp");
leg->AddEntry(graph2,"AK8 softdrop #beta = 2","lp");
leg->AddEntry(graph3,"AK8 softdrop #beta = -1","lp");
leg->AddEntry(graph4,"AK8 softdrop #beta = 1","lp");
leg->AddEntry(graph5,"AK8 Ungroomed ","lp");



leg->Draw();

////////////////////ak10

mg2->Add(graph0t,"lp");//lp
mg2->Add(graph1t,"lp");//lp
mg2->Add(graph2t,"lp");//lp
mg2->Add(graph3t,"lp");//lp
mg2->Add(graph4t,"lp");//lp
mg2->Add(graph5t,"lp");//lp

c2->cd();
mg2->Draw("AP");
leg2 = new TLegend(0.6,0.7,0.89,0.89);
leg2->AddEntry(graph0t,"AK10 prunned","lp");
leg2->AddEntry(graph1t,"AK10 softdrop #beta = 0","lp");
leg2->AddEntry(graph2t,"AK10 softdrop #beta = 2","lp");
leg2->AddEntry(graph3t,"AK10 softdrop #beta = -1","lp");
leg2->AddEntry(graph4t,"AK10 softdrop #beta = 1","lp");
leg2->AddEntry(graph5t,"AK10 Ungroomed ","lp");



leg2->Draw();

///////////////////////ak10

////////////////////////////ak12

mg3->Add(graph0te,"lp");//lp
mg3->Add(graph1te,"lp");//lp
mg3->Add(graph2te,"lp");//lp
mg3->Add(graph3te,"lp");//lp
mg3->Add(graph4te,"lp");//lp
mg3->Add(graph5te,"lp");//lp

c3->cd();
mg3->Draw("AP");
leg3 = new TLegend(0.6,0.7,0.89,0.89);
leg3->AddEntry(graph0te,"AK12 prunned","lp");
leg3->AddEntry(graph1te,"AK12 softdrop #beta = 0","lp");
leg3->AddEntry(graph2te,"AK12 softdrop #beta = 2","lp");
leg3->AddEntry(graph3te,"AK12 softdrop #beta = -1","lp");
leg3->AddEntry(graph4te,"AK12 softdrop #beta = 1","lp");
leg3->AddEntry(graph5te,"AK12 Ungroomed ","lp");



leg3->Draw();


//////////////////////////ak12




}//main programme
