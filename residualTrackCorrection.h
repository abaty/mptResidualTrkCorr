#include "TFile.h"
#include "TH2D.h"
#include <iostream>

TFile * residualTrackCorrFile;
TH2D * leadCorr;
TH2D * subleadCorr;

//Usage: Call residualTrkCorrInit
//Call getTrkResidualCorr() for each trk to get a factor to multiply by

void residualTrkCorrInit(int coneSize = 3, bool isCalo = true, bool useFragJEC = true)
{
//some exception handling for bad inputs
  if(coneSize > 5 || coneSize<2)
  {
    std::cout << "Cone Size Specified not in appropriate range, defaulting to R=0.3" << std::endl;
    coneSize = 3;
  }
  if(!isCalo)
  {
    std::cout << "PF jets selected, defaulting to Cone size of R=0.3" << std::endl;
    coneSize = 3;
  }
  if(!useFragJEC)
  {
    std::cout << "Fragmentation JEC not used, defaulting to Cone size of R=0.3" << std::endl;
  }

//opening files  
  if(isCalo && useFragJEC)
  {
    residualTrackCorrFile = new TFile(Form("residualTrackCorrection/ResidualTrkCorr_AkVs%dCalo_eta16.root",coneSize),"read");
  }
  if(!isCalo)
  {
    residualTrackCorrFile = new TFile(Form("residualTrackCorrection/ResidualTrkCorr_AkVs3PF_NoFFJEC_eta16.root",coneSize),"read");
  }
  if(!useFragJEC && isCalo)
  {
    residualTrackCorrFile = new TFile(Form("residualTrackCorrection/ResidualTrkCorr_AkVs3Calo_NoFFJEC_eta16.root",coneSize),"read");
  }

//getting correction files
  leadCorr = (TH2D*) residualTrackCorrFile->Get("leadCorr");
  subleadCorr = (TH2D*) residualTrackCorrFile->Get("subleadCorr");

  return;
}

double getTrkResidualCorr(double r_lead, double r_sublead, double Aj, double pt)
{
  double correction = 1;
  if(Aj<0 || Aj>1) std::cout << "Error: Aj is not in the range [0,1]" << std::endl;
  if(pt<0.5)
  {
    std::cout << "Residual track Correction not calculated for pt<0.5; returning correction of 1" << std::endl;
    return correction;
  }
  //handling an exception where the track is out of bounds in pt
  if(pt >= 300) pt = 299;

  if(r_lead < 0.2)          correction = 1.0/leadCorr->GetBinContent(leadCorr->FindBin(Aj,pt));
  else if(r_sublead < 0.2)  correction = 1.0/subleadCorr->GetBinContent(subleadCorr->FindBin(Aj,pt));

  return correction;
}
