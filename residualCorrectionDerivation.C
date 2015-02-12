#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TLine.h"
#include "TCut.h"

void residualCorrection()
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  const int najBins = 4;
  const int nptBins = 5;

  const double ajBins[najBins+1] = {0,0.11,0.22,0.33,1};
  const double ptBins[nptBins+1] = {0.5,1,2,4,8,300};

  TFile * f= new TFile("/export/d00/scratch/abaty/trackingEff/closure_ntuples/Correction_Vs5Calo_ntuple_dijet_etaLT16.root","read");
  TTree * nt_track = (TTree*)f->Get("nt_track");
  TTree * nt_particle = (TTree*)f->Get("nt_particle");
  
  TH2D * leadCorr = new TH2D("leadCorr",":A_{j}:p_{t}^{track}",najBins,ajBins,nptBins,ptBins);
  TH2D * subleadCorr = new TH2D("subleadCorr",":A_{j}:p_{t}^{track}",najBins,ajBins,nptBins,ptBins);

  TH2D * leadCorr_gen = new TH2D("leadCorr_gen",":A_{j}:p_{t}^{track}",najBins,ajBins,nptBins,ptBins);
  TH2D * subleadCorr_gen = new TH2D("subleadCorr_gen",":A_{j}:p_{t}^{track}",najBins,ajBins,nptBins,ptBins);

  TCut cutReco_lead = "trackselect && pt>0.5 && r_lead<0.2";
  TCut cutReco_sublead = "trackselect && pt>0.5 && r_sublead<0.2";
  nt_track->Draw("pt:asym>>leadCorr","((1-fake)/eff)*weight"*(cutReco_lead),"colz");
  nt_track->Draw("pt:asym>>subleadCorr","((1-fake)/eff)*weight"*(cutReco_sublead),"colz");

  TCut cutGen_lead = "pt>0.5 && r_lead<0.2";
  TCut cutGen_sublead = "pt>0.5 && r_sublead<0.2";
  nt_particle->Draw("pt:asym>>leadCorr_gen","weight"*(cutGen_lead),"colz");
  nt_particle->Draw("pt:asym>>subleadCorr_gen","weight"*(cutGen_sublead),"colz");

  leadCorr->Divide(leadCorr_gen);
  subleadCorr->Divide(subleadCorr_gen);

  TH2D * diff = (TH2D*) leadCorr->Clone("diff");
  diff->Add(subleadCorr,-1);

  TCanvas * c2 = new TCanvas("c2","c2",500,500);
  c2->SetRightMargin(0.15);
  c2->SetLogy();
  leadCorr->GetYaxis()->SetTitle("p_{t}");
  leadCorr->GetXaxis()->SetTitle("A_{j}");
  leadCorr->Draw("colz");
  c2->SaveAs("plots/leading_correction_AkVs5Calo.png");

  TCanvas * c3 = new TCanvas("c3","c3",500,500);
  c3->SetRightMargin(0.15);
  c3->SetLogy();
  subleadCorr->GetYaxis()->SetTitle("p_{t}");
  subleadCorr->GetXaxis()->SetTitle("A_{j}");
  subleadCorr->Draw("colz");
  c3->SaveAs("plots/subleading_correction_AkVs5Calo.png");

  TCanvas * c4 = new TCanvas("c4","c4",500,500);
  c4->SetRightMargin(0.15);
  c4->SetLogy();
  diff->GetYaxis()->SetTitle("p_{t}");
  diff->GetXaxis()->SetTitle("A_{j}");
  diff->Draw("colz");
  c4->SaveAs("plots/difference_AkVs5Calo.png");

  TFile * outf = new TFile("rootFiles/ResidualTrkCorr_AkVs5Calo_eta16.root","recreate");
  leadCorr->SetDirectory(0);
  subleadCorr->SetDirectory(0);
  leadCorr->Write();
  subleadCorr->Write();
}
