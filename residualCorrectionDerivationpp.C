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
#include "TStyle.h"

void residualCorrection(int cone = 3)
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  const int najBins = 4;
  const int nptBins = 5;

  const double ajBins[najBins+1] = {0,0.11,0.22,0.33,1};
  const double ptBins[nptBins+1] = {0.5,1,2,4,8,300};

  TFile * f= new TFile(Form("/net/hidsk0001/d00/scratch/dgulhan/ntuples_res_trk_pp/track_ntuple_HiForest_Private_PYTHIA_localdb_ppJEC_pthat30_merged_forest_0_ak%dCalo.root",cone),"read");
  TTree * nt_track = (TTree*)f->Get("nt_track");
  TTree * nt_particle = (TTree*)f->Get("nt_particle");
  
  TH2D * leadCorr = new TH2D("leadCorr",":A_{j}:p_{t}^{track}",najBins,ajBins,nptBins,ptBins);
  TH2D * subleadCorr = new TH2D("subleadCorr",":A_{j}:p_{t}^{track}",najBins,ajBins,nptBins,ptBins);

  TH2D * leadCorr_gen = new TH2D("leadCorr_gen",":A_{j}:p_{t}^{track}",najBins,ajBins,nptBins,ptBins);
  TH2D * subleadCorr_gen = new TH2D("subleadCorr_gen",":A_{j}:p_{t}^{track}",najBins,ajBins,nptBins,ptBins);

  //removed "w_vz*weight" factor
  TCut cutReco_lead = "trackselect && pt>0.5 && r_lead<0.2";
  TCut cutReco_sublead = "trackselect && pt>0.5 && r_sublead<0.2";
  nt_track->Draw("pt:-(asym)>>leadCorr","((1-secondary)*(1-fake)/((1+multrec)*eff))"*(cutReco_lead),"colz");
  nt_track->Draw("pt:-(asym)>>subleadCorr","((1-secondary)*(1-fake)/((1+multrec)*eff))"*(cutReco_sublead),"colz");

  TCut cutGen_lead = "pt>0.5 && r_lead<0.2";
  TCut cutGen_sublead = "pt>0.5 && r_sublead<0.2";
  nt_particle->Draw("pt:-(asym)>>leadCorr_gen",(cutGen_lead),"colz");
  nt_particle->Draw("pt:-(asym)>>subleadCorr_gen",(cutGen_sublead),"colz");

  leadCorr->Divide(leadCorr_gen);
  subleadCorr->Divide(subleadCorr_gen);

  TH2D * diff = (TH2D*) leadCorr->Clone("diff");
  diff->Add(subleadCorr,-1);

  TLegend * l1 = new TLegend(0.37,0.5,0.77,0.9);
  l1->AddEntry((TObject*)0,"r_{lead}<0.2 Closure","");
  l1->SetBorderSize(0);  
  l1->SetTextSize(0.04);

  TLegend * l2 = new TLegend(0.37,0.5,0.77,0.9);
  l2->AddEntry((TObject*)0,"r_{sublead}<0.2 Closure","");
  l2->SetBorderSize(0);
  l2->SetTextSize(0.04);

  TLegend * l3 = new TLegend(0.37,0.5,0.77,0.9);
  l3->AddEntry((TObject*)0,"Difference in Closures","");
  l3->SetBorderSize(0);
  l3->SetTextSize(0.04);

  TCanvas * c2 = new TCanvas("c2","c2",500,500);
  c2->SetRightMargin(0.15);
  c2->SetLogy();
  leadCorr->GetYaxis()->SetTitle("p_{t}");
  leadCorr->GetXaxis()->SetTitle("A_{j}");
  leadCorr->Draw("colz");
  l1->Draw("same");
  c2->SaveAs(Form("plots/leading_correction_AkVs%dCalo_pp.pdf",cone));
  c2->SaveAs(Form("plots/leading_correction_AkVs%dCalo_pp.png",cone));

  TCanvas * c3 = new TCanvas("c3","c3",500,500);
  c3->SetRightMargin(0.15);
  c3->SetLogy();
  subleadCorr->GetYaxis()->SetTitle("p_{t}");
  subleadCorr->GetXaxis()->SetTitle("A_{j}");
  subleadCorr->Draw("colz");
  l2->Draw("same");
  c3->SaveAs(Form("plots/subleading_correction_AkVs%dCalo_pp.pdf",cone));
  c3->SaveAs(Form("plots/subleading_correction_AkVs%dCalo_pp.png",cone));

  TCanvas * c4 = new TCanvas("c4","c4",500,500);
  c4->SetRightMargin(0.15);
  c4->SetLogy();
  diff->GetYaxis()->SetTitle("p_{t}");
  diff->GetXaxis()->SetTitle("A_{j}");
  diff->Draw("colz");
  l3->Draw("same");
  c4->SaveAs(Form("plots/difference_AkVs%dCalo_pp.pdf",cone));
  c4->SaveAs(Form("plots/difference_AkVs%dCalo_pp.png",cone));

  TFile * outf = new TFile(Form("rootFiles/ResidualTrkCorr_AkVs%dCalo_pp.root",cone),"recreate");
  leadCorr->SetDirectory(0);
  subleadCorr->SetDirectory(0);
  leadCorr->Write();
  subleadCorr->Write();
  diff->Write();
}
