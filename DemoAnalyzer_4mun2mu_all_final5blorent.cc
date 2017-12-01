// Code continue from code 5b cuma tukar method guna tlorentzvector
// Code ni buang semua tight cuts and election
// If the result is the same as result code 5a (which means i dont fuck up my code hnya sbb buang tight cuts), i try pilih mZa close to Z mass instead of heavier mass

// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <utility>

// user include files, general
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//------ EXTRA HEADER FILES--------------------//
#include "math.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"

// for Root histogramming
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"

// for tracking information
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

// for vertex information 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// for muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"


//for beamspot information
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//for electron informaton
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

// for particle flow information
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

// class declaration

class DemoAnalyzer: public edm::EDAnalyzer {
public:
  explicit DemoAnalyzer(const edm::ParameterSet&);
  ~DemoAnalyzer();

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  bool providesGoodLumisection(const edm::Event& iEvent);

  // Declare Root histograms or tree
  // For a description of their content see below

  TTree *t1;
  TTree *t2;
  TTree *t3;

  TH1D *h_globalmu_size;
  TH1D *h_tmu_size;
  TH1D *h_e_size;

  TH1D *h_nggmu;
  TH1D *h_ngmu_jpsi;
  TH1D *h_ngmu_loose;
  TH1D *h_nge_loose;


  TH1D *h_m1_gmu;
  TH1D *h_m2_gmu;
  TH1D *h_m3_gmu;
  TH1D *h_m_2e_loose;

  TH1D *h_m1_jpsi;
  TH1D *h_m2_jpsi;
  TH1D *h_m3_jpsi;
  TH1D *h_m4_jpsi;
  TH1D *h_m5_jpsi;
  TH1D *h_m6_jpsi;

  TH1D *h_mZ_2mu_loose;

  TH1D *h_mZ_2e_loose;

  TH1D *h_mZ12_4mu_loose;
  TH1D *h_mZ34_4mu_loose;
  TH1D *h_mZ13_4mu_loose;
  TH1D *h_mZ24_4mu_loose;
  TH1D *h_mZ14_4mu_loose;
  TH1D *h_mZ23_4mu_loose;
  TH1D *h_mZa_4mu_loose;
  TH1D *h_mZb_4mu_loose;
  TH1D *h_m1_m4mu_loose;
  TH1D *h_m2_m4mu_loose;
  TH1D *h_m3_m4mu_loose;
  TH1D *h_m4_m4mu_loose;

  TH1D *h_mZ12_4e_loose; 
  TH1D *h_mZ34_4e_loose;
  TH1D *h_mZ13_4e_loose;
  TH1D *h_mZ24_4e_loose;
  TH1D *h_mZ14_4e_loose;
  TH1D *h_mZ23_4e_loose;
  TH1D *h_mZa_4e_loose;
  TH1D *h_mZb_4e_loose;
  TH1D *h_m1_m4e_loose;
  TH1D *h_m2_m4e_loose;
  TH1D *h_m3_m4e_loose;
  TH1D *h_m4_m4e_loose;

  TH1D *h_mZmu_2mu2e_loose;
  TH1D *h_mZe_2mu2e_loose;
  TH1D *h_mZa_2mu2e_loose;
  TH1D *h_mZb_2mu2e_loose;
  TH1D *h_m1_m2mu2e_loose;
  TH1D *h_m2_m2mu2e_loose;
  TH1D *h_m3_m2mu2e_loose;
  TH1D *h_m4_m2mu2e_loose;

  // Control Plot

  // Global Muon
  TH1D *h_p_gmu;
  TH1D *h_pt_gmu_b4;
  TH1D *h_eta_gmu_b4;
  TH1D *h_chi2_gmu;
  TH1D *h_phi_gmu;
  TH1D *h_ndof_gmu;
  TH1D *h_normchi2_gmu;
  TH1D *h_validhits_gmu;
  TH1D *h_pixelhits_gmu;

  TH1D *h_pt_gmu_after;
  TH1D *h_eta_gmu_after;

  // Tracker Muon
  TH1D *h_p_tmu;
  TH1D *h_pt_tmu_b4;
  TH1D *h_eta_tmu_b4;
  TH1D *h_phi_tmu;
  TH1D *h_chi2_tmu;
  TH1D *h_ndof_tmu;
  TH1D *h_normchi2_tmu;
  TH1D *h_validhits_tmu_jpsi; // Jpsi
  TH1D *h_pixelhits_tmu_jpsi;

  TH1D *h_pt_tmu_jpsi_after;
  TH1D *h_eta_tmu_jpsi_after;

  // GlobalTracker Muon
  TH1D *h_goodhit;
  TH1D *h_dxy_mu;
  TH1D *h_mustation;
  TH1D *h_relPFIso_mu;
  TH1D *h_Isotrk_mu;

  TH1D *h_relPFIso_mu_loose_after;
  TH1D *h_validhits_mu;
  TH1D *h_pixelhits_mu;

  TH1D *h_pt_tmu_loose_after_Zto2mu;
  TH1D *h_eta_tmu_loose_after_Zto2mu;

  TH1D *h_pt_tmu_loose_after;
  TH1D *h_eta_tmu_loose_after;
  TH1D *h_pt_tmu_loose_after_2mu2e;
  TH1D *h_eta_tmu_loose_after_2mu2e;

  // Electron
  TH1D *h_p_e;
  TH1D *h_et_e;
  TH1D *h_pt_e_b4;
  TH1D *h_eta_e_b4;
  TH1D *h_phi_e;
  TH1D *h_suclust_eta;
  TH1D *h_suclust_rawE;
  TH1D *h_relPFIso_e;
  TH1D *h_Isotrk_e;
  TH1D *h_relPFIso_e_loose_after;
  TH2D *h_relPFIso_pt_e;

  TH1D *h_dxy_e;

  TH1D *h_pt_e_loose_after_Zto2e;
  TH1D *h_eta_e_loose_after_Zto2e;

  TH1D *h_pt_e_loose_after;
  TH1D *h_eta_e_loose_after;

  TH1D *h_relPFIso_2mu_loose_after;
  TH1D *h_relPFIso_2e_loose_after;

  TH1D *h_pt_e_loose_after_2mu2e;
  TH1D *h_eta_e_loose_after_2mu2e;

  TH1D *h_SIP3d_mu_b4;
  TH1D *h_SIP3d_e_b4;
  TH1D *h_misshite;


  // Declare variables
  int  nGoodGlobalMuon, nGoodMuonJpsi, nGoodMuonLoose;

  double s1, s2, s3, s4, s;
  double dx,dy,dz, rap, pz;

  double mZ12, mZ34, mZ13, mZ24, mZ14, mZ23;
  double mass4mu, pt_4mu, eta_4mu;
  double pt_mu1, pt_mu2, pt_mu3, pt_mu4;
  double eta_mu1, eta_mu2, eta_mu3, eta_mu4;
  double mZa, mZb;

  double sqm1, mZ, eZ12, eZ34, eZ13, eZ24, eZ14, eZ23; 
  double pxZ12, pxZ34, pxZ13, pxZ24,  pxZ14, pxZ23;
  double pyZ12, pyZ34, pyZ13, pyZ24, pyZ14, pyZ23;
  double pzZ12, pzZ34, pzZ13, pzZ24, pzZ14, pzZ23;
  double pZ12, pZ34, pZ13, pZ24, pZ14, pZ23;
  double pTZ12, pTZ34, pTZ13, pTZ24, pTZ14, pTZ23;

  double dZ12, dZ34, dZ13, dZ24, dZ14, dZ23;
  double dZc1, dZc2, dZc3;
  double eZa, pxZa, pyZa, pzZa, pTZa;
  double eZb, pxZb, pyZb, pzZb, pTZb;

  int nGoodElectronLoose;
  double sqme;
  int misshits;

  double mass4e, pt_4e, eta_4e;
  double pt_e1, pt_e2, pt_e3, pt_e4;
  double eta_e1, eta_e2, eta_e3, eta_e4;
  
  double mass2mu2e, pt_2mu2e, eta_2mu2e;;
  double pt_2mu1, pt_2mu2, pt_2e1, pt_2e2;
  double eta_2mu1, eta_2mu2, eta_2e1, eta_2e2;
  
  double goodhit;
  double relPFIso_mu;
  double relPFIso_e;

  double Isotrk_mu;
  double Isotrk_e;

  double IP3d_mu;
  double ErrIP3d_mu;
  double SIP3d_mu;
  double IP3d_e;
  double ErrIP3d_e;
  double SIP3d_e;

  //  int nRun1,nEvt1,nLumi1;
  int nRun,nEvt,nLumi;

  TLorentzVector p4Za, p4Zb, p4H;
  
  // mu12.SetPtEtaPhiE();

  // const int arr = 4;
  // const int type[arr] = {};

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig) {

  // *****************************************************************
  // This is the main analysis routine
  // The goal is to approximately reproduce the dimuon mass spectrum  
  // from MUO-10-004
  // *****************************************************************

  // now do what ever initialization is needed
  edm::Service<TFileService> fs;

  // ************************************
  // book histograms and set axis labels
  // (called once for initialization)
  // ************************************

  // Global Muon size
  h_globalmu_size = fs->make<TH1D>("NGMuons", "GMuon Size", 10, 0., 10.);
  h_globalmu_size->GetXaxis()->SetTitle("Number of GMuons");
  h_globalmu_size->GetYaxis()->SetTitle("Number of Events");

  // Muon size
  h_tmu_size = fs->make<TH1D>("NMuons", "Muon Size", 10, 0., 10.);
  h_tmu_size->GetXaxis()->SetTitle("Number of Muons");
  h_tmu_size->GetYaxis()->SetTitle("Number of Events");

  // Electron size
  h_e_size = fs->make<TH1D>("Nelectrons", "Electron Size", 10, 0., 10.);
  h_e_size->GetXaxis()->SetTitle("Number of Electrons");
  h_e_size->GetYaxis()->SetTitle("Number of Events");

  // No. of Good Global Muon
  h_nggmu = fs->make<TH1D>("NGoodGMuons", "No. of Good Global Muons", 10, 0., 10.);
  h_nggmu->GetXaxis()->SetTitle("Number of GMuons");
  h_nggmu->GetYaxis()->SetTitle("Number of Events");

  // No. of Good Tracker Muon
  h_ngmu_jpsi = fs->make<TH1D>("NGoodMuonsJpsi", "No. of Good Tracker Muons for Jpsi", 10, 0., 10.);
  h_ngmu_jpsi->GetXaxis()->SetTitle("Number of Muons");
  h_ngmu_jpsi->GetYaxis()->SetTitle("Number of Events");

  // No. of Good GlobalTracker Muon with loose cuts 
  h_ngmu_loose = fs->make<TH1D>("NGoodMuonsLoose", "No. of Good GlobalTracker Muons with looser cuts", 10, 0., 10.);
  h_ngmu_loose->GetXaxis()->SetTitle("Number of Muons");
  h_ngmu_loose->GetYaxis()->SetTitle("Number of Events");

  // No. of Good Electron with loose cuts 
  h_nge_loose = fs->make<TH1D>("NGoodElectronLoose", "No. of Good Electron", 10, 0., 10.);
  h_nge_loose->GetXaxis()->SetTitle("Number of Electrons");
  h_nge_loose->GetYaxis()->SetTitle("Number of Events");

  // Dimuon mass spectrum up to 4 GeV (low mass range, rho/omega, phi, psi) in Global Muon
  h_m1_gmu = fs->make<TH1D>("GMmass", "Global Muon mass", 400, 0., 4.);
  h_m1_gmu->GetXaxis()->SetTitle("Invariant Mass for Nmuon=2 (in GeV/c^2)");
  h_m1_gmu->GetYaxis()->SetTitle("Number of Events");

  // Dimuon mass spectrum up to 120 GeV (high mass range: upsilon, Z) in Global Muon
  h_m2_gmu = fs->make<TH1D>("GMmass_extended", "GMmass", 240, 0., 120.);
  h_m2_gmu->GetXaxis()->SetTitle("Invariant Mass for Nmuon=2 (in GeV/c^2)");
  h_m2_gmu->GetYaxis()->SetTitle("Number of Events");

  // Dimuon mass spectrum up to 600 GeV in Global Muon
  h_m3_gmu = fs->make<TH1D>("GMmass_extended_600", "GMmass", 240, 0., 600.);
  h_m3_gmu->GetXaxis()->SetTitle("Invariant Mass for Nmuon=2 (in GeV/c^2)");
  h_m3_gmu->GetYaxis()->SetTitle("Number of Events");

  // Dielectron mass spectrum for loose cut
  h_m_2e_loose = fs->make<TH1D>("Dielectron_mass_loose", "dielectron_mass with isolation", 240, 0., 120.);
  h_m_2e_loose->GetXaxis()->SetTitle("Invariant Mass for Nelectron=2 (in GeV/c^2)");
  h_m_2e_loose->GetYaxis()->SetTitle("Number of Events");

  // Jpsi mass spectrum with different rapidity range in Tracker Muon
  h_m1_jpsi = fs->make<TH1D>("TM_mass1_j_psi", "TM_mass_jpsi_rap_1.2", 120, 0., 20.);
  h_m1_jpsi->GetXaxis()->SetTitle("Invariant Mass for Nmuon=2 with|Y|<1.2(in GeV/c^2)");
  h_m1_jpsi->GetYaxis()->SetTitle("Number of Events");

  h_m2_jpsi = fs->make<TH1D>("TM_mass2_j_psi", "TM_mass_jpsi_rap_1.6", 120, 0., 20.);
  h_m2_jpsi->GetXaxis()->SetTitle("Invariant Mass for Nmuon=2 with |Y|>1.2 and |Y|<1.6 (in GeV/c^2)");
  h_m2_jpsi->GetYaxis()->SetTitle("Number of Events");

  h_m3_jpsi = fs->make<TH1D>("TM_mass3_j_psi", "TM_mass_jpsi_rap_2.4", 120, 0., 20.);
  h_m3_jpsi->GetXaxis()->SetTitle("Invariant Mass for Nmuon=2 with |Y|>1.6 and |Y|<2.4 (in GeV/c^2)");
  h_m3_jpsi->GetYaxis()->SetTitle("Number of Events");

  // Same as above just different mass range and bin
  h_m4_jpsi = fs->make<TH1D>("TM_mass4_j_psi", "TM_mass_jpsi_rap_1.2", 90, 2.6, 3.5);
  h_m4_jpsi->GetXaxis()->SetTitle("Invariant Mass for Nmuon=2 with|Y|<1.2(in GeV/c^2)");
  h_m4_jpsi->GetYaxis()->SetTitle("Number of Events");

  h_m5_jpsi = fs->make<TH1D>("TM_mass5_j_psi", "TM_mass_jpsi_rap_1.6", 90, 2.6, 3.5);
  h_m5_jpsi->GetXaxis()->SetTitle("Invariant Mass for Nmuon=2 with |Y|>1.2 and |Y|<1.6 (in GeV/c^2)");
  h_m5_jpsi->GetYaxis()->SetTitle("Number of Events");

  h_m6_jpsi = fs->make<TH1D>("TM_mass6_j_psi", "TM_mass_jpsi_rap_2.4", 90, 2.6, 3.5);
  h_m6_jpsi->GetXaxis()->SetTitle("Invariant Mass for Nmuon=2 with |Y|>1.6 and |Y|<2.4 (in GeV/c^2)");
  h_m6_jpsi->GetYaxis()->SetTitle("Number of Events");

  // ZTo2mu mass spectrum in GlobalTracker Muon for loose cut
  h_mZ_2mu_loose = fs->make<TH1D>("massZto2muon_loose", "mass of Z to 2 muon (loose)", 120, 40., 120.);
  h_mZ_2mu_loose->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ_2mu_loose->GetYaxis()->SetTitle("Number of Events");

  // ZTo2e mass spectrum in Gsf electron for loose cut
  h_mZ_2e_loose = fs->make<TH1D>("massZto2e_loose", "mass of Z to 2e (loose)", 120, 40., 120.);
  h_mZ_2e_loose->GetXaxis()->SetTitle("Invariant Mass for Nelectron=2 (in GeV/c^2)");
  h_mZ_2e_loose->GetYaxis()->SetTitle("Number of Events");

  // These histograms are for ZZ/ZZ*To4mu reconstruction with different combination and with loose cuts
  // First combination: 1234
  // Mass Z12
  h_mZ12_4mu_loose = fs->make<TH1D>("mZ12_4mu_loose", "mass of Z12", 75, 0., 150.);
  h_mZ12_4mu_loose->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ12_4mu_loose->GetYaxis()->SetTitle("Number of Events");

  // Mass Z34
  h_mZ34_4mu_loose = fs->make<TH1D>("mZ34_4mu_loose", "mass of Z34", 75, 0., 150.);
  h_mZ34_4mu_loose->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ34_4mu_loose->GetYaxis()->SetTitle("Number of Events");

  // Second combination: 1324
  // Mass Z13
  h_mZ13_4mu_loose = fs->make<TH1D>("mZ13_4mu_loose", "mass of Z13", 75, 0., 150.);
  h_mZ13_4mu_loose->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ13_4mu_loose->GetYaxis()->SetTitle("Number of Events");

  // Mass Z24
  h_mZ24_4mu_loose = fs->make<TH1D>("mZ24_4mu_loose", "mass of Z24", 75, 0., 150.);
  h_mZ24_4mu_loose->GetXaxis()->SetTitle("Invariant mass of dimuon");
  h_mZ24_4mu_loose->GetYaxis()->SetTitle("Number of Events");

  // Third combination: 1423
  // Mass Z14
  h_mZ14_4mu_loose = fs->make<TH1D>("mZ14_4mu_loose", "mass of Z14", 75, 0., 150.);
  h_mZ14_4mu_loose->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ14_4mu_loose->GetYaxis()->SetTitle("Number of Events");

  // Mass Z23
  h_mZ23_4mu_loose = fs->make<TH1D>("mZ23_4mu_loose", "mass of Z23", 75, 0., 150.);
  h_mZ23_4mu_loose->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZ23_4mu_loose->GetYaxis()->SetTitle("Number of Events");

  // Mass Za: mass of ZTo2mu with higher mass
  h_mZa_4mu_loose = fs->make<TH1D>("mZa_4mu_loose", "mass Z higher", 120, 0., 120.);
  h_mZa_4mu_loose->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZa_4mu_loose->GetYaxis()->SetTitle("Number of Events");

  // Mass Zb: mass of ZTo2mu with lower mass
  h_mZb_4mu_loose = fs->make<TH1D>("mZb_4mu_loose", "mass Z lower", 120, 0., 120.);
  h_mZb_4mu_loose->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZb_4mu_loose->GetYaxis()->SetTitle("Number of Events");

  // 4muon mass spectrum with loose cut in GlobalTracker Muon paper 7TeV
  h_m1_m4mu_loose = fs->make<TH1D>("mass4mu_loose_7TeV", "mass of 4 muon", 51, 98., 608.);
  h_m1_m4mu_loose->GetXaxis()->SetTitle("Invariant mass of 4muons (in GeV/c^2)");
  h_m1_m4mu_loose->GetYaxis()->SetTitle("Number of Events");

  // 4muon mass spectrum with loose cut in GlobalTracker Muon paper 8TeV
  h_m2_m4mu_loose = fs->make<TH1D>("mass4mu_loose_8TeV", "mass of 4 muon", 74, 70., 810.);
  h_m2_m4mu_loose->GetXaxis()->SetTitle("Invariant mass of 4muons (in GeV/c^2)");
  h_m2_m4mu_loose->GetYaxis()->SetTitle("Number of Events");

  // 4muon mass spectrum with loose cut in GlobalTracker Muon paper 8TeV lower range
  h_m3_m4mu_loose = fs->make<TH1D>("mass4mu_loose_8TeV_low", "mass of 4 muon", 37, 70., 181.);
  h_m3_m4mu_loose->GetXaxis()->SetTitle("Invariant mass of 4muons (in GeV/c^2)");
  h_m3_m4mu_loose->GetYaxis()->SetTitle("Number of Events");

  // 4muon mass spectrum with loose cut in GlobalTracker Muon full range
  h_m4_m4mu_loose = fs->make<TH1D>("mass4mu_loose_full", "mass of 4 muon", 300, 0., 900.);
  h_m4_m4mu_loose->GetXaxis()->SetTitle("Invariant mass of 4muons (in GeV/c^2)");
  h_m4_m4mu_loose->GetYaxis()->SetTitle("Number of Events");

  // These histograms are for 4electrons reconstruction with different combination and loose cut
  // First combination: 1234
  // Mass Z12
  h_mZ12_4e_loose = fs->make<TH1D>("mZ12_4e_loose", "mass of Z12", 75, 0., 150.);
  h_mZ12_4e_loose->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ12_4e_loose->GetYaxis()->SetTitle("Number of Events");

  // Mass Z34
  h_mZ34_4e_loose = fs->make<TH1D>("mZ34_4e_loose", "mass of Z34", 75, 0., 150.);
  h_mZ34_4e_loose->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ34_4e_loose->GetYaxis()->SetTitle("Number of Events");

  // Second combination: 1324
  // Mass Z13
  h_mZ13_4e_loose = fs->make<TH1D>("mZ13_4e_loose", "mass of Z13", 75, 0., 150.);
  h_mZ13_4e_loose->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ13_4e_loose->GetYaxis()->SetTitle("Number of Events");

  // Mass Z24
  h_mZ24_4e_loose = fs->make<TH1D>("mZ24_4e_loose", "mass of Z24", 75, 0., 150.);
  h_mZ24_4e_loose->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ24_4e_loose->GetYaxis()->SetTitle("Number of Events");

  // Third combination: 1423
  // Mass Z14
  h_mZ14_4e_loose = fs->make<TH1D>("mZ14_4e_loose", "mass of Z14", 75, 0., 150.);
  h_mZ14_4e_loose->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ14_4e_loose->GetYaxis()->SetTitle("Number of Events");

  // Mass Z23
  h_mZ23_4e_loose = fs->make<TH1D>("mZ23_4e_loose", "mass of Z23", 75, 0., 150.);
  h_mZ23_4e_loose->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZ23_4e_loose->GetYaxis()->SetTitle("Number of Events");

  // Mass Za: mass of Z with higher mass
  h_mZa_4e_loose = fs->make<TH1D>("mZa_4e_loose", "mass Z higher", 120, 0., 120.);
  h_mZa_4e_loose->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZa_4e_loose->GetYaxis()->SetTitle("Number of Events");

  // Mass Zb: mass of Z with lower mass
  h_mZb_4e_loose = fs->make<TH1D>("mZb_4e_loose", "mass Z lower", 120, 0., 120.);
  h_mZb_4e_loose->GetXaxis()->SetTitle("Invariant mass of dielectron (in GeV/c^2)");
  h_mZb_4e_loose->GetYaxis()->SetTitle("Number of Events");

  // 4electron mass spectrum with loose cut in Gsf Electron paper 7TeV
  h_m1_m4e_loose = fs->make<TH1D>("mass4e_loose_7TeV", "mass of 4 electrons", 51, 98., 608.);
  h_m1_m4e_loose->GetXaxis()->SetTitle("Invariant mass of 4e (in GeV/c^2)");
  h_m1_m4e_loose->GetYaxis()->SetTitle("Number of Events");

  // 4electron mass spectrum with loose cut in Gsf Electron paper 8TeV
  h_m2_m4e_loose = fs->make<TH1D>("mass4e_loose_8TeV", "mass of 4 electrons", 74, 70., 810.);
  h_m2_m4e_loose->GetXaxis()->SetTitle("Invariant mass of 4e (in GeV/c^2)");
  h_m2_m4e_loose->GetYaxis()->SetTitle("Number of Events");

  // 4electron mass spectrum with loose cut in Gsf Electron paper 8TeV lower range
  h_m3_m4e_loose = fs->make<TH1D>("mass4e_loose_8TeV_low", "mass of 4 electrons", 37, 70., 181.);
  h_m3_m4e_loose->GetXaxis()->SetTitle("Invariant mass of 4e (in GeV/c^2)");
  h_m3_m4e_loose->GetYaxis()->SetTitle("Number of Events");

  // 4electron mass spectrum with loose cut in Gsf Electron full range
  h_m4_m4e_loose = fs->make<TH1D>("mass4e_loose_full", "mass of 4 electrons", 300, 0., 900.);
  h_m4_m4e_loose->GetXaxis()->SetTitle("Invariant mass of 4e (in GeV/c^2)");
  h_m4_m4e_loose->GetYaxis()->SetTitle("Number of Events");

 
  // These histograms are for ZZ/ZZ*To2mu2e reconstruction with different combination and with loose cuts
  // Mass of Z to 2mu from 2mu2e
  h_mZmu_2mu2e_loose = fs->make<TH1D>("massZmu_from2mu2e_loose", "mass of Z2mu from 2mu2e", 75, 0., 150.);
  h_mZmu_2mu2e_loose->GetXaxis()->SetTitle("Invariant mass of Z1 (in GeV/c^2)");
  h_mZmu_2mu2e_loose->GetYaxis()->SetTitle("Number of Events");

  // Mass of Z to 2e from 2mu2e
  h_mZe_2mu2e_loose = fs->make<TH1D>("massZe_from2mu2e_loose", "mass of Z2e from 2mu2e", 75, 0., 150.);
  h_mZe_2mu2e_loose->GetXaxis()->SetTitle("Invariant Mass of Z2 (in GeV/c^2)");
  h_mZe_2mu2e_loose->GetYaxis()->SetTitle("Number of Events");

  // Mass Za: mass of Z1 with higher mass
  h_mZa_2mu2e_loose = fs->make<TH1D>("mZa_2mu2e_loose", "mass Z higher", 120, 0., 120.);
  h_mZa_2mu2e_loose->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZa_2mu2e_loose->GetYaxis()->SetTitle("Number of Events");

  // Mass Zb: mass of Z2 with lower mass
  h_mZb_2mu2e_loose = fs->make<TH1D>("mZb_2mu2e_loose", "mass Z lower", 120, 0., 120.);
  h_mZb_2mu2e_loose->GetXaxis()->SetTitle("Invariant mass of dimuon (in GeV/c^2)");
  h_mZb_2mu2e_loose->GetYaxis()->SetTitle("Number of Events");

  // 2muon 2electron mass spectrum with loose cut paper 7TeV
  h_m1_m2mu2e_loose = fs->make<TH1D>("mass2mu2e_loose_7TeV", "mass of 2 muons and 2 electrons", 51, 98., 608.);
  h_m1_m2mu2e_loose->GetXaxis()->SetTitle("Invariant mass of 2mu2e (in GeV/c^2)");
  h_m1_m2mu2e_loose->GetYaxis()->SetTitle("Number of Events");

  // 2muon 2electron mass spectrum with loose cut paper 8TeV
  h_m2_m2mu2e_loose = fs->make<TH1D>("mass2mu2e_loose_8TeV", "mass of 2 muons and 2 electrons", 74, 70., 810.);
  h_m2_m2mu2e_loose->GetXaxis()->SetTitle("Invariant mass of 2mu2e (in GeV/c^2)");
  h_m2_m2mu2e_loose->GetYaxis()->SetTitle("Number of Events");

  // 2muon 2electron mass spectrum with loose cut paper 8TeV lower range
  h_m3_m2mu2e_loose = fs->make<TH1D>("mass2mu2e_loose_8TeV_low", "mass of 2 muons and 2 electrons", 37, 70., 181.);
  h_m3_m2mu2e_loose->GetXaxis()->SetTitle("Invariant mass of 2mu2e (in GeV/c^2)");
  h_m3_m2mu2e_loose->GetYaxis()->SetTitle("Number of Events");

  // 2muons 2electrons mass spectrum with loose cut full range
  h_m4_m2mu2e_loose = fs->make<TH1D>("mass2mu2e_loose_full", "mass of 2 muons and 2 electrons", 300, 0., 900.);
  h_m4_m2mu2e_loose->GetXaxis()->SetTitle("Invariant mass of 2mu2e (in GeV/c^2)");
  h_m4_m2mu2e_loose->GetYaxis()->SetTitle("Number of Events");

   /////////////////////////// CONTROL PLOTS /////////////////////////////////////
  // Below are the histograms for the control plots

  //------------- Global Muon -------------------//
  // Momentum of Global Muon
  h_p_gmu = fs->make<TH1D>("GM_momentum", "GM Momentum", 200, 0., 200.);
  h_p_gmu->GetXaxis()->SetTitle("Momentum (GeV/c)");
  h_p_gmu->GetYaxis()->SetTitle("Number of Events");

  // Transverse momentum of Global Muon
  h_pt_gmu_b4 = fs->make<TH1D>("b4_GM_pT", "GM Transverse Momentum", 200, 0., 200.);
  h_pt_gmu_b4->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_gmu_b4->GetYaxis()->SetTitle("Number of Events");

  // Pseudorapidity of Global Muon
  h_eta_gmu_b4 = fs->make<TH1D>("b4_GM_eta", "GM Eta", 140, -3.5, 3.5);
  h_eta_gmu_b4->GetXaxis()->SetTitle("eta");
  h_eta_gmu_b4->GetYaxis()->SetTitle("Number of Events");

  // Phi of Global Muon
  h_phi_gmu = fs->make<TH1D>("GM_phi", "GM Phi", 314, -3.15, 3.15);
  h_phi_gmu->GetXaxis()->SetTitle("Phi");
  h_phi_gmu->GetYaxis()->SetTitle("Number of Events");

  // Chi square of Global Muon
  h_chi2_gmu = fs->make<TH1D>("GM_chi2", "GM Chi2", 300, 0., 150.);
  h_chi2_gmu->GetXaxis()->SetTitle("Chi2");
  h_chi2_gmu->GetYaxis()->SetTitle("Number of Events");

  // Number of degree of freedom of Global Muon
  h_ndof_gmu = fs->make<TH1D>("GM_ndof", "GM Number degree of freedom", 100, 0., 100.);
  h_ndof_gmu->GetXaxis()->SetTitle("Ndof");
  h_ndof_gmu->GetYaxis()->SetTitle("Number of Events");

  // Normalized chi square of Global Muon
  h_normchi2_gmu = fs->make<TH1D>("GM_normalizedchi2", "GM NormalizedChi2", 200, 0., 20.);
  h_normchi2_gmu->GetXaxis()->SetTitle("NormalizedChi2");
  h_normchi2_gmu->GetYaxis()->SetTitle("Number of Events");

  // Validhits of Global Muon
  h_validhits_gmu = fs->make<TH1D>("GM_validhits", "GM ValidHits", 100, 0., 100.);
  h_validhits_gmu->GetXaxis()->SetTitle("Number of valid hits");
  h_validhits_gmu->GetYaxis()->SetTitle("Number of Events");

  // Pixelhits of Global Muon
  h_pixelhits_gmu = fs->make<TH1D>("GM_pixelhits", "GM Pixelhits", 14, 0., 14.);
  h_pixelhits_gmu->GetXaxis()->SetTitle("Munber of pixel hits");
  h_pixelhits_gmu->GetYaxis()->SetTitle("Number of Events");

  // pT of Global Muon
  h_pt_gmu_after = fs->make<TH1D>("after_GM_pT", "GM Transverse Momentum", 200, 0., 200.);
  h_pt_gmu_after->GetXaxis()->SetTitle("pT of global muons (GeV/c)");
  h_pt_gmu_after->GetYaxis()->SetTitle("Number of Events");

  // eta of Global Muon
  h_eta_gmu_after = fs->make<TH1D>("after_GM_eta", "GM Eta", 140, -3.5, 3.5);
  h_eta_gmu_after->GetXaxis()->SetTitle("eta");
  h_eta_gmu_after->GetYaxis()->SetTitle("Number of Events");

  //--------------------------------------------------//

  //--------------- Tracker Muon ---------------------//
  // Momentum of Tracker Muon 
  h_p_tmu = fs->make<TH1D>("TM_momentum", "TM Momentum", 200, 0., 200.);
  h_p_tmu->GetXaxis()->SetTitle("Momentum (GeV/c)");
  h_p_tmu->GetYaxis()->SetTitle("Number of Events");

  // pT of Tracker Muon
  h_pt_tmu_b4 = fs->make<TH1D>("b4_TM_pt", "TM pT", 200, 0., 200.);
  h_pt_tmu_b4->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_tmu_b4->GetYaxis()->SetTitle("Number of Events");

  // Eta of Tracker Muon
  h_eta_tmu_b4 = fs->make<TH1D>("b4_TM_eta", "TM Eta", 140, -3.5, 3.5);
  h_eta_tmu_b4->GetXaxis()->SetTitle("eta");
  h_eta_tmu_b4->GetYaxis()->SetTitle("Number of Events");

  // Phi of Tracker Muon
  h_phi_tmu = fs->make<TH1D>("TM_phi", "TM Phi", 314, -3.15, 3.15);
  h_phi_tmu->GetXaxis()->SetTitle("Phi");
  h_phi_tmu->GetYaxis()->SetTitle("Number of Events");

  // Chi of Tracker Muon
  h_chi2_tmu = fs->make<TH1D>("TM_chi2", "TM Chi2", 500, 0., 100.);
  h_chi2_tmu->GetXaxis()->SetTitle("Chi2");
  h_chi2_tmu->GetYaxis()->SetTitle("Number of Events");

  // No. degree of freedom of Tracker Muon
  h_ndof_tmu = fs->make<TH1D>("TM_Ndof", "TM Number degree of freedom", 60, 0., 60.);
  h_ndof_tmu->GetXaxis()->SetTitle("Ndof");
  h_ndof_tmu->GetYaxis()->SetTitle("Number of Events");

  // Normalizedchi2 of Tracker Muon
  h_normchi2_tmu = fs->make<TH1D>("TM_NormalizedChi2", "TM NormalizedChi2", 200, 0., 20.);
  h_normchi2_tmu->GetXaxis()->SetTitle("NormalizedChi2");
  h_normchi2_tmu->GetYaxis()->SetTitle("Number of Events");

  // Validhits of Tracker Muon
  h_validhits_tmu_jpsi = fs->make<TH1D>("TM_validhits", "TM ValidHits", 100, 0., 100.); // Jpsi
  h_validhits_tmu_jpsi->GetXaxis()->SetTitle("Number of valid hits");
  h_validhits_tmu_jpsi->GetYaxis()->SetTitle("Number of Events");

  // Pixelhits of Tracker Muon
  h_pixelhits_tmu_jpsi = fs->make<TH1D>("TM_pixelhits", "TM Pixelhits", 14, 0., 14.);
  h_pixelhits_tmu_jpsi->GetXaxis()->SetTitle("Munber of pixel hits");
  h_pixelhits_tmu_jpsi->GetYaxis()->SetTitle("Number of Events");

  // pT of Tracker Muon after cuts in Jpsi
  h_pt_tmu_jpsi_after = fs->make<TH1D>("after_TM_jpsi_pt", "TM pT after cuts in Jpsi", 200, 0., 200.);
  h_pt_tmu_jpsi_after->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_tmu_jpsi_after->GetYaxis()->SetTitle("Number of Events");

  // Eta of Tracker Muon after cuts in Jpsi
  h_eta_tmu_jpsi_after = fs->make<TH1D>("after_TM_jpsi_eta", "TM Eta after cuts in Jpsi", 140, -3.5, 3.5);
  h_eta_tmu_jpsi_after->GetXaxis()->SetTitle("eta");
  h_eta_tmu_jpsi_after->GetYaxis()->SetTitle("Number of Events");

  //--------------------------------------------------------//

  //--------------- Global Tracker Muon (GTM) --------------------// 
  // No. of valid muon hits in GTM
  h_goodhit = fs->make<TH1D>("GTM_goodMuonChamberHit", "GTM No. of goodMuonChamberHit", 40, 0., 40.);
  h_goodhit->GetXaxis()->SetTitle("Number of muon valid hits");
  h_goodhit->GetYaxis()->SetTitle("Number of Events");

  // Transverse impact parameter with respect to primary vertex in GTM
  h_dxy_mu = fs->make<TH1D>("GTM_dxy", "GTM dxy", 100, 0., 1.);
  h_dxy_mu->GetXaxis()->SetTitle("Transverse impact parameter w.r.t. BS");
  h_dxy_mu->GetYaxis()->SetTitle("Number of Events");

  // Muons matching with at least two Muon Stations in GTM
  h_mustation = fs->make<TH1D>("GTM_MuonStations", "GTM MuonStations", 2, 0., 2.);
  h_mustation->GetXaxis()->SetTitle("Muons matching with at least two muon stations");
  h_mustation->GetYaxis()->SetTitle("Number of Events");

  // Relative isolation in GTM
  h_relPFIso_mu = fs->make<TH1D>("GTM_RelPFIso", "GTM Relative Isolation", 100, 0., 5.);
  h_relPFIso_mu->GetXaxis()->SetTitle("Relative Isolation");
  h_relPFIso_mu->GetYaxis()->SetTitle("Number of Events");

  // Track isolation in GTM
  h_Isotrk_mu = fs->make<TH1D>("GTM_Isotrk", "GTM Track Isolation", 100, 0., 5.);
  h_Isotrk_mu->GetXaxis()->SetTitle("Track Isolation");
  h_Isotrk_mu->GetYaxis()->SetTitle("Number of Events");

  // Relative isolation in GTM after the loose and no.of cand. cuts
  h_relPFIso_mu_loose_after = fs->make<TH1D>("GTM_RelPFIso_loose_after", "GTM RelPFIso loose cuts", 100, 0., 5.);
  h_relPFIso_mu_loose_after->GetXaxis()->SetTitle("Relative Isolation");
  h_relPFIso_mu_loose_after->GetYaxis()->SetTitle("Number of Events");

  // Validhits in GTM
  h_validhits_mu = fs->make<TH1D>("GTM_validhits", "GTM ValidHits", 100, 0., 100.);
  h_validhits_mu->GetXaxis()->SetTitle("Number of valid hits");
  h_validhits_mu->GetYaxis()->SetTitle("Number of Events");

  // Pixelhits in GTM
  h_pixelhits_mu = fs->make<TH1D>("GTM_pixelhits", "GTM Pixelhits", 14, 0., 14.);
  h_pixelhits_mu->GetXaxis()->SetTitle("Munber of pixel hits");
  h_pixelhits_mu->GetYaxis()->SetTitle("Number of Events");


  // pT muon after loose cuts for Z to 2mu
  h_pt_tmu_loose_after_Zto2mu = fs->make<TH1D>("after_GTM_loose_pt_Zto2mu", "Muon pT with loose cuts Zto2mu", 200, 0., 200.);
  h_pt_tmu_loose_after_Zto2mu->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_tmu_loose_after_Zto2mu->GetYaxis()->SetTitle("Number of Events");

  // Eta muon after loose cuts for Z to 2mu
  h_eta_tmu_loose_after_Zto2mu = fs->make<TH1D>("after_GTM_loose_eta_Zto2mu", "Muon eta with loose cuts Zto2mu", 140, -3.5, 3.5);
  h_eta_tmu_loose_after_Zto2mu->GetXaxis()->SetTitle("Eta");
  h_eta_tmu_loose_after_Zto2mu->GetYaxis()->SetTitle("Number of Events");


  // pT muon after loose cuts
  h_pt_tmu_loose_after = fs->make<TH1D>("after_GTM_loose_pt", "Muon pT with loose cuts", 200, 0., 200.); // Z, 4mu, 2mu 2e
  h_pt_tmu_loose_after->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_tmu_loose_after->GetYaxis()->SetTitle("Number of Events");

  // Eta muon after loose cuts
  h_eta_tmu_loose_after = fs->make<TH1D>("after_GTM_loose_eta", "Muon eta with loose cuts", 140, -3.5, 3.5);
  h_eta_tmu_loose_after->GetXaxis()->SetTitle("eta");
  h_eta_tmu_loose_after->GetYaxis()->SetTitle("Number of Events");

  // 2mu2e
  // Relative isolation of 2mu for 2mu2e after the loose and no.of cand. cuts
  h_relPFIso_2mu_loose_after = fs->make<TH1D>("after_relPFIso_loose_2mu", "RelPFIso loose cuts", 50, 0., 5.);
  h_relPFIso_2mu_loose_after->GetXaxis()->SetTitle("Relative Isolation");
  h_relPFIso_2mu_loose_after->GetYaxis()->SetTitle("Number of Events");

  // Relative isolation of 2e for 2mu2e after the loose and no.of cand. cuts
  h_relPFIso_2e_loose_after = fs->make<TH1D>("after_relPFIso_loose_2e", "RelPFIso loose cuts", 50, 0., 5.);
  h_relPFIso_2e_loose_after->GetXaxis()->SetTitle("Relative Isolation");
  h_relPFIso_2e_loose_after->GetYaxis()->SetTitle("Number of Events");

  // pT muon after loose cuts for 2mu2e
  h_pt_tmu_loose_after_2mu2e = fs->make<TH1D>("after_GTM_loose_pt_2mu2e", "Muon pT with loose cuts", 200, 0., 200.); 
  h_pt_tmu_loose_after_2mu2e->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_tmu_loose_after_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // Eta muon after loose cuts for 2mu2e
  h_eta_tmu_loose_after_2mu2e = fs->make<TH1D>("after_GTM_loose_eta_2mu2e", "Muon eta with loose cuts", 140, -3.5, 3.5);
  h_eta_tmu_loose_after_2mu2e->GetXaxis()->SetTitle("eta");
  h_eta_tmu_loose_after_2mu2e->GetYaxis()->SetTitle("Number of Events");


  //------------------- Electron Control Plot ------------------------//
  // Electron momentum
  h_p_e = fs->make<TH1D>("Electron_momentum", "Electron momentum", 200, 0., 200.);
  h_p_e->GetXaxis()->SetTitle("Momentum (GeV/c)");
  h_p_e->GetYaxis()->SetTitle("Number of Events");

  // Electron eT
  h_et_e = fs->make<TH1D>("Electron_eT", "Electron eT", 200, 0., 200.);
  h_et_e->GetXaxis()->SetTitle("eT (GeV/c)");
  h_et_e->GetYaxis()->SetTitle("Number of Events");

  // Electron pT before cuts
  h_pt_e_b4 = fs->make<TH1D>("b4_Electron_pT", "Electron pT", 200, 0, 200.);
  h_pt_e_b4->GetXaxis()->SetTitle("pT (GeV/c)");
  h_pt_e_b4->GetYaxis()->SetTitle("Number of Events");

  // Electron eta before cuts
  h_eta_e_b4 = fs->make<TH1D>("b4_Electron_eta", "Electron_eta", 140, -3.5, 3.5);
  h_eta_e_b4->GetXaxis()->SetTitle("eta");
  h_eta_e_b4->GetYaxis()->SetTitle("Number of Events");

  // Electron phi
  h_phi_e = fs->make<TH1D>("Electron_phi", "Electron_phi", 314, -3.17, 3.17);
  h_phi_e->GetXaxis()->SetTitle("Phi");
  h_phi_e->GetYaxis()->SetTitle("Number of Events");

  // Electron supercluster eta
  h_suclust_eta = fs->make<TH1D>("Electron_Super_Cluster_eta", "Super Cluster eta", 140, -3.5, 3.5);
  h_suclust_eta->GetXaxis()->SetTitle("Super Cluster Eta");
  h_suclust_eta->GetYaxis()->SetTitle("Number of Events");

  // Electron super cluster rawenergy
  h_suclust_rawE = fs->make<TH1D>("Electron_Super_Cluster_rawenergy", "Super Cluster rawenergy", 200, 0., 200.);
  h_suclust_rawE->GetXaxis()->SetTitle("Super Cluster Energy");
  h_suclust_rawE->GetYaxis()->SetTitle("Number of Events");

  // Relative isolation electron
  h_relPFIso_e = fs->make<TH1D>("Electron_RelPFIso", "Relative Isolation in Electron", 100, 0., 5.);
  h_relPFIso_e->GetXaxis()->SetTitle("Relative Isolation");
  h_relPFIso_e->GetYaxis()->SetTitle("Number of Events");

  // Track isolation electron
  h_Isotrk_e = fs->make<TH1D>("Electron_Isotrk", "Track Isolation in Electron", 100, 0., 5.);
  h_Isotrk_e->GetXaxis()->SetTitle("Track Isolation");
  h_Isotrk_e->GetYaxis()->SetTitle("Number of Events");

  // Relative isolation electron after
  h_relPFIso_e_loose_after = fs->make<TH1D>("Electron_RelPFIso_loose_after", "Relative Isolation in Electron", 100, 0., 5.);
  h_relPFIso_e_loose_after->GetXaxis()->SetTitle("Relative Isolation");
  h_relPFIso_e_loose_after->GetYaxis()->SetTitle("Number of Events");

  // Just for checking the relation of pT and RelPFIso
  double Isocheck[12] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.5, 10.};
  h_relPFIso_pt_e = fs->make<TH2D>("Electron_RelPFIso_pT", "Relative Isolation in Electron 2D", 11, Isocheck, 100, 0., 50.);
  h_relPFIso_pt_e->GetXaxis()->SetTitle("Relative Isolation");
  h_relPFIso_pt_e->GetYaxis()->SetTitle("Number of Events");

  // Transverse impact parameter with respect to primary vertex in Electron
  h_dxy_e = fs->make<TH1D>("Electron_dxy", "Electron dxy", 100, 0., 1.);
  h_dxy_e->GetXaxis()->SetTitle("Transverse impact parameter w.r.t. primvtx");
  h_dxy_e->GetYaxis()->SetTitle("Number of Events");

  // Electron pT after loose cuts for Zto2e
  h_pt_e_loose_after_Zto2e = fs->make<TH1D>("after_Electron_loose_pT_Zto2e", "Electron pT", 240, 0., 120.);
  h_pt_e_loose_after_Zto2e->GetXaxis()->SetTitle("Electron pT (GeV/c)");
  h_pt_e_loose_after_Zto2e->GetYaxis()->SetTitle("Number of Events");

  // Electron eta after loose cuts for Zto2e
  h_eta_e_loose_after_Zto2e = fs->make<TH1D>("after_Electron_loose_eta_Zto2e", "Electron eta", 140, -3.5, 3.5);
  h_eta_e_loose_after_Zto2e->GetXaxis()->SetTitle("Eta");
  h_eta_e_loose_after_Zto2e->GetYaxis()->SetTitle("Number of Events");


  // Electron pT after loose cuts
  h_pt_e_loose_after = fs->make<TH1D>("after_Electron_loose_pT", "Electron pT", 240, 0., 120.);
  h_pt_e_loose_after->GetXaxis()->SetTitle("Electron pT (GeV/c)");
  h_pt_e_loose_after->GetYaxis()->SetTitle("Number of Events");

  // Electron eta after loose cuts
  h_eta_e_loose_after = fs->make<TH1D>("after_Electron_loose_eta", "Electron eta", 140, -3.5, 3.5);
  h_eta_e_loose_after->GetXaxis()->SetTitle("eta");
  h_eta_e_loose_after->GetYaxis()->SetTitle("Number of Events");

  // 2mu2e
  // Electron pT after loose cuts for 2mu2e
  h_pt_e_loose_after_2mu2e = fs->make<TH1D>("after_Electron_loose_pT_2mu2e", "Electron pT", 240, 0., 120.);
  h_pt_e_loose_after_2mu2e->GetXaxis()->SetTitle("Electron pT (GeV/c)");
  h_pt_e_loose_after_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // Electron eta after loose cuts for 2mu2e
  h_eta_e_loose_after_2mu2e = fs->make<TH1D>("after_Electron_loose_eta_2mu2e", "Electron eta", 140, -3.5, 3.5);
  h_eta_e_loose_after_2mu2e->GetXaxis()->SetTitle("eta");
  h_eta_e_loose_after_2mu2e->GetYaxis()->SetTitle("Number of Events");

  // SIP for muon
  h_SIP3d_mu_b4 = fs->make<TH1D>("SIP3d_mu", "SIP_3D for Muon", 100, 0., 10.);
  h_SIP3d_mu_b4->GetXaxis()->SetTitle("SIP_3D");
  h_SIP3d_mu_b4->GetYaxis()->SetTitle("Number of Events");

  // SIP for electron
  h_SIP3d_e_b4 = fs->make<TH1D>("SIP3d_e", "SIP_3D for Electron", 100, 0., 10.);
  h_SIP3d_e_b4->GetXaxis()->SetTitle("SIP_3D");
  h_SIP3d_e_b4->GetYaxis()->SetTitle("Number of Events");

  h_misshite = fs->make<TH1D>("Electron_misshit" , "Electron Track missing hits", 5, 0., 5.);
  h_misshite->GetXaxis()->SetTitle("gsfTrack Hit type");
  h_misshite->GetYaxis()->SetTitle("Number of Events");
  //  t = fs->make<TTree>("tree", "tree");

}

DemoAnalyzer::~DemoAnalyzer() {
  //do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------//
void DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

// **********************************************
// here each relevant event will get analyzed 
// **********************************************

  // using namespace edm;
  // using namespace reco;
  // using namespace std;

  // nRun1  = iEvent.eventAuxiliary().run();
  // nEvt1  = iEvent.eventAuxiliary().event();
  // nLumi1 = iEvent.eventAuxiliary().luminosityBlock();

  nRun  = iEvent.run();
  nEvt  = (iEvent.id()).event(); // iEvent xd class named event()
  nLumi = iEvent.luminosityBlock();
  
  
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif

  // Event is to be analyzed

  edm::LogInfo("Demo")
  << "Starting to analyze \n"
  << "Event number: " << (iEvent.id()).event()
  << ", Run number: " << iEvent.run()
  << ", Lumisection: " << iEvent.luminosityBlock();

  //------------------Load (relevant) Event information------------------------//
  // INFO: Getting Data From an Event
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookChapter4#GetData
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMGetDataFromEvent#get_ByLabel
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAodDataTable
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable

  // INFO: globalMuons
  // NB: note that when using keyword "globalMuons" getByLabel-function returns 
  // reco::TrackCollection
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel("generalTracks", tracks);

  edm::Handle<reco::TrackCollection> gmuons;
  iEvent.getByLabel("globalMuons", gmuons);

  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel("muons", muons);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);

  edm::Handle<reco::VertexCollection> primvtxHandle;
  iEvent.getByLabel("offlinePrimaryVertices", primvtxHandle);

  edm::Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByLabel("gsfElectrons",electrons);

  reco::BeamSpot beamSpot;
  if ( beamSpotHandle.isValid() )
    {
      beamSpot = *beamSpotHandle;

    } else
    {
      edm::LogInfo("Demo")
	<< "No beam spot available from EventSetup \n";
    }

  reco::VertexCollection primvtx;
  if ( primvtxHandle.isValid() )
    {
      primvtx = *primvtxHandle;

    } else
    {
      edm::LogInfo("Demo")
	<< "No primary vertex available from EventSetup \n";
    }


  // Declare vector that contain a pair of variables u wanna save
  // In this case: size and pT
  std::vector< std::pair<int, double> > vIdPt;
  std::vector< std::pair<int, double> > vIdPtjpsi;
  std::vector< std::pair<int, double> > vIdPtmuloose;
  std::vector< std::pair<int, double> > vIdPteloose;

  /*
  std:: vector<int> listgmuons;
  std::vector<int> listmu4jpsi;
  std::vector<int> listmuloose; // use for Z reco and 2 mu
  std::vector<int> listmutight; // use for 4 mu
  std::vector<int> listeloose; // use for Z and 4e
  std::vector<int> listetight; // use for Z and 4e
  */
  // Initialize variables
  eZ12 = -9999.; eZ34 = -9999.; eZ13 = -9999.; eZ24 = -9999.; eZ14 = -9999.; eZ23 = -9999.; // select largest, init -
  pxZ12 = -9999.; pxZ34 = -9999.; pxZ13 = -9999.; pxZ24 = -9999.; pxZ14 = -9999.; pxZ23 = -9999.;
  pyZ12 = -9999.; pyZ34 = -9999.; pyZ13 = -9999.; pyZ24 = -9999.; pyZ14 = -9999.; pyZ23 = -9999.;
  pzZ12 = -9999.; pzZ34 = -9999.; pzZ13 = -9999.; pzZ24 = -9999.; pzZ14 = -9999.; pzZ23 = -9999.; 
  pZ12 = -9999.; pZ34 = -9999.; pZ13 = -9999.; pZ24 = -9999.; pZ14 = -9999.; pZ23 = -9999.; 
  pTZ12 = -9999.; pTZ34 = -9999.; pTZ13 = -9999.; pTZ24 = -9999.; pTZ14 = -9999.; pTZ23 = -9999.; 

  mZ12 = -9999.; mZ34 = -9999.; mZ13 = -9999.; mZ24 = -9999.; mZ14 = -9999.; mZ23 = -9999.;

  dZ12 = 9999.; dZ34 = 9999.; dZ13 = 9999.; dZ24 = 9999.; dZ14 = 9999.; dZ23 = 9999.; // select smallest, init +
  dZc1 = 9999.; dZc2 = 9999.; dZc3 = 9999.;

  eZa = -9999.; pxZa = -9999.; pyZa = -9999.; pzZa = -9999.; pTZa = -9999.; mZa = -9999.;
  eZb = -9999.; pxZb = -9999.; pyZb = -9999.; pzZb = -9999.; pTZb = -9999.; mZb = -9999.;
  mass4mu = -9999.; pt_4mu = -9999.; eta_4mu = -9999.;
  pt_mu1 = -9999.; pt_mu2 = -9999.; pt_mu3 = -9999.; pt_mu4 = -9999.;
  eta_mu1 = -9999.; eta_mu2 = -9999.; eta_mu3 = -9999.; eta_mu4 = -9999.;

  mass4e = -9999.; pt_4e = -9999.; eta_4e = -9999.;
  pt_e1 = -9999.; pt_e2 = -9999.; pt_e3 = -9999.; pt_e4 = -9999.;
  eta_e1 = -9999.; eta_e2 = -9999.; eta_e3 = -9999.; eta_e4 = -9999.;
  
  s = -9999.;
  s1 = -9999.; s2 = -9999.; s3 = -9999.; s4 = -9999.;
  dx = 9999.; dy = 9999.; dz = 9999.; rap = -9999.; pz = -9999.;

  mass2mu2e = -9999.; pt_2mu2e = -9999.; eta_2mu2e = -9999.;
  pt_2mu1 = -9999.; pt_2mu2 = -9999.; pt_2e1 = -9999.; pt_2e2 = -9999.;
  eta_2mu1 = -9999.; eta_2mu2 = -9999.; eta_2e1 = -9999.; eta_2e2 = -9999.;
  
  goodhit = -9999.;
  relPFIso_mu = -9999.;
  relPFIso_e = -9999.;

  Isotrk_mu = -9999.;
  Isotrk_e = -9999.;

  IP3d_mu = -9999;
  ErrIP3d_mu = -9999;
  SIP3d_mu = -9999.;
  IP3d_e = -9999;
  ErrIP3d_e = -9999;
  SIP3d_e = -9999.;

  p4Za.SetPxPyPzE(0., 0., 0., 0.);
  p4Zb.SetPxPyPzE(0., 0., 0., 0.);
  p4H.SetPxPyPzE(0., 0., 0., 0.);

  // constant value squared muon mass, squared electron mass and Z
  sqm1 = (0.105658) * (0.105658);
  sqme = (0.0005109989) * (0.0005109989);
  mZ = 91.1876;

  h_globalmu_size->Fill(gmuons->size());
  h_tmu_size->Fill(muons->size());
  h_e_size->Fill(electrons->size());

  ///////////////////////////////////////////////////////////////////////////////
  /////////////////////// Global Muon Collection Start //////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  // Loop over global muons size and select good muons

  for (unsigned t = 0; t < gmuons->size(); t++)
    {
      const reco::Track &iMuon = (*gmuons)[t];
      const reco::HitPattern& p = iMuon.hitPattern();

      h_p_gmu->Fill(iMuon.p());
      h_pt_gmu_b4->Fill(iMuon.pt());
      h_eta_gmu_b4->Fill(iMuon.eta());
      h_phi_gmu->Fill(iMuon.phi());

      h_chi2_gmu->Fill(iMuon.chi2());
      h_ndof_gmu->Fill(iMuon.ndof());
      h_normchi2_gmu->Fill(iMuon.normalizedChi2());

      // Counter
      int GM_ValidHits = 0;
      int GM_PixelHits = 0;
   
      // loop over the hits of the track
      for (int i = 0; i < p.numberOfHits(); i++) 
	{
	  uint32_t hit = p.getHitPattern(i);

	  // if the hit is valid and in pixel
	  if (p.validHitFilter(hit) && p.pixelHitFilter(hit)) {GM_PixelHits++;}
	  if (p.validHitFilter(hit)) {GM_ValidHits++;}
	}

      h_validhits_gmu->Fill(GM_ValidHits);
      h_pixelhits_gmu->Fill(GM_PixelHits);

      if (GM_ValidHits >= 12 && GM_PixelHits >= 2 && iMuon.normalizedChi2() < 4.0)
	{
	  // Save a vector that contain 2 variables: gmuon size (.first) and the pT (.second)
	  vIdPt.push_back( std::make_pair(t, iMuon.pt()) );
	}
    }

  // Sort the pT (.second) to decending order (from highest pT to lowest pT)
  std::sort(vIdPt.begin(), vIdPt.end(), [](const std::pair<int, double> &idPt1, const std::pair<int, double> &idPt2) {return (idPt1.second > idPt2.second);});
  //  std::cout << vIdPt.size() << " ";

  // for (unsigned i = 0; i < vIdPt.size(); i++)
  // std::cout << vIdPt.at(i).second << " ";
  // std::cout << std::endl;

  ///////////////////////////////////////////////////////////////////////////////
  //////////////////////// Global Muon Collection End ///////////////////////////
  ///////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////
  ///////////////////////// Reco Muon Collection Start ///////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  // Loop over muons size and select good muons

  for (unsigned u = 0; u < muons->size(); u++)
    {
      const reco::Muon &itMuon = (*muons)[u];

      // math::XYZPoint point(beamSpot.position());
      math::XYZPoint point(primvtx[0].position());

      if (itMuon.isPFMuon() && itMuon.isPFIsolationValid() && (itMuon.globalTrack()).isNonnull())                               // 2
      //if (itMuon.isPFMuon() && itMuon.isPFIsolationValid() && itMuon.isGlobalMuon() && (itMuon.globalTrack()).isNonnull())    // 3
      //if (itMuon.isPFMuon() && itMuon.isPFIsolationValid() && (itMuon.track()).isNonnull())                                   // 4
      //if (itMuon.isPFMuon() && itMuon.isPFIsolationValid() && itMuon.isTrackerMuon() && (itMuon.track()).isNonnull())         // 5
      // if (itMuon.isPFMuon() && itMuon.isPFIsolationValid() && (itMuon.isGlobalMuon() || itMuon.isTrackerMuon()) && (itMuon.track()).isNonnull() && (itMuon.globalTrack()).isNonnull())         // 6a (GTrk) and 6b (TTrk)

	{

	  h_p_tmu->Fill(itMuon.p());
	  h_pt_tmu_b4->Fill(itMuon.pt());
	  h_eta_tmu_b4->Fill(itMuon.eta());
	  h_phi_tmu->Fill(itMuon.phi());

	  // some muons might not have valid track references
	  // if ((itMuon.globalTrack()).isNonnull()) {

	  h_chi2_tmu->Fill((itMuon.globalTrack())->chi2());
	  h_ndof_tmu->Fill((itMuon.globalTrack())->ndof());
	  h_normchi2_tmu->Fill((itMuon.globalTrack())->normalizedChi2());

	  // h_chi2_tmu->Fill((itMuon.globalTrack())->chi2());
	  // h_ndof_tmu->Fill((itMuon.globalTrack())->ndof());
	  // h_normchi2_tmu->Fill((itMuon.globalTrack())->normalizedChi2());

	  //============ Good muons for Jpsi reconstruction start ==============//
	  //========= Use tracker muon only and no relative isolation ==========//

	  const reco::HitPattern& p = (itMuon.globalTrack())->hitPattern();
	  // const reco::HitPattern& p = (itMuon.globalTrack())->hitPattern();

	  int TM_ValidHits = 0;
	  int TM_PixelHits = 0;

	  // loop over the hits of the track
	  for (int i=0; i<p.numberOfHits(); i++)
	    {
	      uint32_t hit = p.getHitPattern(i);

	      // if the hit is valid and in pixel
	      if (p.validHitFilter(hit) && p.pixelHitFilter(hit)) {TM_PixelHits++;}
	      if (p.validHitFilter(hit)) {TM_ValidHits++;}
	    }

	  h_validhits_tmu_jpsi->Fill(TM_ValidHits);
	  h_pixelhits_tmu_jpsi->Fill(TM_PixelHits);

	  if (TM_ValidHits >= 12 && TM_PixelHits >= 2 && (itMuon.globalTrack())->normalizedChi2() < 4.0)
	    // if (TM_ValidHits >= 12 && TM_PixelHits >= 2 && (itMuon.globalTrack())->normalizedChi2() < 4.0)

	    {
	      if (std::abs(itMuon.eta()) < 1.3 && itMuon.pt() > 3.3)
		vIdPtjpsi.push_back( std::make_pair(u, itMuon.pt()) );
	  
	      else if (std::abs(itMuon.eta()) > 1.3 && std::abs(itMuon.eta()) < 2.2 && itMuon.p() > 2.9)
		vIdPtjpsi.push_back( std::make_pair(u, itMuon.pt()) );
	  
	      else if (std::abs(itMuon.eta()) > 2.2  && std::abs(itMuon.eta()) < 2.4 && itMuon.pt() > 2.4)
		vIdPtjpsi.push_back( std::make_pair(u, itMuon.pt()) );

	      // h_pt_tmu_jpsi_after->Fill(itMuon.pt());
	      // h_eta_tmu_jpsi_after->Fill(itMuon.eta());
	    }
	
	  //============= Good muons for Jpsi reconstruction end ===============//


	  //==== Good muons for Z, 4 muons and 2 muon reconstruction start =====//
	  //======= Use Particle Flow (PF) Muon & PF relative isolation ========//

	  // if (itMuon.isGlobalMuon() && (itMuon.globalTrack()).isNonnull() &&
	  // itMuon.isPFMuon() && itMuon.isPFIsolationValid())

	  // atau buat mcm ni:
	  // if (itMuon.isPFMuon() && (itMuon.isGlobalMuon() || itMuon.isTrackerMuon()) && (itMuon.globalTrack()).isNonnull() && itMuon.isPFIsolationValid())
	  // if (itMuon.isPFMuon() && itMuon.isGlobalMuon() && (itMuon.globalTrack()).isNonnull() && itMuon.isPFIsolationValid())

	  // if (itMuon.isPFMuon() && itMuon.isPFIsolationValid())
	
	  // PF Relative isolation for muons
	  relPFIso_mu = ((itMuon.pfIsolationR04()).sumChargedHadronPt +
			 (itMuon.pfIsolationR04()).sumNeutralHadronEt + 
			 (itMuon.pfIsolationR04()).sumPhotonEt) / itMuon.pt(); 

	  h_relPFIso_mu->Fill(relPFIso_mu);
	  
	  // Checking hit pattern info
	  
	  const reco::HitPattern& GTM_p = (itMuon.globalTrack())->hitPattern();

	  goodhit = GTM_p.numberOfValidMuonHits();
	  h_goodhit->Fill(goodhit);

	  h_dxy_mu->Fill((itMuon.globalTrack())->dxy(point));

	  IP3d_mu = sqrt( (itMuon.globalTrack()->dxy(point) * itMuon.globalTrack()->dxy(point)) +
			   (itMuon.globalTrack()->dz(point) * itMuon.globalTrack()->dz(point)) );
	  
	  ErrIP3d_mu = sqrt( (itMuon.globalTrack()->d0Error() * itMuon.globalTrack()->d0Error()) +
			      (itMuon.globalTrack()->dzError() * itMuon.globalTrack()->dzError()) );
	  
	  SIP3d_mu = IP3d_mu / ErrIP3d_mu;

	  h_SIP3d_mu_b4->Fill(SIP3d_mu);

	  int GTM_ValidHits = 0;
	  int GTM_PixelHits = 0;

	  for (int i = 0; i < GTM_p.numberOfHits(); i++) {
	    uint32_t hit = GTM_p.getHitPattern(i);

	    // If the hit is valid and in pixel
	    if (GTM_p.validHitFilter(hit) && GTM_p.pixelHitFilter(hit))
	      {GTM_PixelHits++;}
	    
	      if (GTM_p.validHitFilter(hit)) {GTM_ValidHits++;}
	    }

	  h_validhits_mu->Fill(GTM_ValidHits);
	  h_pixelhits_mu->Fill(GTM_PixelHits);

			 // bool mustation = muon::isGoodMuon(itMuon, muon::TMLastStationTight);
			 // h_mustation->Fill(mustation);

	  // Good muons with loose pt, eta, dxy and relPFIso_mu cuts
	  // if (std::abs(SIP3d_mu) < 4. && std::abs((itMuon.globalTrack())->dxy(point)) < 0.5 && std::abs((itMuon.globalTrack())->dz(point)) < 1. &&
	  if (std::abs(SIP3d_mu) < 4. && std::abs((itMuon.globalTrack())->dxy(point)) < 0.5 && std::abs((itMuon.globalTrack())->dz(point)) < 1. &&

	      relPFIso_mu < 0.4)
	    //Isotrk_mu < 0.7 std::abs((itMuon.globalTrack())->dxy(point)) < 0.5
	    {
	      if (itMuon.pt() > 5. && std::abs(itMuon.eta()) < 2.4)
		{
		  // h_pt_tmu_loose_after->Fill(itMuon.pt());
		  // h_eta_tmu_loose_after->Fill(itMuon.eta());

		  vIdPtmuloose.push_back( std::make_pair(u, itMuon.pt()) );
		  // listmuloose.push_back(u);
		  // muons.at(u) is good, save the index

		}
	    }

	} // end of if(itMuon.isGlobalMuon()..........)

      //=================== Good muons for Z, 4 muons... end ===================//

	// } // end of if((itMuon.globalTrack()).isNonnull())
    } // end of loop over good candidates of muons for all cases

  std::sort(vIdPtjpsi.begin(), vIdPtjpsi.end(), [](const std::pair<int, double> &idPtjpsi1, const std::pair<int, double> &idPtjpsi2) {return (idPtjpsi1.second > idPtjpsi2.second);});
  std::sort(vIdPtmuloose.begin(), vIdPtmuloose.end(), [] (const std::pair<int, double> &idPtmuloose1, const std::pair<int, double> &idPtmuloose2) {return (idPtmuloose1.second > idPtmuloose2.second);});

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////// Muon Collection end /////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////
  ////////////////////// Electron Collection Start //////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  for (unsigned te = 0; te < electrons->size(); te++)
    {
      const reco::GsfElectron &iElectron = (*electrons)[te];

      // math::XYZPoint point(beamSpot.position());
       math::XYZPoint point(primvtx[0].position());

       if (iElectron.passingPflowPreselection()) {
      misshits = ((iElectron.gsfTrack())->trackerExpectedHitsInner()).numberOfHits(); 

      IP3d_e = sqrt ( (iElectron.gsfTrack()->dxy(point) * iElectron.gsfTrack()->dxy(point)) + (iElectron.gsfTrack()->dz(point) * iElectron.gsfTrack()->dz(point)) );
      ErrIP3d_e = sqrt ( (iElectron.gsfTrack()->d0Error() * iElectron.gsfTrack()->d0Error()) + (iElectron.gsfTrack()->dzError() * iElectron.gsfTrack()->dzError()) );
      SIP3d_e = IP3d_e / ErrIP3d_e;

      h_SIP3d_e_b4->Fill(SIP3d_e);

      h_p_e->Fill(iElectron.p());
      h_et_e->Fill(iElectron.et());
      h_pt_e_b4->Fill(iElectron.pt());
      h_eta_e_b4->Fill(iElectron.eta());
      h_phi_e->Fill(iElectron.phi());
      // h_suclust_eta->Fill(std::abs((iElectron.superCluster())->eta()));
      h_suclust_eta->Fill((iElectron.superCluster())->eta());
      h_suclust_rawE->Fill(std::abs((iElectron.superCluster())->rawEnergy()));
      h_misshite->Fill(misshits);
      
      //Relative isolation for electron
      // ni kene tukr name nnti & kita ximplement PF iso for electron
      // relPFIso_e = (iElectron.dr03TkSumPt() + iElectron.dr03EcalRecHitSumEt() + iElectron.dr03HcalTowerSumEt()) / iElectron.pt();

       Isotrk_e = (iElectron.dr03TkSumPt()) / iElectron.pt();

       //       relPFIso_e = ((iElectron.pfIsolationR04()).sumChargedHadronPt +
       //    (iElectron.pfIsolationR04()).sumNeutralHadronEt +
       //		    (iElectron.pfIsolationR04()).sumPhotonEt) / iElectron.pt(); // relPFiso use in paper 7&8TeV

       // Testing for another way for PF iso E
       relPFIso_e = ((iElectron.pfIsolationVariables()).chargedHadronIso +
		      (iElectron.pfIsolationVariables()).neutralHadronIso +
		      (iElectron.pfIsolationVariables()).photonIso) / iElectron.pt();

	
      h_relPFIso_e->Fill(relPFIso_e);
      h_Isotrk_e->Fill(Isotrk_e);

      h_relPFIso_pt_e->Fill(relPFIso_e, iElectron.pt());

      h_dxy_e->Fill((iElectron.gsfTrack())->dxy(point));

      // Good Electron with loose cut (implement style sasha)
      if (iElectron.pt() > 7.)
	{
	  if ((std::abs((iElectron.superCluster())->eta())) < 2.5)
	    {
	      if (misshits <= 1 && std::abs(SIP3d_e) < 4.) //Isotrk_e < 0.7 
		{
		  if (std::abs(iElectron.gsfTrack()->dxy(point)) < 0.5 && std::abs(iElectron.gsfTrack()->dz(point)) < 1.)
		    {
		      if (iElectron.isEB())
			{
			  // if (relPFIso_e < 1.55)
			  if (relPFIso_e < 0.4)
			    {
			      vIdPteloose.push_back( std::make_pair(te, iElectron.pt()) );
			    }
			}
		      else if (iElectron.isEE()) {
			//if (relPFIso_e < 2.4)
			if (relPFIso_e < 0.4)
			  {
			    vIdPteloose.push_back( std::make_pair(te, iElectron.pt()) );
			  }
		      }
		    }
		}
	    }
	}
       }
      /*
      // Good Electron with loose cuts
      if (std::abs(iElectron.eta()) < 2.5) 
	{
	  if (iElectron.et() > 7.) 
	    {
	      if (iElectron.isEB()) // for barrel
		{ 
		  // Relative Isolation
		  if ((iElectron.dr03TkSumPt() / (iElectron.pt())) < 0.12 && (iElectron.dr03EcalRecHitSumEt() / (iElectron.pt())) < 0.09 && (iElectron.dr03HcalTowerSumEt() / (iElectron.pt())) < 0.10 && misshits <= 1 && iElectron.convDist() < 0.02 && iElectron.convDcot() < 0.02) // WP90
		    {
		      // Electron Id
		      if (iElectron.sigmaIetaIeta() < 0.01 && std::abs(iElectron.deltaPhiSuperClusterTrackAtVtx()) < 0.8 && std::abs(iElectron.deltaEtaSuperClusterTrackAtVtx()) < 0.007 && iElectron.hcalOverEcal() < 0.12) // WP90
			{
			  listeloose.push_back(te); // this index is for Z and 4e
			}
		    }
		}

	      else if (iElectron.isEE()) //for endcap
		{
		  if ((iElectron.dr03TkSumPt() / (iElectron.pt())) < 0.05 && (iElectron.dr03EcalRecHitSumEt() / (iElectron.pt())) < 0.06 && (iElectron.dr03HcalTowerSumEt() / (iElectron.pt())) < 0.03 && misshits <= 1 && iElectron.convDist() < 0.02 && iElectron.convDcot() < 0.02) // WP90
		    {
		      if (iElectron.sigmaIetaIeta() < 0.03 && std::abs(iElectron.deltaPhiSuperClusterTrackAtVtx()) < 0.7 && std::abs(iElectron.deltaEtaSuperClusterTrackAtVtx()) < 0.009 && iElectron.hcalOverEcal() < 0.05) // WP90
			{
			  listeloose.push_back(te);
			}
		    }
		}
			  
	      h_et_e_loose_after->Fill(iElectron.et());
	      h_eta_e_loose_after->Fill(iElectron.eta());

	    } // end of if(iElectron.et()>20)
	} // end of if(((std::abs(iElectron.eta()) < 2.5
      */
      
    } // for (unsigned te = 0; te < electrons->size(); te++)
    
  std::sort(vIdPteloose.begin(), vIdPteloose.end(), [](const std::pair<int, double> &idPteloose1, const std::pair<int, double> &idPteloose2) {return (idPteloose1.second > idPteloose2.second);});

  ///////////////////////// Electron Collection end ///////////////////////////////

  // these index (list...) is declare as good candidates in integer
  nGoodGlobalMuon = vIdPt.size(); // good candidate for globalmuon 
  nGoodMuonJpsi = vIdPtjpsi.size(); // good candidate for globalmuon 
  nGoodMuonLoose = vIdPtmuloose.size(); 
  nGoodElectronLoose = vIdPteloose.size(); 

  /*nGoodGlobalMuon = listgmuons.size(); // good candidate for globalmuon 
  nGoodMuonJpsi = listmu4jpsi.size(); // good candidate for jpsi
  nGoodMuonLoose = listmuloose.size(); 
  nGoodMuonTight = listmutight.size();
  nGoodElectronLoose = listeloose.size(); 
  nGoodElectronTight = listetight.size();
  */

  h_nggmu->Fill(nGoodGlobalMuon);
  h_ngmu_jpsi->Fill(nGoodMuonJpsi);
  h_ngmu_loose->Fill(nGoodMuonLoose);
  h_nge_loose->Fill(nGoodElectronLoose);

  /////////////////////// All calculation start here /////////////////////////////

  //===== Dimuon using Global Muon start =====//

  // For the case of nGoodGlobalMuon > 2, the selected two muons are always have the highest pT as we have sorted it before
  if (nGoodGlobalMuon >= 2)
    {
      
      const reco::Track &gmuon1 = (*gmuons)[vIdPt.at(0).first];
      const reco::Track &gmuon2 = (*gmuons)[vIdPt.at(1).first];

      //     const reco::Track &gmuon1 = (*gmuons)[listgmuons.at(0)];
      //     const reco::Track &gmuon2 = (*gmuons)[listgmuons.at(1)];
      
      // The sum of 2 charge global muon should be 0 (neutral)
      if (gmuon1.charge() + gmuon2.charge() == 0)
	{
	  for (unsigned i = 0; i < vIdPt.size(); i++)
	    {
	      // These pT and eta are fill after all the cuts
	      h_pt_gmu_after->Fill(vIdPt.at(i).second); // access directly .second as the second pair is already pT
	      h_eta_gmu_after->Fill(((*gmuons)[vIdPt.at(i).first]).eta());
	    }

	  s1 = sqrt(((gmuon1.p()) * (gmuon1.p()) + sqm1) * ((gmuon2.p()) * (gmuon2.p()) + sqm1));
	  s2 = gmuon1.px() * gmuon2.px() + gmuon1.py() * gmuon2.py() + gmuon1.pz() * gmuon2.pz();
	  s = sqrt(2.0 * (sqm1 + (s1 - s2)));

	  h_m1_gmu->Fill(s);
	  h_m2_gmu->Fill(s);
	  h_m3_gmu->Fill(s);

	}
    }

  //===== Dimuon using Global Muon end =====//


  //===== JpsiTo2Muon start =====//

  if (nGoodMuonJpsi >= 2)
    {
      const reco::Muon &muon1 = (*muons)[vIdPtjpsi.at(0).first];
      const reco::Muon &muon2 = (*muons)[vIdPtjpsi.at(1).first];

      //const reco::Muon &muon1 = (*muons)[listmu4jpsi.at(0)];
      //const reco::Muon &muon2 = (*muons)[listmu4jpsi.at(1)];
      
      // The sum of 2 charges should be 0 (neutral)
      if (muon1.charge() + muon2.charge() == 0)
	{
	  for (unsigned i = 0; i < vIdPtjpsi.size(); i++)
	    {
	      h_pt_tmu_jpsi_after->Fill(vIdPtjpsi.at(i).second);
	      h_eta_tmu_jpsi_after->Fill(((*muons)[vIdPtjpsi.at(i).first]).eta());
	    }

	    // h_pt_tmu_jpsi_after->Fill(itMuon.pt());
	    // h_eta_tmu_jpsi_after->Fill(itMuon.eta());
	  s1 = sqrt(((muon1.p()) * (muon1.p()) + sqm1) * ((muon2.p()) * (muon2.p()) + sqm1)); // E1E2
	  s2 = muon1.px() * muon2.px() + muon1.py() * muon2.py() + muon1.pz() * muon2.pz(); // p1p2
	  s = sqrt(2.0 * (sqm1 + (s1 - s2)));

	  // rapidity
	  s3 = (muon1.p()) * (muon1.p()) + (muon2.p()) * (muon2.p());
	  s4 = sqrt(s3 + 2 * s2 + s * s); // s4 = energy
	  pz = muon2.pz() + muon1.pz();
	  rap = 0.5 * log((s4 + pz) / (s4 - pz));

	  dx = (muon1.globalTrack())->vx() - (muon2.globalTrack())->vx();
	  dy = (muon1.globalTrack())->vy() - (muon2.globalTrack())->vy();
	  dz = (muon1.globalTrack())->vz() - (muon2.globalTrack())->vz();

	  if (dx <= 0.1 && dy <= 0.1 && dz <= 0.3)
	    {
	      if (std::abs(rap)<1.2) { 
		h_m1_jpsi->Fill(s);
 		h_m4_jpsi->Fill(s);
	      }
	      else if (std::abs(rap)<1.6) {
		h_m2_jpsi->Fill(s);
		h_m5_jpsi->Fill(s);
	      }
	      else if (std::abs(rap)<2.4) {
		h_m3_jpsi->Fill(s);
		h_m6_jpsi->Fill(s);
	      }
	    }
	}
    } // end of if (nGoodMuonJpsi == 2)

  //===== JpsiTo2Muon end =====//

  //===== ZTo2Muon using GlobalTrack Muon foor loose cut start =====//

  if (nGoodMuonLoose >= 2)
    {
      const reco::Muon &muon1 = (*muons)[vIdPtmuloose.at(0).first];
      const reco::Muon &muon2 = (*muons)[vIdPtmuloose.at(1).first];

      // const reco::Muon &muon1 = (*muons)[listmuloose.at(0)];
      // const reco::Muon &muon2 = (*muons)[listmuloose.at(1)];

      if (muon1.charge() + muon2.charge() == 0)
	{
	  for (unsigned i = 0; i < vIdPtmuloose.size(); i++)
	    {
	      // These pT and eta are fill after all the cuts
	      h_pt_tmu_loose_after_Zto2mu->Fill(vIdPtmuloose.at(i).second); // access directly .second as the second pair is already pT
	      h_eta_tmu_loose_after_Zto2mu->Fill(((*muons)[vIdPtmuloose.at(i).first]).eta());
	    }

	  s1 = sqrt(((muon1.p()) * (muon1.p()) + sqm1) * ((muon2.p()) * (muon2.p()) + sqm1));
	  s2 = muon1.px() * muon2.px() + muon1.py() * muon2.py() + muon1.pz() * muon2.pz();
	  s = sqrt(2.0 * (sqm1 + (s1 - s2)));

	  h_mZ_2mu_loose->Fill(s);
	}
    }

  //===== ZTo2Muon using GlobalTrack Muon for loose cut end =====//


  //===== ZTo2Electron loose cut start =====//

  if (nGoodElectronLoose >= 2)
    {
      const reco::GsfElectron &elec1 = (*electrons)[vIdPteloose.at(0).first];
      const reco::GsfElectron &elec2 = (*electrons)[vIdPteloose.at(1).first];

      if (elec1.charge() + elec2.charge() == 0)
	{
	  for (unsigned i = 0; i < vIdPteloose.size(); i++)
	    {
	      // These pT and eta are fill after all the cuts
	      h_pt_e_loose_after_Zto2e->Fill(vIdPteloose.at(i).second); // access directly .second as the second pair is already pT
	      h_eta_e_loose_after_Zto2e->Fill((((*electrons)[vIdPteloose.at(i).first]).superCluster())->eta());
	    }

	  s1 = sqrt(((elec1.p()) * (elec1.p()) + sqme) * ((elec2.p()) * (elec2.p()) + sqme));
	  s2 = elec1.px() * elec2.px() + elec1.py() * elec2.py() + elec1.pz() * elec2.pz();
	  s = sqrt(2.0 * (sqme + (s1 - s2)));

	  h_mZ_2e_loose->Fill(s);
	  h_m_2e_loose->Fill(s); // mass of dielectron with isolation
	}
    }

  //===== ZTo2Electron loose cut end =====//


  //===== ZZ/ZZ*To4Muon with loose cuts start =====//

  // Now, for these goodmuons, pair up and calculate mass
  if (nGoodMuonLoose >= 4)
    {
      const reco::Muon &muon1 = (*muons)[vIdPtmuloose.at(0).first];
      const reco::Muon &muon2 = (*muons)[vIdPtmuloose.at(1).first];
      const reco::Muon &muon3 = (*muons)[vIdPtmuloose.at(2).first];
      const reco::Muon &muon4 = (*muons)[vIdPtmuloose.at(3).first];

      if (muon1.charge() + muon2.charge() + muon3.charge() + muon4.charge() == 0)
	{

	  for (unsigned i = 0; i < vIdPtmuloose.size(); i++)
	    {
	      relPFIso_mu = ((((*muons)[vIdPtmuloose.at(i).first]).pfIsolationR04()).sumChargedHadronPt + (((*muons)[vIdPtmuloose.at(i).first]).pfIsolationR04()).sumNeutralHadronEt + (((*muons)[vIdPtmuloose.at(i).first]).pfIsolationR04()).sumPhotonEt) / (((*muons)[vIdPtmuloose.at(i).first]).pt()); // relPFiso use in paper 7&8TeV

	      h_relPFIso_mu_loose_after->Fill(relPFIso_mu);

	      h_pt_tmu_loose_after->Fill(vIdPtmuloose.at(i).second);
	      h_eta_tmu_loose_after->Fill(((*muons)[vIdPtmuloose.at(i).first]).eta());
	    }

	  // First combination: Combine muon 1234
	  if (muon1.charge() + muon2.charge() == 0) // each lepton pair cas = 0
	    {
	      eZ12 = (sqrt(muon1.p() * muon1.p() + sqm1)) +
		     (sqrt(muon2.p() * muon2.p() + sqm1));
	
	      pxZ12 = muon1.px() + muon2.px();
	      pyZ12 = muon1.py() + muon2.py();
	      pzZ12 = muon1.pz() + muon2.pz();

	      if (muon3.charge() + muon4.charge() == 0)
		{
		  eZ34 = (sqrt(muon3.p() * muon3.p() + sqm1)) +
		         (sqrt(muon4.p() * muon4.p() + sqm1));

		  pxZ34 = muon3.px() + muon4.px();
		  pyZ34 = muon3.py() + muon4.py();
		  pzZ34 = muon3.pz() + muon4.pz();

		  // Calculate p4
		  pZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12) + (pzZ12 * pzZ12));
		  pZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34) + (pzZ34 * pzZ34));

		  pTZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12));
		  pTZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34));

		  mZ12 = sqrt((eZ12 * eZ12) - (pZ12 * pZ12));
		  mZ34 = sqrt((eZ34 * eZ34) - (pZ34 * pZ34));

		  if (mZ12 > 0.) h_mZ12_4mu_loose->Fill(mZ12);
		  if (mZ34 > 0.) h_mZ34_4mu_loose->Fill(mZ34);
		}
	    }

	  dZ12 = std::abs( mZ12 - mZ );
	  dZ34 = std::abs( mZ34 - mZ );

	  // take the smallest diff between mass to use for 4muon combination
	  dZc1 = (dZ12 < dZ34) ? dZ12 : dZ34; 
 
	  // Second combination: Combine muon 1324
	  if (muon1.charge() + muon3.charge() == 0)
	    {
	      eZ13 = (sqrt(muon1.p() * muon1.p() + sqm1)) +
		     (sqrt(muon3.p() * muon3.p() + sqm1));

	      pxZ13 = muon1.px() + muon3.px();
	      pyZ13 = muon1.py() + muon3.py();
	      pzZ13 = muon1.pz() + muon3.pz();

	      if (muon2.charge() + muon4.charge() == 0)
		{
		  eZ24 = (sqrt(muon2.p() * muon2.p() + sqm1)) +
		         (sqrt(muon4.p() * muon4.p() + sqm1));

		  pxZ24 = muon2.px() + muon4.px();
		  pyZ24 = muon2.py() + muon4.py();
		  pzZ24 = muon2.pz() + muon4.pz();

		  // Calculate p4
		  pZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13) + (pzZ13 * pzZ13));
		  pZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24) + (pzZ24 * pzZ24));

		  pTZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13));
		  pTZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24));

		  mZ13 = sqrt((eZ13 * eZ13) - (pZ13 * pZ13));
		  mZ24 = sqrt((eZ24 * eZ24) - (pZ24 * pZ24));

		  if (mZ13 > 0.) h_mZ13_4mu_loose->Fill(mZ13);
		  if (mZ24 > 0.) h_mZ24_4mu_loose->Fill(mZ24);
		}
	    }

	  dZ13 = std::abs( mZ13 - mZ );
	  dZ24 = std::abs( mZ24 - mZ );

	  dZc2 = (dZ13 < dZ24) ? dZ13 : dZ24; 

	  // Third combination: Combine muon 1423
	  if (muon1.charge() + muon4.charge() == 0)
	    {
	      eZ14 = (sqrt(muon1.p() * muon1.p() + sqm1)) +
		     (sqrt(muon4.p() * muon4.p() + sqm1));

	      pxZ14 = muon1.px() + muon4.px();
	      pyZ14 = muon1.py() + muon4.py();
	      pzZ14 = muon1.pz() + muon4.pz();

	      if (muon2.charge() + muon3.charge() == 0)
		{
		  eZ23 = sqrt((muon2.p() * muon2.p() + sqm1)) +
		         (sqrt(muon3.p() * muon3.p() + sqm1));
       
		  pxZ23 = muon2.px() + muon3.px();
		  pyZ23 = muon2.py() + muon3.py();
		  pzZ23 = muon2.pz() + muon3.pz();

		  // Calculate p4
		  pZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14) + (pzZ14 * pzZ14));
		  pZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23) + (pzZ23 * pzZ23));

		  pTZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14));
		  pTZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23));

		  mZ14 = sqrt((eZ14 * eZ14) - (pZ14 * pZ14));
		  mZ23 = sqrt((eZ23 * eZ23) - (pZ23 * pZ23));

		  if (mZ14 > 0.) h_mZ14_4mu_loose->Fill(mZ14);
		  if (mZ23 > 0.) h_mZ23_4mu_loose->Fill(mZ23);
		}
	    }

	  dZ14 = std::abs( mZ14 - mZ );
	  dZ23 = std::abs( mZ23 - mZ );

	  dZc3 = (dZ14 < dZ23) ? dZ14 : dZ23;

	  /*
	  if (dZc1 < dZc2 && dZc1 < dZc3)
	    {
	      eZa = (mZ12 > mZ34) ? eZ12 : eZ34;    // eZ12 
	      pxZa = (mZ12 > mZ34) ? pxZ12 : pxZ34;
	      pyZa = (mZ12 > mZ34) ? pyZ12 : pyZ34;
	      pzZa = (mZ12 > mZ34) ? pzZ12 : pzZ34;
	      mZa = (mZ12 > mZ34) ? mZ12 : mZ34;

	      eZb = (mZ12 <= mZ34) ? eZ12 : eZ34;   // eZ34
	      pxZb = (mZ12 <= mZ34) ? pxZ12 : pxZ34;
	      pyZb = (mZ12 <= mZ34) ? pyZ12 : pyZ34;
	      pzZb = (mZ12 <= mZ34) ? pzZ12 : pzZ34;
	      mZb = (mZ12 <= mZ34) ? mZ12 : mZ34;
	    }

	  else if (dZc2 < dZc1 && dZc2 < dZc3)
	    {
	      eZa = (mZ13 > mZ24) ? eZ13 : eZ24;   // eZ13
	      pxZa = (mZ13 > mZ24) ? pxZ13 : pxZ24;
	      pyZa = (mZ13 > mZ24) ? pyZ13 : pyZ24;
	      pzZa = (mZ13 > mZ24) ? pzZ13 : pzZ24;
	      mZa = (mZ13 > mZ24) ? mZ13 : mZ24;

	      eZb = (mZ13 <= mZ24) ? eZ13 : eZ24;  // eZ24
	      pxZb = (mZ13 <= mZ24) ? pxZ13 : pxZ24;
	      pyZb = (mZ13 <= mZ24) ? pyZ13 : pyZ24;
	      pzZb = (mZ13 <= mZ24) ? pzZ13 : pzZ24;
	      mZb = (mZ13 <= mZ24) ? mZ13 : mZ24;
	    }

	  else if (dZc3 < dZc1 && dZc3 < dZc2)
	    {
	      eZa = (mZ14 > mZ23) ? eZ14 : eZ23;   // eZ14
	      pxZa = (mZ14 > mZ23) ? pxZ14 : pxZ23;
	      pyZa = (mZ14 > mZ23) ? pyZ14 : pyZ23;
	      pzZa = (mZ14 > mZ23) ? pzZ14 : pzZ23;
	      mZa = (mZ14 > mZ23) ? mZ14 : mZ23;

	      eZb = (mZ14 <= mZ23) ? eZ14 : eZ23;  // eZ23
	      pxZb = (mZ14 <= mZ23) ? pxZ14 : pxZ23;
	      pyZb = (mZ14 <= mZ23) ? pyZ14 : pyZ23;
	      pzZb = (mZ14 <= mZ23) ? pzZ14 : pzZ23;
	      mZb = (mZ14 <= mZ23) ? mZ14 : mZ23;
	    }
*/

	  
	   bool ptZadaug = false;

	  if (dZc1 < dZc2 && dZc1 < dZc3)
	    {
	      if (dZ12 < dZ34)
		{
		  eZa  = eZ12;     
		  pxZa = pxZ12;
		  pyZa = pyZ12;
		  pzZa = pzZ12;
		  pTZa = pTZ12;
		  mZa  = mZ12;

		  if (muon1.pt() > 20. and muon2.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ34;
		  pxZb = pxZ34;
		  pyZb = pyZ34;
		  pzZb = pzZ34;
		  pTZb = pTZ34;
		  mZb  = mZ34;
		}
	      else
		{
		  eZa  = eZ34;
		  pxZa = pxZ34;
		  pyZa = pyZ34;
		  pzZa = pzZ34;
		  pTZa = pTZ34;
		  mZa  = mZ34;

		  if (muon3.pt() > 20. and muon4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ12;
		  pxZb = pxZ12;
		  pyZb = pyZ12;
		  pzZb = pzZ12;
		  pTZb = pTZ12;
		  mZb  = mZ12;
		}
	    }

	  else if (dZc2 < dZc1 && dZc2 < dZc3)
	    {
	      if (dZ13 < dZ24)
		{
		  eZa  = eZ13;
		  pxZa = pxZ13;
		  pyZa = pyZ13;
		  pzZa = pzZ13;
		  pTZa = pTZ13;
		  mZa  = mZ13;

		  if (muon1.pt() > 20. and muon3.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ24;
		  pxZb = pxZ24;
		  pyZb = pyZ24;
		  pzZb = pzZ24;
		  pTZb = pTZ24;
		  mZb  = mZ24;
		}
	      else
		{
		  eZa  = eZ24;
		  pxZa = pxZ24;
		  pyZa = pyZ24;
		  pzZa = pzZ24;
		  pTZa = pTZ24;
		  mZa  = mZ24;

		  if (muon2.pt() > 20. and muon4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ13;
		  pxZb = pxZ13;
		  pyZb = pyZ13;
		  pzZb = pzZ13;
		  pTZb = pTZ13;
		  mZb  = mZ13;
		}
	    }

	  else if (dZc3 < dZc1 && dZc3 < dZc2)
	    {
	      if (dZ14 < dZ23)
		{
		  eZa  = eZ14;
		  pxZa = pxZ14;
		  pyZa = pyZ14;
		  pzZa = pzZ14;
		  pTZa = pTZ14;
		  mZa  = mZ14;

		  if (muon1.pt() > 20. and muon4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ23;
		  pxZb = pxZ23;
		  pyZb = pyZ23;
		  pzZb = pzZ23;
		  pTZb = pTZ23;
		  mZb  = mZ23;
		}
	      else
		{
		  eZa  = eZ23;
		  pxZa = pxZ23;
		  pyZa = pyZ23;
		  pzZa = pzZ23;
		  pTZa = pTZ23;
		  mZa  = mZ23;

		  if (muon2.pt() > 20. and muon3.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ14;
		  pxZb = pxZ14;
		  pyZb = pyZ14;
		  pzZb = pzZ14;
		  pTZb = pTZ14;
		  mZb  = mZ14;
		}
	    }

	   if (ptZadaug) {
	    if (mZa > 40. && mZa < 120.) {
	      if (mZb > 12. && mZb < 120.) {

		h_mZa_4mu_loose->Fill(mZa);
		h_mZb_4mu_loose->Fill(mZb);

	  //	  pT_Za = sqrt((pxZa * pxZa) + (pyZa * pyZa));

	  // Calculate 4 muon
		p4Za.SetPxPyPzE(pxZa, pyZa, pzZa, eZa);
		p4Zb.SetPxPyPzE(pxZb, pyZb, pzZb, eZb);

		p4H = p4Za + p4Zb;

		mass4mu = p4H.M();
		pt_4mu = p4H.Pt();
		eta_4mu = p4H.Eta();
		
		/*
		mass4mu = sqrt( (mZa * mZa + mZb * mZb) + (2 * (eZa * eZb)) -
				(2 * (pxZa * pxZb + pyZa * pyZb + pzZa * pzZb)) );
		pt_4mu = sqrt( (pxZa * pxZa + pxZb * pxZb + (2 * pxZa * pxZb)) + (pyZa * pyZa + pyZb * pyZb + (2 * pyZa * pyZb)) );
		eta_4mu = -log( ((pxZa + pxZb) / p_4mu) / ((1 + (pzZa + pzZb) / p_4mu)) );
		*/

		pt_mu1 = muon1.pt();
		pt_mu2 = muon2.pt();
		pt_mu3 = muon3.pt();
		pt_mu4 = muon4.pt();

		eta_mu1 = muon1.eta();
		eta_mu2 = muon2.eta();
		eta_mu3 = muon3.eta();
		eta_mu4 = muon4.eta();

		
		if (mass4mu > 70.)
		  {
		    h_m1_m4mu_loose->Fill(mass4mu);
		    h_m2_m4mu_loose->Fill(mass4mu);
		    h_m3_m4mu_loose->Fill(mass4mu);
		    h_m4_m4mu_loose->Fill(mass4mu);

		    t1->Fill();
		  }
	  //h_multiplicity_4mu->Fill();
	  // pt_4mu = sqrt( (pxZa * pxZa + pxZb * pxZb + (2 * pxZa * pxZb)) + (pyZa * pyZa + pyZb * pyZb + (2 * pyZa * pyZb)) );
	  // h_pt_4mu_loose->Fill(pt_4mu);
	      }
	    }
	  }
	  // eta_4mu = -log( ((pxZa + pxZb) / p_4mu) / ((1 + (pzZa + pzZb) / p_4mu)) );
	  // h_eta_4mu_loose->Fill(eta_4mu);

	  // mu->Fill();
	} // end of total charge
    } // end of nGoodMuon

  //===== 4 muons calculation with loose cuts end =====//


  //===== 4 electrons calculation with loose cuts start =====//
  // for (int i = 0; i < vIdPteloose.size(); i++)
    //    h->Fill( (*electrons)[vIdPteloose.at(i).first].eta() );

  // Now, for these goodelectrons, pair up and calculate mass
  if (nGoodElectronLoose >= 4)
    {
      const reco::GsfElectron &elec1 = (*electrons)[vIdPteloose.at(0).first];
      const reco::GsfElectron &elec2 = (*electrons)[vIdPteloose.at(1).first];
      const reco::GsfElectron &elec3 = (*electrons)[vIdPteloose.at(2).first];
      const reco::GsfElectron &elec4 = (*electrons)[vIdPteloose.at(3).first];

      /*
      const reco::GsfElectron &elec1 = (*electrons)[listeloose.at(0)];
      const reco::GsfElectron &elec2 = (*electrons)[listeloose.at(1)];
      const reco::GsfElectron &elec3 = (*electrons)[listeloose.at(2)];
      const reco::GsfElectron &elec4 = (*electrons)[listeloose.at(3)];
      */
      if (elec1.charge() + elec2.charge() + elec3.charge() + elec4.charge() == 0)
	{

	  for (unsigned i = 0; i < vIdPteloose.size(); i++)
	    {
	      // relPFIso_e = ( ((*electrons)[vIdPteloose.at(i).first]).dr03TkSumPt() + ((*electrons)[vIdPteloose.at(i).first]).dr03EcalRecHitSumEt() + ((*electrons)[vIdPteloose.at(i).first]).dr03HcalTowerSumEt()) / (vIdPteloose.at(i).second);

	      relPFIso_e = ((((*electrons)[vIdPteloose.at(i).first]).pfIsolationVariables()).chargedHadronIso +
			    (((*electrons)[vIdPteloose.at(i).first]).pfIsolationVariables()).neutralHadronIso +
			    (((*electrons)[vIdPteloose.at(i).first]).pfIsolationVariables()).photonIso) / (((*electrons)[vIdPteloose.at(i).first]).pt()); 

	      h_relPFIso_e_loose_after->Fill(relPFIso_e);

	      h_pt_e_loose_after->Fill(vIdPteloose.at(i).second);
	      h_eta_e_loose_after->Fill((((*electrons)[vIdPteloose.at(i).first]).superCluster())->eta());
	    }

	  // First combination: Combine elec 1234
	  if (elec1.charge() + elec2.charge() == 0)
	    {
	      eZ12 = (sqrt(elec1.p() * elec1.p() + sqme)) + (sqrt(elec2.p() * elec2.p() + sqme));
	
	      pxZ12 = elec1.px() + elec2.px();
	      pyZ12 = elec1.py() + elec2.py();
	      pzZ12 = elec1.pz() + elec2.pz();

	      if (elec3.charge() + elec4.charge() == 0)
		{
		  eZ34 = (sqrt(elec3.p() * elec3.p() + sqme)) + (sqrt(elec4.p() * elec4.p() + sqme));

		  pxZ34 = elec3.px() + elec4.px();
		  pyZ34 = elec3.py() + elec4.py();
		  pzZ34 = elec3.pz() + elec4.pz();

		  // Calculate the momentum and invariant mass of elec 1 and 2, elec 3 and 4
		  pZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12) + (pzZ12 * pzZ12));
		  pZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34) + (pzZ34 * pzZ34));

		  pTZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12));
		  pTZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34));

		  mZ12 = sqrt((eZ12 * eZ12) - (pZ12 * pZ12));
		  mZ34 = sqrt((eZ34 * eZ34) - (pZ34 * pZ34));

		  if (mZ12 > 0.) h_mZ12_4e_loose->Fill(mZ12);
		  if (mZ34 > 0.) h_mZ34_4e_loose->Fill(mZ34);
		}
	    }

	  dZ12 = std::abs( mZ12 - mZ );
	  dZ34 = std::abs( mZ34 - mZ );

	  dZc1 = (dZ12 < dZ34) ? dZ12 : dZ34; // take the smallest diff between mass to use for 4elec combination

	  // Second combination: Combine elec 1324
	  if (elec1.charge() + elec3.charge() == 0)
	    {

	      eZ13 = (sqrt(elec1.p() * elec1.p() + sqme)) + (sqrt(elec3.p() * elec3.p() + sqme));

	      pxZ13 = elec1.px() + elec3.px();
	      pyZ13 = elec1.py() + elec3.py();
	      pzZ13 = elec1.pz() + elec3.pz();

	      if (elec2.charge() + elec4.charge() == 0)
		{

		  eZ24 = (sqrt(elec2.p() * elec2.p() + sqme)) + (sqrt(elec4.p() * elec4.p() + sqme));

		  pxZ24 = elec2.px() + elec4.px();
		  pyZ24 = elec2.py() + elec4.py();
		  pzZ24 = elec2.pz() + elec4.pz();

		  // Calculate the momentum and invariant mass of elec 1 and 3, elec 2 and 4
		  pZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13) + (pzZ13 * pzZ13));
		  pZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24) + (pzZ24 * pzZ24));

		  pTZ13 = sqrt((pxZ13 * pxZ13) + (pyZ13 * pyZ13) + (pzZ13 * pzZ13));
		  pTZ24 = sqrt((pxZ24 * pxZ24) + (pyZ24 * pyZ24) + (pzZ24 * pzZ24));

		  mZ13 = sqrt((eZ13 * eZ13) - (pZ13 * pZ13));
		  mZ24 = sqrt((eZ24 * eZ24) - (pZ24 * pZ24));

		  if (mZ13 > 0.) h_mZ13_4e_loose->Fill(mZ13);
		  if (mZ24 > 0.) h_mZ24_4e_loose->Fill(mZ24);
		}
	    }

	  dZ13 = std::abs( mZ13 - mZ );
	  dZ24 = std::abs( mZ24 - mZ );

	  dZc2 = (dZ13 < dZ24) ? dZ13 : dZ24; 

	  // Third combination: Combine elec 1423
	  if (elec1.charge() + elec4.charge() == 0)
	    {

	      eZ14 = (sqrt(elec1.p() * elec1.p() + sqme)) + (sqrt(elec4.p() * elec4.p() + sqme));

	      pxZ14 = elec1.px() + elec4.px();
	      pyZ14 = elec1.py() + elec4.py();
	      pzZ14 = elec1.pz() + elec4.pz();

	      if (elec2.charge() + elec3.charge() == 0)
		{

		  eZ23 = sqrt((elec2.p() * elec2.p() + sqme)) + (sqrt(elec3.p() * elec3.p() + sqme));
       
		  pxZ23 = elec2.px() + elec3.px();
		  pyZ23 = elec2.py() + elec3.py();
		  pzZ23 = elec2.pz() + elec3.pz();

		  // Calculate the momentum and invariant mass of elec 1 and 4, elec 2 and 3
		  pZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14) + (pzZ14 * pzZ14));
		  pZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23) + (pzZ23 * pzZ23));

		  pTZ14 = sqrt((pxZ14 * pxZ14) + (pyZ14 * pyZ14) + (pzZ14 * pzZ14));
		  pTZ23 = sqrt((pxZ23 * pxZ23) + (pyZ23 * pyZ23) + (pzZ23 * pzZ23));

		  mZ14 = sqrt((eZ14 * eZ14) - (pZ14 * pZ14));
		  mZ23 = sqrt((eZ23 * eZ23) - (pZ23 * pZ23));

		  if (mZ14 > 0.) h_mZ14_4e_loose->Fill(mZ14);
		  if (mZ23 > 0.) h_mZ23_4e_loose->Fill(mZ23);
		}
	    }

	  dZ14 = std::abs( mZ14 - mZ );
	  dZ23 = std::abs( mZ23 - mZ );

	  dZc3 = (dZ14 < dZ23) ? dZ14 : dZ23; 

	  bool ptZadaug = false;

	  // Now whichever have the smallest diff is considered the best comb 
	  if (dZc1 < dZc2 && dZc1 < dZc3)
	    {
	      if (dZ12 < dZ34)
		{
		  eZa  = eZ12;     
		  pxZa = pxZ12;
		  pyZa = pyZ12;
		  pzZa = pzZ12;
		  pTZa = pTZ12;
		  mZa  = mZ12;

		  if (elec1.pt() > 20. and elec2.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ34;
		  pxZb = pxZ34;
		  pyZb = pyZ34;
		  pzZb = pzZ34;
		  pTZb = pTZ34;
		  mZb  = mZ34;
		}
	      else
		{
		  eZa  = eZ34;  
		  pxZa = pxZ34;
		  pyZa = pyZ34;
		  pzZa = pzZ34;
		  pTZa = pTZ34;
		  mZa  = mZ34;

		  if (elec3.pt() > 20. and elec4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ12;
		  pxZb = pxZ12;
		  pyZb = pyZ12;
		  pzZb = pzZ12;
		  pTZb = pTZ12;
		  mZb  = mZ12;
		}
	    }

	  else if (dZc2 < dZc1 && dZc2 < dZc3)
	    {
	      if (dZ13 < dZ24)
		{
		  eZa  = eZ13;
		  pxZa = pxZ13;
		  pyZa = pyZ13;
		  pzZa = pzZ13;
		  pTZa = pTZ13;
		  mZa  = mZ13;

		  if (elec1.pt() > 20. and elec3.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ24;
		  pxZb = pxZ24;
		  pyZb = pyZ24;
		  pzZb = pzZ24;
		  pTZb = pTZ24;
		  mZb  = mZ24;
		}
	      else
		{
		  eZa  = eZ24;
		  pxZa = pxZ24;
		  pyZa = pyZ24;
		  pzZa = pzZ24;
		  pTZa = pTZ24;
		  mZa  = mZ24;

		  if (elec2.pt() > 20. and elec4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ13;
		  pxZb = pxZ13;
		  pyZb = pyZ13;
		  pzZb = pzZ13;
		  pTZb = pTZ13;
		  mZb  = mZ13;
		}
	    }

	  else if (dZc3 < dZc1 && dZc3 < dZc2)
	    {
	      if (dZ14 < dZ23)
		{
		  eZa  = eZ14;
		  pxZa = pxZ14;
		  pyZa = pyZ14;
		  pzZa = pzZ14;
		  pTZa = pTZ14;
		  mZa  = mZ14;

		  if (elec1.pt() > 20. and elec4.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ23;
		  pxZb = pxZ23;
		  pyZb = pyZ23;
		  pzZb = pzZ23;
		  pTZb = pTZ23;
		  mZb  = mZ23;
		}
	      else
		{
		  eZa  = eZ23;
		  pxZa = pxZ23;
		  pyZa = pyZ23;
		  pzZa = pzZ23;
		  pTZa = pTZ23;
		  mZa  = mZ23;

		  if (elec2.pt() > 20. and elec3.pt() > 10.)
		  ptZadaug = true;

		  eZb  = eZ14;
		  pxZb = pxZ14;
		  pyZb = pyZ14;
		  pzZb = pzZ14;
		  pTZb = pTZ14;
		  mZb  = mZ14;
		}
	    }

	  if (ptZadaug) {
	    if (mZa > 40. && mZa < 120.) {
	      if (mZb > 12. && mZb < 120.) {
		h_mZa_4e_loose->Fill(mZa);
		h_mZb_4e_loose->Fill(mZb);
		
		// Calculate 4 elec
		p4Za.SetPxPyPzE(pxZa, pyZa, pzZa, eZa);
		p4Zb.SetPxPyPzE(pxZb, pyZb, pzZb, eZb);

		p4H = p4Za + p4Zb;

		mass4e = p4H.M();
		pt_4e = p4H.Pt();
		eta_4e = p4H.Eta();

		//mass4e = sqrt( (mZa * mZa + mZb * mZb) + (2 * (eZa * eZb)) - (2 * (pxZa * pxZb + pyZa * pyZb + pzZa * pzZb)) );

		pt_e1 = elec1.pt();
		pt_e2 = elec2.pt();
		pt_e3 = elec3.pt();
		pt_e4 = elec4.pt();

		eta_e1 = elec1.eta();
		eta_e2 = elec2.eta();
		eta_e3 = elec3.eta();
		eta_e4 = elec4.eta();


		if (mass4e > 70.)
		  {
		    h_m1_m4e_loose->Fill(mass4e);
		    h_m2_m4e_loose->Fill(mass4e);
		    h_m3_m4e_loose->Fill(mass4e);
		    h_m4_m4e_loose->Fill(mass4e);

		    t2->Fill();
		  }
		// mu->Fill();
	      }
	    }
	  }
	} // end of total charge
    } // end of nGoodElectron

  //===== 4 electrons calculation with loose cuts end =====//


  //===== 2 muons 2 electrons calculation with loose cuts start =====//

  if (nGoodMuonLoose >= 2 && nGoodElectronLoose >= 2)
    {
      const reco::Muon &muon1 = (*muons)[vIdPtmuloose.at(0).first];
      const reco::Muon &muon2 = (*muons)[vIdPtmuloose.at(1).first];
      const reco::GsfElectron &elec1 = (*electrons)[vIdPteloose.at(0).first];
      const reco::GsfElectron &elec2 = (*electrons)[vIdPteloose.at(1).first];

      /*
      const reco::Muon &muon1 = (*muons)[listmuloose.at(0)];
      const reco::Muon &muon2 = (*muons)[listmuloose.at(1)];
      const reco::GsfElectron &elec1 = (*electrons)[listeloose.at(0)];
      const reco::GsfElectron &elec2 = (*electrons)[listeloose.at(1)];
      */

      if (muon1.charge() + muon2.charge() + elec1.charge() + elec2.charge() == 0)
	{

	  for (unsigned i = 0; i < vIdPtmuloose.size(); i++)
	    {
	      relPFIso_mu = ( (((*muons)[vIdPtmuloose.at(i).first]).pfIsolationR04()).sumChargedHadronPt +
			      (((*muons)[vIdPtmuloose.at(i).first]).pfIsolationR04()).sumNeutralHadronEt +
			      (((*muons)[vIdPtmuloose.at(i).first]).pfIsolationR04()).sumPhotonEt ) / (((*muons)[vIdPtmuloose.at(i).first]).pt());

	      h_relPFIso_2mu_loose_after->Fill(relPFIso_mu);
	      h_pt_tmu_loose_after_2mu2e->Fill(vIdPtmuloose.at(i).second);
	      h_eta_tmu_loose_after_2mu2e->Fill(((*muons)[vIdPtmuloose.at(i).first]).eta());
	    }

	  for (unsigned i = 0; i < vIdPteloose.size(); i++)
	    {
	      // relPFIso_e = ( ((*electrons)[vIdPteloose.at(i).first]).dr03TkSumPt() + ((*electrons)[vIdPteloose.at(i).first]).dr03EcalRecHitSumEt() + ((*electrons)[vIdPteloose.at(i).first]).dr03HcalTowerSumEt()) / (vIdPteloose.at(i).second);


	      relPFIso_e = ((((*electrons)[vIdPteloose.at(i).first]).pfIsolationVariables()).chargedHadronIso +
			    (((*electrons)[vIdPteloose.at(i).first]).pfIsolationVariables()).neutralHadronIso +
			    (((*electrons)[vIdPteloose.at(i).first]).pfIsolationVariables()).photonIso) / (((*electrons)[vIdPteloose.at(i).first]).pt()); 

	      
	      h_relPFIso_2e_loose_after->Fill(relPFIso_e);
	      h_pt_e_loose_after_2mu2e->Fill(vIdPteloose.at(i).second);
	      h_eta_e_loose_after_2mu2e->Fill(((*electrons)[vIdPteloose.at(i).first]).eta());
	    }

	  if (muon1.charge() + muon2.charge() == 0)
	    {
	      eZ12 = (sqrt(muon1.p() * muon1.p() + sqm1)) + (sqrt(muon2.p() * muon2.p() + sqm1));
	
	      pxZ12 = muon1.px() + muon2.px();
	      pyZ12 = muon1.py() + muon2.py();
	      pzZ12 = muon1.pz() + muon2.pz();

	      if (elec1.charge() + elec2.charge() == 0)
		{
		  eZ34 = (sqrt(elec1.p() * elec1.p() + sqme)) + (sqrt(elec2.p() * elec2.p() + sqme));

		  pxZ34 = elec1.px() + elec2.px();
		  pyZ34 = elec1.py() + elec2.py();
		  pzZ34 = elec1.pz() + elec2.pz();

		  // Calculate the momentum and invariant mass of elec 1 and 2, elec 1 and 2
		  pZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12) + (pzZ12 * pzZ12));
		  pZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34) + (pzZ34 * pzZ34));

		  pTZ12 = sqrt((pxZ12 * pxZ12) + (pyZ12 * pyZ12));
		  pTZ34 = sqrt((pxZ34 * pxZ34) + (pyZ34 * pyZ34));

		  mZ12 = sqrt((eZ12 * eZ12) - (pZ12 * pZ12));
		  mZ34 = sqrt((eZ34 * eZ34) - (pZ34 * pZ34));

		  if (mZ12 > 0.) h_mZmu_2mu2e_loose->Fill(mZ12);
		  if (mZ34 > 0.) h_mZe_2mu2e_loose->Fill(mZ34);

		}
	    }

	  dZ12 = std::abs(mZ12 - mZ); // mu
	  dZ34 = std::abs(mZ34 - mZ); // e

	  bool ptZadaug = false;

	  if (dZ12 < dZ34)
	    {
	      eZa  = eZ12;
	      pxZa = pxZ12;
	      pyZa = pyZ12;
	      pzZa = pzZ12;
	      pTZa = pTZ12;
	      mZa  = mZ12;

	      if (muon1.pt() > 20. and muon2.pt() > 10.)
		ptZadaug = true;

	      eZb  = eZ34;
	      pxZb = pxZ34;
	      pyZb = pyZ34;
	      pzZb = pzZ34;
	      pTZb = pTZ34;
	      mZb  = mZ34;
	    }
	  else
	    {
	      eZa  = eZ34;
	      pxZa = pxZ34;
	      pyZa = pyZ34;
	      pzZa = pzZ34;
	      pTZa = pTZ34;
	      mZa  = mZ34;

	      if (elec1.pt() > 20. and elec2.pt() > 10.)
		ptZadaug = true;

	      eZb  = eZ12;
	      pxZb = pxZ12;
	      pyZb = pyZ12;
	      pzZb = pzZ12;
	      pTZb = pTZ12;
	      mZb  = mZ12;
	  }

	  if (ptZadaug) {
	    if (mZa > 40. && mZa < 120.) {
	      if (mZb > 12. && mZb < 120.) {
		h_mZa_2mu2e_loose->Fill(mZa);
		h_mZb_2mu2e_loose->Fill(mZb);

		// Now combine these 2 muons and 2 electrons
 
		// Calculate 4 lepton: 2muons 2electrons

		p4Za.SetPxPyPzE(pxZa, pyZa, pzZa, eZa);
		p4Zb.SetPxPyPzE(pxZb, pyZb, pzZb, eZb);

		p4H = p4Za + p4Zb;

		mass2mu2e = p4H.M();
		pt_2mu2e = p4H.Pt();
		eta_2mu2e = p4H.Eta();

		//mass2mu2e = sqrt( (mZ12 * mZ12 + mZ34 * mZ34) + (2 * (eZ12 * eZ34)) - (2 * (pxZ12 * pxZ34 + pyZ12 * pyZ34 + pzZ12 * pzZ34)) );

		pt_2mu1 = muon1.pt();
		pt_2mu2 = muon2.pt();
		pt_2e1 = elec1.pt();
		pt_2e2 = elec2.pt();

		eta_2mu1 = muon1.eta();
		eta_2mu2 = muon2.eta();
		eta_2e1 = elec1.eta();
		eta_2e2 = elec2.eta();

		if (mass2mu2e > 70.)
		  {
		    h_m1_m2mu2e_loose->Fill(mass2mu2e);
		    h_m2_m2mu2e_loose->Fill(mass2mu2e);
		    h_m3_m2mu2e_loose->Fill(mass2mu2e);
		    h_m4_m2mu2e_loose->Fill(mass2mu2e);

		    t3->Fill();
		  }
		// std::cout << "massmu = " << massmu << std::endl;
	      }
	    }
	  }
	  // mu->Fill();

	} // end of total charge
    } // end of nGoodMuonNElectron

  //===== 2 muons 2 electrons calculation with loose cuts end =====//
  
} // DemoAnalyzer::analyze ends


// ------ method called once each job just before starting event loop ---------//

void DemoAnalyzer::beginJob() {

  t1 = new TTree("tree4mu", "tree4mu");
  t2 = new TTree("tree4e", "tree4e");
  t3 = new TTree("tree2mu2e", "tree2mu2e");

  //  t->Branch("nRun1", &nRun1, "nRun1/I");
  //  t->Branch("nEvt1", &nEvt1, "nEvt1/I");
  //  t->Branch("nLumi1", &nLumi1, "nLumi1/I");
  // tree 4mu
  t1->Branch("nRun", &nRun, "nRun/I");
  t1->Branch("nEvt", &nEvt, "nEvt/I");
  t1->Branch("nLumi", &nLumi, "nLumi/I");
  t1->Branch("mass4mu", &mass4mu, "mass4mu/D");
  t1->Branch("pt_4mu", &pt_4mu, "pt_4mu/D");
  t1->Branch("eta_4mu", &eta_4mu, "eta_4mu/D");
  /*  t1->Branch("pt_mu1", &pt_mu1, "pt_mu1/D");
  t1->Branch("pt_mu2", &pt_mu2, "pt_mu2/D");
  t1->Branch("pt_mu3", &pt_mu3, "pt_mu3/D");
  t1->Branch("pt_mu4", &pt_mu4, "pt_mu4/D");
  t1->Branch("eta_mu1", &eta_mu1, "eta_mu1/D");
  t1->Branch("eta_mu2", &eta_mu2, "eta_mu2/D");
  t1->Branch("eta_mu3", &eta_mu3, "eta_mu3/D");
  t1->Branch("eta_mu4", &eta_mu4, "eta_mu4/D");
  */
  // tree 4e
  t2->Branch("nRun", &nRun, "nRun/I");
  t2->Branch("nEvt", &nEvt, "nEvt/I");
  t2->Branch("nLumi", &nLumi, "nLumi/I");
  t2->Branch("mass4e", &mass4e, "mass4e/D");
  t2->Branch("pt_4e", &pt_4e, "pt_4e/D");
  t2->Branch("eta_4e", &eta_4e, "eta_4e/D");
  /*  t2->Branch("pt_e1", &pt_e1, "pt_e1/D");
  t2->Branch("pt_e2", &pt_e2, "pt_e2/D");
  t2->Branch("pt_e3", &pt_e3, "pt_e3/D");
  t2->Branch("pt_e4", &pt_e4, "pt_e4/D");
  t2->Branch("eta_e1", &eta_e1, "eta_e1/D");
  t2->Branch("eta_e2", &eta_e2, "eta_e2/D");
  t2->Branch("eta_e3", &eta_e3, "eta_e3/D");
  t2->Branch("eta_e4", &eta_e4, "eta_e4/D");
  */
  // tree 2mu 2e
  t3->Branch("nRun", &nRun, "nRun/I");
  t3->Branch("nEvt", &nEvt, "nEvt/I");
  t3->Branch("nLumi", &nLumi, "nLumi/I");
  t3->Branch("mass2mu2e", &mass2mu2e, "mass2mu2e/D");
  t3->Branch("pt_2mu2e", &pt_2mu2e, "pt_2mu2e/D");
  t3->Branch("eta_2mu2e", &eta_2mu2e, "eta_2mu2e/D");
  /*  t3->Branch("pt_2mu1", &pt_2mu1, "pt_2mu1/D");
  t3->Branch("pt_2mu2", &pt_2mu2, "pt_2mu2/D");
  t3->Branch("pt_2e1", &pt_2e1, "pt_2e1/D");
  t3->Branch("pt_2e2", &pt_2e2, "pt_2e2/D");
  t3->Branch("eta_2mu1", &eta_2mu1, "eta_2mu1/D");
  t3->Branch("eta_2mu2", &eta_2mu2, "eta_2mu2/D");
  t3->Branch("eta_2e1", &eta_2e1, "eta_2e1/D");
  t3->Branch("eta_2e2", &eta_2e2, "eta_2e2/D");
  
 
 mu->Branch("nGoodElec", &nGoodMuon, "nGoodMuon/I");
 mu->Branch("mZ12", &mZ12, "mZ12/D");
 mu->Branch("mZ34", &mZ34, "mZ34/D");
 mu->Branch("mZ13", &mZ13, "mZ13/D");
 mu->Branch("mZ24", &mZ24, "mZ24/D");
 mu->Branch("mZ14", &mZ14, "mZ14/D");
 mu->Branch("mZ23", &mZ23, "mZ23/D");
 mu->Branch("mZa", &mZa, "mZa/D");
 mu->Branch("mZb", &mZb, "mZb/D");
 mu->Branch("massmu", &massmu, "massmu/D");
  */
}

// ------------ method called once each job just after ending the event loop  ------------
void DemoAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);

