//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 12 16:51:53 2019 by ROOT version 6.12/04
// from TTree tree/tree
// found on file: ../inputTree/trackTree_MC2017_test_1.root
//////////////////////////////////////////////////////////

#ifndef MatCalc_h
#define MatCalc_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <sstream>
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TLine.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom1.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <vector>
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRint.h"
#include "TString.h"
#include "TSystem.h"
#include "Math/IFunction.h"
#include "Math/Point3D.h"
#include <cmath>


// Header file for the classes stored in the TTree if any.

class MatCalc {
public :
   TChain         *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           n;
   Float_t         z[60];   //[n]
   Float_t         eta[60];   //[n]
   Float_t         phi[60];   //[n]
   Float_t         r[60];   //[n]
   Float_t         pt[60];   //[n]
   Float_t         Bx[60];   //[n]
   Float_t         By[60];   //[n]
   Float_t         Bz[60];   //[n]
   Float_t         xUnc[60];   //[n]
   Float_t         yUnc[60];   //[n]
   Float_t         zUnc[60];   //[n]
   Float_t         etaUnc[60];   //[n]
   Float_t         phiUnc[60];   //[n]
   Float_t         rUnc[60];   //[n]
   Int_t           stereo[60];   //[n]
   Int_t           glued[60];   //[n]
   Int_t           detector[60];   //[n]
   Int_t           layer[60];   //[n]
   Float_t         trackPt;
   Float_t         trackPtErr;
   Float_t         trackEta;
   Float_t         trackPhi;
   Float_t         trackX0;
   Float_t         trackY0;
   Float_t         trackZ0;
   Float_t         trackCharge;
   Float_t         beamSpotX0;
   Float_t         beamSpotY0;
   Float_t         beamSpotZ0;
   Float_t         beamSpotdxdz;
   Float_t         beamSpotdydz;
   Float_t         beamSpotSignmaZ;
   Float_t         beamSpotX0error;
   Float_t         beamSpotY0error;
   Float_t         beamSpotZ0error;

   // List of branches
   TBranch        *b_n;   //!
   TBranch        *b_z;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_r;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_Bx;   //!
   TBranch        *b_By;   //!
   TBranch        *b_Bz;   //!
   TBranch        *b_xUnc;   //!
   TBranch        *b_yUnc;   //!
   TBranch        *b_zUnc;   //!
   TBranch        *b_etaUnc;   //!
   TBranch        *b_phiUnc;   //!
   TBranch        *b_rUnc;   //!
   TBranch        *b_stereo;   //!
   TBranch        *b_glued;   //!
   TBranch        *b_detector;   //!
   TBranch        *b_layer;   //!
   TBranch        *b_trackPt;   //!
   TBranch        *b_trackPtErr;   //!
   TBranch        *b_trackEta;   //!
   TBranch        *b_trackPhi;   //!
   TBranch        *b_trackX0;   //!
   TBranch        *b_trackY0;   //!
   TBranch        *b_trackZ0;   //!
   TBranch        *b_trackCharge;   //!
   TBranch        *b_bsX0;   //!
   TBranch        *b_bsY0;   //!
   TBranch        *b_bsZ0;   //!
   TBranch        *b_bsdxdz;   //!
   TBranch        *b_bsdydz;   //!
   TBranch        *b_bsSigmaZ;   //!
   TBranch        *b_bsX0Error;   //!
   TBranch        *b_bsY0Error;   //!
   TBranch        *b_bsZ0Error;   //!

   bool isMC_;
   bool isFirstRun_;
   float sag=-99;
   int n_pixel=0;

   MatCalc(TChain *tree=0);
   MatCalc(bool isFirstRun, bool isMC, TString dirName);
   virtual ~MatCalc();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void writeoutput();
   void MinBiasAnalysis();

   std::map<std::string, TH1F*> h_1d;
   std::map<std::string, TH1F*> h_1dm;
   std::map<std::string, TH1F*> h_1dpl;
   std::map<std::string, TH2F*> h_2d;
   std::map<std::string, TH2F*> h_2dpl;
   std::map<std::string, TH2F*> h_2dm;
   std::map<std::string, TH2F*> h_2dptpl;
   std::map<std::string, TH2F*> h_2dptm;
   std::map<std::string, TH2F*> h_2dxpl;
   std::map<std::string, TH2F*> h_2dxm;
   std::map<std::string, TH2F*> h_2dxptpl;
   std::map<std::string, TH2F*> h_2dxptm;
   std::map<std::string, TH2F*> h_2dp;
   std::map<std::string, TH2F*> h_2dp2;
   std::map<std::string, TH3F*> h_3dpl;
   std::map<std::string, TH3F*> h_3dm;
   std::map<std::string, TH3F*> h_3dbm;
   std::map<std::string, TH3F*> h_3dbpl;
   // Compute the sagitta given the radii and the phi coordinate of three point in the tracker
   float calculSag(float r1, float r2, float r3, float phi1, float phi2, float phi3){

     float a[2]={0};
     float b[2]={0};

     a[0]=r2*TMath::Cos(phi2)-r1*TMath::Cos(phi1);
     a[1]=r2*TMath::Sin(phi2)-r1*TMath::Sin(phi1);

     b[0]=r3*TMath::Cos(phi3)-r1*TMath::Cos(phi1);
     b[1]=r3*TMath::Sin(phi3)-r1*TMath::Sin(phi1);

     float sin_theta=(a[0]*b[1]-a[1]*b[0])/TMath::Sqrt(b[0]*b[0]+b[1]*b[1])/TMath::Sqrt(a[0]*a[0]+a[1]*a[1]);
     float a_par=TMath::Sqrt(a[0]*a[0]+a[1]*a[1])*TMath::Cos(TMath::ASin(sin_theta));

     if(a[0]*a[1]!=0&&b[0]*b[1]!=0) return sin_theta/(TMath::Sqrt(b[0]*b[0]+b[1]*b[1])-a_par);
     else return -99;
   }

   // Compute the squared distance from point 1 to point 3
   float calcb2(float r1, float r2, float r3, float phi1, float phi2, float phi3){

     float b[2]={0};

     b[0]=r3*TMath::Cos(phi3)-r1*TMath::Cos(phi1);
     b[1]=r3*TMath::Sin(phi3)-r1*TMath::Sin(phi1);

     float cst=(b[0]*TMath::Cos(phi2)+b[1]*TMath::Sin(phi2))/TMath::Sqrt(b[0]*b[0]+b[1]*b[1]);
     return (b[0]*b[0]+b[1]*b[1])*cst;
   }

   // Compute sin given eta
   float calcSin(float eta){

     float exp=TMath::Exp(-eta);

     return 2*exp/(1+exp*exp);
   }

   // Compute cos given eta
   float calcCos(float eta){

     float exp=TMath::Exp(-eta);

     return (1-exp*exp)/(1+exp*exp);
   }

   // Plotting functions
   void plot1D(std::string title, float xval, double weight, std::map<std::string, TH1F*> &allhistos,
       int numbinsx, float xmin, float xmax)
   {
     std::map<std::string, TH1F*>::iterator iter= allhistos.find(title);
     if(iter == allhistos.end()) //no histo for this yet, so make a new one

     {
       TH1F* currentHisto= new TH1F(title.c_str(), title.c_str(), numbinsx, xmin, xmax);
       currentHisto->Fill(xval, weight);
       allhistos.insert(std::pair<std::string, TH1F*> (title,currentHisto) );
     }
     else // exists already, so just fill it

     {
       (*iter).second->Fill(xval, weight);

     }
   }

   void plot2D(std::string title, float xval, float yval, double weight, std::map<string, TH2F*> &allhistos,
       int numbinsx, float xmin, float xmax, int numbinsy, float ymin, float ymax)
   {

     std::map<std::string, TH2F*>::iterator iter= allhistos.find(title);
     if(iter == allhistos.end()) //no histo for this yet, so make a new one

     {
       TH2F* currentHisto= new TH2F(title.c_str(), title.c_str(), numbinsx, xmin, xmax, numbinsy, ymin, ymax);
       currentHisto->Fill(xval, yval, weight);
       allhistos.insert(std::pair<std::string, TH2F*> (title,currentHisto) );
     }
     else // exists already, so just fill it

     {
       (*iter).second->Fill(xval, yval, weight);
     }

     return;

   }

   void plot3D(std::string title, float xval, float yval, float zval, double weight, std::map<std::string, TH3F*> &allhistos,
       int numbinsx, float xmin, float xmax, int numbinsy, float ymin, float ymax, int numbinsz, float zmin, float zmax)
   {

     std::map<std::string, TH3F*>::iterator iter= allhistos.find(title);
     if(iter == allhistos.end()) //no histo for this yet, so make a new one

     {
       TH3F* currentHisto= new TH3F(title.c_str(), title.c_str(), numbinsx, xmin, xmax, numbinsy, ymin, ymax, numbinsz, zmin, zmax);
       currentHisto->Fill(xval, yval, zval, weight);
       allhistos.insert(std::pair<std::string, TH3F*> (title,currentHisto) );
     }
     else // exists already, so just fill it

     {
       (*iter).second->Fill(xval, yval, zval, weight);
     }

     return;

   }

};

#endif

#ifdef MatCalc_cxx
MatCalc::MatCalc(TChain *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../inputTree/trackTree_MC2017_test_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../inputTree/trackTree_MC2017_test_1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("../inputTree/trackTree_MC2017_test_1.root:/hitanalyzer");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}
MatCalc::MatCalc(bool isFirstRun, bool isMC, TString dirName) {
	isMC_ = isMC;
	isFirstRun_ = isFirstRun;
  TSystemDirectory dir(dirName, dirName);
  TList *files = dir.GetListOfFiles();
  fChain = new TChain("hitanalyzer/tree");
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(".root")) {
        //if(!isMC_&&i>2) break; //set the number of files to run
        //cout << fname.Data() << endl;
        TString filename;
        filename = dirName+fname.Data();
        cout << filename<< endl;
        fChain->AddFile(filename);
      }
    }
    Init(fChain);
  }
}

MatCalc::~MatCalc()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MatCalc::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MatCalc::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MatCalc::Init(TChain *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("n", &n, &b_n);
   fChain->SetBranchAddress("z", z, &b_z);
   fChain->SetBranchAddress("eta", eta, &b_eta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("r", r, &b_r);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("Bx", Bx, &b_Bx);
   fChain->SetBranchAddress("By", By, &b_By);
   fChain->SetBranchAddress("Bz", Bz, &b_Bz);
   fChain->SetBranchAddress("xUnc", xUnc, &b_xUnc);
   fChain->SetBranchAddress("yUnc", yUnc, &b_yUnc);
   fChain->SetBranchAddress("zUnc", zUnc, &b_zUnc);
   fChain->SetBranchAddress("etaUnc", etaUnc, &b_etaUnc);
   fChain->SetBranchAddress("phiUnc", phiUnc, &b_phiUnc);
   fChain->SetBranchAddress("rUnc", rUnc, &b_rUnc);
   fChain->SetBranchAddress("stereo", stereo, &b_stereo);
   fChain->SetBranchAddress("glued", glued, &b_glued);
   fChain->SetBranchAddress("detector", detector, &b_detector);
   fChain->SetBranchAddress("layer", layer, &b_layer);
   fChain->SetBranchAddress("trackPt", &trackPt, &b_trackPt);
   fChain->SetBranchAddress("trackPtErr", &trackPtErr, &b_trackPtErr);
   fChain->SetBranchAddress("trackEta", &trackEta, &b_trackEta);
   fChain->SetBranchAddress("trackPhi", &trackPhi, &b_trackPhi);
   fChain->SetBranchAddress("trackX0", &trackX0, &b_trackX0);
   fChain->SetBranchAddress("trackY0", &trackY0, &b_trackY0);
   fChain->SetBranchAddress("trackZ0", &trackZ0, &b_trackZ0);
   fChain->SetBranchAddress("trackCharge", &trackCharge, &b_trackCharge);
   fChain->SetBranchAddress("beamSpotX0", &beamSpotX0, &b_bsX0);
   fChain->SetBranchAddress("beamSpotY0", &beamSpotY0, &b_bsY0);
   fChain->SetBranchAddress("beamSpotZ0", &beamSpotZ0, &b_bsZ0);
   fChain->SetBranchAddress("beamSpotdxdz", &beamSpotdxdz, &b_bsdxdz);
   //Fix This
   //    fChain->SetBranchAddress("beamSpotdxdz", &beamSpotdydz, &b_bsdydz);
   fChain->SetBranchAddress("beamSpotSignmaZ", &beamSpotSignmaZ, &b_bsSigmaZ);
   fChain->SetBranchAddress("beamSpotX0error", &beamSpotX0error, &b_bsX0Error);
   fChain->SetBranchAddress("beamSpotY0error", &beamSpotY0error, &b_bsY0Error);
   fChain->SetBranchAddress("beamSpotZ0error", &beamSpotZ0error, &b_bsZ0Error);
   Notify();
}

Bool_t MatCalc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MatCalc::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MatCalc::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MatCalc_cxx
