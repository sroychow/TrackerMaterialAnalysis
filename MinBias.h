#ifndef __CINT__
#endif
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
#include <cmath>

class MatCalc {
  public:
    MatCalc(bool isFirstRun, bool isMC, TString dirName);
    ~MatCalc(){}

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
    void MinBiasAnalysis();
		void writeoutput();

  private:
    // Define the observables that we need in the analysis
    int n=0;
    int stereo[35]={0};
    float zUnc[35]={0};
    float rUnc[35]={0};
    float phiUnc[35]={0};
    float z[35]={0};
    float eta[35]={0};
    float phi[35]={0};
    float r[35]={0};
    float pt[35]={0};
    float trackEta=0;
    float trackPhi=0;
    float trackZ0=0;
    float trackY0=0;
    float trackX0=0;
    float trackPt=-99;
    float trackPtErr=-99;
    int layer[35]={-99};
    int detector[35]={0};
    float sag=-99;
    int n_pixel=0;
    float bsX0=0.;
    float bsY0=0.;
    float bsZ0=0.;
    float bsdxdz=0.;
    float bsdydz=0.;
    float bsSigmaZ=0.;
    float bsX0Error=0.;
    float bsY0Error=0.;
    float bsZ0Error=0.;
		bool isMC_;
		bool isFirstRun_;
    TChain* Tree;
    TFile* fout;

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
};
