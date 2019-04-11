#define MatCalc_cxx
#include "MatCalc.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MatCalc::Loop()
{
//   In a ROOT session, you can do:
//      root> .L MatCalc.C
//      root> MatCalc t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      std::cout << "i=" << jentry << std::endl;
      // if (Cut(ientry) < 0) continue;
   }
}

//Function to write histos to output
void MatCalc::writeoutput() {
	TFile *fileOUT;
	// Writes 1D histos
	if(isFirstRun_){
		if(isMC_) fileOUT = new TFile("MinBiasMC_sag1D.root","RECREATE");
		else fileOUT = new TFile("MinBiasDATA_sag1D.root","RECREATE");
	}
	else{
		//TFile *fileOUT;
		if(isMC_) fileOUT = new TFile("MinBiasMC_sag1D_new.root","RECREATE");
		else fileOUT = new TFile("MinBiasDATA_sag1D_new.root","RECREATE");
	}
  fileOUT->cd();
  for(auto it1dm: h_1d) {
		it1dm.second->Write();
		delete it1dm.second;
	}

	for(auto it1dm: h_1dm) {
		it1dm.second->Write();
		delete it1dm.second;
	}

	for(auto it1dpl : h_1dpl) {
		it1dpl.second->Write();
		delete it1dpl.second;
	}

   fileOUT->Save();
   fileOUT->Close();
	// Writes 3D histos

	if(!isFirstRun_){
		TFile *file3D;
		if(isMC_) file3D= new TFile("MinBias3D_MC.root","RECREATE");
		else file3D=new TFile("MinBias3D_DATA.root","RECREATE");
        file3D->cd();
	for(auto it3d : h_3dm) {
		it3d.second->Write();
		delete it3d.second;
	}
	for(auto it3d2 : h_3dpl) {
		it3d2.second->Write();
		delete it3d2.second;
	}
        file3D->Save();
        file3D->Close();
	}

}
void MatCalc::MinBiasAnalysis(){

	int nEntries = fChain->GetEntries();

	string position[100]={""};
	string region[100]={""};
	string local_eta[100]={""};
	string triplet="";
	string tripl_name="";
	int i_next=0;
	int i_prev=0;
	int number[100]={0};
	double eta_max[100]={-99};
	double eta_min[100]={-99};
	float sin_track=-99;
	int sel_tracks=0;
	float mean_min=0;
	float mean_pl=0;
	int in_plots=0;

	TFile *fileIN;
	if(!isFirstRun_){
		if(isMC_) fileIN = TFile::Open("MinBiasMC_sag1D.root");
		else fileIN = TFile::Open("MinBiasDATA_sag1D.root");
	}

	// Loop over tracks
	for(int i = 1; i < nEntries; i++) {

	  if(i%100000==0) cout <<"Analyzed entry "     << i <<"/"<< nEntries
		                     <<" selected tracks "   << sel_tracks
												 <<" entering in plots " <<in_plots
												 <<endl;
		if(in_plots > 1000000) break;
		fChain->GetEntry(i);
		n_pixel=0;

//definition from BeamSpot.h
//double x(const double z) const { return x0() + dxdz() * (z - z0()); }
/*
  const BeamSpot::Point BeamSpot::position(const double z) const {
    Point pos(x(z),y(z),z);
    return pos;
  }
*/
		// A bunch of plots
                //std::cout << beamSpotX0 << "\t" << beamSpotY0 << "\t" << beamSpotZ0 << std::endl;
		double bsX0 = beamSpotX0 + beamSpotdxdz * (trackZ0 - beamSpotZ0);
		double bsY0 = beamSpotY0 + beamSpotdydz * (trackZ0 - beamSpotZ0);

                plot1D("beamSpotX0", bsX0, 1, h_1d, 1000, -1., 1.);
                plot1D("beamSpotY0", bsY0, 1, h_1d, 1000, -1., 1.);
                plot1D("beamSpotZ0", trackZ0, 1, h_1d, 2000, -2., 2.);

		// Quality cuts
		if(trackPt>1.5) continue;
		if(trackPt<0.75) continue;
		if(n>100) continue;
		if(n<14) continue;

		trackPtErr=trackPtErr*trackPt*trackPt; // Since Mike defined it wrong in the ntuples :)
		if(trackPtErr>0.01) continue;

		sel_tracks++;

		// Loop over hits j<n
		for(int j = 0; j < n; j++) {

			if(stereo[j]==1) continue; // avoid stereo modules

			// identify layers with strings instead of numbers
			if(detector[j] == 0){
				position[j]="_pixel"+std::to_string(layer[j]);
				region[j]="_barrel";
				number[j]=layer[j];
			}
			if(detector[j] == 1){
				position[j]="_pixeldisk"+std::to_string(layer[j]);
				if(layer[j]>0) position[j]="_pixeldisk+"+std::to_string(layer[j]);
				region[j]="_forward";
				if(layer[j]<0) region[j]="_backward";
				number[j]=TMath::Abs(layer[j]);
			}
			if(detector[j]==2){
				position[j]="_TIB"+std::to_string(layer[j]);
				region[j]="_barrel";
				number[j]=layer[j]+3;
			}
			if(detector[j]==3){
				position[j]="_TOB"+std::to_string(layer[j]);
				region[j]="_barrel";
				number[j]=layer[j]+7;
			}
			if(detector[j]==4){
				position[j]="_TID"+std::to_string(layer[j]);
				if(layer[j]>0) position[j]="_TID+"+std::to_string(layer[j]);
				region[j]="_forward";
				if(layer[j]<0) region[j]="_backward";
				number[j]=TMath::Abs(layer[j])+2;
			}
			if(detector[j]==5){
				position[j]="_TEC"+std::to_string(layer[j]);
				if(layer[j]>0) position[j]="_TEC+"+std::to_string(layer[j]);
				region[j]="_forward";
				if(layer[j]<0) region[j]="_backward";
				number[j]=TMath::Abs(layer[j])+5;
			}
      //std::cout << "Position[j]=" << position[j] << std::endl;
			if(position[j]=="_pixel1"||position[j]=="_pixel2"||position[j]=="_pixel3"||position[j]=="_pixeldisk+1"||position[j]=="_pixeldisk+2"||position[j]=="_pixeldisk-1"||position[j]=="_pixeldisk-2")
			n_pixel++;
		}

		//if(n_pixel<2||n_pixel>3) continue;

		// Loop again to form the triplet strings excluding first and last hits
		for (int j=0; j<n; j++){

			if(stereo[j]==1) continue;
			i_next=0;
			i_prev=0;

      for(int k=0; k<n; k++){

				if(stereo[k]==1) continue;
				if(i_next*i_prev!=0) break;

				if(region[k]==region[j]){
					if(number[j]==number[k]+1) i_prev=k;
					else if(number[j]==number[k]-1) i_next=k;
				}

        //special for first layer
        if(position[j]=="_pixel1"){
          i_prev=500;
          if(position[k]=="_pixe2") i_next=k;
          //std::cout << "Pixel1=" << i_prev << "," << i_next << std::endl;
        }

				if(position[j]=="_TIB1"||position[j]=="_TIB2"||position[j]=="_TIB3"||position[j]=="_TIB4"){
					if(position[k]=="_TID+1"||position[k]=="_TID-1") i_next=k;
				}
				if(position[j]=="_pixeldisk+1"||position[j]=="_pixeldisk-1"){
					if(position[k]=="_pixel2"||position[k]=="_pixel1"||position[k]=="_pixel3") i_prev=k;
				}

				if(position[j]=="_TID+1"||position[j]=="_TID-1"){
					if(position[k]=="_TIB1"||position[k]=="_TIB2"||position[k]=="_TIB3") i_prev=k;
				}
				if(position[j]=="_TEC+1"||position[j]=="_TEC-1"){
					if(position[k]=="_TOB1"||position[k]=="_TOB2"||position[k]=="_TOB3"||position[k]=="_TOB4"||position[k]=="_TOB5") i_prev=k;
				}

				if(position[j]=="_pixeldisk+2"||position[j]=="_pixeldisk-2"){
					if(position[k]=="_TIB1") i_next=k;
				}
				if(position[j]=="_TID+2"||position[j]=="_TID-2"){
					if(position[k]=="_TOB1") i_next=k;
				}
				if(position[j]=="_pixeldisk+1"||position[j]=="_pixeldisk-1"){
					if(position[k]=="_TIB1") i_next=k;
				}
				if(position[j]=="_TID+1"||position[j]=="_TID-1"){
					if(position[k]=="_TOB1") i_next=k;
				}
			}

			if(i_next*i_prev==0) continue;
      if(i_prev == 500) {
        tripl_name="bSpot"+position[j]+position[i_next];
  			triplet=position[j];
      } else {
        tripl_name=position[i_prev]+position[j]+position[i_next];
  			triplet=position[j];
      }
      //std::cout << "Triplet is:" << tripl_name << std::endl;

			if(position[j]=="_pixel2") in_plots++;

			// A bunch of plots
			plot1D("pt", trackPt, 1, h_1d, 1000, 0.5, 2);
			plot1D("local pt"+triplet, pt[j], 1, h_1d, 1000, 0.5, 2);
			plot1D("eta", trackEta, 1, h_1d, 1000, -5, 5);
			plot1D("local eta"+triplet, eta[j], 1, h_1d, 160, -2.0, 2.0);

			// Correct for incidence angle
			if(region[j]=="_barrel") sin_track=calcSin(trackEta);
			else sin_track=TMath::Abs(calcCos(trackEta));
      float rprev = -99., phiprev;
      if(i_prev == 500) {
        ROOT::Math::XYZPoint bs(bsX0, bsY0, trackZ0);
        rprev = bs.R();
        phiprev = bs.Phi();
        //std::cout << "BS:" << rprev << "," << phiprev << std::endl;
      }
      else {
        rprev = r[i_prev];
        phiprev = phi[i_prev];
      }
			float b2=calcb2(rprev,r[j],r[i_next], phiprev,phi[j],phi[i_next])*sin_track;
			sag=calculSag(rprev,r[j],r[i_next], phiprev,phi[j],phi[i_next]);

			plot2D("b2vspt", trackPt*trackPt, sqrt(b2), 1, h_2d, 100, 0.5, 2.25, 500, 0., 50);

			TString sagmin_name="sagmin"+triplet;
			TString sagpl_name="sagpl"+triplet;

			if(!isFirstRun_){
				TH1 *sagmean_min=(TH1*)fileIN->Get(sagmin_name);
				TH1 *sagmean_pl=(TH1*)fileIN->Get(sagpl_name);

				if(!sagmean_min) continue;
				if(!sagmean_pl) continue;

				mean_min=sagmean_min->GetMean();
				mean_pl=sagmean_pl->GetMean();
			}
			if(triplet!=""){
				if(isFirstRun_){
					if(sag<0) plot1D("sagmin"+triplet, pt[j]*sag, 1, h_1dm, 600, -0.01, 0);
					else plot1D("sagpl"+triplet, pt[j]*sag, 1, h_1dpl, 600, 0, 0.01);
				}
				else{
					if(sag<0) plot1D("sagminus_new"+triplet, ((pt[j]*sag)-mean_min)*sqrt(b2), 1, h_1dm, 600, -0.1, 0.1);
					else plot1D("sagplus_new"+triplet, ((pt[j]*sag)-mean_pl)*sqrt(b2), 1, h_1dpl, 600, -0.1, 0.1);

					if(sag<0) plot3D("sag3Dminus"+triplet, eta[j], pt[j]*pt[j], ((pt[j]*sag)-mean_min)*sqrt(b2), 1, h_3dm, 41, -4, 4, 10, 0.5, 2.25, 600, -0.01, 0.01);
					else plot3D("sag3Dplus"+triplet, eta[j], pt[j]*pt[j], ((pt[j]*sag)-mean_pl)*sqrt(b2), 1, h_3dpl, 41, -4, 4, 10, 0.5, 2.25, 600, -0.01, 0.01);
				}

				triplet="";
				tripl_name="";
				sagmin_name="sagminus";
				sagpl_name="sagplus";
			}
		}
		for (int j=0; j<n; j++){
			position[j]="";
			r[j]=0;
			phi[j]=0;
			sag=0;
		}

	}
}
