
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

//------------------------------------------------------------------------------


TH1F* Mll(const char *inputFile, int evNum, const char * histoname)
{
  TChain chain("Delphes");
  chain.Add(inputFile);
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  Int_t i, j;
  TLorentzVector vector;
  TLorentzVector vecLepton[1];
  TLorentzVector vecLepton1, vecLepton2, sum;
  Double_t muon_mass = 0.105658389;
  Double_t electron_mass = 0.000510998902;
  GenParticle *particle;
  Double_t delta_R_ll;
  // Book histograms
  //TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 200, 0.0, 1000.0);
  //TH1 *histJetPTBTag = new TH1F("jet_pt for b-tagged jets", "jet P_{T}", 200, 0.0, 1000.0);
  //TH1 *histNJet = new TH1F("N_jet", "N_jet", 20, 0.0, 20);
  TH1F *histInvMll = new TH1F(histoname, histoname, 200, 0, 10);

  // Event loop
  for(int entry = 0; entry < evNum; ++entry){  
  	if ( entry%1000 == 0){
  		cout << "progress " << entry << endl;
  	}
    treeReader->ReadEntry(entry);
    if( ((branchElectron->GetEntries() == 2) || (branchMuon->GetEntries() == 2)) || ((branchElectron->GetEntries() == 1) && (branchMuon->GetEntries() == 1)) ){    
      if(branchJet->GetEntries() > 0){
        Jet *jet = (Jet*) branchJet->At(0);
       // histJetPT->Fill(jet->PT);
        if (branchJet->GetEntries() == 2){
          for(i = 0; i < branchJet->GetEntriesFast(); ++i){
            if((jet->BTag & (1 << 0))){
            	//vector.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            	if ((branchElectron->GetEntries() == 1) && (branchMuon->GetEntries() == 1)){

            		for(j = 0; j < branchElectron->GetEntriesFast(); ++j){
            			Electron *electron = (Electron*) branchElectron->At(0);
            			particle = (GenParticle*) electron->Particle.GetObject();
            			vecLepton1.SetPtEtaPhiM(electron->PT, electron->Eta, electron->Phi, electron_mass);
            			//cout << particle->PID << endl; 
            			//cout << electron->PT << "  " << j << " " << electron->Eta << " " << electron->Phi <<  endl;
            		}
            		for(j = 0; j < branchMuon->GetEntriesFast(); ++j){
            			Muon *muon = (Muon*) branchMuon->At(0);
            			particle = (GenParticle*) muon->Particle.GetObject();
            			vecLepton2.SetPtEtaPhiM(muon->PT, muon->Eta, muon->Phi, muon_mass);
            			//cout << particle->PID << endl; 
            			//cout << muon->PT << "  " << j << " " << muon->Eta << " " << muon->Phi <<  endl;
            		}
            		sum = vecLepton1 + vecLepton2;
            		delta_R_ll = sqrt ((vecLepton1.Phi() - vecLepton2.Phi())*(vecLepton1.Phi() - vecLepton2.Phi()) + (vecLepton1.Eta() - vecLepton2.Eta())*(vecLepton1.Eta() - vecLepton2.Eta()) );
            		//delta_R_ll = sqrt( (electron->Eta - muon->Eta)*(electron->Eta - muon->Eta) + (electron->Phi - muon->Phi)*(electron->Phi - muon->Phi));
            		histInvMll->Fill(delta_R_ll);
            		//cout << sum.M() << endl;

            	}
            	/*if (branchElectron->GetEntries() == 2){
            		Electron *electron = (Electron*) branchElectron->At(0);
            		particle = (GenParticle*) electron->Particle.GetObject();
            		cout << "two electrons" << endl;
	            	for(j = 0; j < branchElectron->GetEntriesFast(); ++j){
						cout << particle->PID << endl; 
	            		cout << electron->PT << "  " << j << " " << electron->Eta << " " << electron->Phi <<  endl;

	            		//vecLepton[i].SetPtEtaPhiM(Electron.PT, Double_t eta, Double_t phi, Double_t e);
	            	}
            	//}*/

              //histogramm filling and variable selection happen here
             // histJetPTBTag->Fill(jet->PT);
            }

          }
        }
      }
    }
  }
 // histNJet->SetTitle("p p -> t t~ delphes level; N_jet ;");
 // histJetPT->SetTitle("p p -> t t~ delphes level; PT_jet GeV;");
 // histJetPTBTag->SetTitle("p p -> t t~ delphes level b-tagged events; PT_jet GeV;");
  //TCanvas *c1 = new TCanvas("c1","demo bin labels",10,10,900,500);
  //c1->Divide(2,1);
  //c1->cd(1);
  //histInvMll->Draw();
  return histInvMll;
  //c1->cd(2);
  //histJetPTBTag->Draw();
}

void double_higgs_analysis(const char *inputFile_background, const char *inputFile_signal)
{
	TFile *MyFile = new TFile("hh/hh_histogramms_r.root","recreate");
	gFile = MyFile;
	gSystem->Load("libDelphes");
	const char *histname1 = "Mll_background";
	const char *histname2 = "Mll_signal";
	TH1F *hist1, *hist2;
	hist1 = Mll(inputFile_background, 100000, histname1);
	hist2 = Mll(inputFile_signal, 10000, histname2);

 TCanvas *c1 = new TCanvas("c1","demo bin labels",10,10,900,500);
  //c1->Divide(2,1);
  //c1->cd(1);
  hist1->SetTitle("#M_ll; #M_ll [GeV];");
  hist2->SetTitle("#M_ll; #M_ll [GeV];");
  hist1->Draw();
  //c1->cd(2);
  hist2->SetFillStyle(3011);
  hist2->SetFillColor(kRed);
  hist2->Draw("SAMES");
 // hist1->Write("Mll_tt");
 // hist2->Write("Mll_hh");
  MyFile->Write();

}