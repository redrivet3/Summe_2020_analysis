#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <bitset>
#include <TMath.h>
#include "/w/hallb-scifs17exp/clas12/markov/dataQuality/clasqaDB/srcC/include/QADB.h"

#include "clas12reader.h"

using namespace clas12;

double pi = TMath::Pi();
int helicity;
double fCup;
vector<int> sectorE;
vector<double> p4_ele_px;
vector<double> p4_ele_py;
vector<double> p4_ele_pz;
vector<double> p4_ele_vx;
vector<double> p4_ele_vy;
vector<double> p4_ele_vz;
vector<double> p4_ele_E;
vector<double> p4_ele_theta;
vector<double> p4_ele_phi;
vector<double> p4_ele_mo;

// create historgrams for w, q^2, and w vs q^2
TH1F *wHisto = new TH1F("wHisto", "wHisto", 1000, 0, 5);
TH1F *q2Histo = new TH1F("q2Histo", "q2Histo", 1000, 0, 10);
TH2F *wq2Histo = new TH2F("wq2Histo", "wq2Histo",  1000, 0, 5,  1000, 0, 10);
// create histograms of theta angle vs phi angle and theta angle vs momentum
TH2F *thetaphihisto = new TH2F("thetaphihisto","thetaphihisto",1000,5,40,1000,-200,200);
TH2F *thetaphisto = new TH2F("thetaphisto", "thetaphisto", 1000, 5, 40, 1000,0,13);

TH2F *beamphiHisto = new TH2F("beamphiHisto", "beanphiHisto",500, 0, 360, 500, 0, 12);
TH2F *beamphiHisto2 = new TH2F("beamphiHisto2", "beanphiHisto2",500, 0, 360, 500, -10, 30);
// used to create a set of 6 2D histograms to be filled with each sector w value
TH1D *wSector[6];
TH2D *ebeam_sector[6];    
// calculates w from a given set of electrondata and beam energy given
double kin_W(TLorentzVector ele, float Ebeam){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  return fCM.M();
}


// calculates q^2 from a given set of electrondata and beam energy given
double kin_Q2(TLorentzVector ele, float Ebeam){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return -fGamma.M2();
}


// fills out the electron data into a set of vectors to be saved in the root file
void fill_output_vector_electron(TLorentzVector el){
  if(el.E() > 0){
    p4_ele_px.push_back(el.Px());
    p4_ele_py.push_back(el.Py());
    p4_ele_pz.push_back(el.Pz());
    p4_ele_E.push_back(el.E());
  }
}


// not sure what this is for, does not seemed to used in the example
void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	     rp->par()->getPz(),p4.M());
}



//main 
void simpleAnaLC2(){
  float eBeam = 10.646;   //pre determined beam energy
  float mProton = 0.98327; //given mass of the proton in GeV
  auto db=TDatabasePDG::Instance();

  // used to setup the 6 2D histograms to fill w from each sector
  for (int s=0;s<6;s++){
    wSector[s]=new TH1D(Form("wSectorS%d",s + 1),Form("wSectorS%d",s + 1), 300, .7 , 1.1 );
    ebeam_sector[s]=new TH2D(Form("ebeam_sectorS%d",s + 1),Form("ebeam_sectorS%d",s + 1),500, 0, 360, 500, 8, 12);
  }

  TString inputFile;
  TString outputFile = "test2.root"; // name of the resulting root file

  // checks for a hippo file?
  for(Int_t i=1;i<gApplication->Argc();i++){
    TString opt=gApplication->Argv(i);
    if((opt.Contains(".hipo"))){
      inputFile=opt(5,opt.Sizeof());
    }
  }

  // asks for a hippo file if one is nor given
  if(inputFile==TString())  {
    std::cout << " *** please provide a file name..." << std::endl;
    exit(0);
  }
  /////////////////////////////////////
  
  cout<<"Analysing hipo file "<<inputFile<<endl;

  TLorentzVector beam(0, 0, 10.646, 10.646); 
  TLorentzVector el;          //electron vector 
  TLorentzVector pro;         //proton vector
  TLorentzVector elmc;        //unused vector as far as I can tell
  TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());  
  TLorentzVector q;           //vector to be filled with q^2 

  // just a clock setup
  auto start = std::chrono::high_resolution_clock::now();
  gBenchmark->Start("timer");
  int counter=0;
  /////////////////////////////////////
  
  std::cout << " reading file example program (HIPO) "  << __cplusplus << std::endl;

  // initalizes the hippo reader
  hipo::reader  reader;
  cout<<"Analysing hipo file "<<inputFile<<endl;
  reader.open(inputFile);
  hipo::dictionary factory;
  reader.readDictionary(factory);
  TFile *out;

  // creates the output object to later write to the root file
  out = new TFile(outputFile, "RECREATE");
  TTree out_tree("out_tree","out_tree");
  out_tree.Branch("helicity", &helicity);
  out_tree.Branch("fCup", &fCup);
  
  //electrons data 
  out_tree.Branch("sectorE", &sectorE);
  out_tree.Branch("p4_ele_px", &p4_ele_px);
  out_tree.Branch("p4_ele_py", &p4_ele_py);
  out_tree.Branch("p4_ele_pz", &p4_ele_pz);
  out_tree.Branch("p4_ele_E", &p4_ele_E);

  out_tree.Branch("p4_ele_vx", &p4_ele_vx);
  out_tree.Branch("p4_ele_vy", &p4_ele_vy);
  out_tree.Branch("p4_ele_vz", &p4_ele_vz);

  out_tree.Branch("p4_ele_theta", &p4_ele_theta);
  out_tree.Branch("p4_ele_phi", &p4_ele_phi);
  out_tree.Branch("p4_ele_mo", &p4_ele_mo);

  //used to access the data in the hippo file
  hipo::event event;
  hipo::bank PART(factory.getSchema("REC::Particle"));
  hipo::bank PartCalorimeter(factory.getSchema("REC::Calorimeter"));
  hipo::bank EVENT(factory.getSchema("REC::Event"));
  hipo::bank Traj(factory.getSchema("REC::Traj"));
  hipo::bank RConfig(factory.getSchema("RUN::config"));
  // not neccesary for my work 
  QADB * qa = new QADB("/w/hallb-scifs17exp/clas12/markov/dataQuality/clasqaDB/srcC/qaTree.merged.json");
  int proSearch=2;
  // starts the main while loop to access data and manipulated as needed
  while(reader.next()==true){
    //clears out all of the vectors
    p4_ele_px.clear();
    p4_ele_py.clear();
    p4_ele_pz.clear();
    p4_ele_vx.clear();
    p4_ele_vy.clear();
    p4_ele_vz.clear();
    p4_ele_E.clear();
    sectorE.clear();
    p4_ele_theta.clear();
    p4_ele_phi.clear();
    p4_ele_mo.clear();

    counter++; // event counter
   
    // gets the 
    reader.read(event);
    event.getStructure(PART);
    event.getStructure(EVENT);
    event.getStructure(PartCalorimeter);
    event.getStructure(Traj);
    event.getStructure(RConfig);
    //filling out electron data
    int id = PART.getInt("pid",0);
    int charge = PART.getShort("charge",0);
    float startTime = EVENT.getFloat("startTime", 0);
    int status = PART.getInt("status", 0);
    int PCALSize = PartCalorimeter.getSize();
    int partSize = PART.getSize();
    int nEvent = RConfig.getInt("event", 0);
    int runNum = RConfig.getInt("run", 0);

    //prparing to fill out the proton data
    Int_t Np=0;
    int nRows = PART.getRows();
    int pid = PART.getInt("pid",proSearch);
    int chargepro = PART.getShort("charge",proSearch);
    int statuspro = PART.getInt("status", proSearch);
    bool run_pro = false;
    if (run_pro == true){
      for(int i =0; i<partSize;i++){
	//sets up proton data 
	proSearch=i;
	pid = PART.getInt("pid",proSearch);
	chargepro = PART.getShort("charge",proSearch);
	statuspro = PART.getInt("status", proSearch);
	if (pid == 2212 && statuspro<5000 && statuspro>999 && chargepro==1 &&  partSize > 0 && PART.getFloat("vz",0) > -10 && PART.getFloat("vz",0) < 5 && startTime > 0){
	  i=partSize+1;
	}
      }
    }
    // creates some useful variables 
    double eltheta;     //electron theta 
    double protheta;    //proton theta
    float elphi;        //electron phi
    float a0;           //used to calulate beam energy
    float a1;            //used to calulate beam energy
    float ebeam_cal;        //calculate beam energy
    float ebeam_cal_pro;
    float promass = db->GetParticle(2212)->Mass();
    float elmass = db->GetParticle(11)->Mass();
    //starts the check of the data, first 2 ifs not really neccessary 
    if(qa->Query(runNum,nEvent)) {
      if(qa->Golden()) {
	// runs through electron data
	if (id == 11 && status < -1999 && status > -4000 && partSize > 0 && charge == -1 && PART.getFloat("vz",0) > -10 && PART.getFloat("vz",0) < 5 && startTime > 0){
	  el.SetXYZM(PART.getFloat("px",0), PART.getFloat("py",0), PART.getFloat("pz",0), db->GetParticle(11)->Mass());
	  fill_output_vector_electron(el);
	  eltheta = el.Theta() * (180/pi);
	  elphi = el.Phi() * (180/pi);
	  // fill in vectors
	  sectorE.push_back(PartCalorimeter.getInt("sector", 0));                            
	  p4_ele_vx.push_back(PART.getFloat("vx",0));
	  p4_ele_vy.push_back(PART.getFloat("vy",0));
	  p4_ele_vz.push_back(PART.getFloat("vz",0));
	  p4_ele_theta.push_back(eltheta);
	  p4_ele_phi.push_back(elphi);
	  p4_ele_mo.push_back(el.P());
	  fCup = EVENT.getFloat("beamCharge", 0);
	  helicity = EVENT.getInt("helicity", 0);
	  // fills out the histograms
	  wHisto->Fill(kin_W(el, eBeam));
	  q2Histo->Fill(kin_Q2(el, eBeam));
	  wq2Histo->Fill(kin_W(el, eBeam), kin_Q2(el, eBeam));
	  thetaphihisto->Fill(eltheta, elphi);
	  thetaphisto->Fill(eltheta,el.P());

	  if(kin_W(el,eBeam) < 1.1 && kin_W(el,eBeam) > .7 ){
	    wSector[PartCalorimeter.getInt("sector",0)-1]->Fill(kin_W(el, eBeam));
	  }
	  // calculates the beam energy from electron data
       	  ebeam_cal =((promass*el.P())/(promass-el.P()+(el.P()*TMath::Cos(el.Theta()))));

	  if(el.P()>0){
	    if(kin_W(el,eBeam)<.97 && PartCalorimeter.getInt("sector",0) == 1){
	      beamphiHisto->Fill(elphi+180, ebeam_cal);    	
	      ebeam_sector[PartCalorimeter.getInt("sector",0)-1]->Fill(elphi+180,ebeam_cal);
	    }

	    if(kin_W(el,eBeam)<1.2 && PartCalorimeter.getInt("sector",0) == 2){
	      beamphiHisto->Fill(elphi+180, ebeam_cal);    	
	      ebeam_sector[PartCalorimeter.getInt("sector",0)-1]->Fill(elphi+180,ebeam_cal);
	    }

	    if(kin_W(el,eBeam)<1.2 && PartCalorimeter.getInt("sector",0) == 3){
	      beamphiHisto->Fill(elphi+180, ebeam_cal);    	
	      ebeam_sector[PartCalorimeter.getInt("sector",0)-1]->Fill(elphi+180,ebeam_cal);
	    }

	    if(kin_W(el,eBeam)<1.2 && PartCalorimeter.getInt("sector",0) == 4){
	      beamphiHisto->Fill(elphi+180, ebeam_cal);    	
	      ebeam_sector[PartCalorimeter.getInt("sector",0)-1]->Fill(elphi+180,ebeam_cal);
	    }

	    if(kin_W(el,eBeam)<1.2 && PartCalorimeter.getInt("sector",0) == 5){
	      beamphiHisto->Fill(elphi+180, ebeam_cal);    	
	      ebeam_sector[PartCalorimeter.getInt("sector",0)-1]->Fill(elphi+180,ebeam_cal);
	    }

	    if(kin_W(el,eBeam)<1.2 && PartCalorimeter.getInt("sector",0) == 6){
	      beamphiHisto->Fill(elphi+180, ebeam_cal);    	
	      ebeam_sector[PartCalorimeter.getInt("sector",0)-1]->Fill(elphi+180,ebeam_cal);
	    }
	  }




	  if(run_pro == true){
	    if (pid == 2212 && statuspro<5000 && statuspro>999 && chargepro==1 &&  partSize > 0 && PART.getFloat("vz",proSearch) > -10 && PART.getFloat("vz",proSearch) < 5 && startTime > 0){
	      pro.SetXYZM(PART.getFloat("px",proSearch), PART.getFloat("py",proSearch), PART.getFloat("pz",proSearch), db->GetParticle(11)->Mass());
	      a0 = 1 - 1/(TMath::Cos(el.Theta())-TMath::Sin(el.Theta())/TMath::Tan(-1 * pro.Theta()));	     
	      a1 = TMath::Sin(el.Theta())/TMath::Sin(el.Theta() + pro.Theta());
	      ebeam_cal_pro = 2 * (db->GetParticle(2212)->Mass()) * a0/(TMath::Power(a1,2) - TMath::Power(a0,2));
	    }
	    if(el.P()>0){
	      if(kin_W(el,eBeam)<1.2){
		beamphiHisto2->Fill(elphi+180, ebeam_cal_pro);    	
	      }
	    }
	  }
	}
 /////////////////////////////////////////////////////////////////////
      }
      out_tree.Fill();
    }
    else {
     
      cout <<"DEFECTO"<<endl;
      cout <<" file " << qa->GetFilenum() <<endl;
    }
    if (counter%10000 == 0){
      cout << counter <<endl;
    }
    //if(counter==10)break;
  }

  out->Write();
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");
  out->mkdir("Example");
  out->cd("Example");
  wHisto->Write();
  q2Histo->Write();
  wq2Histo->Write();
  for (int i=0;i<6;i++){
    wSector[i]->Write();
    ebeam_sector[i]->Write();
  }
  thetaphihisto->Write();
  thetaphisto->Write();
  beamphiHisto->Write();
  beamphiHisto2->Write();
  out->Close();
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";
}

