#include <iostream>
#include <iomanip>
#include <locale>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "vector"
#include "vector"
using std::vector;
#include "algorithm"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

TTree          *fChain;   //!pointer to the analyzed TTree or TChain               
Int_t           fCurrent; //!current Tree number in a TChain                       


int iDBG=1;




int EMJselect(bool otfile, bool hasPre, const char* inputfilename,const char* outputfilename,
	      float HTcut, float pt1cut, float pt2cut, float pt3cut, float pt4cut, float jetacut,float alphaMaxcut, float maxIPcut, float NemfracCut,float CemfracCut,int ntrk1cut, int NemergingCut,bool blind) {
  // "ntuple.root", "histos.root"
  // suggest cuts 1000., 400.,200.,125.,50.,0.2,0.9,0.9,0,1
  // right now this code hard wires the jet pT cut and requires emerging jets to have at least
  // one track with pT> 1 GeV

  // read the Tree generated by tree1w and fill two histograms
  // note that we use "new" to create the TFile and TTree objects,
  // to keep them alive after leaving this function.
 
  int npass=0;

  TFile *f = new TFile(inputfilename);

  // get histogram of events before trigger
  TH1F *eventCountPreTrigger;

  if(hasPre) {
    if(otfile) eventCountPreTrigger = static_cast<TH1F*>(f->Get("eventCountPreTrigger/eventCountPreTrigger")->Clone());
  } else {
    if(otfile)  eventCountPreTrigger = new TH1F("eventCountPreTrigger","haha",2,0.,2.);
  }



  TTree *tt = (TTree*)f->Get("emJetAnalyzer/emJetTree");

  Int_t nVtx, event, lumi, run, nTrueInt, nTracks;
  Float_t met_pt, met_phi;

  //pv
  float pv_x,pv_y,pv_z;


  // gen particles
  vector<int> *gp_index=new vector<int>;
  vector<int> *gp_pdgId=new vector<int>;
  vector<float> *gp_pt = new vector<float>;
  vector<float> *gp_eta = new vector<float>;
  vector<float> *gp_phi = new vector<float>;

  // jet
  vector<int> *jet_index=new vector<int>;
  vector<int> *jet_source=new vector<int>;
  vector<float> *jet_pt = new vector<float>;
  vector<float> *jet_eta = new vector<float>;
  vector<float> *jet_phi = new vector<float>;
  vector<float> *jet_alphaMax = new vector<float>;
  vector<float> *jet_cef = new vector<float>;
  vector<float> *jet_nef = new vector<float>;
  vector<float> *jet_chf = new vector<float>;
  vector<float> *jet_nhf = new vector<float>;
  //  vector<float> *jet_phf = new vector<float>;
  vector<vector<float> > *track_pt = 0;
  vector<vector<float> > *track_eta = 0;
  vector<vector<int> > *track_source = 0;
  vector<vector<int> > *track_index = 0;
  vector<vector<int> > *track_jet_index = 0;
  vector<vector<int> > *track_vertex_index = 0;
  vector<vector<int> > *track_algo = 0;
  vector<vector<float> > *track_vertex_weight =0;
  vector<vector<float> > *track_ipZ =0;
  vector<vector<float> > *track_ipXY = 0;
  vector<vector<float> > *track_ipXYSig = 0;
  vector<vector<int> > *track_nMissInnerHits = 0;
  vector<vector<int> > *track_nMissInnerPxlLayers = 0;
  vector<vector<int> > *track_nPxlLayers = 0;
  vector<vector<int> > *track_nHits = 0;

  
  

  //get event count pre trigger



  //for ntuple
  // gen particles
  tt->SetBranchAddress("gp_index",&gp_index);
  tt->SetBranchAddress("gp_pdgId",&gp_pdgId);
  tt->SetBranchAddress("gp_pt",&gp_pt);
  tt->SetBranchAddress("gp_eta",&gp_eta);
  tt->SetBranchAddress("gp_phi",&gp_phi);


  // vertices
  tt->SetBranchAddress("nVtx",&nVtx);
  tt->SetBranchAddress("nTrueInt",&nTrueInt);
  tt->SetBranchAddress("nTracks",&nTracks);
  tt->SetBranchAddress("pv_x",&pv_x);
  tt->SetBranchAddress("pv_y",&pv_y);
  tt->SetBranchAddress("pv_z",&pv_z);

  // general event information
  tt->SetBranchAddress("event",&event);
  tt->SetBranchAddress("lumi",&lumi);
  tt->SetBranchAddress("run",&run);
  tt->SetBranchAddress("met_pt",&met_pt);
  tt->SetBranchAddress("met_phi",&met_phi);

  //jets
  tt->SetBranchAddress("jet_index",&jet_index);
  tt->SetBranchAddress("jet_source",&jet_source);
  tt->SetBranchAddress("jet_pt",&jet_pt);
  tt->SetBranchAddress("jet_eta",&jet_eta);
  tt->SetBranchAddress("jet_phi",&jet_phi);
  tt->SetBranchAddress("jet_cef",&jet_cef);
  tt->SetBranchAddress("jet_nef",&jet_nef);
  tt->SetBranchAddress("jet_chf",&jet_chf);
  tt->SetBranchAddress("jet_nhf",&jet_nhf);
  //  tt->SetBranchAddress("jet_phf",&jet_phf);
  tt->SetBranchAddress("jet_alphaMax",&jet_alphaMax);
  tt->SetBranchAddress("track_pt",&track_pt);
  tt->SetBranchAddress("track_eta",&track_eta);
  tt->SetBranchAddress("track_source",&track_source);
  tt->SetBranchAddress("track_index",&track_index);
  tt->SetBranchAddress("track_jet_index",&track_jet_index);
  tt->SetBranchAddress("track_algo",&track_algo);
  tt->SetBranchAddress("track_vertex_index",&track_vertex_index);
  tt->SetBranchAddress("track_vertex_weight",&track_vertex_weight);
  tt->SetBranchAddress("track_ipXY",&track_ipXY);
  tt->SetBranchAddress("track_ipXYSig",&track_ipXYSig);
  tt->SetBranchAddress("track_nMissInnerHits",&track_nMissInnerHits);
  tt->SetBranchAddress("track_nMissInnerPxlLayers",&track_nMissInnerPxlLayers);
  tt->SetBranchAddress("track_nPxlLayers",&track_nPxlLayers);
  tt->SetBranchAddress("track_nHits",&track_nHits);
  tt->SetBranchAddress("track_ipZ",&track_ipZ);

  

  

  // create a histograms
  TH1F *acount,*count,*hjetcut,*hjetchf,*h_nemg,*hnjet,*hpt,*heta,*heta2,*halpha,*H_T,*H_T2,*H_T3,*H_T4,*hbcut_ntrkpt1,*hacut_ntrkpt1,*hbcut_nef,*hacut_nef,*hbcut_cef,*hacut_cef,*hbcut_alphamax,*hacut_alphamax,*hHTnm1,*hnHitsnm1,*hntrk1nm1,*hmaxipnm1,*hpt1nm1,*hpt2nm1,*hpt3nm1,*hpt4nm1,*halphanm1,*hnemnm1,*hpt1,*hpt2,*hpt3,*hpt4,*hipXYEJ,*hipXYnEJ,*htvw,*htvwEJ,*hnmaxipnm1,*hn2maxipnm1,*hjptfrb,*hjptfra1,*hjptfra2,*hjptfrbc,*hjptfra1c,*hjptfra2c,*hjptb,*hjpta,*haMgj,*hHTko,*hpt1ko,*hpt2ko,*hpt3ko,*hpt4ko,*hipXYSigEJ,*hipXYSignEJ,*hmaxipXYEJ,*hmaxipXYnEJ,*hmeanipXYEJ,*hmeanipXYnEJ,*hmass;

  TH2F *aMip,*haMvjpt,*haMvHT,*haMvnvtx;

  if(otfile) {
  acount = new TH1F("acount","counts",20,0.,20.);
  count = new TH1F("count","counts",3,0,3);
  count->SetStats(0);
  count->SetCanExtend(TH1::kAllAxes);

  // 1d
  hjetcut = new TH1F("hjetcut","jetcut counts",20,0.,20.);
  hjetchf = new TH1F("hjetchf","jet charged hadron fr",20,0.,1.2);
  h_nemg = new TH1F("h_nemg","number of emerging jets",20,0.,20.);
  hnjet = new TH1F("hnjet","number of jets",20,0.,20.);
  hpt = new TH1F("hpt","jet pt distribution",200,0.,1000.);
  heta   = new TH1F("heta","jet eta distribution",100,-4.,4.);
  heta2   = new TH1F("heta2","jet eta distribution first 4 jets",100,-4.,4.);
  halpha   = new TH1F("halpha","jet alphaMax distribution",100,0.,1.5);
  haMgj   = new TH1F("haMgj","jet alphaMax distribution, good jets",100,0.,1.5);
  H_T      = new TH1F("H_T"," HT distribution before cut", 100,0.,5000.);
  H_T2      = new TH1F("H_T2"," HT distribution after cut", 100,0.,5000.);
  H_T3      = new TH1F("H_T3"," HT distribution at end", 100,0.,5000.);
  H_T4      = new TH1F("H_T4"," HT distribution test", 100,0.,5000.);
  hpt1 = new TH1F("hpt1"," pT of leading jet",200,0.,1000.);
  hpt2 = new TH1F("hpt2"," pT of second jet",200,0.,1000.);
  hpt3 = new TH1F("hpt3"," pT of third jet",200,0.,1000.);
  hpt4 = new TH1F("hpt4"," pT of fourth jet",200,0.,1000.);
  hbcut_ntrkpt1 = new TH1F("hbcut_ntrkpt1","number tracks pt>1 before cut",20,0.,20.);
  hacut_ntrkpt1 = new TH1F("hacut_ntrkpt1","number tracks pt>1 after cut",20,0.,20.);
  hbcut_nef = new TH1F("hbcut_nef","neutral em fraction before cut",10,0.,1.2);
  hacut_nef = new TH1F("hacut_nef","neutral em fraction after cut",10,0.,1.2);
  hbcut_cef = new TH1F("hbcut_cef","charged em fraction before cut",50,0.,1.2);
  hacut_cef = new TH1F("hacut_cef","charged em fraction after cut",50,0.,1.2);
  hbcut_alphamax = new TH1F("hbcut_alphamax","alphamax before cut",50,0.,1.5);
  hacut_alphamax = new TH1F("hacut_alphamax","alphamax after cut",50,0.,1.5);
  hHTnm1 = new TH1F("hHTnm1","HT n-1",100,0.,5000.);
  hpt1nm1 = new TH1F("hpt1nm1","pt1 n-1",200,0.,1000.);
  hpt2nm1 = new TH1F("hpt2nm1","pt2 n-1",200,0.,1000.);
  hpt3nm1 = new TH1F("hpt3nm1","pt3 n-1",200,0.,1000.);
  hpt4nm1 = new TH1F("hpt4nm1","pt4 n-1",200,0.,1000.);
  halphanm1 = new TH1F("halphanm1","alpha max n-1",200,0.,1.5);
  hmaxipnm1 = new TH1F("hmaxipnm1","ip max n-1",200,0.,10.);
  hnmaxipnm1 = new TH1F("hnmaxipnm1","new 2 ip max n-1",200,0.,10.);
  hn2maxipnm1 = new TH1F("hn2maxipnm1","new 1  ip max n-1",200,0.,10.);
  hnHitsnm1 = new TH1F("hnHitsnm1","number Hits n-1",40,0.,40.);
  hntrk1nm1 = new TH1F("hntrk1nm1","number tracks pt>1 n-1",50,0.,50.);
  hnemnm1 = new TH1F("hnemnm1","N emerging jets n-1",10,0.,10.);
  hipXYEJ = new TH1F("hipXYEJ","impact parameter  tracks of emerging jets",300,-1.,1.);
  hipXYnEJ = new TH1F("hipXYnEJ","impact parameter  tracks of not emerging jets",300,-1.,1.);
  htvw = new TH1F("htvw","track vertex weight ",15,-5.,10.);
  htvwEJ= new TH1F("htvwEJ","track vertex weight Emerging Jets ",15,-5.,10.);
  hipXYSigEJ = new TH1F("hipXYSigEJ","ip sig emerging jets",100,0.,10.);
  hipXYSignEJ = new TH1F("hipXYSignEJ","ip sig not emerging jets",100,0.,10.);
  hmaxipXYEJ = new TH1F("hmaxipXYEJ","max ip emerging jets",1000,0.,10.);
  hmaxipXYnEJ = new TH1F("hmaxipXYnEJ","max ip not emerging jets",1000,0.,10.);
  hmeanipXYEJ = new TH1F("hmeanipXYEJ","mean ip emerging jets",1000,0.,2.);
  hmeanipXYnEJ = new TH1F("hmeanipXYnEJ","mean ip not emerging jets",1000,0.,2.);
  hjptb = new TH1F("hjptb"," pT of basic jet",100,0.,1000.);
  hjpta = new TH1F("hjpta"," pT of emergng jets",100,0.,1000.);
  hjptfrb = new TH1F("hjptfrb"," pT of basic jets passing kine selection and n<4",100,0.,1000.);
  hjptfra1 = new TH1F("hjptfra1"," pT of basic jets passing kine, almost selection and n<4",100,0.,1000.);
  hjptfra2 = new TH1F("hjptfra2"," pT of basic jets passing kine, almost, and emerging selection and n<4",100,0.,1000.);
  hjptfrbc = new TH1F("hjptfrbc"," pT of basic jets passing kine selection",100,0.,1000.);
  hjptfra1c = new TH1F("hjptfra1c"," pT of basic jets passing kine, almost selection",100,0.,1000.);
  hjptfra2c = new TH1F("hjptfra2c"," pT of basic jets passing kine, almost, and emerging selection",100,0.,1000.);
  hHTko      = new TH1F("hHTko"," HT distribution test kine cuts", 100,0.,5000.);
  hpt1ko = new TH1F("hpt1ko"," pT of leading jet kine cuts",200,0.,1000.);
  hpt2ko = new TH1F("hpt2ko"," pT of second jet kine cuts",200,0.,1000.);
  hpt3ko = new TH1F("hpt3ko"," pT of third jet kine cuts",200,0.,1000.);
  hpt4ko = new TH1F("hpt4ko"," pT of fourth jet kine cuts",200,0.,1000.);
  hmass = new TH1F("hmass","mass emerging and non",500,0.,5000.);


  //2d
  aMip = new TH2F("aMip"," alpha Max versus max IP n-1 plot",100,0.,1.,100,0.,10.);
  haMvjpt = new TH2F("haMvjpt"," alpha Max versus jet pT ",100,0.,1.,100,0.,700.);
  haMvHT = new TH2F("haMvHT"," alpha Max versus HT ",100,0.,1.,100,0.,2500.);
  haMvnvtx = new TH2F("haMvnvtx"," alpha Max versus nvtx ",40,0.,1.,100,0.,40.);
  }

  //read all entries and fill the histograms
  Int_t nentries = (Int_t)tt->GetEntries();


  // loop over events
  for (Int_t i=0; i<nentries; i++) {
    if(iDBG>0) std::cout<<"***event "<<event<<std::endl;
 
    if(!hasPre) eventCountPreTrigger->Fill(1.5); 
    
    if(otfile) count->Fill("All",1);  // count number of events
    if(otfile) acount->Fill(0.5);
    tt->GetEntry(i);
    if(iDBG>0) {
      std::cout<<"event number is "<<event<<" number of vertex is "<<nVtx<<std::endl;
      std::cout<<"pv position is "<<pv_x<<","<<pv_y<<","<<pv_z<<std::endl;
    }

    // make some basic plots on all events before any selections




    // gen particles
    // especially find the first dark quark and dark anti quark
    // assume the following particle is d or dbar
    int firstdkq=0;
    int firstadkq=0;
    int firstdq=0;
    int firstadq=0;

    int NNNgp = (*gp_index).size();
    if(iDBG>0) std::cout<<" gen particle id pt eta phi"<<std::endl;
    for(Int_t j=1; j<NNNgp-1; j++) {
      //	if(iDBG>0) std::cout<<"    "<<(*gp_pdgId)[j]<<" "<<(*gp_pt)[j]<<" "<<(*gp_eta)[j]<<" "<<(*gp_phi)[j]<<std::endl;

  
      if(((*gp_pdgId)[j]==4900101)&&(firstdkq==0)
           &&( ((*gp_pdgId)[j+1]==1)||((*gp_pdgId)[j-1]==1)) ) {
	firstdkq=j;
	firstdq=j+1;
	if((*gp_pdgId)[j+1]!=1) firstdq=j-1;
	if(iDBG>0) std::cout<<"    match "<<(*gp_pdgId)[firstdkq]<<" "<<(*gp_pt)[firstdkq]<<" "<<(*gp_eta)[firstdkq]<<" "<<(*gp_phi)[firstdkq]<<std::endl;
	if(iDBG>0) std::cout<<"    match "<<(*gp_pdgId)[firstdq]<<" "<<(*gp_pt)[firstdq]<<" "<<(*gp_eta)[firstdq]<<" "<<(*gp_phi)[firstdq]<<std::endl;
      }
      if(((*gp_pdgId)[j]==-4900101)&&(firstadkq==0)
	 &&( ((*gp_pdgId)[j+1]==-1)||((*gp_pdgId)[j-1]==-1)) ) {
	firstadkq=j;
	firstadq=j+1;
	if((*gp_pdgId)[j+1]!=-1) firstadq=j-1;
	if(iDBG>0) std::cout<<"    match "<<(*gp_pdgId)[firstadkq]<<" "<<(*gp_pt)[firstadkq]<<" "<<(*gp_eta)[firstadkq]<<" "<<(*gp_phi)[firstadkq]<<std::endl;
	if(iDBG>0) std::cout<<"    match "<<(*gp_pdgId)[firstadq]<<" "<<(*gp_pt)[firstadq]<<" "<<(*gp_eta)[firstadq]<<" "<<(*gp_phi)[firstadq]<<std::endl;

      }
    }
    if(iDBG>0) std::cout<<std::endl<<std::endl;
    if(firstdq==0) std::cout<<" first dark quark not found"<<std::endl;
    if(firstdq==0) std::cout<<" first down quark not found"<<std::endl;
    if(firstdq==0) std::cout<<" first anti dark quark not found"<<std::endl;
    if(firstdq==0) std::cout<<" first anti down quark not found"<<std::endl;



    // jets
    vector<int> jet_ntrkpt1((*jet_index).size());
    vector<float> jet_meanip((*jet_index).size());
    vector<float> r0((*jet_index).size());
    vector<float> r1((*jet_index).size());
    vector<int> jntrack((*jet_index).size());
    vector<float> jet_e((*jet_index).size());
    vector<float> jet_theta((*jet_index).size());
    vector<float> jet_px((*jet_index).size());
    vector<float> jet_py((*jet_index).size());
    vector<float> jet_pz((*jet_index).size());
    if(otfile) hnjet->Fill((*jet_index).size()+0.5);
    int NNNjet = (*jet_index).size();
    for(Int_t j=0; j<NNNjet; j++) {
      if(iDBG>0) std::cout<<"jet j = "<<j<<std::endl;
      jet_theta[j]=2.*atan(exp(-(*jet_eta)[j]));
      jet_e[j]=(*jet_pt)[j]/sin(jet_theta[j]);
      jet_px[j]=(*jet_pt)[j]*cos((*jet_phi)[j]);
      jet_py[j]=(*jet_pt)[j]*sin((*jet_phi)[j]);
      jet_pz[j]=(*jet_pt)[j]/tan(jet_theta[j]);
				
      if(otfile) hpt->Fill((*jet_pt)[j]);
      if(otfile) heta->Fill((*jet_eta)[j]);
      if(otfile) hjetchf->Fill((*jet_chf)[j]);
      if(otfile) if(j<4) heta2->Fill((*jet_eta)[j]);
      if(otfile) halpha->Fill((*jet_alphaMax)[j]);
      //      calculate  number of tracks with pt > 1
      jet_ntrkpt1[j]=0;
      jet_meanip[j]=0.;
      if(r0.size()>0) r0[j]=0.;
      if(r1.size()>0) r1[j]=0.;
      vector<float> track_pts = track_pt->at(j);
      vector<int> track_sources = track_source->at(j);
      vector<float> track_vertex_weights = track_vertex_weight->at(j);
      vector<float> track_ipXYs = track_ipXY->at(j);
      vector<float> track_ipXYSigs = track_ipXYSig->at(j);
      vector<float> sort_ip(track_pts.size());
      for(int it=0;it<track_pts.size();it++) sort_ip[it]=0;
      jntrack[j]=0;
      for (unsigned itrack=0; itrack<track_pts.size(); itrack++) {
	if(track_sources[itrack]==0) {
	  sort_ip[jntrack[j]]=fabs(track_ipXYs[itrack]);
	  if(otfile) htvw->Fill(track_vertex_weights[itrack]);
	  if(iDBG>0) std::cout<<"track vertex weight is "<<track_vertex_weights[itrack]<<std::endl;
	  if(track_pts[itrack]>1) jet_ntrkpt1[j]+=1;
	  if(iDBG>0) std::cout<<" track "<<itrack<<" ip "<<track_ipXYs[itrack]<<" mean ip "<<jet_meanip[j]<<std::endl;
	  jet_meanip[j]=jet_meanip[j]+fabs(track_ipXYs[itrack]);
	  jntrack[j]++;
	}
      }
      float atmp = jntrack[j];
      if(jntrack[j]>0) jet_meanip[j]=jet_meanip[j]/atmp;
      std::sort(sort_ip.begin(), sort_ip.end());
      std::reverse(sort_ip.begin(),sort_ip.end());
      if(sort_ip.size()>0) r0[j]=sort_ip[0];
      if(sort_ip.size()>1) r1[j]=sort_ip[1];
      if(iDBG>0) std::cout<<"mean max are "<<jet_meanip[j]<<" "<<r0[j]<<std::endl;
     }  // end of loop over jets



      //now see which jets are emerging
    if(iDBG>0) std::cout<<" in event "<<event<<" number of jets is "<<NNNjet<<std::endl;
    vector<bool> emerging(NNNjet);
    vector<bool> almostemerging(NNNjet);
    vector<bool> basicjet(NNNjet);
      for( int i=0;i<4;i++) {
	  emerging[i]=false;
	  almostemerging[i]=false;
	  basicjet[i]=false;
	}
      int nemerging=0;
      int nalmostemerging=0;
      int iijjkk = 4;
      //if(NNNjet<4) iijjkk=NNNjet;
      //      std::cout<<"iijjkk is "<<iijjkk<<std::endl;
      for(int ij=0;ij<NNNjet;ij++) {
	
        vector<float> track_ipXYs = track_ipXY->at(ij);
        vector<float> track_ipXYSigs = track_ipXYSig->at(ij);
        vector<int> track_sources = track_source->at(ij);
        vector<float> track_vertex_weights = track_vertex_weight->at(ij);
	if(otfile) hjetcut->Fill(0.5);

	if(fabs((*jet_eta)[ij])<jetacut) { // jet eta cut
	    if(otfile) hjetcut->Fill(1.5);

	if(otfile) hbcut_nef->Fill((*jet_nef)[ij]);
	if((*jet_nef)[ij]<NemfracCut) {  // neutral fraction
	    if(otfile) hacut_nef->Fill((*jet_nef)[ij]);
	    if(otfile) hjetcut->Fill(2.5);

	    if(otfile) hbcut_ntrkpt1->Fill(jet_ntrkpt1[ij]);
	    if(jet_ntrkpt1[ij]>ntrk1cut) {  // tracks pt>1
	      if(otfile) hacut_ntrkpt1->Fill(jet_ntrkpt1[ij]);
	      if(otfile) hjetcut->Fill(3.5);

	      if(otfile) hbcut_cef->Fill((*jet_cef)[ij]);
	      if((*jet_cef)[ij]<CemfracCut) {  //charged fraction
	        if(otfile) hacut_cef->Fill((*jet_cef)[ij]);
	        if(otfile) hjetcut->Fill(4.5);
		basicjet[ij]=true;

	        if(otfile) hbcut_alphamax->Fill((*jet_alphaMax)[ij]);
	        if((*jet_alphaMax)[ij]<alphaMaxcut) { // alpha max
	          if(otfile) hacut_alphamax->Fill((*jet_alphaMax)[ij]);
	          if(otfile) hjetcut->Fill(5.5);
		  almostemerging[ij]=true;
		  
		  if(ij<4) nalmostemerging=nalmostemerging+1;
		  if(iDBG>0) {
		if(ij<4) {
		  std::cout<<" an almost emerging jet "<<ij<<std::endl;
		  std::cout<<" with r0 of "<<r0[ij]<<std::endl;
		  std::cout<<" and pt of "<<(*jet_pt)[ij]<<std::endl;
		}
		  }
		if(r0[ij]>maxIPcut) { // max IP cut

	        emerging[ij]=true;
	        if(ij<4) nemerging+=1.;
		if(iDBG>0) {
		if(ij<4) {
		  std::cout<<" an emerging jet "<<ij<<std::endl;
		  std::cout<<" with r0 of "<<r0[ij]<<std::endl;
		  std::cout<<" and pt of "<<(*jet_pt)[ij]<<std::endl;
		}
		}
		// look at tracks in the emerging jets
		if(otfile) hmaxipXYEJ->Fill(r0[ij]);
		if(otfile) hmeanipXYEJ->Fill(jet_meanip[ij]);
		if(jet_meanip[ij]>r0[ij]) std::cout<<"DANGER DANGER"<<std::endl;
                for (unsigned itrack=0; itrack<track_ipXYs.size(); itrack++) {
	          if(track_sources[itrack]==0) {
		    if(otfile) hipXYEJ->Fill(track_ipXYs[itrack]);
		    if(otfile) hipXYSigEJ->Fill(track_ipXYSigs[itrack]);
		    if(otfile) htvwEJ->Fill(track_vertex_weights[itrack]);
	           }
                }
		}
	      }
	    }
	  }
        }}
	if(!emerging[ij]) {
	  if(otfile) hmaxipXYnEJ->Fill(r0[ij]);
	  if(otfile) hmeanipXYnEJ->Fill(jet_meanip[ij]);
                for (unsigned itrack=0; itrack<track_ipXYs.size(); itrack++) {
	          if(track_sources[itrack]==0) {
		    if(otfile) hipXYnEJ->Fill(track_ipXYs[itrack]);
		    if(otfile) hipXYSignEJ->Fill(track_ipXYSigs[itrack]);
	           }
                }

	}
	if(iDBG>0) std::cout<<"event pt alphaM cef nef ntrkpt1 r0 emerging  almost "<<event<<" "<<(*jet_pt)[ij]<<" "<<(*jet_alphaMax)[ij]<<" "<<(*jet_cef)[ij]<<" "<<(*jet_nef)[ij]<<" "<<jet_ntrkpt1[ij]<<" "<<r0[ij]<<" "<<emerging[ij]<<" "<<almostemerging[ij]<<std::endl;
      }
      if(otfile) h_nemg->Fill(nemerging);





      // *************************************************************
    // now start the event selections
      // *************************************************************

    // require at least 4 jets
    bool C4jet=true;
    if(NNNjet<3) C4jet=false;
    // HT
    double HT = (*jet_pt)[0]+(*jet_pt)[1]+(*jet_pt)[2]+(*jet_pt)[3];
    if(otfile) H_T->Fill(HT);
    if(otfile) hpt1->Fill((*jet_pt)[0]);
    if(otfile) hpt2->Fill((*jet_pt)[1]);
    if(otfile) hpt3->Fill((*jet_pt)[2]);
    if(otfile) hpt4->Fill((*jet_pt)[3]);
    bool CHT=true;
    if(HT<HTcut) CHT=false;
    // jet pt
    bool Cpt1=false;
    bool Cpt2=false;
    bool Cpt3=false;
    bool Cpt4=false;
    if(((*jet_pt)[0]>pt1cut)&&(fabs((*jet_eta)[0])<jetacut)) Cpt1=true;
    if(((*jet_pt)[1]>pt2cut)&&(fabs((*jet_eta)[1])<jetacut)) Cpt2=true;
    if(((*jet_pt)[2]>pt3cut)&&(fabs((*jet_eta)[2])<jetacut)) Cpt3=true;
    if(((*jet_pt)[3]>pt4cut)&&(fabs((*jet_eta)[3])<jetacut)) Cpt4=true;
    // number emerging jets
    bool Cnem = true;
    if(nemerging<NemergingCut) Cnem=false;

    //    if(nalmostemerging>=4) Cnem=false;
    bool Canem =true;
    if(nalmostemerging>=4) Canem=false;


    //blind
    if(blind) {
      Cnem=false;
      Canem=false;
    }

    // do some plots
    if(otfile) {

      // kine only plots
      if(C4jet&&CHT&&Cpt1&&Cpt2&&Cpt3&&Cpt4) {
	hHTko->Fill(HT);
	hpt1ko->Fill((*jet_pt)[0]);
	hpt2ko->Fill((*jet_pt)[1]);
	hpt3ko->Fill((*jet_pt)[2]);
	hpt4ko->Fill((*jet_pt)[3]);
      }

      // jet plots
      for(int i=0;i<NNNjet;i++) {
	if(basicjet[i]) {
	  if((*jet_pt)[i]>50 ) {
	    haMgj->Fill((*jet_alphaMax)[i]);
	    if(iDBG>0) {
	    if((*jet_alphaMax)[i]<0.015) {
	      std::cout<<"CHECK"<<std::endl;
	      std::cout<<"alpha max is "<<(*jet_alphaMax)[i]<<std::endl;
	      std::cout<<"jet pt is "<<(*jet_pt)[i]<<std::endl;  
	      std::cout<<"jet eta is "<<(*jet_eta)[i]<<std::endl;  
	      std::cout<<"jet phi is "<<(*jet_phi)[i]<<std::endl;  
	      std::cout<<" number tracks in jet is "<<jntrack[i]<<std::endl;
	      std::cout<<" number true interactions is "<<nTrueInt<<std::endl;
	      std::cout<<" total number tracks is "<<nTracks<<std::endl;
	    }
	    }

	    haMvjpt->Fill((*jet_alphaMax)[i],(*jet_pt)[i]);
	    haMvHT->Fill((*jet_alphaMax)[i],HT);
	    haMvnvtx->Fill((*jet_alphaMax)[i],nVtx);
	    hjptb->Fill((*jet_pt)[i]);
	    if(emerging[i]) {
	      hjpta->Fill((*jet_pt)[i]);
	    }
	  }}
      }

      //N-1 plots
    if(C4jet&&Cpt1&&Cpt2&&Cpt3&&Cpt4&&Cnem&&Canem) hHTnm1->Fill(HT);
    if(C4jet&&CHT&&Cpt2&&Cpt3&&Cpt4&&Cnem&&Canem) hpt1nm1->Fill((*jet_pt)[0]);
    if(C4jet&&CHT&&Cpt1&&Cpt3&&Cpt4&&Cnem&&Canem) hpt2nm1->Fill((*jet_pt)[1]);
    if(C4jet&&CHT&&Cpt1&&Cpt2&&Cpt4&&Cnem&&Canem) hpt3nm1->Fill((*jet_pt)[2]);
    if(C4jet&&CHT&&Cpt1&&Cpt2&&Cpt3&&Cnem&&Canem) hpt4nm1->Fill((*jet_pt)[3]);
    if(C4jet&&CHT&&Cpt1&&Cpt2&&Cpt3&&Cpt4&&Canem) hnemnm1->Fill(nemerging);
    if(C4jet&&CHT&&Cpt1&&Cpt2&&Cpt3&&Cpt4&&Canem) {
      for(int i=0;i<3;i++) {
	if(basicjet[i]) {
	  halphanm1->Fill((*jet_alphaMax)[i]);
	  aMip->Fill((*jet_alphaMax)[i],r0[i]);
	  hntrk1nm1->Fill(jet_ntrkpt1[i]);
	  if(((*jet_alphaMax)[i]<alphaMaxcut)) {
	    hmaxipnm1->Fill(r0[i]);

	    if(iDBG>0) {
            std::cout<<" almost emerging"<<std::endl;
	    if(r0[i]<0.05) std::cout<<"DANGER DANGER"<<std::endl;
            std::cout<<" jet pt is "<<(*jet_pt)[i]
	    	     <<" ntrkpt1 is "<<jet_ntrkpt1[i]
	    	     <<" meanip is "<<jet_meanip[i]
	    	     <<" ip max is "<<r0[i]
	    	     <<" second largest ip is "<<r1[i]
	    	     <<" alpha max is "<<(*jet_alphaMax)[i]
            <<std::endl;
	    }
            vector<float> track_pts = track_pt->at(i);
            vector<int> track_sources = track_source->at(i);
            vector<float> track_vertex_weights = track_vertex_weight->at(i);
            vector<float> track_ipXYs = track_ipXY->at(i);
            vector<float> track_ipXYSigs = track_ipXYSig->at(i);
            vector<int> track_nMissInnerHitss = track_nMissInnerHits->at(i);
            vector<int> track_nMissInnerPxlLayerss = track_nMissInnerPxlLayers->at(i);
            vector<int> track_nPxlLayerss = track_nPxlLayers->at(i);
            vector<int> track_nHitss = track_nHits->at(i);
            vector<float> track_ipZs = track_ipZ->at(i);
            for (unsigned itrack=0; itrack<track_pts.size(); itrack++) {
	      if(track_sources[itrack]==0) {
		if(iDBG>0) {
			std::cout<<"    track pt is "<<track_pts[itrack]
			 <<" ipxy is "<<track_ipXYs[itrack]
			 <<" ipxysig is "<<track_ipXYSigs[itrack]
			 <<" ipZ is "<<track_ipZs[itrack]
			 <<" missinnerhits is "<<track_nMissInnerHitss[itrack]
			 <<" missinnerpxllayers is "<<track_nMissInnerPxlLayerss[itrack]
			 <<" pxllayers is "<<track_nPxlLayerss[itrack]
			 <<" nHits is "<<track_nHitss[itrack]
			 <<std::endl;
		}
		if(otfile) hnHitsnm1->Fill(track_nHitss[itrack]);
	      }
            }
	    


	  }
	}
      }
    }


    if(C4jet&&CHT&&Cpt1&&Cpt2&&Cpt3&&Cpt4&&nalmostemerging>=2) {
    if(otfile) H_T4->Fill(HT);
      for(int i=0;i<3;i++) {
	if(almostemerging[i]) {
	  if(((*jet_alphaMax)[i]<alphaMaxcut)) {
	    hnmaxipnm1->Fill(r0[i]);
	  }
	}
      }
    }

    if(C4jet&&CHT&&Cpt1&&Cpt2&&Cpt3&&Cpt4&&nalmostemerging>=1) {
      for(int i=0;i<3;i++) {
	if(almostemerging[i]) {
	  if(((*jet_alphaMax)[i]<alphaMaxcut)) {
	    hn2maxipnm1->Fill(r0[i]);
	  }
	}
      }
    }



    }


      // make plots for fake rate studes
    if(C4jet&&CHT&&Cpt1&&Cpt2&&Cpt3&&Cpt4&&Canem) {
      for(Int_t j=0; j<NNNjet; j++) {
	if(basicjet[j]) {
	  hjptfrb->Fill((*jet_pt)[j]);
	  if(almostemerging[j]){
	    hjptfra1->Fill((*jet_pt)[j]);
	    if(emerging[j]) {
	      hjptfra2->Fill((*jet_pt)[j]);
	    }
	  }
	}
      }
    }
      // check without Canem
    if(C4jet&&CHT&&Cpt1&&Cpt2&&Cpt3&&Cpt4) {
      for(Int_t j=0; j<NNNjet; j++) {
	if(basicjet[j]) {
	  hjptfrbc->Fill((*jet_pt)[j]);
	  if(almostemerging[j]){
	    hjptfra1c->Fill((*jet_pt)[j]);
	    if(emerging[j]) {
	      hjptfra2c->Fill((*jet_pt)[j]);
	    }
	  }
	}
      }
    }


    // apply cuts sequentially

    if(iDBG>0) std::cout<<"c4jet cht cpt1 cpt2 cpt3 cpt4 cnem "<<C4jet<<" "<<CHT<<" "<<Cpt1<<" "<<Cpt2<<" "<<Cpt3<<" "<<Cpt4<<" "<<Cnem<<std::endl;

    if(C4jet) {
    if(otfile) count->Fill("4 jets",1);
    if(otfile) acount->Fill(1.5);

    // calculate HT and require it greater than some cut value
    if(CHT) {
    if(otfile) count->Fill("HT",1);
    if(otfile) acount->Fill(2.5);
    if(otfile) H_T2->Fill(HT);

    // do pT cuts on jets  
    if(Cpt1) {
    if(otfile) count->Fill("jet pt1",1);
    if(otfile) acount->Fill(3.5);


    if(Cpt2) {
    if(otfile) count->Fill("jet pt2",1);
    if(otfile) acount->Fill(4.5);


    if(Cpt3) {
    if(otfile) count->Fill("jet pt3",1);
    if(otfile) acount->Fill(5.5);


    if(Cpt4) {
    if(otfile) count->Fill("jet pt4",1);
    if(otfile) acount->Fill(6.5);



      // require at least N emerging jets
    if(Cnem) {
      std::cout<<"PASS without almost"<<std::endl;
      std::cout<<"n emerging nealmost emergin is "<<nemerging<<" "<<nalmostemerging<<std::endl;
      std::cout<<Canem<<std::endl;



      if(otfile) count->Fill("emerging",1);
      if(otfile) acount->Fill(7.5);
      if(Canem) {
	std::cout<<"PASS with almost"<<std::endl;

        if(otfile) count->Fill("almostemerging",1);
        if(otfile) acount->Fill(8.5);


          npass+=1;
	  std::cout<<"passing run lumi event filename is "<<run<<" "<<lumi<<" "<<event<<" "<<inputfilename<<std::endl;
	  for(int i=0;i<4;i++) {
	    std::cout<<"  for jet "<<i<<" pt eta nef cfe ntrkpt1 alphamax r0"<<std::endl;
	    std::cout<<"     "<<(*jet_pt)[i]<<" "<<(*jet_eta)[i]<<" "<<(*jet_nef)[i]<<" "<<(*jet_cef)[i]<<" "<<jet_ntrkpt1[i]<<" "<<(*jet_alphaMax)[i]<<" "<<r0[i]<<" "<<std::endl;
	  }
          if(otfile) {
	    H_T3->Fill(HT);   
	    float mass;
	    for(int i5=0;i5<4;i5++) {
	    for(int i6=i5+1;i6<4;i6++) {
	      if((emerging[i5]&&!emerging[i6])||(!emerging[i5]&&emerging[i6])) {
	      mass = sqrt(
			  pow((jet_e[i5]+jet_e[i6]),2) -
			  pow((jet_px[i5]+jet_px[i6]),2) -
			  pow((jet_py[i5]+jet_py[i6]),2) -
			  pow((jet_pz[i5]+jet_pz[i6]),2)

                );
	      hmass->Fill(mass);
	    }}}

	      
	  }

    

    std::cout<<"npass  event is "<<npass<<" "<<event<<std::endl;
    std::cout<<"nemerging nalmostemerging "<<nemerging<<" "<<nalmostemerging<<std::endl;

    }}}}}}}}

  }  // end of loop over events

  if(otfile) {
    TFile myfile(outputfilename,"RECREATE");
    count->LabelsDeflate();
    count->LabelsOption("v");
  //  count->LabelsOption("a");

    eventCountPreTrigger->Write();
    acount->Write();
    count->Write();
    hjetcut->Write();
    hpt->Write();
    hnjet->Write();
    heta->Write();
    heta2->Write();
    halpha->Write();
    haMgj->Write();
    H_T->Write();
    H_T2->Write();
    H_T3->Write();
    H_T4->Write();
    hpt1->Write();
    hpt2->Write();
    hpt3->Write();
    hpt4->Write();
    h_nemg->Write();
    hjetchf->Write();
    hbcut_ntrkpt1->Write();
    hacut_ntrkpt1->Write();
    hbcut_nef->Write();
    hacut_nef->Write();
    hbcut_cef->Write();
    hacut_cef->Write();
    hbcut_alphamax->Write();
    hacut_alphamax->Write();
    hHTnm1->Write();
    hpt1nm1->Write();
    hpt2nm1->Write();
    hpt3nm1->Write();
    hpt4nm1->Write();
    halphanm1->Write();
    hmaxipnm1->Write();
    hnmaxipnm1->Write();
    hn2maxipnm1->Write();
    hnHitsnm1->Write();
    hntrk1nm1->Write();
    hnemnm1->Write();
    hipXYEJ->Write();
    hipXYnEJ->Write();
    htvw->Write();
    htvwEJ->Write();
    hipXYSigEJ->Write();
    hipXYSignEJ->Write();
    hmaxipXYEJ->Write();
    hmaxipXYnEJ->Write();
    hmeanipXYEJ->Write();
    hmeanipXYnEJ->Write();
    hjptb->Write();
    hjpta->Write();
    hjptfrb->Write();
    hjptfra1->Write();
    hjptfra2->Write();
    hjptfrbc->Write();
    hjptfra1c->Write();
    hjptfra2c->Write();
    hHTko->Write();
    hpt1ko->Write();
    hpt2ko->Write();
    hpt3ko->Write();
    hpt4ko->Write();
    hmass->Write();

    //2d
    aMip->Write();
    haMvjpt->Write();
    haMvHT->Write();
    haMvnvtx->Write();

    myfile.Close();
  }

  tt->ResetBranchAddresses();
  
  delete jet_index;
  delete jet_source;
  delete jet_pt;
  delete jet_eta;
  delete jet_phi;
  delete jet_alphaMax;
  delete jet_cef;
  delete jet_nef;
  delete jet_chf;
  //  delete jet_phf;
  delete track_pt;
  delete track_eta;
  delete track_source;
  delete track_index;
  delete track_jet_index;
  delete track_vertex_index;
  delete track_algo;
  delete track_vertex_weight;
  delete track_ipZ;
  delete track_ipXY;
  delete track_ipXYSig;
  


  f->Close();
  


  return npass;
}
