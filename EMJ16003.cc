#include <iostream>
#include <iomanip>
#include <locale>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "vector"
using std::vector;
#include "algorithm"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//TTree          *fChain;   //!pointer to the analyzed TTree or TChain
//Int_t           fCurrent; //!current Tree number in a TChain

float CalcMedian(std::vector<float> vec) {
    if(vec.empty()) return 0;
    else {
        std::sort(vec.begin(), vec.end());
        if(vec.size() % 2 == 0)
            return (vec[vec.size()/2-1] + vec[vec.size()/2])/2;
        else
            return vec[vec.size()/2];
    }
}



int EMJ16003(bool otfile, bool hasPre, const char* inputfilename,const char* outputfilename) {


 
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

    /*
      vector<int> *vertex_index=new vector<int>;
      vector<int> *vertex_source=new vector<int>;
      vector<float> *vertex_z = new vector<float>;
    */

    //get event count pre trigger



    //for ntuple
    tt->SetBranchAddress("nVtx",&nVtx);
    tt->SetBranchAddress("nTrueInt",&nTrueInt);
    tt->SetBranchAddress("nTracks",&nTracks);
    tt->SetBranchAddress("event",&event);
    tt->SetBranchAddress("lumi",&lumi);
    tt->SetBranchAddress("run",&run);
    tt->SetBranchAddress("met_pt",&met_pt);
    tt->SetBranchAddress("met_phi",&met_phi);
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

    /*
      tt->SetBranchAddress("vertex_index",&vertex_index);
      tt->SetBranchAddress("vertex_source",&vertex_source);
      tt->SetBranchAddress("vertex_z",&vertex_z);
    */

    // create a histograms
    TH1F *acount,*count,*hjetcut,*hjetchf,*h_nemg,*hnjet,*hpt,*heta,*heta2,*halpha,*H_T,*H_T2,*H_T3,*H_T4,*hbcut_ntrkpt1,*hacut_ntrkpt1,*hbcut_nef,*hacut_nef,*hbcut_cef,*hacut_cef,*hbcut_alphamax,*hacut_alphamax,*hHTnm1,*hnHitsnm1,*hntrk1nm1,*hmaxipnm1,*hpt1nm1,*hpt2nm1,*hpt3nm1,*hpt4nm1,*halphanm1,*hnemnm1,*hpt1,*hpt2,*hpt3,*hpt4,*hipXYEJ,*hipXYnEJ,*htvw,*htvwEJ,*hnmaxipnm1,*hn2maxipnm1,*hjptfrb,*hjptfra1,*hjptfra2,*hjptfrbc,*hjptfra1c,*hjptfra2c,*hjptb,*hjpta,*haMgj,*hHTko,*hpt1ko,*hpt2ko,*hpt3ko,*hpt4ko,*hipXYSigEJ,*hipXYSignEJ,*hmaxipXYEJ,*hmaxipXYnEJ,*hmeanipXYEJ,*hmeanipXYnEJ,*hmass,*hmedipXYEJ,*hmedipXYnEJ,*hmeanipXYSigEJ,*hmeanipXYSignEJ,*hmedipXYSigEJ,*hmedipXYSignEJ,*hlogmeanipXYSigEJ,*hlogmeanipXYSignEJ,*hlogmedipXYSigEJ,*hlogmedipXYSignEJ;

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
        hmedipXYEJ = new TH1F("hmedipXYEJ","median ip emerging jets",1000,0.,2.);
        hmedipXYnEJ = new TH1F("hmedipXYnEJ","median ip not emerging jets",1000,0.,2.);
        hmeanipXYSigEJ = new TH1F("hmeanipXYSigEJ","mean ip_{sig} emerging jets",1000,0.,100.);
        hmeanipXYSignEJ = new TH1F("hmeanipXYSignEJ","mean ip_{sig} not emerging jets",1000,0.,100.);
        hmedipXYSigEJ = new TH1F("hmedipXYSigEJ","median ip_{sig} emerging jets",1000,0.,100.);
        hmedipXYSignEJ = new TH1F("hmedipXYSignEJ","median ip_{sig} not emerging jets",1000,0.,100.);
        hlogmeanipXYSigEJ = new TH1F("hlogmeanipXYSigEJ","mean log_{10} ip_{sig} emerging jets",1000,-1.,4.);
        hlogmeanipXYSignEJ = new TH1F("hlogmeanipXYSignEJ","mean log_{10} ip_{sig} not emerging jets",1000,-1.,4.);
        hlogmedipXYSigEJ = new TH1F("hlogmedipXYSigEJ","median log_{10} ip_{sig} emerging jets",1000,-1.,4.);
        hlogmedipXYSignEJ = new TH1F("hlogmedipXYSignEJ","median log_{10} ip_{sig} not emerging jets",1000,-1.,4.);
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

    std::cout<<"!!!!!!! entering EMJ16003"<<std::endl;


    // loop over events
    for (Int_t i=0; i<nentries; i++) {
        //    std::cout<<"***event "<<event<<std::endl;
 
        if(!hasPre) eventCountPreTrigger->Fill(1.5); 
    
        if(otfile) count->Fill("All",1);  // count number of events
        if(otfile) acount->Fill(0.5);
        tt->GetEntry(i);
        //    std::cout<<"event number is "<<event<<" number of vertex is "<<nVtx<<std::endl;



        // calculate some variables about jets
        vector<int> jet_ntrkpt1((*jet_index).size());
        vector<float> jet_meanip((*jet_index).size());
        vector<float> jet_meanipsig((*jet_index).size());
        vector<float> jet_logmeanipsig((*jet_index).size());
        vector<float> jet_medip((*jet_index).size());
        vector<float> jet_medipsig((*jet_index).size());
        vector<float> jet_logmedipsig((*jet_index).size());
        vector<float> jet_medtheta2D((*jet_index).size());
        vector<float> jet_logmedtheta2D((*jet_index).size());
        vector<float> r0((*jet_index).size());
        vector<float> r0sig((*jet_index).size());
        vector<float> r1((*jet_index).size());
        vector<int> jntrack((*jet_index).size());
        vector<int> jntrkip1mm((*jet_index).size());
        vector<float> jet_e((*jet_index).size());
        vector<float> jet_theta((*jet_index).size());
        vector<float> jet_px((*jet_index).size());
        vector<float> jet_py((*jet_index).size());
        vector<float> jet_pz((*jet_index).size());
        double HT = 0.;
        if(otfile) hnjet->Fill((*jet_index).size()+0.5);
        int NNNjet = (*jet_index).size();
        for(Int_t j=0; j<NNNjet; j++) {
            //      std::cout<<"jet j = "<<j<<std::endl;
            jet_theta[j]=2.*atan(exp(-(*jet_eta)[j]));
            jet_e[j]=(*jet_pt)[j]/sin(jet_theta[j]);
            jet_px[j]=(*jet_pt)[j]*cos((*jet_phi)[j]);
            jet_py[j]=(*jet_pt)[j]*sin((*jet_phi)[j]);
            jet_pz[j]=(*jet_pt)[j]/tan(jet_theta[j]);
            if(((*jet_pt)[j]>40)&&(fabs((*jet_eta)[j])<3.0)) {//UMD uses pf jet (|eta|<2.4); Princeton uses Calo jets (|eta|<3.0)
                HT=HT+(*jet_pt)[j];
            }
				
            if(otfile) hpt->Fill((*jet_pt)[j]);
            if(otfile) heta->Fill((*jet_eta)[j]);
            if(otfile) hjetchf->Fill((*jet_chf)[j]);
            if(otfile) if(j<4) heta2->Fill((*jet_eta)[j]);
            if(otfile) halpha->Fill((*jet_alphaMax)[j]);
            //      calculate  number of tracks with pt > 1
            jet_ntrkpt1[j]=0;
            jet_meanip[j]=0.;
            jet_logmeanipsig[j]=0.;
            jet_medip[j]=0.;
            jet_medipsig[j]=0.;
            jet_logmedipsig[j]=0.;
            jet_medtheta2D[j]=999.;//FIXME: dummy value
            jet_logmedtheta2D[j]=999.;//FIXME: dummy value

            r0[j]=0.;
            r1[j]=0.;
            r0sig[j]=0.;
            vector<float> track_pts = track_pt->at(j);
            vector<int> track_sources = track_source->at(j);
            vector<float> track_vertex_weights = track_vertex_weight->at(j);
            vector<float> track_ipXYs = track_ipXY->at(j);
            vector<float> track_ipXYSigs = track_ipXYSig->at(j);
            vector<float> sort_ip(track_pts.size());
            vector<float> sort_ipsig(track_pts.size());
            for(uint it=0;it<track_pts.size();it++) sort_ip[it]=0;
            for(uint it=0;it<track_pts.size();it++) sort_ipsig[it]=0;
            vector<float> jet_trkip;
            vector<float> jet_trkipsig;
            vector<float> jet_theta2D;

            jntrack[j]=0;
            for (unsigned itrack=0; itrack<track_pts.size(); itrack++) {
                if((track_sources[itrack]==0)&&(track_pts[itrack]>1)) {
                    if(fabs(track_ipXYs[itrack]<1)) {
                        jntrkip1mm[j]=jntrkip1mm[j]+1;
                    }
                    sort_ip[jntrack[j]]=fabs(track_ipXYs[itrack]);
                    sort_ipsig[jntrack[j]]=fabs(track_ipXYSigs[itrack]);
                    if(otfile) htvw->Fill(track_vertex_weights[itrack]);
                    //	  std::cout<<"track vertex weight is "<<track_vertex_weights[itrack]<<std::endl;
                    if(track_pts[itrack]>1) jet_ntrkpt1[j]+=1;
                    //	  std::cout<<" track "<<itrack<<" ip "<<track_ipXYs[itrack]<<" mean ip "<<jet_meanip[j]<<std::endl;
                    jet_meanip[j]=jet_meanip[j]+fabs(track_ipXYs[itrack]);
                    jet_meanipsig[j]=jet_meanipsig[j]+fabs(track_ipXYSigs[itrack]);
                    jet_logmeanipsig[j]=jet_logmeanipsig[j]+fabs(track_ipXYSigs[itrack]);
                    jet_trkip.push_back(fabs(track_ipXYs[itrack]));
                    jet_trkipsig.push_back(fabs(track_ipXYSigs[itrack]));
                    //jet_theta2D.push_back(track_theta2D[itrack]);
                    jntrack[j]++;
                }
            }
            float atmp = jntrack[j];
            if(jntrack[j]>0) {
                // mean
                jet_meanip[j]=jet_meanip[j]/atmp;
                jet_logmeanipsig[j]=jet_logmeanipsig[j]/atmp;
                jet_logmeanipsig[j]=log10(jet_logmeanipsig[j]);

                // median
                jet_medip[j] = CalcMedian(jet_trkip);
                jet_medipsig[j] = CalcMedian(jet_trkipsig);
                jet_logmedipsig[j]=log10(jet_medipsig[j]);
                //jet_medtheta2D[j] = CalcMedian(jet_theta2D); //FIXME
                //jet_logmedtheta2D[j] = log10(jet_medtheta2D[j]); //FIXME
            }

            std::sort(sort_ip.begin(), sort_ip.end(), std::greater<float>());
            if(sort_ip.size()>0) r0[j]=sort_ip[0];
            if(sort_ip.size()>1) r1[j]=sort_ip[1];
            std::sort(sort_ipsig.begin(), sort_ipsig.end(), std::greater<float>());
            if(sort_ipsig.size()>0) r0sig[j]=sort_ipsig[0];
            //std::cout << "ipSigMax = " << r0sig[j] << std::endl;
            //      std::cout<<"mean max are "<<jet_meanip[j]<<" "<<r0[j]<<std::endl;
        }  // end of loop over jets
        //std::cout<<"HT = " << HT << std::endl;

        //now see which jets are emerging
        //    std::cout<<" in event "<<event<<" number of jets is "<<NNNjet<<std::endl;
        vector<bool> tagged(NNNjet);
        vector<bool> emerging(NNNjet);
        vector<bool> almostemerging(NNNjet);
        vector<bool> basicjet(NNNjet);
        vector<bool> trigjet(NNNjet);
        vector<bool> strigjet(NNNjet);
        for( int i=0;i<NNNjet;i++) {
            tagged[i]=false;
            emerging[i]=false;
            almostemerging[i]=false;
            basicjet[i]=false;
            trigjet[i]=false;
            strigjet[i]=false;
	}
        int ntagged=0;
        int nemerging=0;
        int nalmostemerging=0;
        int nkine=0;
        int ntrigjet=0;
        int nstrigjet=0;
        //int iijjkk = 4;
        //if(NNNjet<4) iijjkk=NNNjet;
        //      std::cout<<"iijjkk is "<<iijjkk<<std::endl;
        for(int ij=0;ij<NNNjet;ij++) {
            vector<float> track_ipXYs = track_ipXY->at(ij);
            vector<float> track_ipXYSigs = track_ipXYSig->at(ij);
            vector<int> track_sources = track_source->at(ij);
            vector<float> track_vertex_weights = track_vertex_weight->at(ij);
            if(otfile) hjetcut->Fill(0.5);

            if((fabs((*jet_eta)[ij])<2.0)) { // Trig/kinematic eta cut
                if(otfile) hjetcut->Fill(1.5);

                // Not in EXO-16-003
                if(otfile) hbcut_nef->Fill((*jet_nef)[ij]);
                if((*jet_nef)[ij]<0.9) {  // neutral fraction
                    if(otfile) hacut_nef->Fill((*jet_nef)[ij]);
                    if(otfile) hjetcut->Fill(2.5);

                    if(otfile) hbcut_ntrkpt1->Fill(jet_ntrkpt1[ij]);
                    if(jet_ntrkpt1[ij]>1) {  // tracks pt>1
                        if(otfile) hacut_ntrkpt1->Fill(jet_ntrkpt1[ij]);
                        if(otfile) hjetcut->Fill(3.5);

                        if(otfile) hbcut_cef->Fill((*jet_cef)[ij]);
                        if((*jet_cef)[ij]<0.9) {  //charged fraction
                            if(otfile) hacut_cef->Fill((*jet_cef)[ij]);
                            if(otfile) hjetcut->Fill(4.5);
                            basicjet[ij]=true;
                        }
                    }
                }

                // check for trigger jets
                if((*jet_pt)[ij]>40) {
                    if(jntrkip1mm[ij]<3) {
                        trigjet[ij]=true;
                        ntrigjet=ntrigjet+1;
                        if(r0sig[ij]>5) {
                            strigjet[ij]=true;
                            nstrigjet=nstrigjet+1;
                        }
                    }
                }

                if( (*jet_pt)[ij]>60 ) {  // jet pt kinematic cut
                    nkine+=1;
                    if(otfile) hjetcut->Fill(5.5);
                    if(otfile) hbcut_alphamax->Fill((*jet_alphaMax)[ij]);
                    if((*jet_alphaMax)[ij]<0.05) { // alpha max 5%
                        if(otfile) hacut_alphamax->Fill((*jet_alphaMax)[ij]);
                        if(otfile) hjetcut->Fill(6.5);
                        almostemerging[ij]=true;
                        nalmostemerging=nalmostemerging+1;

                        if(jet_logmedipsig[ij]>1.5) { // log median ipsig
                            emerging[ij]=true;
                            nemerging+=1.;

                            if (jet_logmedtheta2D[ij]>-1.6) { // FIXME: need to calculate this variable; dummy value used
                                tagged[ij]=true;
                                ntagged+=1;

                                // look at tracks in the emerging jets
                                if(otfile) hmaxipXYEJ->Fill(r0[ij]);
                                if(otfile) hmeanipXYEJ->Fill(jet_meanip[ij]);
                                if(otfile) hmedipXYEJ->Fill(jet_medip[ij]);
                                if(otfile) hmeanipXYSigEJ->Fill(jet_meanipsig[ij]);
                                if(otfile) hmedipXYSigEJ->Fill(jet_medipsig[ij]);
                                if(otfile) hlogmeanipXYSigEJ->Fill(jet_logmeanipsig[ij]);
                                if(otfile) hlogmedipXYSigEJ->Fill(jet_logmedipsig[ij]);

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
            }
            if(!emerging[ij]) {
                if(otfile) hmaxipXYnEJ->Fill(r0[ij]);
                if(otfile) hmeanipXYnEJ->Fill(jet_meanip[ij]);
                if(otfile) hmedipXYnEJ->Fill(jet_medip[ij]);
                if(otfile) hmeanipXYSignEJ->Fill(jet_meanipsig[ij]);
                if(otfile) hmedipXYSignEJ->Fill(jet_medipsig[ij]);
                if(otfile) hlogmeanipXYSignEJ->Fill(jet_logmeanipsig[ij]);
                if(otfile) hlogmedipXYSignEJ->Fill(jet_logmedipsig[ij]);
                for (unsigned itrack=0; itrack<track_ipXYs.size(); itrack++) {
                    if(track_sources[itrack]==0) {
                        if(otfile) hipXYnEJ->Fill(track_ipXYs[itrack]);
                        if(otfile) hipXYSignEJ->Fill(track_ipXYSigs[itrack]);
                    }
                }

            }
            //	std::cout<<"event pt alphaM cef nef ntrkpt1 r0 emerging  almost "<<event<<" "<<(*jet_pt)[ij]<<" "<<(*jet_alphaMax)[ij]<<" "<<(*jet_cef)[ij]<<" "<<(*jet_nef)[ij]<<" "<<jet_ntrkpt1[ij]<<" "<<r0[ij]<<" "<<emerging[ij]<<" "<<almostemerging[ij]<<std::endl;
        }
        if(otfile) h_nemg->Fill(nemerging);



        // *************************************************************
        // now start the event selections
        // *************************************************************


        // trig1 requires at least 2 jets (pt>40GeV; jntrkip1mm<3)
        bool Ctrig1=false;
        if((HT>500)&&(ntrigjet>1)) {
            Ctrig1=true;
        }

        // trig2
        bool Ctrig2=false;
        if((HT>350)&&(nstrigjet>1)) {
            Ctrig2=true;
        }

        // HT

        if(otfile) H_T->Fill(HT);
        if(otfile) hpt1->Fill((*jet_pt)[0]);
        if(otfile) hpt2->Fill((*jet_pt)[1]);
        if(otfile && NNNjet>2) hpt3->Fill((*jet_pt)[2]);
        if(otfile && NNNjet>3) hpt4->Fill((*jet_pt)[3]);
        bool CHT=false;
        if((Ctrig1)&&(HT>650)) CHT=true;
        if((Ctrig2)&&(HT>450)) CHT=true;

        // Kinematic cuts
        bool Ckine = false;
        if(nkine>1) Ckine=true;

        // number emerging jets
//         bool Cnem = false;
//         if(nemerging>1) Cnem=true;

        // number tagged jets
        bool Cntag = false;
        if (ntagged>1) Cntag=true;


        // apply cuts sequentially

        //    std::cout<<"c4jet cht cpt1 cpt2 cpt3 cpt4 cnem "<<C4jet<<" "<<CHT<<" "<<Cpt1<<" "<<Cpt2<<" "<<Cpt3<<" "<<Cpt4<<" "<<Cnem<<std::endl;

        if(Ctrig1||Ctrig2) {
            if(otfile) count->Fill("trigger",1);
            if(otfile) acount->Fill(1.5);

            // calculate HT and require it greater than some cut value
            if(CHT) {
                if(otfile) count->Fill("HT",1);
                if(otfile) acount->Fill(2.5);
                if(otfile) H_T2->Fill(HT);

                if(Ckine) {
                    if(otfile) count->Fill("Kinematic",1);
                    if(otfile) acount->Fill(3.5);
                    if(otfile) H_T3->Fill(HT);

                    // require at least 2 tagged jets
                    if(Cntag) {
                        if(otfile) count->Fill("Ntag",1);
                        if(otfile) {
                            acount->Fill(4.5);
                            H_T4->Fill(HT);
                            npass+=1;
                            if(ntagged>2) acount->Fill(5.5);
                            if(ntagged>3) acount->Fill(6.5);
                        }

                        std::cout<<"passing run lumi event filename is "<<run<<" "<<lumi<<" "<<event<<" "<<inputfilename<<std::endl;
                        for(int i=0;i<NNNjet;i++) {
                            std::cout<<"  for jet "<<i<<" pt eta nef cfe ntrkpt1 alphamax r0"<<std::endl;
                            std::cout<<"     "<<(*jet_pt)[i]<<" "<<(*jet_eta)[i]<<" "<<(*jet_nef)[i]<<" "<<(*jet_cef)[i]<<" "<<jet_ntrkpt1[i]<<" "<<(*jet_alphaMax)[i]<<" "<<r0[i]<<" "<<std::endl;
                        }
                        if(otfile) {
                            float mass;
                            for(int i5=0;i5<NNNjet;i5++) {
                                for(int i6=i5+1;i6<NNNjet;i6++) {
                                    if((emerging[i5]&&!emerging[i6])||(!emerging[i5]&&emerging[i6])) {
                                        mass = sqrt(
                                                    pow((jet_e[i5]+jet_e[i6]),2) -
                                                    pow((jet_px[i5]+jet_px[i6]),2) -
                                                    pow((jet_py[i5]+jet_py[i6]),2) -
                                                    pow((jet_pz[i5]+jet_pz[i6]),2)
                                                    );
                                        hmass->Fill(mass);
                                    }
                                }
                            }
                        }
                        std::cout<<"npass,  event "<<npass<<" "<<event<<std::endl;
                        std::cout<<"ntagged, nemerging, nalmostemerging "<<ntagged<<" "<<nemerging<<" "<<nalmostemerging<<std::endl;
                    }
                }
            }
        }

    }  // end of loop over events

    if(otfile) {
        TFile myfile(outputfilename,"RECREATE");
        //count->LabelsDeflate();
        //count->LabelsOption("v");
        //count->LabelsOption("a");

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
        hmedipXYEJ->Write();
        hmedipXYnEJ->Write();
        hmeanipXYSigEJ->Write();
        hmeanipXYSignEJ->Write();
        hmedipXYSigEJ->Write();
        hmedipXYSignEJ->Write();
        hlogmeanipXYSigEJ->Write();
        hlogmeanipXYSignEJ->Write();
        hlogmedipXYSigEJ->Write();
        hlogmedipXYSignEJ->Write();
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
