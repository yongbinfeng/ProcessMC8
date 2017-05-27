
#include <iostream>
#include <string>
#include <map>

void QCDhists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* binnames,std::string* aaname,std::string ohname, bool hasPre,bool norm, bool blind, bool b16003, float themass) ;



int main(int argc, char *argv[])
{ 
    int imode=atoi(argv[1]);
    int iblind=atoi(argv[2]);
    int i16003=atoi(argv[3]);
    float TheMass=atof(argv[4]);


    std::cout<<" doing cuts for hypothesis signal mass of "<<TheMass<<std::endl;

    bool b16003 = false;
    if(i16003>0) b16003=true;

    bool blind=false;
    if(iblind!=0) blind=true;

    bool hasPre=true;


    if(imode==0) {
        std::cout<<"doing background"<<std::endl;
    } else if(imode==1) {
        std::cout<<"doing modelA"<<std::endl;
    } else if(imode==2) {
        std::cout<<"doing modelB"<<std::endl;
    } else if(imode==3) {
        std::cout<<"doing quick QCD"<<std::endl;
    } else if(imode==4) {
        std::cout<<"doing debug sample"<<std::endl;
    } else if(imode==5) {
        std::cout<<"doing Wjet data sample"<<std::endl;
        hasPre=false;
    } else if(imode==6) {
        std::cout<<"doing Wjet MC sample"<<std::endl;
        hasPre=false;
    } else if(imode==7) {
        std::cout<<"doing DATA"<<std::endl;
        hasPre=true;
        blind=true;
    } else if(imode==8) {
        std::cout<<"doing QCD74"<<std::endl;
    } else {
        std::cout<<"invalid choice"<<std::endl;
    }


    float goalintlum=20; // fb-1                                                                                        
    //float goalintlum=0.07956; // fb-1                                                                                        

    // for background 

    const int nbin=3; // 1000-1500,1500-2000,200toInf
    float xsec[nbin]={1064000,121500,25420}; // fb 
    int nfiles[nbin]={47, 49, 14};
    std::string binnames[nbin]={"QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"};
    std::string aaname[nbin]={"configs/QCD_HT1000to1500_80X_config.txt", "configs/QCD_HT1500to2000_80X_config.txt", "configs/QCD_HT2000toInf_80X_config.txt"};

    /*
    const int nbin=5; // 500-700,700-1000,1000-1500,1500-2000,200toInf
    float xsec[nbin]={29370000,6524000,1064000,121500,25420}; // fb 
    int nfiles[nbin]={229, 170, 47, 49, 14};
    std::string binnames[nbin]={"QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"};
    std::string aaname[nbin]={"configs/QCD_HT500to700_80X_config.txt", "configs/QCD_HT700to1000_80X_config.txt", "configs/QCD_HT1000to1500_80X_config.txt", "configs/QCD_HT1500to2000_80X_config.txt", "configs/QCD_HT2000toInf_80X_config.txt"};
    */

    // quick background
    const int qnbin=1; // 500-700,700-1000,1000-1500,1500-2000,200toInf
    float qxsec[qnbin]={25420}; // fb 
    int qnfiles[qnbin]={3};
    std::string qbinnames[qnbin]={"QCD_HT2000toInf"};
    std::string aaname_q[qnbin]={"configs/QCD_HT2000toInf_80X_config.txt"};

    // for signal models A.  mediat mass is 1000
    const int anbin=1; 
    float axsec[anbin]={18.45}; // fb 
    int anfiles[anbin]={31}; 
    //int anfiles[anbin]={5}; 
    std::string abinnames[anbin]={"modelA"};
    std::string aaname_a[anbin]={"configs/ModelA_80X_config.txt"};

    // for signal models B.  mediat mass is 1000
    const int bnbin=1; 
    float bxsec[bnbin]={18.45}; // fb 
    int bnfiles[bnbin]={30}; 
    //int bnfiles[bnbin]={5}; 
    std::string bbinnames[bnbin]={"modelB"};
    std::string aaname_b[bnbin]={"configs/ModelB_80X_config.txt"};


    // for debugging
    const int dnbin=1; 
    float dxsec[dnbin]={18.45}; // fb 
    int dnfiles[dnbin]={1}; 
    std::string dbinnames[dnbin]={"tmpStore"};
    std::string aaname_d[dnbin]={"configs/debug_config.txt"};

    // Wjets data sample
    const int wnbin=1;
    float wxsec[wnbin]={11811000}; // fb 
    //int wnfiles[wnbin]={345};
    int wnfiles[wnbin]={150};
    std::string wbinnames[wnbin]={"WSkim"};
    std::string aaname_w[wnbin]={"configs/WSkim_80X_config.txt"};

    // Wjets MC sample
    const int wmcnbin=1; 
    float wmcxsec[wmcnbin]={11811000}; // fb 
    int wmcnfiles[wmcnbin]={898};
    std::string wmcbinnames[wmcnbin]={"WMCSkim"};
    std::string aaname_wmc[wmcnbin]={"configs/WMCSkim_80X_config.txt"};

    // DATA
    const int datanbin=1; 
    float dataxsec[datanbin]={11811000}; // fb 
    int datanfiles[datanbin]={19};
    std::string databinnames[datanbin]={"DATA"};
    std::string aaname_data[dnbin]={"configs/Data_80X_config.txt"};


    //QCD74
    const int q74nbin=5; // 500-700,700-1000,1000-1500,1500-2000,200toInf
    float q74xsec[q74nbin]={29370000,6524000,1064000,121500,25420}; // fb 
    int q74nfiles[q74nbin]={229, 170, 47, 49, 14};
    std::string q74binnames[q74nbin]={"QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"};
    std::string aaname_q74[q74nbin]={"configs/QCD_HT500to700_80X_config.txt", "configs/QCD_HT700to1000_80X_config.txt", "configs/QCD_HT1000to1500_80X_config.txt", "configs/QCD_HT1500to2000_80X_config.txt", "configs/QCD_HT2000toInf_80X_config.txt"};




    if(imode==0) {
      QCDhists(goalintlum,nbin,xsec,nfiles,binnames,aaname,"SumHistsQCD.root",hasPre,true,blind,b16003,TheMass);
    } else if (imode==1) {
        QCDhists(goalintlum,anbin,axsec,anfiles,abinnames,aaname_a,"SumHistsModelA.root",hasPre,true,blind,b16003,TheMass);
    } else if (imode==2) {
        QCDhists(goalintlum,bnbin,bxsec,bnfiles,bbinnames,aaname_b,"SumHistsModelB.root",hasPre,true,blind,b16003,TheMass);
    } else if (imode==3) {
        QCDhists(goalintlum,qnbin,qxsec,qnfiles,qbinnames,aaname_q,"SumHistsQQCD.root",hasPre,true,blind,b16003,TheMass);
    } else if (imode==4) {
        QCDhists(goalintlum,dnbin,dxsec,dnfiles,dbinnames,aaname_d,"SumHistsDebug.root",hasPre,true,blind,b16003,TheMass);
    } else if (imode==5) {
        QCDhists(goalintlum,wnbin,wxsec,wnfiles,wbinnames,aaname_w,"SumHistsWSkim.root",hasPre,false,blind,b16003,TheMass);
    } else if (imode==6) {
        QCDhists(goalintlum,wmcnbin,wmcxsec,wmcnfiles,wmcbinnames,aaname_wmc,"SumHistsWMCSkim.root",hasPre,true,blind,b16003,TheMass);
    } else if (imode==7) {
        QCDhists(goalintlum,datanbin,dataxsec,datanfiles,databinnames,aaname_data,"SumHistsDATA.root",hasPre,false,blind,b16003,TheMass);
    } else if (imode==8) {
      QCDhists(goalintlum,q74nbin,q74xsec,q74nfiles,q74binnames,aaname_q74,"SumHistsQCD74.root",hasPre,true,blind,b16003,TheMass);
    }
}
