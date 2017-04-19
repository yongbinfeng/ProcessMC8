
#include <iostream>
#include <string>
#include <map>

void QCDhists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* binnames,std::string aaname,std::string ohname, int dooptk, int doopta,bool hasPre,bool norm, bool blind, bool b16003) ;



int main(int argc, char *argv[])
{ 
    int dooptk =*(argv[1])-'0';
    int doopta =*(argv[2])-'0';
    int imode=*(argv[3])-'0';
    int iblind=*(argv[4])-'0';
    int i16003=*(argv[5])-'0';

    bool donorm=true;

    bool b16003 = false;
    if(i16003>0) b16003=true;

    bool blind=false;
    if(iblind!=0) blind=true;

    bool hasPre=true;

    if(dooptk==0) {
        std::cout<<"not doing kinematic optimization"<<std::endl;
    } else if(dooptk==1) {
        std::cout<<"doing kinematic optimization"<<std::endl;
    } else {
        std::cout<<"invalid choice"<<std::endl;
    }

    if(doopta==0) {
        std::cout<<"not doing optimization on an emerging jet variable"<<std::endl;
    } else if(doopta==1) {
        std::cout<<"doing alphamax optimization"<<std::endl;
    } else if(doopta==2) {
        std::cout<<"doing maxIP optimization"<<std::endl;
    } else {
        std::cout<<"invalid choice"<<std::endl;
    }

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
	donorm=false;
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

    std::string aaname = "/data/users/eno/outputQCD/";  // area containing subdirectors with YHS's ntuples

    // for background 
    //const int nbin=2; // 500-700,700-1000,1000-1500,1500-2000,200toInf
    //float xsec[nbin]={29370000,6524000}; // fb 
    //int nfiles[nbin]={1,1};                                                                                
    //std::string binnames[nbin]={"QCD_HT1500to2000","QCD_HT2000toInf"};

    const int nbin=5; // 500-700,700-1000,1000-1500,1500-2000,200toInf
    float xsec[nbin]={29370000,6524000,1064000,121500,25420}; // fb 
    int nfiles[nbin]={138,133,50,40,23};
    std::string binnames[nbin]={"QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"};

    // quick background
    const int qnbin=1; // 500-700,700-1000,1000-1500,1500-2000,200toInf
    float qxsec[qnbin]={25420}; // fb 
    int qnfiles[qnbin]={3};
    std::string qbinnames[qnbin]={"QCD_HT2000toInf"};


    // for signal models A.  mediat mass is 1000
    const int anbin=1; 
    float axsec[anbin]={18.45}; // fb 
    int anfiles[anbin]={50}; 
    //int anfiles[anbin]={5}; 
    std::string abinnames[anbin]={"modelA"};

    // for signal models B.  mediat mass is 1000
    const int bnbin=1; 
    float bxsec[bnbin]={18.45}; // fb 
    //int bnfiles[bnbin]={50}; 
    int bnfiles[bnbin]={5}; 
    std::string bbinnames[bnbin]={"modelB"};


    // for debugging
    const int dnbin=1; 
    float dxsec[dnbin]={18.45}; // fb 
    int dnfiles[dnbin]={1}; 
    if(imode==4) aaname = "/home/eno/em8/ProcessMC3/";
    //std::string dbinnames[dnbin]={"debugmodelB80"};
    std::string dbinnames[dnbin]={"debug"};
    //std::string dbinnames[dnbin]={"debugmodelB74"};


    // Wjets data sample
    const int wnbin=1; 
    float wxsec[wnbin]={11811000}; // fb 
    //int wnfiles[wnbin]={345};
    int wnfiles[wnbin]={150};
    std::string wbinnames[wnbin]={"WSkim"};

    // Wjets MC sample
    const int wmcnbin=1; 
    float wmcxsec[wmcnbin]={11811000}; // fb 
    int wmcnfiles[wmcnbin]={898};
    std::string wmcbinnames[wmcnbin]={"WMCSkim"};


    // DATA
    const int datanbin=1; 
    float dataxsec[datanbin]={11811000}; // fb 
    int datanfiles[datanbin]={19};
    std::string databinnames[datanbin]={"DATA"};


    //QCD74

    const int q74nbin=5; // 500-700,700-1000,1000-1500,1500-2000,200toInf
    float q74xsec[q74nbin]={29370000,6524000,1064000,121500,25420}; // fb 
    int q74nfiles[q74nbin]={183,153,65,44,26};
    std::string q74binnames[q74nbin]={"QCD74_HT500to700","QCD74_HT700to1000","QCD74_HT1000to1500","QCD74_HT1500to2000","QCD74_HT2000toInf"};




    if(imode==0) {
        QCDhists(goalintlum,nbin,xsec,nfiles,binnames,aaname,"SumHistsQCD.root",dooptk,doopta,hasPre,donorm,blind,b16003);
    } else if (imode==1) {
        QCDhists(goalintlum,anbin,axsec,anfiles,abinnames,aaname,"SumHistsModelA.root",dooptk,doopta,hasPre,donorm,blind,b16003);
    } else if (imode==2) {
        QCDhists(goalintlum,bnbin,bxsec,bnfiles,bbinnames,aaname,"SumHistsModelB.root",dooptk,doopta,hasPre,donorm,blind,b16003);
    } else if (imode==3) {
        QCDhists(goalintlum,qnbin,qxsec,qnfiles,qbinnames,aaname,"SumHistsQQCD.root",dooptk,doopta,hasPre,donorm,blind,b16003);
    } else if (imode==4) {
        QCDhists(goalintlum,dnbin,dxsec,dnfiles,dbinnames,aaname,"SumHistsDebug.root",0,0,hasPre,donorm,blind,b16003);
    } else if (imode==5) {
        QCDhists(goalintlum,wnbin,wxsec,wnfiles,wbinnames,aaname,"SumHistsWSkim.root",0,0,hasPre,donorm,blind,b16003);
    } else if (imode==6) {
        QCDhists(goalintlum,wmcnbin,wmcxsec,wmcnfiles,wmcbinnames,aaname,"SumHistsWMCSkim.root",0,0,hasPre,donorm,blind,b16003);
    } else if (imode==7) {
        QCDhists(goalintlum,datanbin,dataxsec,datanfiles,databinnames,aaname,"SumHistsDATA.root",0,0,hasPre,donorm,blind,b16003);
    } else if (imode==8) {
        QCDhists(goalintlum,q74nbin,q74xsec,q74nfiles,q74binnames,aaname,"SumHistsQCD74.root",dooptk,doopta,hasPre,donorm,blind,b16003);
    }
}
