/*
 * Description: Electron on Proton beam, tagged lambda and K+ final states 
 *              Please see README for instructions
 * ================================================================
 * Time-stamp: "2023-01-11 12:59:28 trottar"
 * ================================================================
 *
 * Author:  Kijun Park and Richard L. Trotta III <trotta@cua.edu>
 *
 * Copyright (c) trottar
 */

// General definition in this header
#include "EIC_mesonMC.h"
//#include "functions/tim_hobbs/timhobbs.h"
#include "functions/sf.h"
#include <time.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <TROOT.h>

// init DIS cteq pdf
void initcteqpdf();

int getIndexOfLargestElement(Double_t arr[], int FromIndex, int ToIndex);

// k-Lambda splitting function
double fykL(double y, double kT, double L, int type);

// Define form factor function
double f2p(double x);
double f2kp(double p, double x, double th);
double f2kptmono(double p, double x, double th);

// Define the cross-section
double F2N(double x, double q2, int inucl);
// Define lambda SF
double F2L(double x, double q2);
// new definition for collider frame to call CTEQ function
double cdissigma(double x, double y, double q2, double nu, double ep, int inucl );
double cdissigma_n(double x, double y_D, double q2, double nu, double eprime);
double cdissigma_p(double x, double y_D, double q2, double nu, double eprime);

// Define function calls for beam smearing
double sigma_th(double pInc, double mInc, double NormEmit, double betaSt);

int mainx(double xMin,double xMax, double Q2Min,double Q2Max, double rnum, const int nevts, const double pbeam, const double kbeam, bool smear){

  NEvts = nevts;
  PBeam = pbeam;
  kBeam = kbeam;

  Bool_t iran;
  
  if(smear){
    iran=kTRUE;      // TRUE==include incident beam emittance
  }else{
    iran=kFALSE;      // NO incident beam emittance
  }

  // Define the DIS PDF from CTEQ directory:  cteq-tbls/ctq66m/ctq66.00.pds
  initcteqpdf();
  
  TRandom3 ran3;// Random number generator [45ns/call]
  ran3.SetSeed(rnum);// Sets random generator seed which is used initialize the random number generator

  // Define particle properties
  double emass =  mElectron, prmass = MProton, MSpectator = MLambda,
    spmass = MSpectator, mKaon = mKp, kmass = mKaon, photonmass = 0.;
  int e_particle_charge = -1, pr_particle_charge = 1, sp_particle_charge = 0, k_particle_charge = 1;
  int e_particle_id = 11, pr_particle_id =  2212, sp_particle_id = 3122, k_particle_id = 321, photon_particle_id = 22;

  // Number of incident and final state particles of reaction
  int NumPtls = 5, inucl = 1;

  // Incident electron information
  double e_momx_incident,e_momy_incident,e_momz_incident,e_energy_incident ;

  // Incident proton information
  double p_momx_incident,p_momy_incident,p_momz_incident,p_energy_incident ;

  // vertex for the initial electron
  double_t ve0X_Lab,ve0Y_Lab,ve0Z_Lab ;
  
  // vertex for the initial proton
  double_t vp0X_Lab,vp0Y_Lab,vp0Z_Lab ;

  // TDIS scattered electron information
  double pex_Res,pey_Res,pez_Res,EeE_Res;
  double e_momx_lab,e_momy_lab,e_momz_lab,e_energy_lab ;
  double_t vex_Lab,vey_Lab,vez_Lab ;
  double_t vex_Res,vey_Res,vez_Res ;
 
  // TDIS detected kaon information
  double pkx_Res,pky_Res,pkz_Res,EkE_Res;
  double k_momx_lab,k_momy_lab,k_momz_lab,k_energy_lab ;
  double_t vkx_Lab,vky_Lab,vkz_Lab ;
  double_t vkx_Res,vky_Res,vkz_Res ;
 
  // Virtual photon information
  double pphotonx_Res,pphotony_Res,pphotonz_Res,EphotonE_Res;
  double pphotonx_Lab,pphotony_Lab,pphotonz_Lab,EphotonE_Lab ;
  double_t vphotonx_Lab,vphotony_Lab,vphotonz_Lab ;
  double_t vphotonx_Res,vphotony_Res,vphotonz_Res ;

  // TDIS spectator lambda information
  double lambda_momx_lab,lambda_momy_lab,lambda_momz_lab,lambda_energy_lab ;
  double_t vlambx_Lab,vlamby_Lab,vlambz_Lab ;

  // missing particle (X) information
  double pXx_Res,pXy_Res,pXz_Res,EXE_Res ;
  double_t vXx_Res,vXy_Res,vXz_Res ;
  double X_momx_lab,X_momy_lab,X_momz_lab,X_energy_lab ;
  double_t vXx_Lab,vXy_Lab,vXz_Lab ;

  TLorentzVector kIncident_Vertex, PIncident_Vertex, qVirtual_Vertex, kScattered_Vertex;

  // Define daughter particle Vertex
  TLorentzVector pSpectator_Vertex,  PX_Vertex, PX_Vertex_Rest;
  TLorentzVector pScatterLambda_Vertex;
  TLorentzVector pScatterKaon_Vertex;

  // TDIS-SBS xbj means xbj_D (xbj of nucleon should be xbj/2 ? Asking Dasuni for xbj definition
  double Px_p2,Py_p2,Pz_p2,P_p2,theta_p2,phi_p2,E_p2, p2_pt, p2_z;

  double xVkR, xSkR;
  double yy,rho_PS, W2;
  double F2k,kTk,F2kK;
  
  // Kaon structure function selection
  //  typ = 1 ! s-exp For factor (L=,dis=0, FLAG=0) | J = 0 + 1/2
  //  typ = 2 ! t-exp For factor (L=,dis=0, FLAG=0) | J = 0 + 1/2
  //  typ = 3 ! ZEUS parameterization with F2L ("lambda" SF)
  double typ = 3;

  //WE SET THE RENORMALIZATION CUT-OFF PARAMETER LAMBDA "L" = ... IN UNITS OF GeV
  // From 3Var_th.f
  //const double      L = 1.180      !COV DIPOLE: HSS NORM
  //const double      L = 1.710      !IMF DIPOLE: HSS NORM
  //const double      L = 1.630       !HSS +
  // const double      L = 1.560;      //!HSS CENT. VALUE
  //const double      L = 1.48D0      !HSS -

  // Kaon
  // From ssbar_conv/CS-gfor.f
  double      L = 1.003;  //!FOR LAMBDA    (L = [1.003 +/- 0.008] GeV);
  //const double          L = 1.011 !UPPER BOUND --- STAT.  | -- LAMBDA PRODUCTION
  //const double          L = 0.995 !LOWER BOUND --- STAT.  |
  //const double      L = 1.170  !FOR SIGMA+    (L = [1.17 +/- 0.01] GeV)
  //const double      L = 1.240  !FOR SIGMA*+   (L = [1.24 +/- 0.02] GeV)
  //***********************************************************************
  //HERE WE PLACE A GLOBAL FLAG FOR THE CHOICE OF THE WAVEFUNCTION SUPPRESSION FACTOR
  int type = 1; // DIPOLE FORM FACTOR
  // int type = 2; // EXPONENTIAL FORM FACTOR
  // int type = 3; // COV. DIPOLE FORM FACTOR
  // int type = 4; // DIPOLE -- s-channel Lambda exchange??
  
  double weight_tdis;
  
  // char tTitle[80], tName[18], rName[65];
  char tTitle[80], tName[18], rName[67];
  
  sprintf(tTitle,"p(e,e'K\u039B)X Event Generation %3.1f GeV/c x%4.1f GeV/c",kBeam, PBeam);
  sprintf(rName,"/w/eic-scshelf2104/users/singhav/EIC_mesonsf_generator/OUTPUTS/kaon_lambda_v2/k_lambda_crossing_%.3f-%.1fon%.1f_x%.4f-%.4f_q%.1f-%.1f.root", CrossingTheta, kBeam, PBeam,xMin,xMax,Q2Min,Q2Max);

  TFile fRoot(rName,"Recreate", tTitle);
  sprintf(tName,"k_lambda");

  TTree *tree = new TTree("Evnts",tTitle);
  TTree *proctree = new TTree("Process",tTitle);

  typedef struct{
    Double_t s_e, s_q,  Q2, xBj, nu, W, p_RT, tempVar,
      pDrest, x_D, y_D, Yplus,  tSpectator, tPrime,
      TwoPdotk, TwoPdotq, MX2,
      alphaS, pPerpS, pPerpZ;
  } Invariants;
  static Invariants invts;

  Double_t Jacob;

  //  ************************************************
  // Store physics variable for TDIS into ROOTs
  //  ************************************************
  double tk, yk, y, fk, k_int, xk;
  double  xL, tmin;
  double xsec_inclusive_dis;
  double xsec_tagged_dis;
  double f2N,f2L;   
  Double_t qMag, pDotq;
	
  double Q2, xbj, TDIS_t, TDIS_znq, TDIS_Mx2, y_inelasticity;

  const Int_t bufsize=32000;

  // invariants
  tree->Branch("invts",&invts,"s_e/D:s_q/D:Q2/D:xBj/D:nu/D:W/D:p_RT/D:tempVar/D:pDrest/D:x_D/D:y_D/D:Yplus/D:tSpectator/D:tPrime/D:TwoPdotk/D:TwoPdotq/D:MX2/D:alphaS/D:pPerpS/D:pPerpZ/D");
  // incoming particle electron, deuteron(proton)
  tree->Branch("e_Inc.", &kIncident_Vertex,bufsize,1);
  tree->Branch("P_Inc.",  &PIncident_Vertex, bufsize, 1);
  // outgoing particles (electron, virtual photon, proton1, proton2, kaon)
  tree->Branch("e_Scat.", &kScattered_Vertex, bufsize, 1);
  tree->Branch("q_Vir.", &qVirtual_Vertex, bufsize, 1);
  tree->Branch("lamb_scat.", &pScatterLambda_Vertex, bufsize, 1);
  tree->Branch("k.", &pScatterKaon_Vertex, bufsize, 1);

  // Store physics varaibles into ROOTs
  proctree->Branch("xk", &xk, "xk/D");
  proctree->Branch("yk", &yk, "yk/D");
  proctree->Branch("tk", &tk, "tk/D");
  proctree->Branch("fk", &fk, "fk/D");
  proctree->Branch("xL", &xL, "xL/D");
  proctree->Branch("tmin", &tmin, "tmin/D");
  proctree->Branch("k_int", &k_int, "k_int/D");
	
  proctree->Branch("f2N", &f2N, "f2N/D");
  proctree->Branch("f2L", &f2L, "f2L/D");
	
  // TDIS kinematic variables
  proctree->Branch("Q2", &Q2, "Q2/D");
  proctree->Branch("xbj", &xbj, "xbj/D");
  proctree->Branch("TDIS_t", &TDIS_t, "TDIS_t/D");
  proctree->Branch("TDIS_znq", &TDIS_znq, "TDIS_znq/D");
  proctree->Branch("TDIS_Mx2", &TDIS_Mx2, "TDIS_Mx2/D");
  proctree->Branch("y_inelasticity", &y_inelasticity, "y_inelasticity/D");

  proctree->Branch("xsec_inclusive_dis", &xsec_inclusive_dis, "xsec_inclusive_dis/D");
  proctree->Branch("xsec_tagged_dis", &xsec_tagged_dis, "xsec_tagged_dis/D");

  // Incident electron
  proctree->Branch("e_momx_incident", &e_momx_incident, "e_momx_incident/D");
  proctree->Branch("e_momy_incident", &e_momy_incident, "e_momy_incident/D");
  proctree->Branch("e_momz_incident", &e_momz_incident, "e_momz_incident/D");
  proctree->Branch("e_energy_incident", &e_energy_incident, "e_energy_incident/D");

  // Incident proton
  proctree->Branch("p_momx_incident", &p_momx_incident, "p_momx_incident/D");
  proctree->Branch("p_momy_incident", &p_momy_incident, "p_momy_incident/D");
  proctree->Branch("p_momz_incident", &p_momz_incident, "p_momz_incident/D");
  proctree->Branch("p_energy_incident", &p_energy_incident, "p_energy_incident/D");
  
  // TDIS scattered electron information
  proctree->Branch("e_momx_lab", &e_momx_lab, "e_momx_lab/D");
  proctree->Branch("e_momy_lab", &e_momy_lab, "e_momy_lab/D");
  proctree->Branch("e_momz_lab", &e_momz_lab, "e_momz_lab/D");
  proctree->Branch("e_energy_lab", &e_energy_lab, "e_energy_lab/D");

  // TDIS detected kaon information
  proctree->Branch("k_momx_lab", &k_momx_lab, "k_momx_lab/D");
  proctree->Branch("k_momy_lab", &k_momy_lab, "k_momy_lab/D");
  proctree->Branch("k_momz_lab", &k_momz_lab, "k_momz_lab/D");
  proctree->Branch("k_energy_lab", &k_energy_lab, "k_energy_lab/D");

  // TDIS spectator lambda information
  proctree->Branch("lambda_momx_lab", &lambda_momx_lab, "lambda_momx_lab/D");
  proctree->Branch("lambda_momy_lab", &lambda_momy_lab, "lambda_momy_lab/D");
  proctree->Branch("lambda_momz_lab", &lambda_momz_lab, "lambda_momz_lab/D");
  proctree->Branch("lambda_energy_lab", &lambda_energy_lab, "lambda_energy_lab/D");

  // missing particle (X) information
  proctree->Branch("X_momx_lab", &X_momx_lab, "X_momx_lab/D");
  proctree->Branch("X_momy_lab", &X_momy_lab, "X_momy_lab/D");
  proctree->Branch("X_momz_lab", &X_momz_lab, "X_momz_lab/D");
  proctree->Branch("X_energy_lab", &X_energy_lab, "X_energy_lab/D");
	
  double sig_eThx=0.0, sig_eThy=0.0, sig_iThx=0.0, sig_iThy=0.0;

  double MIon = MProton;
  
  // Global Phase Space Factor;
  //s_e - M_ion^2 - m_electron^2
  double uu, uv, uw, ux, uy, norm;
  
  double PIon = PBeam*ZBeam;

  double s_e0 = 2.*kBeam*(sqrt(PIon*PIon+MIon*MIon)+PIon*cos(CrossingTheta))+ MIon*MIon + mElectron*mElectron;
  
  printf("Your kinematics: [xBj_min:xBj_max] = [%9.6f:%9.6f] \n", xMin,xMax);
  printf("Your kinematics: [Q2_min:Q2_max] = [%9.6f:%9.6f] \n", Q2Min,Q2Max);
  printf("Incident Ion Mass %9.5f GeV \n", MIon);
  printf("Incident Electron, Ion Momenta: %8.4f, %8.2f GeV/c | s_0 = %10.4f GeV^2 \n", kBeam, PIon, s_e0);   
	
  if (iran) {
    sig_eThx = sigma_th(kBeam, mElectron, eEpsNX, eBetaStarX);
    sig_eThy = sigma_th(kBeam, mElectron, eEpsNY, eBetaStarY);
    sig_iThx = sigma_th(PBeam, MIon, iEpsNX, iBetaStarX);
    sig_iThy = sigma_th(PBeam, MIon, iEpsNY, iBetaStarY);
  }
  
  double kBeamMC, kBeamMCx, kBeamMCy, kBeamMCz;
  double PBeamMC, PBeamMCx, PBeamMCy, PBeamMCz;
  double e_energy_rest, e_momentum_rest, e_cos_theta_rest, e_phi_rest, e_cos_phi_rest;


  proctree->Branch("e_energy_rest", &e_energy_rest, "e_energy_rest/D");
  proctree->Branch("e_momentum_rest", &e_momentum_rest, "e_momentum_rest/D");
  proctree->Branch("e_phi_rest", &e_phi_rest, "e_phi_rest/D");
  proctree->Branch("e_cos_phi_rest", &e_cos_phi_rest, "e_cos_phi_rest/D");
  proctree->Branch("e_cos_theta_rest", &e_cos_theta_rest, "e_cos_theta_rest/D");
	
  TVector3       UnitXLab, UnitYLab, UnitZLab;
  
  UnitXLab.SetXYZ(1.0,0.0,0.0);
  UnitYLab.SetXYZ(0.0,1.0,0.0);
  UnitZLab.SetXYZ(0.0,0.0,1.0);
  
  TVector3       UnitXRest, UnitYRest, UnitZRest;
  TVector3       UnitXqCM,  UnitYqCM,  UnitZqCM;
  TVector3       BoostCM, BoostRest, BoostRest0;
  TVector3       kScat3vec, pS_3vec;
  TVector3       kScat3vec0, pS_3vec0;
  
  TLorentzVector kIncident_Rest, PIncident_Rest, qVirtual_Rest, p_DT;
  TLorentzVector kIncident_Rest0, PIncident_Rest0, qVirtual_Rest0, p_DT0, p_ST0, p_ST;
  TLorentzVector kScattered_Rest, pSpectator_Rest, PX_Vector_Rest;
  TLorentzVector kScattered_Rest0, pSpectator_Rest0, PX_Vector_Rest0;
  TLorentzVector kIncident_0, PIncident_0; // unsmeared
  
  kIncident_0.SetXYZM(0.0,0.0,-kBeam,mElectron);

  if(smear){
    PIncident_0.SetXYZM(PIon*sin(CrossingTheta)*cos(CrossingPhi),PIon*sin(CrossingTheta)*sin(CrossingPhi),PIon*cos(CrossingTheta),MIon);
  }else{
    PIncident_0.SetXYZM(PIon*-sin(CrossingTheta),0.,PIon*cos(CrossingTheta),MIon); // New IP6 specs
  }
  
  TTree *itree = new TTree("Meta",tTitle);
  itree->Branch("e0.", &kIncident_0,bufsize,1);
  itree->Branch("P0.",  &PIncident_0, bufsize, 1);
  
  typedef struct {
    Int_t nEvts;
    Float_t PhSpFct;
  } MonteCarlo;
  static	MonteCarlo mc;
  mc.nEvts = NEvts;
  mc.PhSpFct   = (Q2Max-Q2Min)*(log(xMax)-log(xMin))*4.*pi*pow((pSMax),3)/3.;
  itree->Branch("MC", &mc, "nEvts/I:PhSpFct/F");

  typedef struct {
    Double_t p2_pt, p2_z,
      phi_p2, theta_p2,
      Pz_p2, P_p2, E_p2,
      Px_p2, Py_p2,Pz_p2_1,
      Pz_p2_2,Pz_p2_3,Pz_p2_4; 
  } P2;
  static P2 p2;

  itree->Branch("P2", &p2, "p2_pt/D:p2_z/D:phi_p2/D:theta_p2/D:Pz_p2/D:P_p2/D:E_p2/D:Px_p2/D:Py_p2/D:Pz_p2_1/D:Pz_p2_2/D:Pz_p2_3/D:Pz_p2_4/D");
  itree->Branch("Jacob",&Jacob,"Jacob/D");

  double pS_rest, csThRecoil, phiRecoil;

  //name of output lund file for GEANT4 use
  ofstream OUT_lund (Form("/w/eic-scshelf2104/users/singhav/EIC_mesonsf_generator/OUTPUTS/kaon_lambda_v2/k_lambda_crossing_%.3f_%.1fon%.1f_x%.4f-%.4f_q%.1f-%.1f_lund.dat", CrossingTheta, kBeam, PBeam,xMin,xMax,Q2Min,Q2Max), ios::trunc);
  ofstream OUT_pythia (Form("/w/eic-scshelf2104/users/singhav/EIC_mesonsf_generator/OUTPUTS/kaon_lambda_v2/k_lambda_crossing_%.3f_%.1fon%.1f_x%.4f-%.4f_q%.1f-%.1f_pythia.dat", CrossingTheta, kBeam, PBeam,xMin,xMax,Q2Min,Q2Max), ios::trunc);
  ofstream OUT_hepmc (Form("/w/eic-scshelf2104/users/singhav/EIC_mesonsf_generator/OUTPUTS/kaon_lambda_v2/k_lambda_crossing_%.3f_%.1fon%.1f_x%.4f-%.4f_q%.1f-%.1f_hepmc.dat", CrossingTheta, kBeam, PBeam,xMin,xMax,Q2Min,Q2Max), ios::trunc);

  // // Pythia output header
  // OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "SIMPLE Event FILE" << endl;
  // OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "============================================" << endl;
  // OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "I, ievent, nParticles" << endl;
  // OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "============================================" << endl;
  // OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "I" << "\t" << "K(I,1)" << "\t" << "K(I,2)" << "\t" << "K(I,3)" << "\t" << "K(I,4)" << "\t" << "K(I,5)" << "\t" << "P(I,1)" << "\t" << "P(I,2)" << "\t" << "P(I,3)" << "\t" << "P(I,4)" << "\t" << "P(I,5)" << "\t" << "V(I,1)" << "\t" << "V(I,2)" << "\t" << "V(I,3)" << "\t" << endl;
  // OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "============================================" << endl;

   // Adjusted LUND header
  OUT_lund << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "Adjusted LUND Type" << endl;
  OUT_lund << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "============================================" << endl;
  OUT_lund << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "Event Index" << endl;
  OUT_lund << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "============================================" << endl;
  OUT_lund << setiosflags(ios::left)  << setiosflags(ios::fixed)  <<"                 "  <<  "Number of Particles" << " \t "  << "xBj" << " \t " << "Q2"  << " \t " << "s_e"  << " \t " << "1.0" << " \t " << "xpi" << " \t" << "ypi" << " \t"  << "tpi"  << " \t"  <<  " \t" << "xsec_inclusive_dis" << " \t" << "xsec_tagged_dis" << endl;
  OUT_lund << setiosflags(ios::left)  << setiosflags(ios::fixed) << "\t" << "Particle Index" << " \t " << "Particle Charge" << " \t " << "Particle States (0-initial, 1-final)" << " \t " << "Particle ID" << " \t " << "Index of Parent" <<  " \t "<< "Index of First Daughter" << " \t " << "x-momentum" << " \t " << "y-momentum" << " \t " << "z-momentum" << " \t " << "Particle energy" << " \t " << "Particle Mass" << " \t " << "x-vertex"  << " \t " << "y-momentum" << " \t " << "z-momentum" << endl;
  OUT_lund << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "============================================" << endl;

  // Pythia output header
  OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "SIMPLE Event FILE" << endl;
  OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "============================================" << endl;
  OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "I, ievent, nParticles" << endl;
  OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "============================================" << endl;
  OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "I" << "\t" << "K(I,1)" << "\t" << "K(I,2)" << "\t" << "K(I,3)" << "\t" << "K(I,4)" << "\t" << "K(I,5)" << "\t" << "P(I,1)" << "\t" << "P(I,2)" << "\t" << "P(I,3)" << "\t" << "P(I,4)" << "\t" << "P(I,5)" << "\t" << "V(I,1)" << "\t" << "V(I,2)" << "\t" << "V(I,3)" << "\t" << endl;
  OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "============================================" << endl;

  // HEPMC output header
  OUT_hepmc << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "HepMC::Version 3.02.02" << endl;
  OUT_hepmc << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "HepMC::Asciiv3-START_EVENT_LISTING" << endl;
  

  // Get LorentzVector for the pScattered Proton for TDIS in rest frame
  TLorentzVector pScatterLambda_Rest;

  double plamb_Lab;

  // *****************************************************************************
  // Generate missing Hadron(Kaon) for TDIS in rest frame
  // initial neutron 4momentum(-P1) - TDIS scattered proton 4momentum(P2)
  // *****************************************************************************
  double E_k,Px_k,Py_k,Pz_k;

  TLorentzVector pScatterKaon_Rest;

  TVector3 pScatterKaon_V3, pScatterLambda_V3;

  double P_k;

  double fkfac = 1;

  double minE = 0.0;
  
  //for progress bar
  double progress=0.0;

  Int_t MEvts=1;
  
  // Start of event generator
  // for (int iEvt=0; iEvt<=NEvts; iEvt++) {
  while(MEvts<NEvts+1){

    // Progress bar
    if(MEvts%1000==0) {	    
      int barWidth = 70;
      progress = ((double)MEvts/(double)NEvts);	    
      // cout<<iEvt<<"/"<<NEvts<<endl;
      // cout << progress << endl;
      std::cout << "[";
      double pos = barWidth * progress;
      for (double i = 0.; i < barWidth; ++i) {
	if (i < pos) std::cout << "=";
	else if (i == pos) std::cout << ">";
	else std::cout << " ";
      }
      std::cout << "] " << int(progress * 100.0) << " %\r";
      std::cout.flush();
    }	 
    
    Jacob = 1.0;
	  
    if (iran) {  // taking into account initial beam smearing
      kBeamMC = kBeam*ran3.Gaus(1.0,eDkOverk);
      kBeamMCx= kBeamMC*ran3.Gaus(0.0,sig_eThx);
      kBeamMCy= kBeamMC*ran3.Gaus(0.0,sig_eThy);
      kBeamMCz=-kBeamMC;  //  Angles are really Tangents
      PBeamMC = PIon*ran3.Gaus(1.0,iDPoverP);
      PBeamMCx= PBeamMC*ran3.Gaus(0.0,sig_iThx); // IP1 configuration
      PBeamMCy= PBeamMC*ran3.Gaus(0.0,sig_iThy);
      PBeamMCz= PBeamMC;	    
    } else { // no beam smearing
      kBeamMC = kBeam;
      kBeamMCx= 0.0;
      kBeamMCy= 0.0;
      kBeamMCz= -kBeam;
      PBeamMC = PIon;
      PBeamMCx= 0.0;
      PBeamMCy= 0.0;
      PBeamMCz= PBeamMC;
    }
    
    //  definition of beam particle vertex
    ve0X_Lab = 0;
    ve0Y_Lab = 0;
    ve0Z_Lab = 0;
    vp0X_Lab = 0;
    vp0Y_Lab = 0;
    vp0Z_Lab = 0;	 	
		
    // these are with beam smearing
    kIncident_Vertex.SetXYZM(kBeamMCx, kBeamMCy, kBeamMCz, mElectron);// Set 4-momentum of incident electron beam
    PIncident_Vertex.SetXYZM(PBeamMCx, PBeamMCy, PBeamMCz, MIon);// Set 4-momentum of incident ion beam
    
    //  Crossing angle ehric with smearing
    PIncident_Vertex.RotateY(CrossingTheta);
    PIncident_Vertex.RotateZ(CrossingPhi);
		
    // Calculation of invariant variable from two beam particles (s_e, s_q) with smearing
    invts.TwoPdotk = 2.*PIncident_Vertex.Dot(kIncident_Vertex);
    invts.s_e = MIon*MIon + mElectron*mElectron + invts.TwoPdotk;
    invts.TwoPdotq = 2.*PIncident_Vertex.Dot(qVirtual_Vertex); 
    invts.s_q = MIon*MIon + invts.TwoPdotq;   
	  
    // Boosting vectors from electron+Ion CM frame
    // Boosting vector from Ion rest frame
    BoostRest = PIncident_Vertex.BoostVector(); 		  
	  
    PIncident_Rest = PIncident_Vertex; 

    PIncident_Rest.Boost(-BoostRest); // should result in (0.,0.,0.,MIon)
    kIncident_Rest = kIncident_Vertex;
    kIncident_Rest.Boost(-BoostRest);

    // **** for debugging purpose
    // for the dubugging purpose: OK, CHECKED
    /*
      cout << " ---------------------------------------------------------------------------------------- "  << endl;
      cout << " ---------------------------------------------------------------------------------------- "  << endl;
      cout << "(A) kIncident_Vertex (Lab)= " << kIncident_Vertex.X() << ", " << kIncident_Vertex.Y() << ", " << kIncident_Vertex.Z() << endl;
      cout << "(B) kIncident_Rest (rest)= " << kIncident_Rest.X() << ", " << kIncident_Rest.Y() << ", " << kIncident_Rest.Z() << endl;
    //*/
		
    // We generated xBj and Q2 binning...
    // Q is not calculated by Ei- Ef
    // Ef is calculated by given Q2 from random number
    // Generate Q2, ln(xBj) uniformly    
    uu   = ran3.Uniform();
    invts.Q2   = Q2Max*uu + Q2Min*(1.-uu);
    uv   = ran3.Uniform();
    invts.xBj  = pow(xMin,1.-uv)*pow(xMax,uv);
    invts.x_D  = invts.xBj*(MProton/MIon);
    invts.y_D  = invts.Q2/(invts.x_D*invts.TwoPdotk);
    invts.Yplus = 1 + ((1-invts.y_D)*(1-invts.y_D));

    // **** for debugging purpose
    // for the dubugging purpose: OK, CHECKED
    if (invts.y_D>=(1.0-2.*mElectron*MIon/invts.TwoPdotk) ) {
      // Unphysical kinematics
      continue;
    }
    if (invts.y_D>=1./(1.+invts.x_D*MIon*MIon/invts.TwoPdotk) ) {
      // Unphysical kinematics
      continue;
    }
	  
    Jacob *=invts.xBj;   // Jacobian  dx/dlnx

    // scattered electron, phi angle in rest frame 
    e_phi_rest = pi*(2.*ran3.Uniform()-1.0);
    e_cos_phi_rest = cos(e_phi_rest);
    //  mElectron-->0 approximation
    e_energy_rest = invts.TwoPdotk*(1.-invts.y_D)/(2.*MIon);

    // **** for debugging purpose
    // for the dubugging purpose: OK, CHECKED
    if (e_energy_rest<mElectron) {
      // should never happen
      printf("illegal Rest frame scattered electron energy =%10.6f \n",e_energy_rest);
      continue;
    }
    
    e_momentum_rest = sqrt(e_energy_rest*e_energy_rest-mElectron*mElectron);// scattered electron momentum in rest frame
    e_cos_theta_rest = (2.*e_energy_rest*kIncident_Rest.E() - invts.Q2 - 2.*mElectron*mElectron) /(2.*e_momentum_rest*kIncident_Rest.P());   // scattered electron \cosine\theta in rest frame	                                                 	  

    // **** for debugging purpose
    // for the dubugging purpose: OK, CHECKED
    if (e_cos_theta_rest*e_cos_theta_rest>1.0) {
      // should never happen
      printf("illegal Rest frame cos(the) = %6.2f \n", e_cos_theta_rest);
      printf(" (k_Rest, k'_Rest, Q2, xBj) = (%8.3f,  %8.4f, %6.2f, %5.3f) \n",kIncident_Rest.E(),e_momentum_rest,invts.Q2,invts.xBj);
      printf(" (2k.P, invts.s_e, y_D, x_D) = (%6.2f, %6.2f,%10.6f, %8.4f) \n",invts.TwoPdotk, invts.s_e, invts.y_D, invts.x_D);
      continue;
    }
	  
    // Definition of unit vector in rest frame (norm vector)
    UnitZRest  = -kIncident_Rest.Vect();
    norm       = 1./UnitZRest.Mag();
    UnitZRest *= norm;
    UnitYRest  = UnitZRest.Cross(UnitXLab);
    norm       = 1./UnitYRest.Mag();
    UnitYRest *= norm;
    UnitXRest  = UnitYRest.Cross(UnitZRest);
		
    // 3-vector of scattered electron in rest frame with smearing
    kScat3vec = e_momentum_rest*(sin(acos(e_cos_theta_rest))*(cos(e_phi_rest)*UnitXRest+sin(e_phi_rest)*UnitYRest)-e_cos_theta_rest*UnitZRest);
    kScattered_Rest.SetVectM(kScat3vec, mElectron);
    qVirtual_Rest = kIncident_Rest-kScattered_Rest;
    kScattered_Vertex = kScattered_Rest;		

    // *** Scattered electron in REST frame
    // For debugging purpose: checked
    /*
      cout << kScat3vec(0) << "  " << kScat3vec(1) << "  " << kScat3vec(2)  << " sqrt(x*x+y*y+z*z) =  " << e_momentum_rest << endl;
    */
	  
    // Proton momentum in rest frame with and without smearing
    invts.pDrest = sqrt(PIncident_Rest(0)*PIncident_Rest(0)+PIncident_Rest(1)*PIncident_Rest(1)+PIncident_Rest(2)*PIncident_Rest(2));

    // ***
    // for the debugging purpose: OK, CHECKED
    /*
      cout << "electron rest frame =" << kScat3vec.X() << ", " << kScat3vec.Y() << ", " << kScat3vec.Z() << endl;
    */
	  
    // Back to Lab frame
    kScattered_Vertex.Boost(BoostRest);
    qVirtual_Vertex = kIncident_Vertex-kScattered_Vertex;
	  
    // Hadronic Unit vectors relative to q
    UnitZqCM = -qVirtual_Rest.Vect();
    norm     = 1./UnitZqCM.Mag();
    UnitZqCM*= norm;
    UnitYqCM = -kScat3vec.Cross(kIncident_Rest.Vect());
    norm     = 1./UnitYqCM.Mag();
    UnitYqCM*=norm;
    UnitXqCM = UnitYqCM.Cross(UnitZqCM);

    qMag = qVirtual_Rest.P();
    pDotq = BoostRest(0)*qVirtual_Rest(0)+BoostRest(1)*qVirtual_Rest(1)+BoostRest(2)*qVirtual_Rest(2);
    p_DT = PIncident_Vertex - (pDotq/qMag/qMag)*qVirtual_Rest;
    p_ST = pSpectator_Vertex -(pDotq/qMag/qMag)*qVirtual_Rest;		

    // Null for proton beam
    pSpectator_Rest  = PIncident_Rest;
    
    //  definition are moved at the beginning of code
    Q2 = invts.Q2;
    xbj = invts.Q2/(2*pSpectator_Rest.Dot(qVirtual_Rest));
    TDIS_znq = p2.p2_z*pSpectator_Rest.Dot(qVirtual_Rest);
    TDIS_Mx2 = (qVirtual_Rest + pSpectator_Rest).Mag2();
    y_inelasticity   = (pSpectator_Rest.Dot(qVirtual_Rest))/(pSpectator_Rest.Dot(kIncident_Rest));

    //***
    // For debugging purpose:
    /*
      cout << "(3) kinematics: TDIS.xBj = " << xbj <<  ", invts.xBj= " << invts.xBj <<  ", invts.Q2= " << invts.Q2<< endl; 
    */
    
    // Not needed!
    // neutron config = proton + kaon(-):  this proton called additional speactator in TDIS concept
    // randomize in pt and z
    p2_pt = gRandom->Uniform(0.005*PBeam); // .5% of incoming ion beam momentum, this is the limit of the transverse momentum of  recoil particle...
    p2_z = gRandom->Uniform(1.);
    phi_p2 = gRandom->Uniform(360.0*D2R);
    Px_p2 = p2_pt*cos(phi_p2);
    Py_p2 = p2_pt*sin(phi_p2);
    Pz_p2 = ( -1.0*TDIS_znq*qVirtual_Rest.Z() + sqrt( TDIS_znq*qVirtual_Rest.Z()*TDIS_znq*qVirtual_Rest.Z() + invts.Q2*(qVirtual_Rest.E()*qVirtual_Rest.E())*(MProton*MProton + p2_pt*p2_pt)- invts.Q2*TDIS_znq*TDIS_znq) )/invts.Q2/1000.;  // unit match with "p2_pt" (GeV)
    P_p2  = sqrt (Pz_p2*Pz_p2 + p2_pt*p2_pt);
    theta_p2 = acos (Pz_p2/P_p2);
    E_p2 = sqrt ( P_p2*P_p2 + MProton*MProton);
    
    // Define scattered proton kinematics
    // Generate pSpectator Spherically Symmetric in rest frame     
    uw       = ran3.Uniform();
    pS_rest   = (pSMax)*pow(uw,1./3.); // uniform in 3p^2 dp = d(p^3), pSMax=0.3
    ux       = ran3.Uniform();
    csThRecoil = (2.*ux-1.);
    uy       = ran3.Uniform();
    phiRecoil  = pi*(2.*uy-1.);
    
    pScatterLambda_V3  = pS_rest*sin(acos(csThRecoil))*(cos(phiRecoil)*UnitXqCM
							+ sin(phiRecoil)*UnitYqCM);
    pScatterLambda_V3 += pS_rest*csThRecoil*UnitZqCM;
    
    pScatterLambda_Rest.SetVectM(pScatterLambda_V3,MSpectator);

    Jacob     *= 1./(2.*pScatterLambda_Rest.E());

    PX_Vertex_Rest = kIncident_Rest+PIncident_Rest-(kScattered_Rest+pScatterLambda_Rest);

    invts.MX2 = PX_Vertex.M2();
    invts.alphaS = ABeam*(pS_rest*csThRecoil+pSpectator_Rest.E())/MIon;
    invts.pPerpS = pS_rest*sqrt(1.-csThRecoil*csThRecoil);

    // alpha cut for select event with minimizing coherrent effects
    if(pow(invts.alphaS-1,2)<0.0001){
      continue;
    }
	  
    // for the debugging purpose: SEEMS NOT CRAZY NUMBER....
    // if (TMath::Sqrt(TDIS_Mx2) > PBeam/2){
    //   cout << "---->TDIS missing mass =" << TMath::Sqrt(TDIS_Mx2)  << endl;
    // }
    
    E_k  = pSpectator_Rest.E() - pScatterLambda_Rest.E();
    Px_k = pSpectator_Rest.X() - pScatterLambda_Rest.X(); 
    Py_k = pSpectator_Rest.Y() - pScatterLambda_Rest.Y();
    Pz_k = pSpectator_Rest.Z() - pScatterLambda_Rest.Z();
      
    // pScatterKaon_Rest.SetXYZM(Px_k,Py_k,Pz_k,mKaon);
    pScatterKaon_Rest = pSpectator_Rest - (pScatterLambda_Rest+PX_Vertex_Rest);
    // pScatterKaon_Rest = pSpectator_Rest - pScatterLambda_Rest;
      
    pScatterKaon_V3.SetXYZ(Px_k,Py_k,Pz_k);

    // Back to Lab frame
    pScatterLambda_Vertex = pScatterLambda_Rest;	  
    pScatterLambda_Vertex.Boost(BoostRest);
    pSpectator_Vertex = pScatterLambda_Vertex;
    pScatterKaon_Vertex = pScatterKaon_Rest;
    pScatterKaon_Vertex.Boost(BoostRest);

    P_k = pScatterKaon_V3.Mag();
      
    lambda_momx_lab = pScatterLambda_Vertex.X();
    lambda_momy_lab = pScatterLambda_Vertex.Y();
    lambda_momz_lab = pScatterLambda_Vertex.Z();
    lambda_energy_lab = sqrt(lambda_momx_lab*lambda_momx_lab+lambda_momy_lab*lambda_momy_lab+lambda_momz_lab*lambda_momz_lab+MLambda*MLambda);
      
    vlambx_Lab = 0.0;
    vlamby_Lab = 0.0;
    vlambz_Lab = 0.0;

    plamb_Lab = sqrt(lambda_momx_lab*lambda_momx_lab+lambda_momy_lab*lambda_momy_lab+lambda_momz_lab*lambda_momz_lab);
    
    pphotonx_Lab = qVirtual_Vertex.X();
    pphotony_Lab = qVirtual_Vertex.Y();
    pphotonz_Lab = qVirtual_Vertex.Z();
    EphotonE_Lab = sqrt(pphotonx_Lab*pphotonx_Lab+pphotony_Lab*pphotony_Lab+pphotonz_Lab*pphotonz_Lab+photonmass*photonmass);
    
    vphotonx_Lab = 0.0;
    vphotony_Lab = 0.0;
    vphotonz_Lab = 0.0;

    k_momx_lab = pScatterKaon_Vertex.X();
    k_momy_lab = pScatterKaon_Vertex.Y();
    k_momz_lab = pScatterKaon_Vertex.Z();
    k_energy_lab = sqrt(k_momx_lab*k_momx_lab+k_momy_lab*k_momy_lab+k_momz_lab*k_momz_lab+mKaon*mKaon);

    vkx_Lab = 0.0;
    vky_Lab = 0.0;
    vkz_Lab = 0.0;
    
    // For debugging purpose
    if (plamb_Lab > (PBeam/ABeam) ) { 
       // Unphysical kinematics  // if TDIS spectator momentum larger than 50% ion momentum
       //printf("impossible of TDIS spectator momentum= %6.2f \n", ppr_Lab);
       continue;
    }

    // for debuggin purpose
    // cout  << "(7L) kaon Vertex, px= " << k_momx_lab << ", py= " << k_momy_lab <<
    // ", pz= " << k_momz_lab << ", Ek= " << k_energy_lab << endl;
	       
    // Definition of kaon variable
    xL = (pSpectator_Vertex.Dot(kIncident_Vertex))/(PIncident_Vertex.Dot(kIncident_Vertex));
    //xL = gRandom->Uniform(1.);
    xk = xbj/(1 - xL);
    tk = (E_k*E_k) - (pScatterKaon_V3.Mag()*pScatterKaon_V3.Mag());
    yk = pScatterKaon_Rest.Dot(qVirtual_Rest)/(pScatterKaon_Rest.Dot(kIncident_Rest));
    
    // *****
    // for debugging purpose
    /*
      cout  << "(7R) kaon REST, px= " << pScatterKaon_Rest.X() << ", py= " << pScatterKaon_Rest.Y()<<
      ", pz= " << pScatterKaon_Rest.Z() << ", E= " << pScatterKaon_Rest.E() << ", M= " << pScatterKaon_Rest.M() << endl;
    */
	  
    if (xbj > 0.001 && xbj < 1.0){

      /*
	L = renormalization cut-off parameter... ex) typ = 1 ! 1.56GeV
	dis = charge state configuration (0=charge exchange, 1=neutral exchange)
	FLAG = combination of splitting function/PDF (0:\pi, 1:\rho, 2:\pi\Delta)
      */
      //  typ = 1 ! s-exp For factor (L=,dis=0, FLAG=0) | J = 0 + 1/2
      if (typ == 1){
	    fk = fkfac*f2kp(P_k, xbj, theta_p2/D2R);
      }

      //  typ = 2 ! t-exp For factor (L=,dis=0, FLAG=0) | J = 0 + 1/2
      if (typ == 2){
	    fk = fkfac*f2kptmono(P_k, xbj, theta_p2/D2R);
      }

      //  typ = 3 ! ZEUS parameterization with F2L ("lambda" SF)      
      if (typ == 3){
	    fk = fkfac*f2kZEUS(xbj,invts.Q2);
      }

      yy = (P_k/MSpectator)*cos(theta_p2)+(1/MSpectator)*(MSpectator-sqrt(MSpectator*MSpectator+P_k*P_k));
      kTk = P_k * sin(theta_p2);
      k_int = fykL(yy, kTk, L, type);
	      		

    }
    // ***
    // For debugging purpose
    /*
      cout << "(8) " <<"fk= " << fk <<  ", 2nd proton Pz=" << p2_pt  << ", " << Pz_p2  << ", P_k= "
      << P_k << ", xbj= "<< xbj << ", theta_p2(deg)= "<< theta_p2 << endl; 
    */

    invts.tSpectator = MIon*MIon+MSpectator*MSpectator - 2.*pSpectator_Vertex.Dot(PIncident_Vertex);
    TDIS_t = invts.tSpectator;
    
    tmin = (((1-xL)/xL)*((1-xL)/xL))*(MSpectator*MSpectator-xL*MIon*MIon);// tmin is the minimumn kinematically allowed value of |t| for a given xL and thus positive

    //TDIS_t = -((xL*p_energy_incident*pScatterLambda_Vertex.Theta()*xL*p_energy_incident*pScatterLambda_Vertex.Theta())/xL)-tmin;
		
    if (TDIS_t < tmin && TDIS_t > -tmin){
      continue;
    }

    //cout << "---->t compare (invts vs TDIS) =" << invts.tSpectator << " " << TDIS_t << endl;
    // invts.tPrime     = 2.*pSpectator_Vertex.Dot(PIncident_Vertex) - MIon*MIon; // invts.tPrime = -MSpectator??
    invts.tPrime     = invts.tSpectator - tmin;
    p_ST = pSpectator_Vertex -(pDotq/qMag/qMag)*qVirtual_Rest;		
    invts.nu = invts.Q2 / (2 * MIon * invts.xBj) ;
    invts.W = MIon*MIon + invts.Q2/invts.xBj -invts.Q2;
    W2 = invts.W*invts.W;

        //  **********************************************************
    //
    //    Cross-section assignment
    //
    //  **********************************************************		
    // From the DIS event generator

    if (xbj > 0.001 && xbj < 1.0){

      f2N   = F2N(invts.xBj, invts.Q2, inucl);  // F2N is define at G4SBSDIS.h
      xsec_inclusive_dis    = cdissigma_p(xbj,invts.y_D,invts.Q2,invts.nu,kScattered_Rest.E());
      // xsec_inclusive_dis    = cdissigma_p(xbj,y_inelasticity,invts.Q2,invts.nu,kScattered_Rest.E());
      xsec_tagged_dis   = xsec_inclusive_dis * (fk/f2N); 
    }

    // ***
    // For debugging purpose
    /*
      if(fk>0.){
      cout << " what is f2N= "  << f2N  << ", fk=" << fk << ", momentum = " << P_k << ", theta = " << theta_p2/D2R << endl;
      }
      cout << "(10) " <<"xsec_inclusive_dis=  " << xsec_inclusive_dis  <<  ", xsec_tagged_dis= " << xsec_tagged_dis
      << ",  f2N= " <<  f2N << endl;
      cout << "kIncident_Rest.E()= "<< kIncident_Rest.E() << ", acos( kScattered_Rest.Z()/kScattered_Rest.Mag())= "
      <<kScattered_Rest.Z() <<", " << kScat3vec.Mag() << ", kScattered_Rest.E()= " <<kScattered_Rest.E()
      << ", nanobarn= "<< nanobarn<< endl;
      cout << "(11) " <<"invts.nu=  " << qVirtual_Rest.E()/MIon << ", " << invts.nu<< ",  kIncident_Rest.E()= "
      << kIncident_Rest.E() << ",  incoming electron energy= "<< invts.nu*2. + kScattered_Rest.E() << endl; 
    */
	  
    // qVirtual_Rest.E(): the virtual 4- momentum transfer between electron and deuteron in rest frame
    // invts.nu : the virtual 4- momentum transfer between electron and target nucleon in rest frame
	  
    //  invts.nu*2. + kScattered_Rest.E() : initial electron energy in rest frame		

    e_momx_incident = kIncident_Vertex.X();
    e_momy_incident = kIncident_Vertex.Y();
    e_momz_incident = kIncident_Vertex.Z();
    
    // p_energy_incident = PIncident_Vertex.E();
    e_energy_incident = sqrt(e_momx_incident*e_momx_incident+e_momy_incident*e_momy_incident+e_momz_incident*e_momz_incident+mElectron*mElectron);

    p_momx_incident = PIncident_Vertex.X();
    p_momy_incident = PIncident_Vertex.Y();
    p_momz_incident = PIncident_Vertex.Z();
    
    // p_energy_incident = PIncident_Vertex.E();
    p_energy_incident = sqrt(p_momx_incident*p_momx_incident+p_momy_incident*p_momy_incident+p_momz_incident*p_momz_incident+MProton*MProton);

    PX_Vertex = kIncident_Vertex+PIncident_Vertex-kScattered_Vertex-pSpectator_Vertex;
	  
    // store the electron 4-vector in rest frame
    pex_Res = kScat3vec.X() ;
    pey_Res = kScat3vec.Y() ;
    pez_Res = kScat3vec.Z() ;
    EeE_Res = sqrt(kScat3vec.Mag2()+mElectron*mElectron);
    vex_Res = 0.0;
    vey_Res = 0.0;
    vez_Res = 0.0;
    
    // store the electron 4-vector in Lab frame
    e_momx_lab = kScattered_Vertex.X() ;
    e_momy_lab = kScattered_Vertex.Y() ;
    e_momz_lab = kScattered_Vertex.Z() ;
    
    e_energy_lab = sqrt(e_momx_lab*e_momx_lab+e_momy_lab*e_momy_lab+e_momz_lab*e_momz_lab+mElectron*mElectron);
    vex_Lab = 0.0;
    vey_Lab = 0.0;
    vez_Lab = 0.0;

    // HERE!!!
    /*
    if ((kScattered_Vertex.Theta() < 1.5) || (kScattered_Vertex.Theta() > 1.8)){
      continue;
    }
    */

    /************************************************************************************
     * 
     * GEANT4 input
     * Legacy LUND format (https://gemc.jlab.org/gemc/html/documentation/generator/lund.html)
     * 
     * **********************************************************************************/

    OUT_lund << setiosflags(ios::left)  << setiosflags(ios::fixed)  <<"                 "  <<  NumPtls << " \t " <<  scientific  << invts.xBj << " \t " << invts.Q2  << " \t " << invts.s_e  << " \t " << "1.0" << " \t " << xk << " \t" << yk << " \t"  << tk  << " \t"  <<  " \t" << xsec_inclusive_dis << " \t" << xsec_tagged_dis << endl;
      
    // incoming particles
    OUT_lund << setiosflags(ios::left) << setiosflags(ios::fixed) << "\t" << "1" << " \t " << e_particle_charge << " \t " << "0" << " \t " << e_particle_id << " \t " << "0" <<  " \t "<< "1" << " \t "<< scientific << kBeamMCx << " \t " << kBeamMCy << " \t " << kBeamMCz << " \t " << kBeamMC << " \t " << emass << " \t " << ve0X_Lab  << " \t " << ve0Y_Lab << " \t " << ve0Z_Lab << endl; 

    OUT_lund << setiosflags(ios::left) << setiosflags(ios::fixed) <<  "\t" << "2" << " \t " << pr_particle_charge << " \t " << "0" << " \t " << pr_particle_id << " \t " << "0" <<  " \t "<< "1" << " \t "<< scientific << PBeamMCx+PBeamMC*sin(CrossingTheta)*cos(CrossingPhi) << " \t " << PBeamMCy+PBeamMC*sin(CrossingTheta)*sin(CrossingPhi) << " \t " << PBeamMC*cos(CrossingTheta) << " \t " << PBeamMC << " \t " << MProton << " \t " << vp0X_Lab  << " \t " << vp0Y_Lab << " \t " << vp0Z_Lab << endl; 

    // daughter particles
    // the scattered electron
    OUT_lund << setiosflags(ios::left) << setiosflags(ios::fixed) << "\t" << "3" << " \t " << e_particle_charge << " \t " << "1" << " \t " << e_particle_id << " \t " << "0" <<  " \t "<< "1" <<  " \t "<< scientific << e_momx_lab << " \t " << e_momy_lab << " \t " << e_momz_lab << " \t " << e_energy_lab << " \t " << emass << " \t " << vex_Lab  << " \t " << vey_Lab << " \t " << vez_Lab << endl; 
    // final state kaon
    OUT_lund << setiosflags(ios::left) << setiosflags(ios::fixed) <<  "\t" << "4" << " \t " << k_particle_charge << " \t " << "1" << " \t " << k_particle_id << " \t " << "0" <<  " \t "<< "1" <<  " \t "<< scientific << k_momx_lab << " \t " << k_momy_lab << " \t " << k_momz_lab << " \t " << k_energy_lab << " \t " << mKaon << " \t " << vkx_Lab  << " \t " << vky_Lab << " \t " << vkz_Lab << endl; 
    // final state lambda (TDIS)
    OUT_lund << setiosflags(ios::left) << setiosflags(ios::fixed) <<  "\t" << "5" << " \t " << sp_particle_charge << " \t " << "1" << " \t " << sp_particle_id << " \t " << "0" <<  " \t "<< "1" <<  " \t "<< scientific << lambda_momx_lab << " \t " << lambda_momy_lab << " \t " << lambda_momz_lab << " \t " << lambda_energy_lab << " \t " << spmass << " \t " << vlambx_Lab  << " \t " << vlamby_Lab << " \t " << vlambz_Lab << endl;  

    /************************************************************************************
     * 
     * Fun4All inputs (https://wiki.bnl.gov/eicug/images/6/6b/Fun4all_format_help.pdf)
     * Pythia6 format (https://eic.github.io/software/pythia6.html#output-file-structure)
     * Example (https://github.com/eic/eicsmeardetectors/blob/master/tests/simple_gen.txt)
     * 
     * **********************************************************************************/
    OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "0" << "\t" << MEvts << "\t" << "1" << "\t" << endl;
    OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "============================================" << endl;
      
    // incoming particles
    OUT_pythia << setiosflags(ios::left) << setiosflags(ios::fixed) << "1" << " \t " << "21" << " \t " << e_particle_id << " \t " << "0" << " \t " << "3" <<  " \t "<< "4" << " \t " << e_momx_incident << " \t " << e_momy_incident << " \t " << e_momz_incident << " \t " << e_energy_incident << " \t " << emass << " \t " << ve0X_Lab  << " \t " << ve0Y_Lab << " \t " << ve0Z_Lab << endl; 

    // Crossing angle 
    //OUT_pythia << setiosflags(ios::left) << setiosflags(ios::fixed) << "2" << " \t " << "21" << " \t " << pr_particle_id << " \t " << "0" << " \t " << "5" <<  " \t "<< "6" << " \t "<< PBeamMCx+PBeamMC*sin(CrossingTheta)*cos(CrossingPhi) << " \t " << PBeamMCy+PBeamMC*sin(CrossingTheta)*sin(CrossingPhi) << " \t " << PBeamMCz*cos(CrossingTheta) << " \t " << PBeamMC << " \t " << MProton << " \t " << vp0X_Lab  << " \t " << vp0Y_Lab << " \t " << vp0Z_Lab << endl; 
    OUT_pythia << setiosflags(ios::left) << setiosflags(ios::fixed) << "2" << " \t " << "21" << " \t " << pr_particle_id << " \t " << "0" << " \t " << "5" <<  " \t "<< "6" << " \t "<< p_momx_incident+PBeamMC*sin(CrossingTheta)*cos(CrossingPhi) << " \t " << p_momy_incident+PBeamMC*sin(CrossingTheta)*sin(CrossingPhi) << " \t " << p_momz_incident*cos(CrossingTheta) << " \t " << p_energy_incident << " \t " << MProton << " \t " << vp0X_Lab  << " \t " << vp0Y_Lab << " \t " << vp0Z_Lab << endl; 

    // daughter particles
    // Virtual photon
    OUT_pythia << setiosflags(ios::left) << setiosflags(ios::fixed) << "3" << " \t " << "21" << " \t " << photon_particle_id << " \t " << "1" << " \t " << "0" <<  " \t "<< "0" <<  " \t " << pphotonx_Lab << " \t " << pphotony_Lab << " \t " << pphotonz_Lab << " \t " << EphotonE_Lab << " \t " << photonmass << " \t " << vphotonx_Lab  << " \t " << vphotony_Lab << " \t " << vphotonz_Lab << endl; 
    // the scattered electron
    OUT_pythia << setiosflags(ios::left) << setiosflags(ios::fixed) << "4" << " \t " << "1" << " \t " << e_particle_id << " \t " << "1" << " \t " << "0" <<  " \t "<< "0" <<  " \t " << e_momx_lab << " \t " << e_momy_lab << " \t " << e_momz_lab << " \t " << e_energy_lab << " \t " << emass << " \t " << vex_Lab  << " \t " << vey_Lab << " \t " << vez_Lab << endl; 
        // final state neutron (TDIS)
    OUT_pythia << setiosflags(ios::left) << setiosflags(ios::fixed) << "5" << " \t " << "1" << " \t " << sp_particle_id << " \t " << "2" << " \t " << "0" <<  " \t "<< "0" <<  " \t " << lambda_momx_lab << " \t " << lambda_momy_lab << " \t " << lambda_momz_lab << " \t " << lambda_energy_lab << " \t " << spmass << " \t " << vlambx_Lab  << " \t " << vlamby_Lab << " \t " << vlambz_Lab << endl;  
    // final state pion
    OUT_pythia << setiosflags(ios::left) << setiosflags(ios::fixed) << "6" << " \t " << "1" << " \t " << k_particle_id << " \t " << "2" << " \t " << "0" <<  " \t "<< "0" <<  " \t " << k_momx_lab << " \t " << k_momy_lab << " \t " << k_momz_lab << " \t " << k_energy_lab << " \t " << mKaon << " \t " << vkx_Lab  << " \t " << vky_Lab << " \t " << vkz_Lab << endl; 
    OUT_pythia << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "=============== Event finished ===============" << endl;

    /************************************************************************************
     * 
     * HEPMC Format for ATHENA simulations (https://arxiv.org/pdf/1912.08005.pdf)
     * 
     * **********************************************************************************/
    OUT_hepmc << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "E" << " " << MEvts << " " << "1" << " " << NumPtls << endl;
    OUT_hepmc << setiosflags(ios::left)  << setiosflags(ios::fixed)  << "U" << " " << "GEV" << " " << "MM" << endl;
      
    // incoming particles
    OUT_hepmc << setiosflags(ios::left) << setiosflags(ios::fixed) << "P" << " " << "1"  << " " << "0" << " " << e_particle_id << " " << e_momx_incident << " " << e_momy_incident << " " << e_momz_incident << " " << e_energy_incident << " " << emass << " " << "4" << endl; 

    // Crossing angle 
    OUT_hepmc << setiosflags(ios::left) << setiosflags(ios::fixed) << "P" << " " << "2" << " " << "0" << " " << pr_particle_id << " " << p_momx_incident+PBeamMC*sin(CrossingTheta)*cos(CrossingPhi) << " " << p_momy_incident+PBeamMC*sin(CrossingTheta)*sin(CrossingPhi) << " " << p_momz_incident*cos(CrossingTheta) << " " << p_energy_incident << " " << MProton <<  " " << "4" << endl;

    // Vertex line
    OUT_hepmc << "V" << " " << "-1" << " " << "0" << " " << "[1,2]" << endl;

    // daughter particles
    // the scattered electron
    OUT_hepmc << setiosflags(ios::left) << setiosflags(ios::fixed) << "P" << " " << "3"  << " " << "-1" << " " << e_particle_id << " " << e_momx_lab << " " << e_momy_lab << " " << e_momz_lab << " " << e_energy_lab << " " << emass <<  " " << "1" << endl;
    // final state neutron (TDIS)
    OUT_hepmc << setiosflags(ios::left) << setiosflags(ios::fixed) << "P" << " " << "4"  << " " << "-1" << " " << sp_particle_id << " " << lambda_momx_lab << " " << lambda_momy_lab << " " << lambda_momz_lab << " " << lambda_energy_lab << " " << spmass <<  " " << "1" << endl;    
    // final state pion
    OUT_hepmc << setiosflags(ios::left) << setiosflags(ios::fixed) << "P" << " " << "5"  << " " << "-1" << " " << k_particle_id << " " << k_momx_lab << " " << k_momy_lab << " " << k_momz_lab << " " << k_energy_lab << " " << kmass <<  " " << "1" << endl;


    tree->Fill();
    proctree->Fill();
    itree->Fill();

    MEvts++;

  }
  
  fRoot.Write();
  fRoot.Close();
  printf("Total of %d events out of %d Trials \n",MEvts,NEvts);
  
  OUT_lund.close();
  OUT_pythia.close();

  return MEvts;
}
