/*
 * Description: Structure function definitions
 * ================================================================
 * Time-stamp: "2020-09-23 14:36:05 trottar"
 * ================================================================
 *
 * Author:  Kijun Park and Richard L. Trotta III <trotta@cua.edu>
 *
 * Copyright (c) trottar
 */


// subroutine to calculate the f2p as a function xbj
double f2p( double x ){
  double f2 = 0.0;
  double p0 = 0.1736;
  double p1 = 4.537;
  double p2 = -48.66;
  double p3 = 236.8;
  double p4 = -665.8;
  double p5 = 1094;
  double p6 = -973.9;
  double p7 = 363.2;

  f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5) + p6*pow(x,6) + p7*pow(x,7);
  
  return f2;
}

double F2N(double x, double Q2, int nucl){
  // nucl =2; // for neutron;

  // Define the DIS PDF from CTEQ directory:  cteq-tbls/ctq66m/ctq66.00.pds
  /* initcteqpdf(); */
  
  double qu = cteq_pdf_evolvepdf(__dis_pdf, 1, x, sqrt(Q2) );
  double qd = cteq_pdf_evolvepdf(__dis_pdf, 2, x, sqrt(Q2) );
  double qubar = cteq_pdf_evolvepdf(__dis_pdf, -1, x, sqrt(Q2) );
  double qdbar = cteq_pdf_evolvepdf(__dis_pdf, -2, x, sqrt(Q2) );

  double quv = qu-qubar;
  double qdv = qd-qdbar;

  double qs = cteq_pdf_evolvepdf(__dis_pdf, 3, x, sqrt(Q2) );

  double F2 = 0.0; 
  double e_u =  2.0/3.0;
  double e_d = -1.0/3.0;

  if( nucl == 1 ){
    F2 += x*( e_u*e_u*quv + e_d*e_d*qdv ); 
  }
  if( nucl == 2){
    F2 += x*( e_u*e_u*qdv + e_d*e_d*quv ); 
  }
  // Sea quarks
  F2  += x*(2.0*e_u*e_u*qubar + 2.0*e_d*e_d*(qdbar + qs));
  return F2;

}


// subroutine to calculate the f2pi as function of recoiled nucleon momentum, xbj, theta
// This is user parametrization by fit the Wally's codes 3Var_x.f() with integration of finite momentum range
// typ = 6 ! t-dependent expon FORM FACTOR
// dis = 1 !NEUTRAL EXCHANGE
// FLAG = 0  --- THE PION CONTRIBUTION      | J = 0 + 1/2
double f2pipaulvillas(double p, double x, double th){
  
  double p0, p1, p2, p3, p4, p5;
  int xflag = 0;
  double fk = 0.0;
  double fth = 0.0;

  if (p > 0.05 && p <= 0.1){
    p0 = 0.37229E-03;
    p1 = 0.69442E-02;
    p2 = -0.16933;
    p3 = 0.31326;
    p4 = 6.6281;
    p5 = -26.467;
    //    if (x < 0.0555 || x > 0.083) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = -0.954 + 66.5*p -1632.4*p*p + 14573.*p*p*p; 
  }

  if (p > 0.1 && p <= 0.2){
    p0 = 0.37229E-03;
    p1 = 0.69442E-02;
    p2 = -0.16933;
    p3 = 0.31326;
    p4 = 6.6281;
    p5 = -26.467;
    //    if (x < 0.0555 || x > 0.16) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = 0.464 -15.4*p + 126.5*p*p;
  }
  
  if (p > 0.2 && p <= 0.3){
    p0 = 0.42299E-03;
    p1 = 0.33955E-01;
    p2 = -0.64850;
    p3 = 4.3910;
    p4 = -13.446;
    p5 = 15.984;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    //     if (x < 0.0555 || x > 0.226) xflag =1;
    fk = -1.133 + 8.5354*p;
  }

  if (p > 0.3 && p <= 0.5){
    p0 = 0.89681E-03;
    p1 = 0.60363E-01;
    p2 = -0.88363;
    p3 = 4.7711;
    p4 = -11.958;
    p5 = 11.768;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    //     if (x < 0.0555 || x > 0.281) xflag = 1;
    fk = -1.345 + 9.47*p -7.91*p*p;
  }

  if (p < 0.05 || p > 0.5){
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    p3 = 0.0;
    p4 = 0.0;
    p5 = 0.0;
  }

  if (th < 1.8 || th > 74){
    fth = 0.0;}
  else{
    fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; }

  // Comment out for EIC
  /* if( xflag == 1 || x < 0.0555 || x > 0.3){ // Hall-A fixed target exp. limit */
  if( xflag == 1){
    return 0.0;}
  else{
    double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5);  
    double f2temp = f2*fk*fth;
    return f2temp;
  }
  
  
 }

// subroutine to calculate the f2pi as function of recoiled nucleon momentum, xbj, theta
// This is user parametrization by fit the Wally's codes 3Var_x.f() with integration of finite momentum range
// typ = 5 ! t-dependent expon FORM FACTOR
// dis = 1 !NEUTRAL EXCHANGE
// FLAG = 0  --- THE PION CONTRIBUTION      | J = 0 + 1/2
double f2pitexp(double p, double x, double th){
  
  double p0, p1, p2, p3, p4, p5;
  int xflag = 0;
  double fk = 0.0;
  double fth = 0.0;

  if (p > 0.05 && p <= 0.1){
    p0 = 0.13151E-04;
    p1 = 0.24655E-02;
    p2 = -0.11011;
    p3 = 1.6664;
    p4 = -10.813;
    p5 = 25.662;
    //    if (x < 0.0555 || x > 0.083) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = -0.954 + 66.5*p -1632.4*p*p + 14573.*p*p*p; 
  }

  if (p > 0.1 && p <= 0.2){
    p0 = 0.21109E-03;
    p1 = 0.24422E-01;
    p2 = -0.65408;
    p3 = 6.1526;
    p4 = -25.771;
    p5 = 41.141;
    //    if (x < 0.0555 || x > 0.16) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = 0.464 -15.4*p + 126.5*p*p;
  }
  
  if (p > 0.2 && p <= 0.3){
    p0 = 0.53755E-03;
    p1 = 0.42349E-01;
    p2 = -0.80465;
    p3 = 5.4029;
    p4 = -16.362;
    p5 = 19.186;
    //     if (x < 0.0555 || x > 0.226) xflag =1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = -1.133 + 8.5354*p;
  }

  if (p > 0.3 && p <= 0.5){
    p0 = 0.13141E-02;
    p1 = 0.89010E-01;
    p2 = -1.3003;
    p3 = 7.0225;
    p4 = -17.631;
    p5 = 17.393;
    //     if (x < 0.0555 || x > 0.281) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = -1.345 + 9.47*p -7.91*p*p;
  }

  if (p < 0.05 || p > 0.5){
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    p3 = 0.0;
    p4 = 0.0;
    p5 = 0.0;
  }

  if (th < 1.8 || th > 74){
    fth = 0.0;}
  else{
    fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; }

  // Comment out for EIC
  /* if( xflag == 1 || x < 0.0555 || x > 0.3){ */
  if( xflag == 1){
    return 0.0;}
  else{
    double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5);  
    double f2temp = f2*fk*fth;
    return f2temp;
  }
  
  
}

// subroutine to calculate the f2pi as function of recoiled nucleon momentum, xbj, theta
// This is user parametrization by fit the Wally's codes 3Var_x.f() with integration of finite momentum range
// typ = 3 ! COV DIP FORM FACTOR
// dis = 1 !NEUTRAL EXCHANGE
// FLAG = 0  --- THE PION CONTRIBUTION      | J = 0 + 1/2
double f2picov(double p, double x, double th){
  
  double p0, p1, p2, p3, p4, p5;
  int xflag = 0;
  double fk = 0.0;
  double fth = 0.0;

  if (p > 0.05 && p <= 0.1){
    p0 = 0.30433e-04;
    p1 = 0.60118e-03;
    p2 = -0.30863e-01;
    p3 = 0.74473e-01;
    p4 = 4.2508;
    p5 = -28.398;
    //    if (x < 0.0555 || x > 0.083) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = -0.954 + 66.5*p -1632.4*p*p + 14573.*p*p*p; 
  }

  if (p > 0.1 && p <= 0.2){
    p0 = 0.43685e-03;
    p1 = 0.83395e-02;
    p2 = -0.20040;
    p3 = 0.36447;
    p4 = 7.8242;
    p5 = -31.101;
    //    if (x < 0.0555 || x > 0.16) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = 0.464 -15.4*p + 126.5*p*p;
  }
  
  if (p > 0.2 && p <= 0.3){
    p0 = 0.60436e-03;
    p1 = 0.45603e-01;
    p2 = -0.86055;
    p3 = 5.6867;
    p4 = -16.812;
    p5 = 19.095;
    //     if (x < 0.0555 || x > 0.226) xflag =1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = -1.133 + 8.5354*p;
  }

  if (p > 0.3 && p <= 0.5){
    p0 = 0.16068e-02;
    p1 = 0.10964;
    p2 = -1.5980;
    p3 = 8.6327;
    p4 = -21.714;
    p5 = 21.475;
    //     if (x < 0.0555 || x > 0.281) xflag = 1;
    if ( x > 0.3) xflag = 1;
    fk = -1.345 + 9.47*p -7.91*p*p;
  }

  if (p < 0.05 || p > 0.5){
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    p3 = 0.0;
    p4 = 0.0;
    p5 = 0.0;
  }

  if (th < 1.8 || th > 74){
    fth = 0.0;}
  else{
    fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; }

  // Comment out for EIC
  /* if( xflag == 1 || x < 0.0555 || x > 0.3){ */
  if( xflag == 1){
    return 0.0;}
  else{
    double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5);  
    double f2temp = f2*fk*fth;
    return f2temp;
  }
  
  
}

// subroutine to calculate the f2pi as function of recoiled nucleon momentum, xbj, theta (Timothy J. Hobbs)
// what is "th" = theta angle
//
// This is user parametrization by fit the Wally's codes 3Var_x.f() with integration of finite momentum range
// typ = 2 !EXPONENTIAL FORM FACTOR (s-depdendent exp)
// dis = 1 !NEUTRAL EXCHANGE
// FLAG = 0  --- THE PION CONTRIBUTION      | J = 0 + 1/2
double f2pi(double p, double x, double th){
  
  double p0, p1, p2, p3, p4, p5;
  int xflag = 0;
  double fk = 0.0;
  double fth = 0.0;

  if (p > 0.05 && p <= 0.1){
    p0 = 3.656e-5;
    p1 = -0.000402;
    p2 = -0.008886;
    p3 = 0.07359;
    p4 = 1.079;
    p5 = -8.953;
    //if (x < 0.0555 || x > 0.083) xflag = 1;
    //    fk = -0.954 + 66.5*p -1632.4*p*p + 14573*p*p*p; // original value from Dasuni
    fk = -0.2473 + 14.79*p -283.55*p*p + 1766.9*p*p*p;   // fitting 3Var_k.f  x1000
  }

  if (p > 0.1 && p <= 0.2){
    p0 = 0.000287;
    p1 = 0.009397;
    p2 = -0.2632;
    p3 = 2.029;
    p4 = -5.878;
    p5 = 4.664;
    //if (x < 0.0555 || x > 0.16) xflag = 1;
    fk = 0.464 -15.4*p + 126.5*p*p;
  }
  
  if (p > 0.2 && p <= 0.3){
    p0 = 0.0003662;
    p1 = 0.02811;
    p2 = -0.4566;
    p3 = 2.428;
    p4 = -5.107;
    p5 = 3.222;
    // if (x < 0.0555 || x > 0.226) xflag =1;
    fk = -1.133 + 8.5354*p;
  }

  if (p > 0.3 && p <= 0.5){
    p0 = 0.0009412;
    p1 = 0.01366;
    p2 = -0.1744;
    p3 = 0.3864;
    p4 = 0.6615;
    p5 = -2.113;
    // if (x < 0.0555 || x > 0.281) xflag = 1;
    fk = -1.345 + 9.47*p -7.91*p*p;
  }

  if (p < 0.05 || p > 0.5){
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    p3 = 0.0;
    p4 = 0.0;
    p5 = 0.0;
  }

  if (th < 1.8 || th > 74){
    fth = 0.0;}
  else{
    fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; }
  
  //  if( xflag == 1 || x < 0.0555 || x > 0.3){
  if( xflag == 1 ){      // EIC collider no limit
    return 0.0;}
  else{
    double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5);  
    double f2temp = f2*fk*fth;
    return f2temp;
  }
  
  
}

// using the ZEUS parameterization with f2
double f2piZEUS(double x, double Q2, int nucl){

  return F2N(x,Q2,nucl)*0.361; // see fig 18 HERA paper for reference
}


double F2L(double x, double Q2){

  // Define the DIS PDF from CTEQ directory:  cteq-tbls/ctq66m/ctq66.00.pds
  /* initcteqpdf(); */
  
  double qu = cteq_pdf_evolvepdf(__dis_pdf, 1, x, sqrt(Q2) );
  double qd = cteq_pdf_evolvepdf(__dis_pdf, 2, x, sqrt(Q2) );
  double qs = cteq_pdf_evolvepdf(__dis_pdf, 3, x, sqrt(Q2) );
  double qubar = cteq_pdf_evolvepdf(__dis_pdf, -1, x, sqrt(Q2) );
  double qdbar = cteq_pdf_evolvepdf(__dis_pdf, -2, x, sqrt(Q2) );
  double qsbar = cteq_pdf_evolvepdf(__dis_pdf, -3, x, sqrt(Q2) );

  double quv = qu-qubar;
  double qdv = qd-qdbar;
  double qsv = qs-qsbar;


  double F2 = 0.0; 
  double e_u =  2.0/3.0;
  double e_d = -1.0/3.0;
  double e_s = -1.0/3.0;

  // Lambda quark composition
  F2 += x*( e_u*e_u*quv + e_s*e_s*qsv ); 
  
  // Sea quarks
  F2  += x*(2.0*e_u*e_u*qubar + 2.0*e_s*e_s*(qsbar + qd));
  return F2;

}

// using the "ZEUS" parameterization with f2 of lambda
double f2kZEUS(double x, double Q2){

  return F2L(x,Q2)*0.361;
}

// Testing K. Park  Nov. 08 2016
// subroutine to calculate the f2kp as function of recoiled nucleon momentum, xbj, theta
// This is user parametrization by fit the Wally's codes 3KVar_x.f() with integration of finite momentum range
// typ = 2 ! s-exp For factor
// dis = 0 ! Charge EXCHANGE
// FLAG = 0  --- THE Kaon CONTRIBUTION      | J = 0 + 1/2
double f2kp(double p, double x, double th){
  
  double p0, p1, p2, p3, p4, p5;
  double fk = 0.0;
  double fth = 0.0;

  if (p > 0.05 && p <= 0.1){
      p0 = 0.83443E-07;
      p1 = -0.13408E-05;
      p2 = 0.22470E-03;
      p3 = -0.75778E-02;
      p4 = 0.85898E-01;
      p5 = -0.32083;
      fk = -0.954 + 66.5*p -1632.4*p*p + 14573.*p*p*p;
      // cout << "--> momentum region 1: "  << endl;  
  }

  if (p > 0.1 && p <= 0.2){
      p0 = 0.32072E-05;
      p1 = 0.50358E-03;
      p2 = -0.76921E-02;
      p3 = 0.23637E-01;
      p4 = 0.87552E-01;
      p5 = -0.39142;
      fk = 0.464 -15.4*p + 126.5*p*p;
      // cout << "--> momentum region 2: "  << endl;
  }
  
  if (p > 0.2 && p <= 0.3){
      p0 = 0.26275E-04;
      p1 = 0.30344E-02;
      p2 = -0.38997E-01;
      p3 = 0.15174;
      p4 = -0.16257;
      p5 = -0.80353E-01;
      fk = -1.133 + 8.5354*p;
      // cout << "--> momentum region 3: "  << endl;
  }

  if (p > 0.3 && p <= 0.5){
    p0 = 0.18012E-03;
    p1 = 0.18799E-01;
    p2 = -0.20900;
    p3 = 0.88713;
    p4 = -1.8633;
    p5 = 1.7046;
    fk = -1.345 + 9.47*p -7.91*p*p;
    // cout << "--> momentum region 4: "  << endl;
  }


  if (p > 0.5 && p <= 1.){
    p0 = 0.12115E-02;
    p1 = 0.72634E-01;
    p2 = -0.60468;
    p3 = 2.4920;
    p4 = -6.1261;
    p5 = 6.3972;
    fk = -1.345 + 9.47*p -7.91*p*p; // same as p=0.3-0.5 GeV/c 
    // cout << "--> momentum region 5: "  << endl;
  }

  
  if (p > 1.0 && p <= 5.){
    p0 = 0.50862E-03;
    p1 = 0.30568E-01;
    p2 = -0.23605;
    p3 = 0.98762;
    p4 = -2.4529;
    p5 =  2.5330;
    fk = -1.345 + 9.47*p -7.91*p*p; // same as p=0.3-0.5 GeV/c 
    // cout << "--> momentum region 6: "  << endl;
  }



  if (p > 5.0 && p <= 100.){
    p0 = 0.20262E-06;
    p1 = 0.25285E-04;
    p2 = -0.45309E-03;
    p3 = 0.28817E-02;
    p4 = -0.80264E-02;
    p5 =  0.83156E-02;
    fk = -1.345 + 9.47*p -7.91*p*p; // same as p=0.3-0.5 GeV/c 
    // cout << "--> momentum region 7: "  << endl;
  }
  
  if (p < 0.05){
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    p3 = 0.0;
    p4 = 0.0;
    p5 = 0.0;
    // cout << "--> momentum region 0: "  << endl;
  }
  

  if (th < 1.8 || th > 74){
    fth = 0.0;}
  else{
    fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; }
  
  double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5);  
  double f2temp = f2*fk*fth;
  return f2temp;
  
  
}


// typ = 3 ! t-exp For factor
// dis = 0 ! Charge EXCHANGE
// FLAG = 0  --- THE Kaon CONTRIBUTION      | J = 0 + 1/2

double f2kptmono(double p, double x, double th){
  
  double p0, p1, p2, p3, p4, p5;
  double fk = 0.0;
  double fth = 0.0;

  if (p > 0.05 && p <= 0.1){
    p0 = 0.20488E-05;
    p1 = 0.63463E-03;
    p2 = -0.27479E-01; 
    p3 = 0.41698;
    p4 = -2.7294;
    p5 = 6.5210;
    fk = -0.954 + 66.5*p -1632.4*p*p + 14573.*p*p*p; 
  }

  if (p > 0.1 && p <= 0.2){
      p0 = 0.41053E-04;
      p1 = 0.55437E-02;
      p2 = -0.14278;
      p3 = 1.3096;
      p4 = -5.3818;
      p5 = 8.4883;
      fk = 0.464 -15.4*p + 126.5*p*p;
  }
  
  if (p > 0.2 && p <= 0.3){
      p0 = 0.16415E-03;
      p1 = 0.13833E-01;
      p2 = -0.25650;
      p3 = 1.6955;
      p4 = -5.0825;
      p5 = 5.9353;
      fk = -1.133 + 8.5354*p;
  }

  if (p > 0.3 && p <= 0.5){
      p0 = 0.69820E-03;
      p1 = 0.48810E-01;
      p2 = -0.70672;
      p3 = 3.8132;
      p4 = -9.6119;
      p5 = 9.5429;
      fk = -1.345 + 9.47*p -7.91*p*p;
  }

  if (p < 0.05 || p > 0.5){
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    p3 = 0.0;
    p4 = 0.0;
    p5 = 0.0;
  }

  if (th < 1.8 || th > 74){
    fth = 0.0;}
  else{
    fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; }

  double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5);  
  double f2temp = f2*fk*fth;
  return f2temp;  
  
 }


// ==============================================================================
/*  
 *
 *   START OF SPLITTING FUNCTION DEFINITIONS
 *
 */
// ==============================================================================

// ******************************************************************************
/* FUNCTION GIVING f(y) FOR N-PION VERTEX, WHERE
 *    y IS THE IMF MOMENTUM OF THE INTERMEDIATE MESON.
 *
 *  WRITTEN: W. Melnitchouk (1999)
 *  MODIFIED: T. HOBBS (2013)
 */
// ******************************************************************************
double fypiN(double y,double kT,double L,int typ){

  double ss,kT2,SpiN;
  double fypiN,t;
  double pi,mN,mpi=0,mP=0,g_piNN,gg=0,FF=0,sM;

  pi = 4*atan(1.0);
  mN  = 0.93891897;  //!masses in GeV!!
  //!!*** WE USE THE COUPLINGS INFERRED FROM HAIDENBAUER ET AL. ***
  g_piNN = sqrt(14.40 * 4*pi);                 //! g_{pi NN}
  gg = 2.0 * pow(g_piNN,2.) / (16.0 * pi*pi);   //!2 ISOSPIN FACTOR

  if(y<=0.0 || y>=1.0){
    fypiN = 0.;
    return fypiN;
  }
  kT2 = kT*kT;
  SpiN = (kT2 + mpi*mpi)/y + (kT2 + mP*mP)/(1.0-y);

  if(typ==0)       FF = ((L*L + mN*mN) / (L*L + SpiN));             //! monopole
  else if(typ==1)  FF = pow(((L*L + mN*mN) / (L*L + SpiN)),2.) ;     //! dipole
  else if(typ==2)  FF = exp( (mN*mN - SpiN)/(L*L) );                //! expon
  else if(typ==3) {
    t = (- kT2 - mN*mN*y*y) /(1.0-y);
    FF = pow(((L*L - mpi*mpi) / (L*L - t)),2.);                      //! cov dip
  }
  else if(typ==4){                       
    sM = (kT2 + (1.0+y)*mpi*mpi)/y +
      (kT2 + y*mP*mP)/(1.0-y) + mN*mN;
    FF = (pow(L,4.) + pow(mP,4.))/(pow(L,4.) + sM*sM); // ! DIPOLE -- s-channel Lambda exchange
  }

  // integrand
  ss = ( kT2 + pow((mP - (1.0-y) * mN),2.)) / (1.0-y) / pow(( (1.0-y)*(SpiN - mN*mN) ),2.) * FF*FF;

  fypiN =  gg * (1.0-y) / y * ss ;
  return fypiN;

}



// ******************************************************************************
/*  Function giving numerical value of f(y) for N-Lambda-K vertex
 *  y is l.c. momentum fraction on baryon.
 *  Written: W. Melnitchouk (1999)
 *  Modified: T. Hobbs (2012)
 ****************************************************************************
 */
double fykL(double y,double kT,double L,int typ){

  double ss,ss0,ikT,kT2,kTmax,kTint,SkL;
  double fykL,t;
  double pi,mN,mk=0,mL=0,g_kLN,gg=0,FF=0,sM;

  pi = 4*atan(1.0);
  mN  = 0.93891897;  //!masses in GeV!!

  mL = 1.1157; // Mass of LAMBDA^0
  mk = 0.4937; // Mass of K^+
  //!!*** WE USE THE COUPLINGS INFERRED FROM MUELLER-GROELING ET AL. ***
  g_kLN = sqrt(15.56 * 4*pi);    //! Mueller-Groeling - PS
  //g_kLN = sqrt(14.40 * 4*pi)    //! g_{pi NN}
  //g_kLN = sqrt(9. * 4*pi)       //! Navarra I
  //g_kLN = sqrt(3.61 * 4*pi)     //! Navarra II
  gg = g_kLN*g_kLN / (16. * pi*pi);
    
  ss0 = 0.;
  ikT=1;
  kTmax = 10.;
  kTint = kTmax/1000.;

  for(kT=kTint;kT<=kTmax;kT++){
    kT2 = kT*kT;
    SkL = (kT2 + mk*mk)/(1.-y) + (kT2 + mL*mL)/y;

    if(typ==0)       FF = ((L*L + mN*mN) / (L*L + SkL));             //! monopole
    else if(typ==1)  FF = pow(((L*L + mN*mN) / (L*L + kT)),2.) ;     //! dipole
    else if(typ==2)  FF = exp( (mN*mN - kT)/(L*L) );                //! expon
    else if(typ==3) {
      t = (- kT2 - mN*mN*y*y) /(1.0-y);
      FF = pow(((L*L - mk*mk) / (L*L - t)),2.);                      //! cov dip
    }
    else if(typ==4){                       
      sM = (kT2 + (1.0+y)*mk*mk)/y +
	(kT2 + y*mL*mL)/(1.0-y) + mN*mN;
      FF = (pow(L,4.) + pow(mL,4.))/(pow(L,4.) + sM*sM); // ! DIPOLE -- s-channel Lambda exchange
    }
    ss = ( kT2 + pow((mL - y*mN),2.)) / y / pow(( (1.0-y)*(SkL - mN*mN) ),2.) * FF*FF * (2*kT);

    if(ikT/2*2!=ikT) ss0 = ss0 + 4*ss;
    else if(ikT/2*2==ikT) ss0 = ss0 + 2*ss;
    else ikT = ikT +1; 
  }

  fykL =  gg * (1.0-y) / y * (kTint/3) * ss ;
  return fykL;

} 
