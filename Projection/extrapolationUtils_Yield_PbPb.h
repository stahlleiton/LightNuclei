#ifndef extrapolationUtils_Yield_PbPb_h
#define extrapolationUtils_Yield_PbPb_h


// Custom headers
#include "../Macros/fitUtils.h"


double yield_PtRapCent_PbPb(const double& pT, const double& rap, const double& cnt, const double& sqrtSNN, const double& mass, const int& charge)
{
  // check input parameters
  if (cnt<0 || pT<0) throw runtime_error("[ERROR] Negative centrality or pT value for yield!");
  if (cnt>1) throw runtime_error("[ERROR] Invalid centrality value for yield!");
  if (mass<=0.05) throw runtime_error("[ERROR] Invalid mass value for yield!");
  if (sqrtSNN<=0) throw runtime_error("[ERROR] Invalid beam energy value for yield!");
  if (charge==0) throw runtime_error("[ERROR] Can not describe neutral particles for yield!");

  // extract the mass dependent parameters
  double T_kin(0), beta(0), n(0), dNdy(0), beamEneRef(0);
  if (sqrtSNN>=3.5) {
    beamEneRef = 5.02;
    if (mass>0.05) {
      // source: Pion, Kaon and Proton spectra in PbPb 5.02 TeV, https://arxiv.org/pdf/1910.07678.pdf
      T_kin = 0.1138*cnt*cnt - 0.0152*cnt + 0.0919;
      beta  = -0.1332*cnt*cnt - 0.0353*cnt + 0.907;
      n     = 11.694*cnt*cnt*cnt*cnt - 13.632*cnt*cnt*cnt + 6.4735*cnt*cnt - 0.8855*cnt + 0.7612;
    }
    if (mass>0.9) {
      // source: Proton spectra in PbPb 5.02 TeV, https://arxiv.org/pdf/1910.07678.pdf
      dNdy = 1161.1*cnt*cnt*cnt*cnt*cnt*cnt - 3474.9*cnt*cnt*cnt*cnt*cnt + 4135.2*cnt*cnt*cnt*cnt - 2558.2*cnt*cnt*cnt + 1001*cnt*cnt - 339.25*cnt + 82.409;
    }
    if (mass>1.5) {
      // source: Deuteron spectra in PbPb 5.02 TeV, https://alice-figure.web.cern.ch/node/11817  ,  https://alice-figure.web.cern.ch/node/11380
      if (charge > 0) {
	T_kin = -0.17*cnt*cnt*cnt + 0.2455*cnt*cnt - 0.0096*cnt + 0.0802;
	beta  = -0.1255*cnt*cnt*cnt - 0.183*cnt*cnt - 0.0182*cnt + 0.8834;
	n     = 8.9736*cnt*cnt*cnt - 6.0804*cnt*cnt + 1.805*cnt + 0.5053;
	dNdy  = -0.0453*cnt*cnt*cnt*cnt + 0.0356*cnt*cnt*cnt + 0.1951*cnt*cnt - 0.3127*cnt + 0.1277;
      }
      else {
	T_kin = 0.3101*cnt*cnt*cnt - 0.3837*cnt*cnt + 0.1583*cnt + 0.075;
	beta  = -0.3768*cnt*cnt*cnt + 0.2109*cnt*cnt - 0.1215*cnt + 0.8876;
	n     = 14.314*cnt*cnt*cnt - 10.08*cnt*cnt + 2.8313*cnt + 0.4492;
	dNdy  = -0.0942*cnt*cnt*cnt*cnt + 0.1313*cnt*cnt*cnt + 0.1192*cnt*cnt - 0.2808*cnt + 0.1222;
      }
    }
    if (mass>2.5) {
      // source: Triton spectra in PbPb 5.02 TeV, ALI−PREL−329128 , doi:10.1088/1742-6596/1602/1/012022
      const auto cent = std::min(cnt, 0.94);
      T_kin = -0.1157*cent + 0.1136;
      beta  = 0.75;
      n     = 4.5616*cent*cent - 1.3877*cent + 0.4536;
      dNdy  = 0.00020844*cent*cent - 0.00039958*cent + 0.00019402;
    }
    if (std::abs(charge)==2) {
      if (mass>2.5) {
	// source: Helium3 spectra in PbPb 5.02 TeV, https://alice-figure.web.cern.ch/node/11477 , https://alice-figure.web.cern.ch/node/11523
	if (charge > 0) {
	  T_kin = 0.2829*cnt + 0.047;
	  beta  = -0.151061*cnt +0.894127;
	  n     = 7.3192*cnt*cnt - 0.1284*cnt + 0.6548;
	  dNdy  = -0.00034001*cnt + 0.00024267;
	}
	else {
	  T_kin = 0.0433*cnt + 0.0915;
	  beta  = -0.104430*cnt +0.869943;
	  n     = 13.596*cnt*cnt - 3.9116*cnt + 0.6269;
	  dNdy  = -0.00027404*cnt + 0.00020721;
	}
      }
      if (mass>3.0) {
	// source: Helium4 spectra in PbPb 5.02 TeV, ALI−PREL−327164 , doi:10.1088/1742-6596/1602/1/012022
	if (charge > 0) dNdy *= 4.9287878e-07/0.00017429285;
	else            dNdy *= 1.0560697e-06/0.00017298111;
      }
    }
  }
  else {
    beamEneRef = 2.76;
    if (mass>0.05) {
      // source: Pion, Kaon and Proton spectra in PbPb 2.76 TeV, 
      T_kin = 0.0731*cnt*cnt + 0.0025*cnt + 0.0959;
      beta  = -0.4301*cnt*cnt + 0.0338*cnt + 0.6463;
      n     = 1.787*cnt*cnt*cnt - 0.0372*cnt*cnt + 0.1798*cnt + 0.7076;
      dNdy  = 1367.4*cnt*cnt*cnt*cnt - 3743.1*cnt*cnt*cnt + 4501.1*cnt*cnt - 2895.1*cnt + 801.15;
    }
    else if (mass>0.9) {
      // source: Proton spectra in PbPb 2.76 TeV, 
      dNdy  = 87.663*cnt*cnt*cnt*cnt - 217.86*cnt*cnt*cnt + 232.95*cnt*cnt - 137.73*cnt + 37.211;
    }
    if (mass>1.5) {
      // source: Deuteron and Helium3 spectra in PbPb 2.76 TeV, Phys. Rev. C 93 (2016) 024917
      T_kin = -0.1861*cnt*cnt + 0.1907*cnt + 0.068;
      beta  = -0.153*cnt*cnt - 0.1056*cnt + 0.8746;
      n     = 3.4038*cnt*cnt - 1.0073*cnt + 0.8364;
      dNdy  = 0.1542*cnt*cnt - 0.2599*cnt + 0.1112;
    }
    if (mass>2.5) {
      // source: Helium3 spectra in PbPb 2.76 TeV, Phys. Rev. C 93 (2016) 024917
      dNdy *= 0.0014/150./0.0036;
    }
    if (mass>3.0) {
      // source: Helium4 spectra in PbPb 2.76 TeV, https://arxiv.org/pdf/1710.07531.pdf
      if (charge > 0) dNdy *= 7.83501324724194E-07/0.000299792795642062;
      else            dNdy *= 1.07922645539243E-06/0.000261347616124723; 
    }
  }

  // compute pT spectrum
  static TF1* fPt(0);
  if (!fPt) fPt = new TF1("fPtBlastWaveYield", blastWaveYield, 0., 100., 5);
  const auto res_pT = evalTF(*fPt, {1., mass, T_kin, beta, n}, pT);

  // compute y spectrum using data from -0.5 < y < 0.5
  // Based on: 10.1103/PhysRevLett.116.222302 (charged-particle at mid-rapidity)
  const auto width_ref = 2*0.23994808*TMath::Log(beamEneRef/0.000938);
  const auto yInt0p5 = std::sqrt(2*TMath::Pi()) * width_ref * TMath::Erf(std::sqrt(2)/width_ref/4); // from dNch/dy vs sNN
  const auto beamEneScale = TMath::Power(sqrtSNN, 2*0.155)/TMath::Power(beamEneRef, 2*0.155); // from dNch/deta vs sNN
  
  const auto width = 2*0.23994808*TMath::Log(sqrtSNN/0.000938);
  const auto res_y = beamEneScale * dNdy * TMath::Gaus(rap, 0., width) / yInt0p5;

  return res_y * res_pT;
};

double yield_PtRapCent_PbPb(const double* x, const double* p)
{
  // input variables
  const auto& t = p[3];
  const auto& pT  = t==0 ? x[0] : (t==1 ? p[4] : (t==2 ? x[0] : (t==3 ? x[0] : (t==4 ? x[0] : (t==5 ? p[4] : p[4])))));
  const auto& rap = t==0 ? x[1] : (t==1 ? x[0] : (t==2 ? p[4] : (t==3 ? x[1] : (t==4 ? p[4] : (t==5 ? x[0] : p[5])))));
  const auto& cnt = t==0 ? x[2] : (t==1 ? x[1] : (t==2 ? x[1] : (t==3 ? p[4] : (t==4 ? p[5] : (t==5 ? p[5] : x[0])))));
  if (t<0 || t>6) throw runtime_error("[ERROR] Wrong type value for yield3D!");
  
  // input parameters
  const auto& sqrtSNN = p[0];
  const auto& mass    = p[1];
  const auto& charge  = p[2];

  return yield_PtRapCent_PbPb(pT, rap, cnt, sqrtSNN, mass, charge);
};

double yield_PtEtaCent_PbPb(const double* x, const double* p)
{
  // input variables
  const auto& pT = x[0];
  const auto& eta = x[1];
  const auto& cnt = x[2];
  
  // input parameters
  const auto& sqrtSNN = p[0];
  const auto& mass    = p[1];
  const auto& charge  = p[2];

  const auto rap = computeRapidity(pT, eta, mass);

  return yield_PtRapCent_PbPb(pT, rap, cnt, sqrtSNN, mass, charge);
};

double spectra_PtRapCent_PbPb(const double * x, const double* p)
{
  const auto& pT = (p[3]>4 || p[3]==1) ? p[4] : x[0];
  return yield_PtRapCent_PbPb(x, p) / (2*TMath::Pi()*pT);
};

double dist_Inc2D_PbPb(TF1& f1D, const UChar_t& id, const double* x, const double* p)
{
  // input parameters
  const auto& scale   = p[0];
  const auto& type    = p[1]; // integrate over: (0) pT,    (1) rapidity,    (2) centrality,    (3) |rapidity|
  const auto& sqrtSNN = p[2];
  const auto& mass    = p[3];
  const auto& charge  = p[4];
  const auto& var1Min = p[5];
  const auto& var1Max = p[6];
  
  // check input parameters
  if (type<0 || type>3) throw runtime_error("[ERROR] Wrong type value for dist 2D!");
  if (scale<=0) throw runtime_error("[ERROR] Wrong scale value for dist 2D!");
  if (var1Min>=var1Max) throw runtime_error("[ERROR] Wrong range values for dist 2D!");
  
  double res(0);
  if (type<3) {
    res += evalTF(f1D, {sqrtSNN, mass, charge, type+4, x[0], x[1]}, var1Min, var1Max, id);
  }
  else if (type==3 && var1Min==0) {
    res += evalTF(f1D, {sqrtSNN, mass, charge, 5, x[0], x[1]}, -var1Max, var1Max, id);
  }
  else if (type==3 && var1Min>0) {
    res += evalTF(f1D, {sqrtSNN, mass, charge, 5, x[0], x[1]}, -var1Max, -var1Min, id);
    res += evalTF(f1D, {sqrtSNN, mass, charge, 5, x[0], x[1]}, var1Min, var1Max, id);
  }

  return scale * res;
};

double yield_Inc2D_PbPb(const double* x, const double* p)
{
  // input parameters
  const auto& var1Max = p[6];
  
  static TF1* f1D(0);
  if (!f1D) f1D = new TF1("f1DYield", yield_PtRapCent_PbPb, -var1Max, var1Max, 6);

  return dist_Inc2D_PbPb(*f1D, 0, x, p);
};

double yield_Diff2D_PbPb(const double* x, const double* p)
{
  // input parameters
  const auto& type    = p[1]; // integrate over: (0) pT,    (1) rapidity,    (2) centrality,    (3) |rapidity|
  const auto& var1Min = p[5];
  const auto& var1Max = p[6];
  
  auto binWidth = (var1Max - var1Min);
  if (type==3) binWidth *= 2;

  return yield_Inc2D_PbPb(x, p) / binWidth;
};

double spectra_Inc2D_PbPb(const double * x, const double* p)
{
  // input parameters
  const auto& var1Max = p[6];
  
  static TF1* f1D(0);
  if (!f1D) f1D = new TF1("f1DSpectra", spectra_PtRapCent_PbPb, -var1Max, var1Max, 6);

  return dist_Inc2D_PbPb(*f1D, 1, x, p);
};

double spectra_Diff2D_PbPb(const double * x, const double* p)
{
  // input parameters
  const auto& type    = p[1]; // integrate over: (0) pT,    (1) rapidity,    (2) centrality,    (3) |rapidity|
  const auto& var1Min = p[5];
  const auto& var1Max = p[6];
  
  auto binWidth = (var1Max - var1Min);
  if (type==3) binWidth *= 2;

  return spectra_Inc2D_PbPb(x, p) / binWidth;
};

double dist_Inc1D_PbPb(TF2& f2D, const UChar_t& id, const double* x, const double* p)
{
  // input parameters
  const auto& scale   = p[0];
  const auto& type    = p[1]; // function of:   (0) pT,     (1) rapidity,     (2) centrality,     (3) pT (integrate |y|),     (4) centrality (integrate |y|)
  const auto& sqrtSNN = p[2];
  const auto& mass    = p[3];
  const auto& charge  = p[4];
  const auto& var1Min = p[5];
  const auto& var1Max = p[6];
  const auto& var2Min = p[7];
  const auto& var2Max = p[8];
  
  // check input parameters
  if (type<0 || type>4) throw runtime_error("[ERROR] Wrong type value for dist 1D!");
  if (scale<=0) throw runtime_error("[ERROR] Wrong scale value for dist 1D!");
  if (var1Min>=var1Max) throw runtime_error("[ERROR] Wrong var1 range values for dist 1D!");
  if (var2Min>=var2Max) throw runtime_error("[ERROR] Wrong var2 range values for dist 1D!");
  
  double res(0);
  if (type<3) {
    res += evalTF(f2D, {sqrtSNN, mass, charge, type+1, x[0]}, var1Min, var1Max, var2Min, var2Max, id);
  }
  else if (type==3 && var1Min==0) {
    res += evalTF(f2D, {sqrtSNN, mass, charge, 1, x[0]}, -var1Max, var1Max, var2Min, var2Max, id);
  }
  else if (type==3 && var1Min>0) {
    res += evalTF(f2D, {sqrtSNN, mass, charge, 1, x[0]}, -var1Max, -var1Min, var2Min, var2Max, id);
    res += evalTF(f2D, {sqrtSNN, mass, charge, 1, x[0]}, var1Min, var1Max, var2Min, var2Max, id);
  }
  else if (type==4 && var2Min==0) {
    res += evalTF(f2D, {sqrtSNN, mass, charge, 3, x[0]}, var1Min, var1Max, -var2Max, var2Max, id);
  }
  else if (type==4 && var2Min>0) {
    res += evalTF(f2D, {sqrtSNN, mass, charge, 3, x[0]}, var1Min, var1Max, -var2Max, -var2Min, id);
    res += evalTF(f2D, {sqrtSNN, mass, charge, 3, x[0]}, var1Min, var1Max, var2Min, var2Max, id);
  }

  return scale * res;
};

double yield_Inc1D_PbPb(const double * x, const double* p)
{
  // input parameters
  const auto& var1Max = p[6];
  const auto& var2Max = p[8];
  
  static TF2* f2D(0);
  if (!f2D) f2D = new TF2("f2DYield", yield_PtRapCent_PbPb, -var1Max, var1Max, -var2Max, var2Max, 5);

  return dist_Inc1D_PbPb(*f2D, 0, x, p);
};

double yield_Diff1D_PbPb(const double * x, const double* p)
{
  // input parameters
  const auto& type    = p[1]; // function of:   (0) pT,     (1) rapidity,     (2) centrality,     (3) pT (integrate |y|),     (4) centrality (integrate |y|)
  const auto& var1Min = p[5];
  const auto& var1Max = p[6];
  const auto& var2Min = p[7];
  const auto& var2Max = p[8];
  
  auto binWidth = (var2Max - var2Min) * (var1Max - var1Min);
  if (type>2) binWidth *= 2;

  return yield_Inc1D_PbPb(x, p) / binWidth;
};

double spectra_Inc1D_PbPb(const double * x, const double* p)
{
  // input parameters
  const auto& var1Max = p[6];
  const auto& var2Max = p[8];
  
  static TF2* f2D(0);
  if (!f2D) f2D = new TF2("f2DSpectra", spectra_PtRapCent_PbPb, -var1Max, var1Max, -var2Max, var2Max, 5);

  return dist_Inc1D_PbPb(*f2D, 1, x, p);
};

double spectra_Diff1D_PbPb(const double * x, const double* p)
{
  // input parameters
  const auto& type    = p[1]; // function of:   (0) pT,     (1) rapidity,     (2) centrality,     (3) pT (integrate |y|),     (4) centrality (integrate |y|)
  const auto& var1Min = p[5];
  const auto& var1Max = p[6];
  const auto& var2Min = p[7];
  const auto& var2Max = p[8];
  
  auto binWidth = (var2Max - var2Min) * (var1Max - var1Min);
  if (type>2) binWidth *= 2;

  return spectra_Inc1D_PbPb(x, p) / binWidth;
};

double dist_Inc_PbPb(TF3& f3D, const UChar_t& id, const std::array<double, 6> range, const double& sqrtSNN, const double& mass, const double& charge, const bool& useAbsRap)
{
  // check input parameters
  if (range[0]>=range[1]) throw runtime_error("[ERROR] Wrong var1 range values for dist inclusive!");
  if (range[2]>=range[3]) throw runtime_error("[ERROR] Wrong var2 range values for dist inclusive!");
  if (range[4]>=range[5]) throw runtime_error("[ERROR] Wrong var3 range values for dist inclusive!");
  
  double res(0);
  if (useAbsRap && range[2]==0) {
    res += evalTF(f3D, {sqrtSNN, mass, charge, 0}, range[0], range[1], -range[3], range[3], range[4], range[5], id);
  }
  else if (useAbsRap && range[2]>0) {
    res += evalTF(f3D, {sqrtSNN, mass, charge, 0}, range[0], range[1], -range[3], -range[2], range[4], range[5], id);
    res += evalTF(f3D, {sqrtSNN, mass, charge, 0}, range[0], range[1], range[2], range[3], range[4], range[5], id);
  }
  else {
    res += evalTF(f3D, {sqrtSNN, mass, charge, 0}, range[0], range[1], range[2], range[3], range[4], range[5], id);
  }

  return res;
};

double yield_Inc_PbPb(const std::array<double, 6> range, const double& sqrtSNN, const double& mass, const int& charge, const bool& useAbsRap=true)
{
  static TF3* f3D(0);
  if (!f3D) f3D = new TF3("f3DYield", yield_PtRapCent_PbPb, range[0], range[1], -range[3], range[3], range[4], range[5], 4);

  return dist_Inc_PbPb(*f3D, 0, range, sqrtSNN, mass, charge, useAbsRap);
};

double yield_Diff_PbPb(const std::array<double, 6> range, const double& sqrtSNN, const double& mass, const int& charge, const bool& useAbsRap=true)
{
  auto binWidth = (range[5] - range[4]) * (range[3] - range[2]) * (range[1] - range[0]);
  if (useAbsRap) binWidth *= 2;

  return yield_Inc_PbPb(range, sqrtSNN, mass, charge, useAbsRap) / binWidth;
};

double spectra_Inc_PbPb(const std::array<double, 6> range, const double& sqrtSNN, const double& mass, const int& charge, const bool& useAbsRap=true)
{
  static TF3* f3D(0);
  if (!f3D) f3D = new TF3("f3DSpectra", spectra_PtRapCent_PbPb, range[0], range[1], -range[3], range[3], range[4], range[5], 4);

  return dist_Inc_PbPb(*f3D, 1, range, sqrtSNN, mass, charge, useAbsRap);
};

double spectra_Diff_PbPb(const std::array<double, 6> range, const double& sqrtSNN, const double& mass, const int& charge, const bool& useAbsRap=true)
{
  auto binWidth = (range[5] - range[4]) * (range[3] - range[2]) * (range[1] - range[0]);
  if (useAbsRap) binWidth *= 2;

  return spectra_Inc_PbPb(range, sqrtSNN, mass, charge, useAbsRap) / binWidth;
};


#endif // #ifndef extrapolationUtils_Yield_PbPb_h
