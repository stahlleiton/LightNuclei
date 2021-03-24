#ifndef draw_Projections_C
#define draw_Projections_C

// Custom headers
#include "../Macros/tdrstyle.C"
#include "../Macros/CMS_lumi.C"
#include "../Macros/utility_LightNuclei.h"
#include "extrapolationUtils.h"
#include "projectedData_Binning.h"
#include "utils.h"


// Global constants
const std::map<std::string, float> PROJ_NMB({{"PbPb_5p5TeV", 25.E9}, {"pp_14TeV", 27E12}});

// Functions
void fillProjectedDataGraphs(GRAPH& input, const DATA_CONTBIN& data, GRAPH gSigniCount, GRAPH gSigniFit);
void plot_totalYield_PbPb(GRAPH& input);
void plot_v2_PbPb(GRAPH& input);
void plot_ratio_pp(GRAPH& input);


void draw_Projections()
{
  // Extract significance
  GRAPH g_SIGNI_Fit, g_SIGNI_Count;
  auto m_SIGNI_Count = extractGraph(g_SIGNI_Count, "SIGNICOUNT", "dEdxCut", binning_ProjectedData);
  auto m_SIGNI_Fit = extractGraph(g_SIGNI_Fit, "SIGNIFIT", "dEdxCut", binning_ProjectedData);
  if (!m_SIGNI_Count.empty() && !m_SIGNI_Fit.empty()) return;

  // Fill graphs with projected data
  GRAPH g_PROJ;
  fillProjectedDataGraphs(g_PROJ, binning_ProjectedData, g_SIGNI_Count, g_SIGNI_Fit);

  // Plot ALICE data graphs
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetFillColor(kWhite);
  
  plot_totalYield_PbPb(g_PROJ);
  plot_v2_PbPb(g_PROJ);
  plot_ratio_pp(g_PROJ);
};


void fillProjectedDataGraphs(GRAPH& input, const DATA_CONTBIN& data, GRAPH gSigniCountM, GRAPH gSigniFitM)
{
  // Fill graphs using projected data bins
  for (const auto& col : data) {
    const auto sqrtS = getBeamEnergy(col.first);
    for (const auto& chg : col.second) {
      for (const auto& par : chg.second) {
	for (const auto& obs : par.second) {
	  for (const auto& r : obs.second) {
	    auto bin = r.first;
	    for (const auto& c : r.second) {
	      if (c.second.empty()) continue;
	      bin.add(c.first);
	      const auto& xVar = c.second.begin()->first;
	      const auto& bins = c.second.begin()->second;
	      if (bins.size()<2) continue;
	      std::cout << "0" << std::endl;
	      const auto& gSigniC = gSigniCountM[col.first][obs.first][r.first][c.first][chg.first][par.first];
	      const auto& gSigniF = gSigniFitM[col.first][obs.first][r.first][c.first][chg.first][par.first];
	      std::vector<std::array<double, 4> > points;
	      for (size_t i=1; i<bins.size(); i++) {
		bin[xVar] = {bins[i-1], bins[i]};
		const auto xVal = (bins[i] + bins[i-1])/2.;
		const auto xErr = roundVal(bins[i] - bins[i-1], 4)/2.;
		// Get projection value
		const auto charge = (chg.first=="Anti" ? -1 : 1)*CHARGE.at(par.first);
		const auto range = std::array<double, 6>({bin.at("pT").first, bin.at("pT").second, bin.at("absRap").first, bin.at("absRap").second, bin.at("cent").first/100., bin.at("cent").second/100.});
		double yVal = getProjection(col.first, obs.first, range, MASS.at(par.first), charge, false, true);
		if (obs.first=="yield") yVal *= PROJ_NMB.at(col.first);
		// Get projected uncertainty
		const auto signiC = (gSigniC.GetN()>0 ? gSigniC.GetY()[i-1] : 0);
		const auto signiF = (gSigniF.GetN()>0 ? gSigniF.GetY()[i-1] : 0);;
		auto signi = std::abs(signiC);//signiC>0. ? signiC : signiF;
		std::cout << " >> signi: " << signi << " , " << signiC << " , " << signiF << std::endl;
		signi *= std::sqrt(PROJ_NMB.at(col.first)/1.E9);
		if (col.first.rfind("pp_",0)!=0 && signi<=2.) continue;
		const auto yErr = (obs.first=="v2" ? 1 : yVal) * (signi>0 ? 1/signi : 0);
		std::cout << "[INFO] Deriving extrapolation for: " << col.first << " , " << chg.first+par.first << " in bin " << bin.str() << " >> " << yVal << " , " << yErr;
		std::cout << " >> signi: " << signi << " , " << signiC << " , " << signiF << std::endl;
		if (yVal==0) throw std::runtime_error("[ERROR] Null y value!");
		points.push_back({xVal, yVal, xErr, yErr});
	      }
	      auto& graph = input[col.first][obs.first][r.first][c.first][chg.first][par.first];
	      graph.Set(points.size());
	      graph.SetTitle(xVar.c_str());
	      for (size_t i=0; i<points.size(); i++) {
		const auto& point = points[i];
		graph.SetPoint(i, point[0], point[1]);
		graph.SetPointError(i, point[2], point[2], point[3], point[3]);
	      }
	    }
	  }
	}
      }
    }
  }
};


void formatGraph(TGraphAsymmErrors& graph, const std::string& varX, const std::string& varY, const std::array<double, 4>& range, const int& color, const int& marker, const float& mSize=2.)
{
  graph.SetMarkerSize(mSize);
  graph.SetFillStyle(0);
  graph.SetMarkerStyle(marker);
  graph.SetMarkerColor(color);
  graph.SetLineColor(color);
  graph.GetXaxis()->CenterTitle(true);
  graph.GetXaxis()->SetLabelSize(0.04);
  graph.GetXaxis()->SetTitleOffset(0.81);
  graph.GetXaxis()->SetTitleSize(0.063);
  graph.GetXaxis()->SetTitleFont(42);
  const auto& xVar = VAR.at(varX);
  graph.GetXaxis()->SetTitle((xVar.first+(xVar.second!="" ? " ("+xVar.second+")" : "")).c_str());
  graph.GetXaxis()->SetLimits(range[0], range[1]);
  graph.GetYaxis()->SetTitle(varY.c_str());
  graph.GetYaxis()->CenterTitle(true);
  graph.GetYaxis()->SetLabelSize(0.040);
  graph.GetYaxis()->SetTitleSize(0.063);
  graph.GetYaxis()->SetTitleOffset(0.85);
  graph.GetYaxis()->SetTitleFont(42);
  graph.GetYaxis()->SetRangeUser(range[2], range[3]);
};


void plot_totalYield_PbPb(GRAPH& input)
{
  for (auto& r : input["PbPb_5p5TeV"]["yield"]) {
    for (auto& c : r.second) {
      for (auto& chg : c.second) {
	TCanvas canvas("c", "c", 1000, 800); canvas.cd();
	  
	size_t j(0);
	std::map<std::string, TGraphAsymmErrors> wGraphs;
	for (auto& gg : chg.second) {
	  auto& graph = wGraphs[gg.first];
	  graph = gg.second;
	  formatGraph(graph, graph.GetTitle(), "Expected yield", {0.0, 3.5, 1., 1.E16}, COLOR[j++], 20);
	}	  
	if (wGraphs.empty()) continue;
      
	wGraphs.begin()->second.Draw("AP2");
	for (auto& g : wGraphs) g.second.Draw("SAMEP2");
	  
	TLegend leg(0.20, 0.60, 0.35, 0.78);
	for (auto& g : wGraphs) {
	  const auto& pLbl = LABEL.at(chg.first).at(g.first);
	  leg.AddEntry(&g.second, ("#bf{"+pLbl+"}").c_str(), "l")->SetTextSize(0.035);
	}
	leg.Draw("SAME");

	TLatex tex; tex.SetNDC();
	tex.SetTextSize(0.040); tex.SetTextFont(61); tex.DrawLatex(0.2, 0.86, "CMS Upgrade projection");
	const auto& centBin = c.first.at("cent");
	tex.SetTextSize(0.040); tex.SetTextFont(62); tex.DrawLatex(0.2, 0.81, Form("#bf{Pb-Pb, #sqrt{s_{NN}} = 5.5 TeV (%.0f-%.0f%%), 3 nb^{-1}}", centBin.first, centBin.second));
    
	CMS_lumi(&canvas, 33, "", "", false, -1, false);

	canvas.SetLogy();
	canvas.Update();
	  
	gSystem->mkdir("Plot/Projection", kTRUE);
	canvas.SaveAs(("Plot/Projection/Projection_yield_PbPb_5p5TeV_"+chg.first+"_"+r.first.str(1)+"_"+c.first.str(1)+".png").c_str());
      }
    }
  }
};


void plot_v2_PbPb(GRAPH& input)
{
  for (auto& c : input["PbPb_5p5TeV"]["v2"]) {

    const auto nBins = c.second.size();
    
    TCanvas canvas("c", "c", nBins>1? 1500 : 1000, nBins>1 ? 500 : 1000); canvas.cd();
    
    std::map<BINF, TPad*> pads;
    std::map<BINF, std::unique_ptr<TLegend> > legs;
    std::map<BINF, std::map<std::string, TGraphAsymmErrors> > wGraphsM;
    
    size_t iPad(0);
    for (auto& r : c.second) {
      if ((r.first.at("absRap").second-r.first.at("absRap").first)>1) continue;
      for (auto& chg : r.second) {

	canvas.cd();
	auto& pad = pads[r.first];
	if (nBins==1) pad = new TPad("", "", 0, 0, 1, 1);
	else if (iPad==0) pad = new TPad("", "", 0, 0, 0.36, 1);
	else if (iPad==1) pad = new TPad("", "", 0.36, 0, 0.68, 1);
	else if (iPad==2) pad = new TPad("", "", 0.68, 0, 1, 1);
	pad->SetFixedAspectRatio(1);
	if (iPad>0) pad->SetLeftMargin(0.0);
	if (nBins>1 && iPad<2) pad->SetRightMargin(0.0);
	pad->SetFillStyle(4000);
	pad->SetFrameFillStyle(4000);
	pad->Draw();
	pad->cd();

	size_t j(0);
	auto& wGraphs = wGraphsM[r.first];
	for (auto& gg : chg.second) {
	  auto& graph = wGraphs[gg.first];
	  graph = gg.second;
	  formatGraph(graph, graph.GetTitle(), iPad==0?"v_{2}":"", {0.1, 11, 0, 0.8}, COLOR[j], MARKER[j]);
	  j += 1;
	}
	if (wGraphs.empty()) continue;
	
	wGraphs.begin()->second.Draw("AP");
	for (auto& g : wGraphs) g.second.Draw("SAMEP");

	legs[r.first].reset(new TLegend(iPad==0?0.2:0.05, 0.65, iPad==0?0.35:0.2, 0.78));
	auto& leg = *legs[r.first];
	for (auto& g : wGraphs) {
	  const auto& pLbl = LABEL.at(chg.first).at(g.first);
	  leg.AddEntry(&g.second, ("#bf{"+pLbl+"}").c_str(), "p")->SetTextSize(0.035);
	}
	leg.Draw("SAME");

	TLatex tex; tex.SetNDC();
	const auto& cntBin = c.first.at("cent");
	const auto& rapBin = r.first.at("absRap");
	const std::string cntLbl = Form("%.0f-%.0f%%", cntBin.first, cntBin.second);
	const std::string rapLbl = rapBin.first>0 ? Form("%.0f<|y|<%.0f", rapBin.first, rapBin.second) : Form("|y|<%.0f", rapBin.second);
	tex.SetTextSize(0.045); tex.SetTextFont(62); tex.DrawLatex(iPad==0?0.2:0.05, 0.84, Form("#bf{%s, %s}", cntLbl.c_str(), rapLbl.c_str()));
    
	CMS_lumi(pad, 33, "", "", false, -1, false);

	pad->Update();
	iPad += 1;
      }
    }
    
    canvas.cd();
    TLatex tex; tex.SetNDC();
    tex.SetTextSize(0.060); tex.SetTextFont(61); tex.DrawLatex(0.06, 0.94, "CMS"); tex.SetTextFont(62);
    tex.SetTextSize(0.060); tex.SetTextFont(52); tex.DrawLatex(0.11, 0.94, "Phase-2 Simulation"); tex.SetTextFont(62);
    tex.SetTextSize(0.060); tex.SetTextFont(62); tex.DrawLatex(0.81, 0.94, "#bf{PbPb 3 nb^{-1} (5.5 TeV)}");
    
    gSystem->mkdir("Plot/Projection", kTRUE);
    canvas.SaveAs(("Plot/Projection/Projection_v2_PbPb_5p5TeV_"+c.first.str(1)+"_"+".png").c_str());
  }
};


double modelFunc(const double* x, const double* p)
{
  // source: https://doi.org/10.1016/j.physletb.2019.03.033
  // input variables
  const auto& dNdeta = x[0];

  // input parameters
  const auto& type = p[0]; // (0) dOverP  ,  (1) He3OverP [two-body coa.]   ,  (2) H3OverHe3 [two-body coa.]   ,   (3) He3OverP [three-body coa.]   ,   (4) H3OverHe3 [three-body coa.]
  const auto scale = (p[1]>0 ? p[1] : 1.);

  // compute model parameters
  const auto Tk = 80.6 + 83.*TMath::Power(1. + 2.33*dNdeta/67.3,  -1./2.33);
  const auto R = 197.3 * TMath::Power(4.18125*dNdeta, 1./3.) / TMath::Sqrt(938.272*Tk);

  double res(0);
  if      (p[0]==0) res = 4.E-3/TMath::Power(1. + 1.6*1.6/R/R, 1.5);
  else if (p[0]==1) res = 7.1E-6/TMath::Power(1. + 1.15*1.15/R/R, 1.5)/TMath::Power(1. + 1.6*1.6/R/R, 1.5);
  else if (p[0]==2) res = TMath::Power(1. + 1.15*1.15/R/R, 1.5)/TMath::Power(1. + 1.039*1.039/R/R, 1.5);
  else if (p[0]==3) res = 7.1E-6/TMath::Power(1. + 1.24*1.24/R/R, 3.);
  else if (p[0]==4) res = TMath::Power(1. + 1.24*1.24/R/R, 3.)/TMath::Power(1. + 1.12*1.12/R/R, 3.);

  return scale * res;
};


void plot_ratio_pp(GRAPH& input)
{
  const auto& nMB = PROJ_NMB.at("pp_14TeV");
  for (auto& r : input["pp_14TeV"]["yield"]) {
    for (auto& c : r.second) {
      for (auto& chg : c.second) {
	TCanvas canvas("c", "c", 2000, 2000); canvas.cd();

	// extract light nuclei yields
	if (chg.second.at("Triton").GetN()==0) throw std::runtime_error("[ERROR] Bins empty for "+chg.first+"Triton");
	if (chg.second.at("Helium3").GetN()==0) throw std::runtime_error("[ERROR] Bins empty for "+chg.first+"Helium3");
	if (chg.second.at("Deuteron").GetN()==0) throw std::runtime_error("[ERROR] Bins empty for "+chg.first+"Deuteron");
	if (chg.second.at("Helium4").GetN()==0) throw std::runtime_error("[ERROR] Bins empty for "+chg.first+"Helium4");
	const auto Nd   = chg.second.at("Deuteron").GetY()[0];
	const auto NHe3 = chg.second.at("Helium3").GetY()[0];
	const auto NH3  = chg.second.at("Triton").GetY()[0];
	const auto NHe4 = chg.second.at("Helium4").GetY()[0];

	// create model function
	const auto dNdeta = 6.01*TMath::Power(14./7., 0.206); // using as reference the averaged dN/deta at 7 TeV
	const auto Np = nMB*0.0223*dNdeta;
	TF1 model("model", modelFunc, 1, 300, 2);

	// draw first pad (top)
	canvas.cd();
	TPad pad1("pad1", "", 0.0, 0.5, 0.5, 1);
	pad1.SetBottomMargin(0.0);
	pad1.SetRightMargin(0.0);
	pad1.SetLeftMargin(0.15);
	pad1.Draw(); pad1.cd();
	const auto f = 1000.; // scale factor

	TGraphAsymmErrors g_dOverP(1);
	g_dOverP.SetPoint(0, dNdeta, Nd*f/Np);
	const auto unc_dOverP = ratioUnc(Nd*f, Np, chg.second.at("Deuteron").GetEYlow()[0]*f, 0);
	g_dOverP.SetPointError(0, 0, 0, unc_dOverP, unc_dOverP);
	std::cout << "dOverP UNCERTAINTY: " << unc_dOverP/(Nd*f/Np) << std::endl;
	formatGraph(g_dOverP, "dNdeta", "d/p (#times10^{-3})", {1, 300, 0.1, 6}, kBlack, 20, 3.);
	g_dOverP.Draw("AP");

	model.SetParameters(0, f);
	const auto m_dOverP = model.DrawCopy("SAME");
	m_dOverP->SetLineColor(kBlack);
	m_dOverP->SetLineWidth(3);

	TGraphAsymmErrors d_dOverP_0p9TeV(1), d_dOverP_2p76TeV(1), d_dOverP_7TeV(1);
	// source: https://www.hepdata.net/record/ins1625294
	if (chg.first=="Anti") {
	  d_dOverP_0p9TeV.SetPoint(0, 3.81, 0.00139*f);
	  d_dOverP_2p76TeV.SetPoint(0, 4.88, 0.001309*f);
	  d_dOverP_7TeV.SetPoint(0, 6.01, 0.001563*f);
	  d_dOverP_0p9TeV.SetPointError(0, 0, 0, sumError(0.00015, 0.00011)*f, sumError(0.00015, 0.00011)*f);
	  d_dOverP_2p76TeV.SetPointError(0, 0, 0, sumError(0.00016, 4.7e-05)*f, sumError(0.00016, 4.7e-05)*f);
	  d_dOverP_7TeV.SetPointError(0, 0, 0, sumError(0.00017, 1.3e-05)*f, sumError(0.00017, 1.3e-05)*f);
	}
	else {
	  d_dOverP_0p9TeV.SetPoint(0, 3.81, 0.00138*f);
	  d_dOverP_2p76TeV.SetPoint(0, 4.88, 0.001482*f);
	  d_dOverP_7TeV.SetPoint(0, 6.01, 0.001626*f);
	  d_dOverP_0p9TeV.SetPointError(0, 0, 0, sumError(0.00015, 0.00014)*f, sumError(0.00015, 0.00014)*f);
	  d_dOverP_2p76TeV.SetPointError(0, 0, 0, sumError(0.00014, 3.6e-05)*f, sumError(0.00014, 3.6e-05)*f);
	  d_dOverP_7TeV.SetPointError(0, 0, 0, sumError(0.00017, 1.3e-05)*f, sumError(0.00017, 1.3e-05)*f);
	}
	formatGraph(d_dOverP_0p9TeV, "dNdeta", "", {1, 300, 1.E-8, 9.E-6}, kViolet, 23, 3.);
	formatGraph(d_dOverP_2p76TeV, "dNdeta", "", {1, 300, 1.E-8, 9.E-6}, kBlue, 24, 3.);
	formatGraph(d_dOverP_7TeV, "dNdeta", "", {1, 300, 1.E-8, 9.E-6}, kRed, 22, 3.);
	d_dOverP_0p9TeV.Draw("SAMEP");
	d_dOverP_2p76TeV.Draw("SAMEP");
	d_dOverP_7TeV.Draw("SAMEP");

	CMS_lumi(&pad1, 33, "", "", false, -1, false);
	
	pad1.SetLogx();
	pad1.Update();

	// draw second pad
	canvas.cd();
	TPad pad2("pad2", "", 0.0, 0, 0.5, 0.5);
	pad2.SetTopMargin(0.0);
	pad2.SetRightMargin(0.0);
	pad2.SetLeftMargin(0.15);
	pad2.Draw(); pad2.cd();

	TGraphAsymmErrors g_He3OverP(1);
	g_He3OverP.SetPoint(0, dNdeta, NHe3/Np);
	const auto unc_He3OverP = ratioUnc(NHe3, Np, chg.second.at("Helium3").GetEYlow()[0], 0);
	g_He3OverP.SetPointError(0, 0, 0, unc_He3OverP, unc_He3OverP);
	std::cout << "He3OverP UNCERTAINTY: " << unc_He3OverP/(NHe3/Np) << std::endl;
	formatGraph(g_He3OverP, "dNdeta", "{}^{3}He/p", {1, 300, 1.E-8, 9.E-6}, kBlack, 20, 3.);
	g_He3OverP.Draw("AP");

	model.SetParameters(1, 1.);
	const auto m_He3OverP_2Body = model.DrawCopy("SAME");
	m_He3OverP_2Body->SetLineColor(kRed);
	m_He3OverP_2Body->SetLineWidth(3);
	m_He3OverP_2Body->SetLineStyle(7);
	model.SetParameters(3, 1.);
	const auto m_He3OverP_3Body = model.DrawCopy("SAME");
	m_He3OverP_3Body->SetLineColor(kBlue);
	m_He3OverP_3Body->SetLineWidth(3);

	TGraphAsymmErrors d_He3OverP(1);
	// source: https://journals.aps.org/prc/pdf/10.1103/PhysRevC.97.024615
	if (chg.first=="Anti") {
	  d_He3OverP.SetPoint(0, 6.01, 0.000000895468749999998);
	  d_He3OverP.SetPointError(0, 0, 0, 0.00000051926806406426, 0.00000051926806406426);
	}
	else {
	  d_He3OverP.SetPoint(0, 6.01, 0.000000885445544554456);
	  d_He3OverP.SetPointError(0, 0, 0, 0.000000512023910359952, 0.000000512023910359952);
	}
	formatGraph(d_He3OverP, "dNdeta", "{}^{3}He/p", {1, 300, 1.E-8, 9.E-6}, kRed, 22, 3.);
	d_He3OverP.Draw("SAMEP");

	CMS_lumi(&pad2, 33, "", "", false, -1, false);
	
	pad2.SetLogx();
	pad2.SetLogy();
	pad2.Update();

	// draw third pad
	canvas.cd();
	TPad pad3("pad3", "", 0.5, 0.5, 1.0, 1);
	pad3.SetBottomMargin(0.0);
	pad3.SetLeftMargin(0.0);
	pad3.SetRightMargin(0.15);
	pad3.Draw(); pad3.cd();

	TGraphAsymmErrors g_H3OverHe3(1);
	g_H3OverHe3.SetPoint(0, dNdeta, NH3/NHe3);
	const auto unc_H3OverHe3 = ratioUnc(NH3, NHe3, chg.second.at("Triton").GetEYlow()[0], chg.second.at("Helium3").GetEYlow()[0]);
	std::cout << "H3OverHe3 UNCERTAINTY: " << unc_H3OverHe3/(NH3/NHe3) << std::endl;
	g_H3OverHe3.SetPointError(0, 0, 0, unc_H3OverHe3, unc_H3OverHe3);
	formatGraph(g_H3OverHe3, "dNdeta", "{}^{3}H/{}^{3}He", {1, 300, 0.11, 2.5}, kBlack, 20, 3.);
	g_H3OverHe3.Draw("APY+");

	model.SetParameters(2, 1.);
	const auto m_H3OverHe3_2Body = model.DrawCopy("SAME");
	m_H3OverHe3_2Body->SetLineColor(kRed);
	m_H3OverHe3_2Body->SetLineWidth(3);
	m_H3OverHe3_2Body->SetLineStyle(7);
	model.SetParameters(4, 1.);
	const auto m_H3OverHe3_3Body = model.DrawCopy("SAME");
	m_H3OverHe3_3Body->SetLineColor(kBlue);
	m_H3OverHe3_3Body->SetLineWidth(3);

	CMS_lumi(&pad3, 33, "", "", false, -1, false);

	pad3.SetLogx();
	pad3.Update();

	// draw fourth pad
	canvas.cd();
	TPad pad4("pad4", "", 0.5, 0, 1.0, 0.5);
	pad4.SetTopMargin(0.0);
	pad4.SetLeftMargin(0.0);
	pad4.SetRightMargin(0.15);
	pad4.Draw(); pad4.cd();

	TGraphAsymmErrors g_He4OverHe3(1);
	g_He4OverHe3.SetPoint(0, dNdeta, NHe4/NHe3);
	const auto unc_He4OverHe3 = ratioUnc(NHe4, NHe3, chg.second.at("Helium4").GetEYlow()[0], chg.second.at("Helium3").GetEYlow()[0]);
	std::cout << "He4OverHe3 UNCERTAINTY: " << unc_He4OverHe3/(NHe4/NHe3) << std::endl;
	g_He4OverHe3.SetPointError(0, 0, 0, unc_He4OverHe3, unc_He4OverHe3);
	formatGraph(g_He4OverHe3, "dNdeta", "{}^{4}He/{}^{3}He", {1, 300, 1.E-5, 9.E-3}, kBlack, 20, 3.);
	g_He4OverHe3.Draw("APY+");

	CMS_lumi(&pad4, 33, "", "", false, -1, false);

	pad4.SetLogx();
	pad4.SetLogy();
	pad4.Update();

	// add legend
	pad1.cd();
	TLegend legData(0.2, 0.65, 0.7, 0.85);
	legData.SetHeader("ALICE published");
	legData.AddEntry(&d_dOverP_0p9TeV, "p+p @ 0.9 TeV", "p");
	legData.AddEntry(&d_dOverP_2p76TeV, "p+p @ 2.76 TeV", "p");
	legData.AddEntry(&d_dOverP_7TeV, "p+p @ 7 TeV", "p");
	legData.Draw("SAME");
	
	TLegend legProj(0.2, 0.50, 0.6, 0.65);
	legProj.SetHeader("CMS projection");
	legProj.AddEntry(&g_dOverP, "p+p @ 14 TeV", "p");
	legProj.Draw("SAME");
	
	pad2.cd();
	TLegend legModel(0.5, 0.3, 0.99, 0.5);
	legModel.SetHeader("PLB 792 (2019) 132-137");
	legModel.AddEntry(m_dOverP, "COAL.", "l");
	legModel.AddEntry(m_He3OverP_2Body, "Two-body COAL.", "l");
	legModel.AddEntry(m_He3OverP_3Body, "Three-body COAL.", "l");
	legModel.Draw("SAME");
    
	canvas.cd();
	TLatex tex; tex.SetNDC();
	tex.SetTextSize(0.030); tex.SetTextFont(61); tex.DrawLatex(0.08, 0.968, "CMS"); tex.SetTextFont(62);
	tex.SetTextSize(0.030); tex.SetTextFont(52); tex.DrawLatex(0.15, 0.968, "Phase-2 Simulation"); tex.SetTextFont(62);
	tex.SetTextSize(0.030); tex.SetTextFont(62); tex.DrawLatex(0.69, 0.968, "#bf{pp 300 pb^{-1} (14 TeV)}");
	  
	gSystem->mkdir("Plot/Projection", kTRUE);
	canvas.SaveAs(("Plot/Projection/Projection_ratio_pp_14TeV_"+chg.first+"_"+r.first.str(1)+"_"+c.first.str(1)+".png").c_str());
      }
    }
  }
};


#endif // #ifndef draw_Projections_C
