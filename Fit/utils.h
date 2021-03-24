#ifndef LightNuclei_Fit_utils_h
#define LightNuclei_Fit_utils_h

// Custom headers
#include "../Macros/tdrstyle.C"
#include "../Macros/CMS_lumi.C"
#include "../Macros/bin.h"
#include "../Macros/utility_LightNuclei.h"
#include "../Macros/fitUtils.h"
// ROOT headers
#include "TSystem.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLegendEntry.h"


// Type definition
typedef std::map< std::string , std::map< BINF , std::map< BINF , std::unique_ptr<TF1> > > > FITFUNC1;
typedef std::map< std::string , std::map< std::string , FITFUNC1 > > FITFUNC2;
typedef std::map< std::string , FITFUNC2 > FITFUNC;
typedef std::map< std::string , std::map< BINF , std::map< BINF , TGraphAsymmErrors > > > GRAPH1;
typedef std::map< std::string , std::map< std::string , GRAPH1 > > GRAPH2;
typedef std::map< std::string , GRAPH2 > GRAPH;
typedef std::map<std::string, GRAPH> GRAPHM;


// Functions
void processDataGraphs(GRAPH& input, const std::set<std::string>& obsS={})
{
  // Fill antiparticle graphs using ratio data
  for (const auto& col : input) {
    if (col.second.find("All")==col.second.end()) continue;
    for (const auto& par : col.second.at("All")) {
      if (par.second.find("ratio")==par.second.end()) continue;
      for (const auto& obs : col.second.at("").at(par.first)) {
	for (const auto& r : obs.second) {
	  for (const auto& rr : par.second.at("ratio")) {
	    if (!rr.first.contain(r.first)) continue;
	    for (const auto& c : r.second) {
	      for (const auto& cc : rr.second) {
		if (!c.first.contain(cc.first)) continue;
		auto& graph = input[col.first]["Anti"][par.first][obs.first][r.first][c.first];
		graph = c.second;
		for (int i=0; i<graph.GetN(); i++) {
		  for (int j=0; j<cc.second.GetN(); j++) {
		    if (!BINF(cc.second, j).contain(BINF(graph, i))) continue;
		    const auto yVal = graph.GetY()[i]*cc.second.GetY()[j];
		    const auto yErrLo = sumError(graph.GetEYlow ()[i]*cc.second.GetY()[j], graph.GetY()[i]*cc.second.GetEYlow ()[j]);
		    const auto yErrHi = sumError(graph.GetEYhigh()[i]*cc.second.GetY()[j], graph.GetY()[i]*cc.second.GetEYhigh()[j]);
		    graph.SetPointY(i, yVal);
		    graph.SetPointError(i, graph.GetEXlow()[i], graph.GetEXhigh()[i], yErrLo, yErrHi);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  // Fill yield and spectra graphs
  for (const auto& col : input) {
    for (const auto& chg : col.second) {
      for (const auto& par : chg.second) {
	for (const auto& obs : par.second) {
	  if (obs.first!="yield" && obs.first!="spectra") continue;
	  for (const auto& r : obs.second) {
	    for (const auto& c : r.second) {
	      auto& graph = input[col.first][chg.first][par.first][obs.first=="yield"?"spectra":"yield"][r.first][c.first];
	      graph = c.second;
	      for (int i=0; i<graph.GetN(); i++) {
		auto f = 2.*TMath::Pi()*graph.GetX()[i];
		if (obs.first=="yield") f = 1./f;
		graph.SetPointY(i, graph.GetY()[i]*f);
		graph.SetPointError(i, graph.GetEXlow()[i], graph.GetEXhigh()[i], graph.GetEYlow()[i]*f, graph.GetEYhigh()[i]*f);
	      }
	    }
	  }
	}
      }
    }
  }
  if (!obsS.empty()) {
    auto tmp = input;
    for (const auto& col : tmp) {
      for (const auto& chg : col.second) {
	for (const auto& par : chg.second) {
	  for (const auto& obs : par.second) {
	    if (obsS.find(obs.first)==obsS.end()) {
	      input[col.first][chg.first][par.first].erase(obs.first);
	    }
	  }
	}
      }
    }
  }
};


void fillDataGraphs(GRAPH& input, const std::set<std::string>& parS, const DATA_CONT& data, const std::set<std::string>& obsS={})
{ 
  // Fill graphs using data
  const bool excVn = (obsS.find("v2")==obsS.end() && obsS.find("v3")==obsS.end());
  for (const auto& col : data) {
    for (const auto& chg : col.second) {
      for (const auto& par : chg.second) {
	if (!parS.empty() && parS.find(par.first)==parS.end()) continue;
	for (const auto& obs : par.second) {
	  if (excVn && (obs.first=="v2" || obs.first=="v3")) continue;
	  for (const auto& r : obs.second) {
	    for (const auto& c : r.second) {
	      auto& graph = input[col.first][chg.first][par.first][obs.first][r.first][c.first];
	      graph.Set(c.second.size());
	      graph.SetTitle(c.second.begin()->first.begin()->first.c_str());
	      size_t i(0);
	      for (const auto& b : c.second) {
		const auto& xVar = b.first.begin()->second;
		const auto xVal = (xVar.second + xVar.first)/2.;
		const auto xErr = roundVal(xVar.second - xVar.first, 4)/2.;
		const auto yVal = b.second[0];
		const auto fc1 = b.second[2]<=0. ? 0. : yVal;
		const auto fc2 = b.second[4]<=0. ? 0. : yVal;
		auto yErrLo = sumError(b.second[2] - fc1, b.second[4] - fc2);
		auto yErrHi = sumError(b.second[1] - fc1, b.second[3] - fc2);
		if (yErrLo==0.) { yErrLo = 0.05*yVal; yErrHi = 0.05*yVal; }
		graph.SetPoint(i, xVal, yVal);
		graph.SetPointError(i++, xErr, xErr, yErrLo, yErrHi);
	      }
	    }
	  }
	}
      }
    }
  }
  processDataGraphs(input, obsS);
};


void fillDataGraphs(GRAPH& input, const std::set<std::string>& colS, const std::set<std::string>& parS, const std::set<std::string>& obsS)
{ 
  // Fill graphs using data
  const bool excVn = (obsS.find("v2")==obsS.end() && obsS.find("v3")==obsS.end());
  for (const auto& col : colS) {
    for (const auto& chg : std::set<std::string>({"", "Anti", "All"})) {
      for (const auto& par : std::set<std::string>({"Deuteron", "Helium3", "Helium4", "Triton"})) {
	if (!parS.empty() && parS.find(par)==parS.end()) continue;
	for (const auto& obs : std::set<std::string>({"yield", "spectra", "ratio", "v2", "v3"})) {
	  if (excVn && (obs=="v2" || obs=="v3")) continue;
	  const auto hist = getHist(col, chg+par, obs, false);
	  if (hist.GetEntries()==0) continue;
	  const auto& xAxis = hist.GetXaxis();
	  const auto& yAxis = hist.GetYaxis();
	  const auto& zAxis = hist.GetZaxis();
	  for (int iz=2; iz<zAxis->GetNbins(); iz++) {
	    BINF r(zAxis->GetTitle(), zAxis->GetBinLowEdge(iz), zAxis->GetBinUpEdge(iz));
	    for (int iy=2; iy<yAxis->GetNbins(); iy++) {
	      BINF c(yAxis->GetTitle(), yAxis->GetBinLowEdge(iy), yAxis->GetBinUpEdge(iy));
	      std::vector<std::array<double, 4>> points;
	      for (int ix=2; ix<xAxis->GetNbins(); ix++) {
		if (hist.GetBinContent(ix, iy, iz)==0. && hist.GetBinError(ix, iy, iz)==0.) continue;
		if (ix>2 && hist.GetBinContent(ix, iy, iz)==hist.GetBinContent(ix-1, iy, iz) && hist.GetBinError(ix, iy, iz)==hist.GetBinError(ix-1, iy, iz)) {
		  points.back()[1] = xAxis->GetBinUpEdge(ix);
		}
		else {
		  points.push_back({xAxis->GetBinLowEdge(ix), xAxis->GetBinUpEdge(ix), hist.GetBinContent(ix, iy, iz), hist.GetBinError(ix, iy, iz)});
		}
	      }
	      auto& graph = input[col][chg][par][obs][r][c];
	      graph.Set(points.size());
	      graph.SetTitle(xAxis->GetTitle());
	      for (size_t i=0; i<points.size(); i++) {
		const auto& point = points[i];
		const auto xVal = (point[1] + point[0])/2.;
		const auto xErr = roundVal(point[1] - point[0], 4)/2.;
		const auto yVal = point[2];
		auto yErr = point[3];
		if (yErr==0.) { yErr = 0.05*yVal; }
		graph.SetPoint(i, xVal, yVal);
		graph.SetPointError(i, xErr, xErr, yErr, yErr);
	      }
	    }
	  }
	}
      }
    }
  }
  processDataGraphs(input, obsS);
};


void plotGraphs(GRAPH& input, const std::string& type, FITFUNC& fitM)
{
  for (auto& col : input) {
    if (col.first.rfind("PbPb_", 0)!=0) continue;
    for (auto& chg : col.second) {
      for (auto& par : chg.second) {
	for (auto& obs : par.second) {
	  for (auto& r : obs.second) {
	    TCanvas c("c", "c", 1000, 800); c.cd();
	    std::map<std::string, std::map< BINF , TGraphAsymmErrors > > graphs;
	    std::string colPP = (col.first=="PbPb_5p5TeV" ? "pp_14TeV" : (col.first=="PbPb_5p02TeV" ? "pp_13TeV" : (col.first=="PbPb_2p76TeV" ? "pp_7TeV" : "")));
	    if (colPP=="") throw std::runtime_error("Invalid pp energy");
	    if (input.at(colPP).at(chg.first).find(par.first)==input.at(colPP).at(chg.first).end()) colPP = "pp_7TeV";
	    const auto enPbPb = col.first.substr(col.first.find("_")+1).replace(1,1,".");
	    const auto enPP = colPP.substr(colPP.find("_")+1);
	    std::string rapLbl;
	    for (const auto& rr : r.first) rapLbl += Form("%g #leq %s < %g %s, ", rr.second.first, VAR.at(rr.first).first.c_str(), rr.second.second, VAR.at(rr.first).second.c_str());
	    if (rapLbl.rfind(", ")!=std::string::npos) rapLbl = rapLbl.substr(0, rapLbl.rfind(", "));
	    
	    graphs["PbPb"] = r.second;
	    bool isYield = (obs.first=="yield" || obs.first=="spectra");
	    const bool hasPP = (input.find(colPP)!=input.end() && input.at(colPP).find(chg.first)!=input.at(colPP).end() && input.at(colPP).at(chg.first).find(par.first)!=input.at(colPP).at(chg.first).end());
	    if (hasPP && (isYield || obs.first=="significance")) graphs["pp"] = input.at(colPP).at(chg.first).at(par.first).at(obs.first).at(r.first);
	    
	    std::map<std::string, std::vector<TGraphAsymmErrors> > wGraphs;
	    for (auto& gg : graphs) {
	      wGraphs[gg.first].resize(gg.second.size());
	      size_t j(0), iG(0);
	      for (auto& g : gg.second) {
		auto& graph = wGraphs[gg.first][iG++];
		graph.Set(g.second.GetN());
		const auto& xName = g.second.GetTitle();
		double factor = 1.0;
		if (gg.first=="PbPb" && isYield) factor = std::pow(2., gg.second.size()-1-j);
		for (int i=0; i<graph.GetN(); i++) {
		  graph.SetPoint(i, g.second.GetX()[i], g.second.GetY()[i]*factor);
		  graph.SetPointError(i, g.second.GetEXlow()[i], g.second.GetEXhigh()[i], g.second.GetEYlow()[i]*factor, g.second.GetEYhigh()[i]*factor);
		}
		graph.SetMarkerSize(1.5);
		graph.SetFillStyle(0);
		graph.SetMarkerStyle(gg.first=="PbPb" ? 20 : 24);
		graph.SetMarkerColor(gg.first=="PbPb" ? COLOR[j] : kBlack);
		graph.SetLineColor(gg.first=="PbPb" ? COLOR[j++] : kBlack);
		graph.GetXaxis()->CenterTitle(false);
		graph.GetXaxis()->SetLabelSize(0.04);
		graph.GetXaxis()->SetTitleOffset(0.78);
		graph.GetXaxis()->SetTitleSize(0.063);
		graph.GetXaxis()->SetTitleFont(42);
		const auto& xVar = VAR.at(xName);
		graph.GetXaxis()->SetTitle((xVar.first+" ("+xVar.second+")").c_str());
		if (par.first=="Deuteron") graph.GetXaxis()->SetLimits(0.47, 7.0);
		else if (par.first=="Helium3") graph.GetXaxis()->SetLimits(0.47, 10.0);
		else if (par.first=="Helium4") graph.GetXaxis()->SetLimits(0.47, 10.0);
		else if (par.first=="Triton") graph.GetXaxis()->SetLimits(0.47, 5.0);
		std::string vlbl = "";
		if (obs.first=="yield") vlbl = "1/N_{ev}d^{2}N/(dp_{T}dy) (GeV/c)^{-1}";
		else if (obs.first=="spectra") vlbl = "1/(2#pip_{T}N_{ev})d^{2}N/(dp_{T}dy) (GeV/c)^{-2}";
		else if (obs.first=="rapYield") vlbl = "dN/dy";
		else if (obs.first=="v2") vlbl = "v_{2} {SP, |#Delta#eta| > 2}";
		else if (obs.first=="v3") vlbl = "v_{3} {SP, |#Delta#eta| > 2}";
		else if (obs.first=="ratio") vlbl = Form("%s/%s", LABEL.at("Anti").at(par.first).c_str(), LABEL.at("").at(par.first).c_str());
		else if (obs.first=="significance") vlbl = "Significance";
		graph.GetYaxis()->SetTitle(vlbl.c_str());
		graph.GetYaxis()->CenterTitle(false);
		if (obs.first=="ratio") graph.GetYaxis()->CenterTitle(true);
		graph.GetYaxis()->SetLabelSize(0.040);
		graph.GetYaxis()->SetTitleSize(0.063);
		graph.GetYaxis()->SetTitleOffset(1.1);
		graph.GetYaxis()->SetTitleFont(42);
		std::pair<double, double> range;
		if (par.first=="Deuteron") {
		  if (obs.first=="yield" && type.rfind("BKG_",0)==0) range = std::make_pair(2.E-10, 900000.);
		  else if (obs.first=="yield") range = std::make_pair(2.E-9, 900.);
		  else if (obs.first=="spectra") range = std::make_pair(2.E-11, 90.);
		  else if (obs.first=="v2") range = std::make_pair(-0.1, 0.7);
		  else if (obs.first=="v3") range = std::make_pair(-0.1, 0.7);
		  else if (obs.first=="ratio") range = std::make_pair(0.3, 1.7);
		  else if (obs.first=="significance") range = std::make_pair(-1800., 2500.);
		}
		else if (par.first=="Helium3") {
		  if (obs.first=="yield") range = std::make_pair(2.E-12, 1E-2);
		  else if (obs.first=="spectra") range = std::make_pair(2.E-14, 1E-3);
		  else if (obs.first=="v2") range = std::make_pair(-0.1, 1.0);
		  else if (obs.first=="ratio") range = std::make_pair(0.3, 1.7);
		  else if (obs.first=="significance") range = std::make_pair(-1800., 2500.);
		}
		else if (par.first=="Helium4") {
		  if (obs.first=="yield") range = std::make_pair(2.E-10, 1E-5);
		  else if (obs.first=="spectra") range = std::make_pair(2.E-11, 1E-6);
		  else if (obs.first=="v2") range = std::make_pair(-0.1, 0.7);
		  else if (obs.first=="ratio") range = std::make_pair(0.3, 1.7);
		  else if (obs.first=="significance") range = std::make_pair(-1800., 2500.);
		}
		else if (par.first=="Triton") {
		  if (obs.first=="yield") range = std::make_pair(2.E-9, 1E-2);
		  else if (obs.first=="spectra") range = std::make_pair(2.E-10, 1E-3);
		  else if (obs.first=="v2") range = std::make_pair(-0.1, 0.7);
		  else if (obs.first=="ratio") range = std::make_pair(0.3, 1.7);
		  else if (obs.first=="significance") range = std::make_pair(-1800., 2500.);
		}
		graph.GetYaxis()->SetRangeUser(range.first, range.second);
	      }
	    }

	    if (wGraphs.empty() || wGraphs.begin()->second.empty()) continue;
	    wGraphs.begin()->second.begin()->Draw("AP2");
	    for (auto& gg : wGraphs) for (auto& g : gg.second) g.Draw("SAMEP2");

	    TLine line1(0.47, 0.0, 7.0, 0.0), line2(0.47, 1.0, 7.0, 1.0);
	    line1.SetLineStyle(7); line2.SetLineStyle(7);
	    if (obs.first=="v2" || obs.first=="v3") line1.Draw("SAME");
	    else if (obs.first=="ratio") line2.Draw("SAME");

	    TLatex tex; tex.SetNDC();
	    double yPos = 0;
	    if (type=="ALICE") {
	      tex.SetTextSize(0.060); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "ALICE"); tex.SetTextFont(62);
	    }
	    else {
	      tex.SetTextSize(0.060); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
	      tex.SetTextSize(0.046); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62);
	      yPos = 0.05;
	    }
	    std::vector<std::unique_ptr<TLegend> > legV;
	    const auto& pLbl = LABEL.at(chg.first).at(par.first);
	    if (isYield || obs.first=="significance") {
	      if (graphs["PbPb"].size()>6) {
		legV.push_back(std::unique_ptr<TLegend>(new TLegend(0.7, 0.15, 0.95, 0.40)));
		legV.push_back(std::unique_ptr<TLegend>(new TLegend(0.5, 0.19, 0.75, 0.44)));
		auto parN = par.first; parN[0] = std::tolower(parN[0]); parN = chg.first + parN + "s"; parN[0] = std::toupper(parN[0]);
		if (type.rfind("BKG_",0)==0) parN = "HYDJET";
		legV[1]->SetHeader((parN+", Pb-Pb #sqrt{s_{NN}} = "+enPbPb).c_str());
		dynamic_cast<TLegendEntry*>(legV[1]->GetListOfPrimitives()->First())->SetTextSize(0.03);
		for (size_t i=0; i<graphs["PbPb"].size(); i++) {
		  auto& gr = *std::next(graphs["PbPb"].begin(), i);
		  const auto& bin = *gr.first.begin();
		  const auto fmt = (bin.first=="cent" ? "%.0f-%.0f" : "%.1f-%.1f")+VAR.at(bin.first).second+"% %s";
		  const auto eLbl = Form(fmt.c_str(), bin.second.first, bin.second.second, (isYield ? Form("(%.0fx)", std::pow(2., graphs["PbPb"].size()-1-i)) : ""));
		  (i%2!=0 ? legV[0] : legV[1])->AddEntry(&wGraphs["PbPb"][i], eLbl, "pf")->SetTextSize(0.03);
		}
		if (hasPP) for (const auto& gr : wGraphs["pp"]) legV[0]->AddEntry(&gr, ("pp, #sqrt{s}="+enPP+(isYield?" (1x)":"")).c_str(), "pf")->SetTextSize(0.025);
		legV[1]->Draw("SAME");
		legV[0]->Draw("SAME");
	      }
	      else {
		legV.push_back(std::unique_ptr<TLegend>(new TLegend(0.68, 0.15, 0.88, 0.45)));
		legV[0]->SetHeader(("Pb-Pb #sqrt{s_{NN}} = "+enPbPb).c_str());
		dynamic_cast<TLegendEntry*>(legV[0]->GetListOfPrimitives()->First())->SetTextSize(0.03);
		for (size_t i=0; i<graphs["PbPb"].size(); i++) {
		  auto& gr = *std::next(graphs["PbPb"].begin(), i);
		  const auto& bin = *gr.first.begin();
		  const auto fmt = (bin.first=="cent" ? "%.0f-%.0f" : "%.1f-%.1f")+VAR.at(bin.first).second+"% %s";
		  const auto eLbl = Form(fmt.c_str(), bin.second.first, bin.second.second, (isYield ? Form("(%.0fx)", std::pow(2., graphs["PbPb"].size()-1-i)) : ""));
		  legV[0]->AddEntry(&wGraphs["PbPb"][i], eLbl, "pf")->SetTextSize(0.025);
		}
		if (hasPP) for (const auto& gr : wGraphs["pp"]) legV[0]->AddEntry(&gr, ("pp, #sqrt{s}="+enPP+(isYield?" (1x)":"")).c_str(), "pf")->SetTextSize(0.025);
		legV[0]->Draw("SAME");
		if (type.rfind("BKG_",0)==0) { tex.SetTextSize(0.045); tex.DrawLatex(0.78, 0.50, "HYDJET"); }
		else { tex.SetTextSize(0.045); tex.DrawLatex(0.78, 0.50, pLbl.c_str()); }
	      }
	      tex.SetTextSize(0.035); tex.DrawLatex(0.78, 0.79-yPos, rapLbl.c_str());
	      std::string projLbl;
	      if (type.find("_NoDeDx")!=std::string::npos) projLbl = "Phase-2 reco.";
	      else if (type.find("_WithDeDx")!=std::string::npos) projLbl = "Phase-2 reco. with dE/dx cut";
	      if (projLbl!="") tex.SetTextSize(0.045); tex.DrawLatex(0.20, 0.84, projLbl.c_str());
	    }
	    else if (obs.first=="v2" || obs.first=="v3") {
	      legV.push_back(std::unique_ptr<TLegend>(new TLegend(0.2, 0.50, 0.35, 0.82)));
	      for (size_t i=0; i<graphs["PbPb"].size(); i++) {
		auto& gr = *std::next(graphs["PbPb"].begin(), i);
		const auto& bin = *gr.first.begin();
		const auto fmt = (bin.first=="cent" ? "%.0f-%.0f" : "%.1f-%.1f")+VAR.at(bin.first).second+"%";
		legV[0]->AddEntry(&wGraphs["PbPb"][i], Form(fmt.c_str(), bin.second.first, bin.second.second), "pf")->SetTextSize(0.025);
	      }
	      legV[0]->Draw("SAME");
	      tex.SetTextSize(0.035); tex.DrawLatex(0.78, 0.79-yPos, rapLbl.c_str());
	      tex.SetTextSize(0.045); tex.DrawLatex(0.2, 0.84, ("Pb-Pb #sqrt{s_{NN}} = "+enPbPb).c_str());
	      tex.SetTextSize(0.045); tex.DrawLatex(0.4, 0.74, pLbl.c_str());
	    }
	    else if (obs.first=="ratio") {
	      legV.push_back(std::unique_ptr<TLegend>(new TLegend(0.7, 0.15, 0.88, 0.45)));
	      legV[0]->SetHeader(("Pb-Pb #sqrt{s_{NN}} = "+enPbPb).c_str());
	      dynamic_cast<TLegendEntry*>(legV[0]->GetListOfPrimitives()->First())->SetTextSize(0.03);
	      for (size_t i=0; i<graphs["PbPb"].size(); i++) {
		auto& gr = *std::next(graphs["PbPb"].begin(), i);
		const auto& bin = *gr.first.begin();
		const auto fmt = (bin.first=="cent" ? "%.0f-%.0f" : "%.1f-%.1f")+VAR.at(bin.first).second+"%";
		legV[0]->AddEntry(&wGraphs["PbPb"][i], Form(fmt.c_str(), bin.second.first, bin.second.second), "pf")->SetTextSize(0.025);
	      }
	      legV[0]->Draw("SAME");
	      tex.SetTextSize(0.035); tex.DrawLatex(0.78, 0.79-yPos, rapLbl.c_str());
	    }

	    std::vector<TF1> fitV;
	    if (!fitM[col.first][chg.first][par.first][obs.first][r.first].empty()) {
	      const auto& fM = fitM.at(col.first).at(chg.first).at(par.first).at(obs.first).at(r.first);
	      size_t i(0);
	      for (const auto& f : fM) {
		if (!f.second) continue;
		fitV.emplace_back(*f.second);
		if (isYield) {
		  setFunction(fitV.back());
		  fitV.back().SetParameters(f.second->GetParameters());
		  fitV.back().SetParameter(0, f.second->GetParameter(0)*std::pow(2., graphs["PbPb"].size()-1-i));
		  i += 1;
		}
	      }
	    }
	    if (!fitM[colPP][chg.first][par.first][obs.first][r.first].empty()) {
	      const auto& fitPP = fitM.at(colPP).at(chg.first).at(par.first).at(obs.first).at(r.first);
	      for (const auto& f : fitPP) if (f.second) fitV.push_back(*f.second);
	    }
	    for (auto& f : fitV) {
	      f.SetLineColor(kBlack);
	      f.Draw("SAMEL");
	    }
    
	    CMS_lumi(&c, 33, "", "", false, -1, false);

	    if (isYield) c.SetLogy();
	    c.Update();

	    const auto dir = "Plot/"+type.substr(0, type.find("_"));
	    gSystem->mkdir(dir.c_str(), kTRUE);
	    c.SaveAs((dir+"/"+type+"_"+obs.first+"_"+chg.first+par.first+"_"+col.first+"_"+r.first.str(1)+".png").c_str());
	  }
	}
      }
    }
  }
};


void plotGraphs(GRAPH& input, const std::string& type)
{
  FITFUNC fM; 
  plotGraphs(input, type, fM);
};


#endif // #ifndef LightNuclei_Fit_utils_h
