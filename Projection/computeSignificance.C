#ifndef computeSignificance_C
#define computeSignificance_C

// Custom headers
#include "../Macros/PIDSelector.h"
#include "../Macros/PIDCutParam.h"
#include "../Macros/TreeReader.h"
#include "utils.h"
// ROOT headers
#include "TGraph.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TVectorF.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TH3F.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooKeysPdf.h"
#include "RooGaussian.h"
#include "RooPlot.h"


// Type definition
typedef std::map< std::string , std::map< std::string , std::map< std::string , std::unique_ptr<RooDataSet> > > > DSCONT;


const std::string DIR = "/Users/andre/Analysis/MTD/TREE/";
const std::map<std::string, std::map<std::string, std::string > > FILES =
  {
   {"BKG", {{"PbPb_5p5TeV", DIR+"particlesForBkg_PbPb.root"},
	    {"pp_14TeV", DIR+"particlesForBkg_pp.root"}}
   },
   {"SIG", {{"PbPb_5p5TeV", DIR+"particlesForSignal_PbPb.root"},
	    {"pp_14TeV", DIR+"particlesForSignal_pp.root"}}
   }
  };
const std::map<std::string, std::map<std::string, std::vector<std::string> > > TREE =
  {
   {"PbPb_5p5TeV", {{"Deuteron", {"Deuteron_FlatPt7", "Deuteron_FlatP20"}}, {"Helium3", {"Helium3_FlatPt8", "Helium3_FlatP20"}}, {"Helium4", {"Helium4_FlatP20"}}, {"Triton", {"Triton_FlatP20"}}, {"MinBias", {"MinBias"}}}},
   {"pp_14TeV", {{"Deuteron", {"Deuteron_FlatPt7", "Deuteron_FlatP20"}}, {"Helium3", {"Helium3_FlatPt8", "Helium3_FlatP20"}}, {"Helium4", {"Helium4_FlatP20"}}, {"Triton", {"Triton_FlatP20"}}, {"MinBias", {"MinBiasPU200"}}}}
  };
const std::map<std::string, std::map<std::string, std::string> > CUT =
  {
   {"PbPb_5p5TeV",
    {
     {"Deuteron", "(int(cutF)&3)==3 && pT*cosh(eta)<6"},
     {"Triton",   "(abs(eta)<1.6 ? pT*cosh(eta)<6 : pT>2) && (int(cutF)&3)==3 && pT*cosh(eta)<8"},
     // OPTIMAL TO REMOVE ALL BKG
     {"Helium3", "dEdx>=(abs(eta)<1.6 ? (pT*cosh(eta)<2.5 ? 7.8 : 6.2) : (pT*cosh(eta)>23 ? 4.9 : 3.9)) && (int(cutF)&3)==3"},
     {"Helium4", "dEdx>=(abs(eta)<1.6 ? (pT*cosh(eta)<2.5 ? 8.3 : 6.2) : (pT<1.6 ? 3.8 : 3.3)) && (int(cutF)&3)==3"},
     {"AntiHelium3", "dEdx>=(abs(eta)<1.6 ? 5.9 : 3.8) && (int(cutF)&3)==3"},
     {"AntiHelium4", "dEdx>=(abs(eta)<1.6 ? 5.9 : 3.8) && (int(cutF)&3)==3"}
    }
   },
   {"pp_14TeV",
    {
     {"Deuteron", "abs(eta)<1.6 && (int(cutF)&3)==3 && pT*cosh(eta)<6"},
     {"Triton",   "abs(eta)<1.6 && (int(cutF)&3)==3 && pT*cosh(eta)<3"},
     // OPTIMAL TO REMOVE ALL BKG
     {"Helium3", "dEdx>=(abs(eta)<1.6 ? (pT*cosh(eta)<2.6 ? 8.5 : 6.2) : (pT<0.6 ? 3.6 : 3.0)) && (int(cutF)&3)==3"},
     {"Helium4", "dEdx>=(abs(eta)<1.6 ? (pT*cosh(eta)<2.6 ? 7.7 : 0.0) : (pT<0.5 ? 4.3 : 3.7)) && (int(cutF)&3)==3"},
     {"AntiHelium3", "dEdx>=3.8 && (int(cutF)&3)==3"},
     {"AntiHelium4", "dEdx>=(abs(eta)<1.6 ? 6.9 : 3.4) && (int(cutF)&3)==3"}
    }
   }
  };

const std::map<std::string, int> CATID = {{"MinBias", 0}, {"Deuteron", 1}, {"Triton", 2}, {"Helium3", 3}, {"Helium4", 4}};
const std::map<std::string, UChar_t> CUTLABEL({{"RECO", 1}, {"dEdxCut", 3}, {"MTDCut", 7}});


// Functions
void iniSelector(std::map<std::string, PIDSelector>& selM, std::map<std::string, float>& varM, const DATA_CONTBIN& data, const float& thr);
void extractDS(DSCONT& ds, const std::string& type, const DATA_CONTBIN& data, std::map<std::string, float>& varM, std::map<std::string, PIDSelector>& selM, const DSCONT& dsSig={});
void computeSignificanceWithinSignalMassSqrRange(GRAPH& gM, const std::string& sel, const DSCONT& dsSigM, const DSCONT& dsBkgM, const DATA_CONTBIN& data, BOOLCONT missing);
void computeSignificanceFromFit(GRAPH& gM, const std::string& sel, const DSCONT& dsSigM, const DSCONT& dsBkgM, const DATA_CONTBIN& data, BOOLCONT missing);


void computeSignificance(const std::string& sel="dEdxCut")
{
  // Suppress Messages for RooFit
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Integration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  // Extract significance
  GRAPH g_SIGNI_Fit, g_SIGNI_Count;
  auto m_SIGNI_Count = extractGraph(g_SIGNI_Count, "SIGNICOUNT", sel, binning_ProjectedData);
  auto m_SIGNI_Fit = extractGraph(g_SIGNI_Fit, "SIGNIFIT", sel, binning_ProjectedData);

  if (!m_SIGNI_Count.empty() || !m_SIGNI_Fit.empty()) {
    // Initialize PID selector
    std::map<std::string, float> varMap;
    std::map<std::string, PIDSelector> selMap;
    iniSelector(selMap, varMap, binning_ProjectedData, 1);

    // Extract signal and background yields from simulation
    DSCONT d_SIG, d_BKG;
    extractDS(d_SIG, "SIG", binning_ProjectedData, varMap, selMap);
    extractDS(d_BKG, "BKG", binning_ProjectedData, varMap, selMap, d_SIG);
    return;

    // Compute significance
    computeSignificanceWithinSignalMassSqrRange(g_SIGNI_Count, sel, d_SIG, d_BKG, binning_ProjectedData, m_SIGNI_Count);
    computeSignificanceFromFit(g_SIGNI_Fit, sel, d_SIG, d_BKG, binning_ProjectedData, m_SIGNI_Fit);
  }
};


void iniSelector(std::map<std::string, PIDSelector>& selM, std::map<std::string, float>& varM, const DATA_CONTBIN& data, const float& thr)
{
  std::set<std::string> pars;
  for (const auto& col : data) {
    for (const auto& chg : col.second) {
      for (const auto& par : chg.second) pars.insert(par.first);
    }
  }
  for (const auto& cut : std::map<std::string, std::pair<std::string, float> >({{"dEdx", {"dEdx", thr}}, {"MTD", {"invBeta", thr}}, {"MTDTight", {"invBeta", 2}}})) {
    for (const auto& etaReg : std::vector<std::string>({"BTL", "ETL"})) {
      if (FUNCMAP_.find(etaReg)==FUNCMAP_.end()) throw std::runtime_error("[ERROR] Region "+etaReg+" not in FUNCMAP_!");
      const BINF bin = { {"absEta", (etaReg=="BTL" ? std::make_pair(0.0, 1.6) : std::make_pair(1.6, 3.0))} };
      for (const auto& par : pars) {
	if (FUNCMAP_.at(etaReg).find(par)==FUNCMAP_.at(etaReg).end()) throw std::runtime_error("[ERROR] Particule "+par+" not in FUNCMAP_!");
	if (FUNCMAP_.at(etaReg).at(par).find(cut.second.first)==FUNCMAP_.at(etaReg).at(par).end()) throw std::runtime_error("[ERROR] Variable "+cut.second.first+" not in FUNCMAP_!");
	const auto& func = FUNCMAP_.at(etaReg).at(par).at(cut.second.first);
	selM[cut.first].setThreshold(cut.second.second, bin, par, "FNC");
	selM[cut.first].setFunction(func.at("Mean"), func.at("Error"), bin, par, cut.second.first, "p", {0.7, 100.0});
      }
    }
    for (const auto& n : std::vector<std::string>({"absEta", "p", cut.second.first})) {
      selM[cut.first].setVariable(n, &varM[n]);
    }
  }
};


bool getDS(DSCONT& ds, TFile& file, const std::array<std::string, 3>& info)
{
  if (!file.IsOpen() || file.IsZombie()) return false;
  for (const auto& chg : {"", "Anti"}) {
    const auto name = info[0]+"_"+info[1]+"_"+chg+info[2];
    const auto& input_ds = dynamic_cast<RooDataSet*>(file.Get(name.c_str()));
    if (!input_ds) {
      std::cout << "[INFO] RooDataSet " << info[0]<< " , " << info[1] << " , " << chg << " , " << info[2] << " not found in " << file.GetName() << std::endl;
      return false;
    }
    ds[info[1]][chg][info[2]].reset(input_ds);
  }
  std::cout << "[INFO] Found all " << info[0] << " RooDataSets in " << file.GetName() << std::endl;
  return true;
};


void extractDS(DSCONT& ds, const std::string& type, const DATA_CONTBIN& data, std::map<std::string, float>& varM, std::map<std::string, PIDSelector>& selM, const DSCONT& dsSig)
{
  // Define RooRealVars
  RooRealVar pTV ("pT",        "p_{T}",   -1.0,     1000.0, "");
  RooRealVar etaV("eta",        "#eta", -100.0,      100.0, "");
  RooRealVar rapV("rap",           "y", -100.0,      100.0, "");
  RooRealVar cent("cent", "centrality",  -1.0,       300.0, "");
  RooRealVar cutF("cutF", "cut flag" ,   -1.0,      1000.0, "");
  RooRealVar wght("weight", "weight" ,   -1.0,       1.E16, "");
  RooRealVar invBeta("invBeta", "TOF 1/#beta",  -1.0,     10.0, "");
  RooRealVar dEdx   ("dEdx",          "dE/dx",  -1.0,   1000.0, "MeV/cm");
  RooRealVar mSqr("mSqr", "TOF m^{2}", -100.0,  1.E6, "GeV^{2}/c^{4}");

  const bool isBkg = (type=="BKG" || type=="MINBIAS");
  const auto typeLbl = (isBkg ? "BKG" : type);
  // Create RooDataSets
  for (const auto& col : data) {
    const auto sqrtS = getBeamEnergy(col.first);
    std::set<std::string> missing;
    for (const auto& par : col.second.at("")) {
      if (type=="MINBIAS") { missing.insert("MINBIAS"); break; }
      if (par.second.find("yield")==par.second.end()) continue;
      // Extract RooDataSets from file
      const auto fileName = "Output/FILE_CUTWEIGHTED/"+typeLbl+"_"+par.first+"_"+col.first+"_CUTWEIGHTED.root";
      TFile inputFile(fileName.c_str(), "READ");
      const auto found = getDS(ds, inputFile, {typeLbl, col.first, par.first});
      inputFile.Close();
      if (found) continue;
      missing.insert(par.first);
      // Initialize RooDataSets
      for (const auto& chg : {"", "Anti"}) {
	const auto name = typeLbl+"_"+col.first+"_"+chg+par.first;
	ds[col.first][chg][par.first].reset(new RooDataSet(name.c_str(), "", {pTV, etaV, cent, mSqr, dEdx, cutF, wght}));
      }
    }
    if (missing.empty()) continue;
    // Loop over particles
    for (const auto& col : FILES.at(typeLbl)) {
      if (data.find(col.first)==data.end()) { std::cout << "[WARNING] System "+col.first+" not present in projection data!" << std::endl; continue; }
      if (data.at(col.first).find("")==data.at(col.first).end()) throw std::runtime_error("[ERROR] Charge not present in "+col.first+" in projection data!");
      std::set<std::string> missing2;
      for (const auto& par : (isBkg ? std::set<std::string>({"MinBias"}) : missing)) {
	// Extract RooDataSets from file
	const auto fileName = "Output/FILE/"+typeLbl+"_"+par+"_"+col.first+".root";
	TFile inputFile(fileName.c_str(), "READ");
	const auto found = getDS(ds, inputFile, {typeLbl, col.first, par});
	inputFile.Close();
	if (found) continue;
	missing2.insert(par);
	// Initialize RooDataSets
	for (const auto& chg : {"", "Anti"}) {
	  const auto name = typeLbl+"_"+col.first+"_"+chg+par;
	  ds[col.first][chg][par].reset(new RooDataSet(name.c_str(), "", {pTV, etaV, cent, invBeta, dEdx, cutF}));
	}
      }
      if (missing2.empty()) continue;
      for (const auto& par : missing2) {
	const int catIdx = std::pow(2, 5+CATID.at(par));
	// Fill RooDataSets using tree
	for (const auto& tree : TREE.at(col.first).at(par)) {
	  TreeReader reader(col.second.c_str(), tree, false);
	  const auto& nentries = reader.getEntries();
	  // Loop over the events
	  std::cout << "[INFO] Starting to process " << nentries << " " << par << " nentries for " << col.first << std::endl;
	  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
	    // Get the entry in the trees
	    reader.setEntry(jentry);
	    if(!(jentry % (isBkg?1000000:100000))) std::cout << "Processed " << jentry << " " << par << " events out of " << nentries << " for " << col.first << std::endl;
	    // Extract information in detector coordinates
	    const auto& eta = reader.getFloat("eta");
	    const auto& p = reader.getFloat("p");
	    const auto pT = p/std::cosh(eta);
	    // Apply selection
	    if (isBkg && reader.getFloat("beta")<=0.) continue;
	    if ((std::abs(eta)<1.4) ? pT <= 0.8 : p <= 0.7) continue;
	    if (std::abs(eta)>=3.0) continue;
	    // Set RooRealVars
	    pTV.setVal(pT);
	    etaV.setVal(eta);
	    cent.setVal(reader.getInt("cent")/2.);
	    const auto& betaVal = reader.getFloat("beta");
	    invBeta.setVal(betaVal>0. ? 1.0/betaVal : -1.);
	    dEdx.setVal(reader.getFloat("dEdx0"));
	    cutF.setVal(catIdx);
	    //
	    // Fill entry
	    const auto chg = (reader.getInt("charge")>0 ? "" : "Anti");
	    ds[col.first][chg][par]->addFast({pTV, etaV, cent, invBeta, dEdx, cutF});
	  }
	}
      }
      // Store RooDataSets in file
      for (const auto& par : missing2) {
	gSystem->mkdir("Output/FILE/", kTRUE);
	const auto fileName = "Output/FILE/"+typeLbl+"_"+par+"_"+col.first+".root";
	TFile outputFile(fileName.c_str(), "RECREATE");
	std::cout << "[INFO] Proceed to store RooDataSets in " << fileName << std::endl;
	if (outputFile.IsOpen() && !outputFile.IsZombie()) {
	  for (const auto& chg : ds.at(col.first)) {
	    const auto& d = chg.second.at(par);
	    d->Write(d->GetName());
	  }
	  outputFile.Write();
	  outputFile.Close();
	  std::cout << "[INFO] Stored RooDataSets in " << fileName << std::endl;
	}
      }
    }
    // Set weights
    float nMinBias(0);
    if (isBkg) {
      TFile file(FILES.at("BKG").at(col.first).c_str(), "READ");
      const auto lbl = "Entries_"+TREE.at(col.first).at("MinBias")[0];
      if (!file.IsOpen() || file.IsZombie() || !file.Get(lbl.c_str())) continue;
      const auto& v = *dynamic_cast<TVectorF*>(file.Get(lbl.c_str()));
      nMinBias = v[lbl.rfind("PU")!=std::string::npos ? 1 : 0];
      std::cout << "MINBIAS for " << col.first << " , " << lbl << " : " << v[0] << " , " << v[1] << " >> " << nMinBias << std::endl;
    }
    std::cout << "[INFO] Setting weights to RooDataSet for " << type << std::endl;
    for (auto& chg : ds.at(col.first)) {
      for (auto& par : (isBkg ? std::set<std::string>({"MinBias"}) : missing)) {
	auto& ds = chg.second.at(par);
	if (type=="SIG") {
	  // Step 1) Fill histogram
	  std::cout << "[INFO] Filling histogram for " << type << " , " << col.first << " , " << chg.first+par << std::endl;
	  TH3F hist("h", "h", 100, 0., 20., 31, -3.1, 3.1, 50, 0., 1.);
	  const auto binVolume = hist.GetXaxis()->GetBinWidth(0) * hist.GetYaxis()->GetBinWidth(0) * hist.GetZaxis()->GetBinWidth(0);
	  for (int j=0; j<ds->numEntries(); j++) {
	    const auto& set = *ds->get(j);
	    TLorentzVector p4; p4.SetPtEtaPhiM(set.getRealValue("pT")*CHARGE.at(par), set.getRealValue("eta"), 0, MASS.at(par));
	    hist.Fill(p4.Pt(), p4.Rapidity(), set.getRealValue("cent")/100.);
	  }
	  // Step 2) Divide projected yield by histogram content
	  std::cout << "[INFO] Setting weight for " << type << " , " << col.first << " , " << chg.first+par << std::endl;
	  const auto charge = (chg.first=="Anti" ? -1 : 1)*CHARGE.at(par);
	  std::unique_ptr<TF3> projYield3D; std::unique_ptr<TF2> projYield2D;
	  if (col.first.rfind("PbPb_", 0)==0) projYield3D.reset(getProjection3D(col.first, "yield", MASS.at(par), charge));
	  else if (col.first.rfind("pp_", 0)==0) projYield2D.reset(getProjection2D(col.first, "yield", MASS.at(par), charge));
	  const auto tmpDS = new RooDataSet(ds->GetName(), ds->GetTitle(), {pTV, etaV, cent, invBeta, dEdx, cutF, wght}, RooFit::WeightVar(wght));
	  for (int j=0; j<ds->numEntries(); j++) {
	    const auto& set = *ds->get(j);
	    TLorentzVector p4; p4.SetPtEtaPhiM(set.getRealValue("pT")*CHARGE.at(par), set.getRealValue("eta"), 0, MASS.at(par));
	    const auto pT = p4.Pt();
	    const auto rap = p4.Rapidity();
	    const auto cent = col.first.rfind("PbPb_", 0)==0 ? set.getRealValue("cent")/100. : 0.5;;
	    // Get projected yield (project to 1B minbias events)
	    double proY(0);
	    if (col.first.rfind("PbPb_", 0)==0) {
	      proY = projYield3D->Eval(hist.GetXaxis()->GetBinCenter(hist.GetXaxis()->FindBin(pT)),
				       hist.GetYaxis()->GetBinCenter(hist.GetYaxis()->FindBin(rap)),
				       hist.GetZaxis()->GetBinCenter(hist.GetZaxis()->FindBin(cent)));
	    }
	    else if (col.first.rfind("pp_", 0)==0) {
	      proY = projYield2D->Eval(hist.GetXaxis()->GetBinCenter(hist.GetXaxis()->FindBin(pT)),
				       hist.GetYaxis()->GetBinCenter(hist.GetYaxis()->FindBin(rap)));
	    }
	    const auto meaY = hist.GetBinContent(hist.FindBin(pT, rap, cent));
	    if (std::isnan(meaY)) throw std::runtime_error("[ERROR] NAN expected entries!");
	    if (std::isnan(binVolume)) throw std::runtime_error("[ERROR] NAN binVolumen!");
	    if (std::isnan(proY)) throw std::runtime_error(Form("[ERROR] NAN projected entries: %g , %g , %g >> %g , %g , %g!", pT, rap, cent,
								hist.GetXaxis()->GetBinCenter(hist.GetXaxis()->FindBin(pT)), hist.GetYaxis()->GetBinCenter(hist.GetYaxis()->FindBin(rap)), hist.GetZaxis()->GetBinCenter(hist.GetZaxis()->FindBin(cent))));
	    const auto weight = (meaY>0 ? (proY * binVolume * 1.E9 / meaY) : 0);
	    // Fill weight
	    tmpDS->addFast(set, weight);
	  }
	  ds.reset(tmpDS);
	}
	else if (isBkg) {
	  const auto tmpDS = new RooDataSet(ds->GetName(), ds->GetTitle(), {pTV, etaV, cent, invBeta, dEdx, cutF, wght}, RooFit::WeightVar(wght));
	  // Project to 1B minbias events
	  for (int j=0; j<ds->numEntries(); j++) {
	    tmpDS->addFast(*ds->get(j), 1.E9/nMinBias);
	  }
	  // Merge HYDJET + light nuclei samples to form full background sample
	  std::cout << "[INFO] Background events before merging with light nuclei: " << tmpDS->numEntries() << " , " << tmpDS->sumEntries() << std::endl;
	  for (const auto& c : CATID) {
	    if (c.first=="MinBias") continue;
	    const auto& sigDS = dsSig.at(col.first).at(chg.first).at(c.first);
	    std::cout << "APPENDING: " << c.first << " , " << sigDS->numEntries() << " , " << sigDS->sumEntries() << std::endl;
	    tmpDS->append(*sigDS);
	  }
	  std::cout << "[INFO] Background events after merging with light nuclei: " << tmpDS->numEntries() << " , " << tmpDS->sumEntries() << std::endl;
	  ds.reset(tmpDS);
	}
      }
    }
    if (type=="MINBIAS") continue;
    // Set cut selections
    std::cout << "[INFO] Setting cut selections to RooDataSet for " << type << std::endl;
    for (auto& par : missing) {
      const int catIdx = std::pow(2, 5+CATID.at(par));
      for (auto& chg : ds.at(col.first)) {
	auto& ds = chg.second.at(par);
	auto& dp = chg.second.at(isBkg ? "MinBias" : par);
	std::cout << "[INFO] Setting cut selections for " << type << " , " << col.first << " , " << chg.first+par << " , " << ds->GetName() << std::endl;
	const auto tmpDS = new RooDataSet(ds->GetName(), ds->GetTitle(), {pTV, etaV, rapV, cent, mSqr, cutF, dEdx, wght}, RooFit::WeightVar(wght));
	for (int j=0; j<dp->numEntries(); j++) {
	  const auto& set = *dp->get(j);
	  // Exclude signal events in background sample
	  const auto& cat = set.getRealValue("cutF", -1);
	  const auto isSignal = (int(cat) & catIdx)>0;
	  if (isBkg && isSignal) continue;
	  const auto& invBeta  = set.getRealValue("invBeta", -1.);
	  if (isBkg && invBeta<1) continue;
	  // Variables for selector
	  const auto& eta = set.getRealValue("eta", -1.);
	  const auto& pT = set.getRealValue("pT", -1.);
	  varM["absEta"] = std::abs(eta);
	  varM["p"] = pT*std::cosh(eta);
	  varM["invBeta"] = set.getRealValue("invBeta", -1.);
	  varM["dEdx"] = set.getRealValue("dEdx", -1.);
	  // Set cut variable
	  int cutV(0);
	  if (varM["invBeta"]>0.) cutV += 1; // MTD reco
	  if (selM.at("dEdx").passSelector(selM.at("dEdx").findBin(), par, +1)) cutV += 2; // dEdx cut
	  //if (varM["invBeta"]>0. && selM.at("MTD").passSelector(selM.at("MTD").findBin(), par)) cutV += 4; // MTD cut
	  if (varM["invBeta"]>0. && selM.at("MTDTight").passSelector(selM.at("MTDTight").findBin(), par)) cutV += 8; // MTD tight cut
	  // Set RooRealVars using particle coordinates
	  TLorentzVector p4; p4.SetPtEtaPhiM(CHARGE.at(par)*pT, eta, 0, MASS.at(par));
	  pTV.setVal(p4.Pt());
	  etaV.setVal(p4.Eta());
	  rapV.setVal(p4.Rapidity());
	  cent.setVal(set.getRealValue("cent"));
	  mSqr.setVal(varM["invBeta"]>0. ? p4.P()*p4.P()*(varM["invBeta"]*varM["invBeta"] - 1.) : -99.);
	  cutF.setVal(cutV);
	  dEdx.setVal(varM["dEdx"]);
	  // Fill RooDataSet
	  tmpDS->addFast({pTV, etaV, rapV, cent, mSqr, cutF, dEdx, wght}, dp->weight());
	}
	ds.reset(tmpDS);
      }
      // Store RooDataSets in file
      gSystem->mkdir("Output/FILE_CUTWEIGHTED/", kTRUE);
      const auto fileName = "Output/FILE_CUTWEIGHTED/"+typeLbl+"_"+par+"_"+col.first+"_CUTWEIGHTED.root";
      TFile outputFile(fileName.c_str(), "RECREATE");
      if (outputFile.IsOpen() && !outputFile.IsZombie()) {
	for (const auto& chg : ds.at(col.first)) {
	  const auto& d = chg.second.at(par);
	  d->Write(d->GetName());
	}
	outputFile.Write();
	outputFile.Close();
	std::cout << "[INFO] Stored RooDataSets in " << fileName << std::endl;
      }
    }
    for (auto& chg : ds.at(col.first)) chg.second["MinBias"].reset();
  }
};


void computeSignificanceWithinSignalMassSqrRange(GRAPH& gM, const std::string& sel, const DSCONT& dsSigM, const DSCONT& dsBkgM, const DATA_CONTBIN& data, BOOLCONT missing)
{
  std::cout << "[INFO] Start to compute significance" << std::endl;
  for (const auto& col : data) {
    for (const auto& chg : col.second) {
      for (const auto& par : chg.second) {
	const auto ch = chg.first!="All" ? chg.first : "";
	auto dsSig = *dsSigM.at(col.first).at(ch).at(par.first);
	auto dsBkg = *dsBkgM.at(col.first).at(ch).at(par.first);
	if (chg.first=="All") {
	  dsSig.append(*dsSigM.at(col.first).at("Anti").at(par.first));
	  dsBkg.append(*dsBkgM.at(col.first).at("Anti").at(par.first));
	}
	if (!dsSig.get()->find("mSqr")) throw std::runtime_error("[ERROR] mSqr variable not found");
	auto var = *dynamic_cast<RooRealVar*>(dsSig.get()->find("mSqr"));
	var.setRange(-10., 40.); var.setBins((var.getMax()-var.getMin())*100/2.);
	for (const auto& obs : par.second) {
	  if (!missing[col.first][chg.first][par.first][obs.first]) continue;
	  const auto pN = CUT.at(col.first).find(chg.first+par.first)!=CUT.at(col.first).end() ? chg.first+par.first : par.first;
	  const auto cut = CUT.at(col.first).at(pN);
	  std::cout << "[INFO] Applying cut: " << cut << std::endl;
	  auto dsSig_sel = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsSig.reduce(cut.c_str())));
	  auto dsBkg_sel = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsBkg.reduce(cut.c_str())));
	  for (const auto& r : obs.second) {
	    auto dsSig_cutR = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsSig_sel->reduce(r.first.dsCut().c_str())));
	    auto dsBkg_cutR = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsBkg_sel->reduce(r.first.dsCut().c_str())));
	    auto bin = r.first;
	    for (const auto& c : r.second) {
	      if (c.second.empty()) continue;
	      auto dsSig_cutC = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsSig_sel->reduce(c.first.dsCut().c_str())));
	      auto dsBkg_cutC = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsBkg_sel->reduce(c.first.dsCut().c_str())));
	      bin.add(c.first);
	      const auto& xVar = c.second.begin()->first;
	      const auto& bins = c.second.begin()->second;
	      if (bins.size()<2) continue;
	      auto& graph = gM[col.first][obs.first][r.first][c.first][chg.first][par.first];
	      graph.Set(bins.size());
	      graph.SetTitle(xVar.c_str());
	      for (size_t i=1; i<bins.size(); i++) {
		bin[xVar] = {bins[i-1], bins[i]};
		auto dsProjSig = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsSig_cutC->reduce(bin.dsCut().c_str())));
		auto dsProjBkg = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsBkg_cutC->reduce(bin.dsCut().c_str())));
		// Compute significance as S/sqrt(S+B)
		auto dsProjSigCut = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsProjSig->reduce("((int(cutF) & 8)==8)")));
		std::cout << col.first << " , " << chg.first+par.first << " , " << bin.str() << " >> " <<  dsProjSigCut->numEntries() << std::endl;
		double signi(0.);
		if (dsProjSigCut->numEntries()>1) {
		  // Derive signal correction factor (correcting for lack of phase space in MC generation at high p)
		  auto binForCorr = bin;
		  binForCorr["pT"] = {std::max(CHARGE.at(par.first)*0.9f, bin["pT"].first), std::max(CHARGE.at(par.first)*1.0f, bin["pT"].second)};
		  binForCorr["absRap"] = {std::min(2.9f, bin["absRap"].first), std::min(3.0f, bin["absRap"].second)};
		  const auto range = std::array<double, 6>({binForCorr["pT"].first, binForCorr["pT"].second, binForCorr["absRap"].first, binForCorr["absRap"].second, bin["cent"].first/100., bin["cent"].second/100.});
		  const auto charge = (chg.first=="Anti" ? -1 : 1)*CHARGE.at(par.first);
		  const auto genEntries = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsSig.reduce(binForCorr.dsCut().c_str())))->sumEntries();
		  auto modelEntries = getProjection(col.first, "yield", range, MASS.at(par.first), charge, false, true);
		  if (chg.first=="All") modelEntries += getProjection(col.first, "yield", range, MASS.at(par.first), -charge, false, true);
		  modelEntries *= 1.E9;
		  const auto corrSig = modelEntries/genEntries;
		  // Compute the mean and RMS of signal m^2 spectrum
		  const auto mean = dsProjSigCut->mean(var);
		  const auto rms = std::sqrt(dsProjSigCut->moment(var, 2, mean));
		  // Compute signal and background entries within signal m^2 region (+- 1-sigma)
		  const std::string mRange = Form("mSqr>=%g && mSqr<=%g", mean-rms, mean+rms);
		  const auto dsBkgSigRange = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsProjBkg->reduce(mRange.c_str())));
		  const auto nBkg = dsBkgSigRange->sumEntries();
		  if (par.first.rfind("Helium",0)==0 && nBkg>0) {
		    for (int i=0; i<dsBkgSigRange->numEntries(); i++) {
		      dsBkgSigRange->get(i); dsBkgSigRange->Print("v");
		    }
		  }
		  const auto dsSigSigRange = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsProjSig->reduce(mRange.c_str())));
		  auto nSig = dsSigSigRange->sumEntries();
		  if (nBkg==0) {
		    for (size_t i=6; i<=50; i++) {
		      const std::string mRange = Form("mSqr>=%g && mSqr<=%g", mean-i*0.2*rms, mean+i*0.2*rms);
		      const auto nSig2 = corrSig * std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsProjSig->reduce(mRange.c_str())))->sumEntries();
		      if (std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsProjBkg->reduce(mRange.c_str())))->numEntries()>0) break;
		      nSig = nSig2;
		    }
		  }
		  nSig *= corrSig;
		  signi = nSig/std::sqrt(nSig + nBkg);
		  if (nBkg>0 && mean < 3*rms) signi *= -1; // set negative when mean is consistent with 0
		  std::cout << "ALL NSIG: " << dsProjSig->sumEntries() << " , NBKG: " << dsProjBkg->sumEntries() << " , " << genEntries << " , " << modelEntries << std::endl;
		  std::cout << "MEAN: " << mean << " , RMS: " << rms << " , NSIG: " << nSig << " , NBKG: " << nBkg << " , SIGNI: " << signi << " , CORR FACTOR: " << corrSig << " ( " << binForCorr.str() << " ) " << std::endl;
		}
		// Store result in graph
		const auto xVal = (bin.at(xVar).second + bin.at(xVar).first)/2.;
		const auto xErr = roundVal(bin.at(xVar).second - bin.at(xVar).first, 4)/2.;
		graph.SetPoint(i-1, xVal, signi);
		graph.SetPointError(i-1, xErr, xErr, 1., 1.);
	      }
	    }
	  }
	  // Store graphs in file
	  gSystem->mkdir("Output/GRAPH_SIGNICOUNT/", kTRUE);
	  const auto fileName = "Output/GRAPH_SIGNICOUNT/GRAPH_SIGNICOUNT_"+chg.first+par.first+"_"+col.first+"_"+obs.first+"_"+sel+".root";
	  TFile outputFile(fileName.c_str(), "RECREATE");
	  if (outputFile.IsOpen() && !outputFile.IsZombie()) {
	    for (const auto& r : gM.at(col.first).at(obs.first)) {
	      for (const auto& c : r.second) {
		auto& graph = c.second.at(chg.first).at(par.first);
		const auto name = "graph_"+col.first+"_"+chg.first+par.first+"_"+obs.first+"_"+r.first.str(1)+"_"+c.first.str(1);
		graph.Write(name.c_str());
	      }
	    }
	    outputFile.Write();
	    outputFile.Close();
	    std::cout << "[INFO] Stored graphs in " << fileName << std::endl;
	  }
	}
      }
    }
  }
};


void computeSignificanceFromFit(GRAPH& gM, const std::string& sel, const DSCONT& dsSigM, const DSCONT& dsBkgM, const DATA_CONTBIN& data, BOOLCONT missing)
{
  std::cout << "15" << std::endl;
  for (const auto& col : data) {
    for (const auto& chg : col.second) {
      for (const auto& par : chg.second) {
	const auto ch = chg.first!="All" ? chg.first : "";
	auto dsSig = *dsSigM.at(col.first).at(ch).at(par.first);
	auto dsBkg = *dsBkgM.at(col.first).at(ch).at(par.first);
	if (chg.first=="All") {
	  dsSig.append(*dsSigM.at(col.first).at("Anti").at(par.first));
	  dsBkg.append(*dsBkgM.at(col.first).at("Anti").at(par.first));
	}
	if (!dsSig.get()->find("mSqr")) throw std::runtime_error("[ERROR] mSqr variable not found");
	auto var = *dynamic_cast<RooRealVar*>(dsSig.get()->find("mSqr"));
	var.setRange(-10.0, 60.); var.setBins((var.getMax()-var.getMin())*100/2.);
	var.setRange("DrawRange", -1., 25.); var.setBins((var.getMax("DrawRange")-var.getMin("DrawRange"))*100/2.);
	for (const auto& obs : par.second) {
	  if (!missing[col.first][chg.first][par.first][obs.first]) continue;
	  const auto pN = CUT.at(col.first).find(chg.first+par.first)!=CUT.at(col.first).end() ? chg.first+par.first : par.first;
	  const auto cut = CUT.at(col.first).at(pN);
	  std::cout << "[INFO] Applying cut: " << cut << std::endl;
	  auto dsSig_sel = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsSig.reduce(cut.c_str())));
	  auto dsBkg_sel = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsBkg.reduce(cut.c_str())));
	  for (const auto& r : obs.second) {
	    auto dsSig_cutR = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsSig_sel->reduce(r.first.dsCut().c_str())));
	    auto dsBkg_cutR = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsBkg_sel->reduce(r.first.dsCut().c_str())));
	    auto bin = r.first;
	    for (const auto& c : r.second) {
	      if (c.second.empty()) continue;
	      auto dsSig_cutC = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsSig_sel->reduce(c.first.dsCut().c_str())));
	      auto dsBkg_cutC = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsBkg_sel->reduce(c.first.dsCut().c_str())));
	      bin.add(c.first);
	      const auto& xVar = c.second.begin()->first;
	      const auto& bins = c.second.begin()->second;
	      if (bins.size()<2) continue;
	      auto& graph = gM[col.first][obs.first][r.first][c.first][chg.first][par.first];
	      graph.Set(bins.size());
	      graph.SetTitle(xVar.c_str());
	      for (size_t i=1; i<bins.size(); i++) {
		bin[xVar] = {bins[i-1], bins[i]};
		auto dsProjSig = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsSig_cutC->reduce(bin.dsCut().c_str())));
		auto dsProjBkg = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsBkg_cutC->reduce(bin.dsCut().c_str())));
		auto dsProjSigCut = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsProjSig->reduce("((int(cutF) & 8)==8)")));
		// Compute significance as S/unc(S) from fit
		double signi(0.);
		std::cout << col.first << " , " << chg.first+par.first << " , " << bin.str() << " >> " <<  dsProjSigCut->numEntries() << std::endl;
		if (dsProjSigCut->numEntries()>5) {
		  // Fit only if mean-5*RMS is <= 0.
		  const auto meanVal = dsProjSigCut->mean(var);
		  const auto rmsVal = std::sqrt(dsProjSigCut->moment(var, 2, meanVal));
		  const auto meanBkgVal = dsProjBkg->mean(var);
		  const auto rmsBkgVal = std::sqrt(dsProjBkg->moment(var, 2, meanBkgVal));
		  const auto diff = (meanVal-meanBkgVal)/std::sqrt(rmsVal*rmsVal + rmsBkgVal*rmsBkgVal);
		  std::cout << par.first << " , " << bin.str() << std::endl;
		  std::cout << "MEAN VALUE: SIG " << meanVal << " , " << rmsVal << " , BKG: " << meanBkgVal << " , " << rmsBkgVal << " CROSS: " << diff << std::endl;
		  if (par.first=="Deuteron" ? diff>1. : diff>8.) continue;
		  // Define datasets to fit
		  const std::string fitRange = Form("mSqr >= %g && mSqr <= %g", var.getMin(), var.getMax());
		  auto dsSigFit = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsProjSig->reduce(fitRange.c_str())));
		  auto dsBkgFit = std::unique_ptr<RooDataSet>(dynamic_cast<RooDataSet*>(dsProjBkg->reduce(fitRange.c_str())));
		  // Fill Phase-2 projected dataset
		  auto dsDataFit =*dsSigFit;
		  dsDataFit.append(*dsBkgFit);
		  if (dsDataFit.numEntries()<40) { std::cout << Form("[ERROR] Total entries %d too small!", dsDataFit.numEntries()) << std::endl; continue; }
		  const auto nSig = dsSigFit->sumEntries();
		  const auto nBkg = dsBkgFit->sumEntries();
		  const auto nTot = dsDataFit.sumEntries();
		  // Create model PDF
		  RooRealVar mean("mean", "", meanVal*0.5, meanVal*1.5, ""); mean.setVal(meanVal); mean.setConstant();
		  RooRealVar sigma("sigma", "", rmsVal*0.9, rmsVal*1.1, ""); sigma.setVal(rmsVal); sigma.setConstant();
		  RooGaussian pdfSig("PDF_SIG", "", var, mean, sigma);
		  {
		    auto frame = std::unique_ptr<RooPlot>(dynamic_cast<RooPlot*>(var.frame(RooFit::Range(0., 20.))));
		    dsProjSig->plotOn(frame.get(), RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.2));
		    pdfSig.plotOn(frame.get(), RooFit::Normalization(dsProjSig->sumEntries(), RooAbsReal::NumEvent), RooFit::LineColor(kBlue), RooFit::LineStyle(2));
		    TCanvas canv("", "", 1000, 1000);
		    frame->Draw();
		    // Save plot
		    gSystem->mkdir("Plot/DSSIG/", kTRUE);
		    canv.SaveAs(Form("Plot/DSSIG/FIT_%s_%s.png", sel.c_str(), (col.first+"_"+chg.first+par.first+"_"+obs.first+"_"+bin.str(1)).c_str()));
		  }
		  RooKeysPdf pdfBkg("PDF_BKG", "", var, *dsProjBkg, RooKeysPdf::NoMirror, 0.5);
		  RooRealVar nVarSig("N_SIG", "", 0, 2*nTot, "");
		  RooRealVar nVarBkg("N_BKG", "", 0, 2*nTot, "");
		  RooAddPdf pdfTot("PDF_TOT", "", {pdfSig, pdfBkg}, {nVarSig, nVarBkg});
		  // Fit projected dataset
		  std::unique_ptr<RooFitResult> fitResult;
		  for (size_t iFit=0; iFit<3; iFit++) {
		    // Initialize fit variables
		    nVarSig.setVal(nSig*1.1);
		    nVarBkg.setVal(nBkg); nVarBkg.setConstant();
		    // Perform fit
		    std::vector<RooCmdArg> cmdL = {RooFit::Extended(1), RooFit::SumW2Error(0), RooFit::Strategy(2-iFit), RooFit::Minimizer("Minuit2"),
						   RooFit::NumCPU(24), RooFit::Save(1), RooFit::PrintLevel(-1), RooFit::Verbose(0)};
		    RooLinkedList fitConf; for (auto& cmd : cmdL) { fitConf.Add(dynamic_cast<TObject*>(&cmd)); };
		    fitResult.reset(pdfTot.fitTo(dsDataFit, fitConf));
		    const auto& edm = fitResult->edm();
		    bool fitFailed(edm==0. || isnan(edm) || edm>0.01);
		    if (!fitFailed) break;
		    fitResult->Print("v");
		    fitResult.reset();
		  }
		  // Print fit result
		  if (!fitResult) std::cout << "[WARNING] Fit failed!" << std::endl;
		  if (fitResult) fitResult->Print("v");
		  // Plot fit result
		  auto frame = std::unique_ptr<RooPlot>(dynamic_cast<RooPlot*>(var.frame(RooFit::Range("DrawRange"))));
		  dsDataFit.plotOn(frame.get(), RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.2));
		  pdfTot.plotOn(frame.get(), RooFit::Normalization(nTot, RooAbsReal::NumEvent), RooFit::LineColor(kBlack), RooFit::LineStyle(2));
		  pdfBkg.plotOn(frame.get(), RooFit::Normalization(nBkg, RooAbsReal::NumEvent), RooFit::LineColor(kRed), RooFit::LineStyle(2));
		  pdfSig.plotOn(frame.get(), RooFit::Normalization(nSig, RooAbsReal::NumEvent), RooFit::LineColor(kBlue), RooFit::LineStyle(2));
		  TCanvas canv("", "", 1000, 1000);
		  frame->Draw();
		  TLatex tex; tex.SetNDC(1); tex.SetTextSize(0.040); tex.SetTextFont(42); 
		  tex.DrawLatex(0.5, 0.87, Form("Signal: %.0f #pm %.0f", nVarSig.getVal(), nVarSig.getError()));
		  tex.DrawLatex(0.5, 0.82, Form("Background: %.0f #pm %.0f", nVarBkg.getVal(), nVarBkg.getError()));
		  if (!fitResult) tex.DrawLatex(0.5, 0.77, "FIT FAILED");
		  frame->SetMinimum(5.);
		  canv.SetLogy();
		  // Save plot
		  gSystem->mkdir("Plot/DSFIT/", kTRUE);
		  canv.SaveAs(Form("Plot/DSFIT/FIT_%s_%s.png", sel.c_str(), (col.first+"_"+chg.first+par.first+"_"+obs.first+"_"+bin.str(1)).c_str()));
		  // Derive significance
		  signi = (fitResult ? nVarSig.getVal()/nVarSig.getError() : 0.);
		  std::cout << "SIGNI: " << signi << std::endl;
		}
		// Store result in graph
		const auto xVal = (bin.at(xVar).second + bin.at(xVar).first)/2.;
		const auto xErr = roundVal(bin.at(xVar).second - bin.at(xVar).first, 4)/2.;
		graph.SetPoint(i-1, xVal, signi);
		graph.SetPointError(i-1, xErr, xErr, 1., 1.);
	      }
	    }
	  }
	  // Store graphs in file
	  gSystem->mkdir("Output/GRAPH_SIGNIFIT/", kTRUE);
	  const auto fileName = "Output/GRAPH_SIGNIFIT/GRAPH_SIGNIFIT_"+chg.first+par.first+"_"+col.first+"_"+obs.first+"_"+sel+".root";
	  TFile outputFile(fileName.c_str(), "RECREATE");
	  if (outputFile.IsOpen() && !outputFile.IsZombie()) {
	    for (const auto& r : gM.at(col.first).at(obs.first)) {
	      for (const auto& c : r.second) {
		auto& graph = c.second.at(chg.first).at(par.first);
		const auto name = "graph_"+col.first+"_"+chg.first+par.first+"_"+obs.first+"_"+r.first.str(1)+"_"+c.first.str(1);
		graph.Write(name.c_str());
	      }
	    }
	    outputFile.Write();
	    outputFile.Close();
	    std::cout << "[INFO] Stored graphs in " << fileName << std::endl;
	  }
	}
      }
    }
  }
};


#endif // #ifndef computeSignificance_C
