#ifndef PIDSelector_h
#define PIDSelector_h

#include "TF1.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodBase.h"

#include <iostream>
#include <string>
#include <map>


class BIN {
 public:
  typedef std::map<std::string, std::pair<float, float> > DATA;
  BIN() {};
  BIN(const DATA& d) { data_ = d; };
  BIN(const BIN&  d) { data_ = d.data_; };
  ~BIN() {};

  bool operator <(const BIN& a) const { return data_  < a.data_; };
  bool operator==(const BIN& a) const { return data_ == a.data_; };

  bool empty() const { return data_.empty(); };

  bool contain(const std::map<std::string, float*>& varMap) const
  {
    for (const auto& d : data_) {
      if (varMap.find(d.first)==varMap.end()) return false;
      const auto& val = *varMap.at(d.first);
      if (d.second.first==d.second.second ? val!=d.second.first : (val < d.second.first || val >= d.second.second)) return false;
    }
    return true;
  };

  static std::string str(const DATA& data)
  {
    std::string s;
    for (const auto& d : data) s += Form("%s: [%g, %g], ", d.first.c_str(), d.second.first, d.second.second);
    if (s.rfind(", ")!=std::string::npos) s = s.substr(0, s.rfind(", "));
    return s;
  };

  std::string str() const
  {
    return str(data_);
  };

  static std::string label(const DATA& data)
  {
    std::string s;
    for (const auto& d : data) s += Form("%s_%g_%g_", d.first.c_str(), d.second.first, d.second.second);
    if (s.rfind("_")!=std::string::npos) s = s.substr(0, s.rfind("_"));
    return s;
  };

  std::string label() const
  {
    return label(data_);
  };

 private:
  DATA data_;
};


class PIDSelector {
 public:
  PIDSelector() {};
  ~PIDSelector()
  {
    for (const auto& v : mvaReaderMap_) if (v.second) delete v.second;
  };
  
  // Setters
  void setVariable(const std::string& name, int* var, const bool& isSpectator=false)
  {
    (isSpectator ? specIntMap_ : varIntMap_).emplace(name, var);
  };
  
  void setVariable(const std::string& name, float* var, const bool& isSpectator=false)
  {
    (isSpectator ? specFloatMap_ : varFloatMap_).emplace(name, var);
  };

  void setThreshold(const float& value, const BIN& bin, const std::string& par, const std::string& var)
  {
    cutMap_[bin][par].emplace(var, value);
  };

  void setMVA(const std::string& method, const std::string& file, const BIN& bin)
  {
    if (cutMap_.find(bin)==cutMap_.end()) throw std::logic_error("[ERROR] setMVA: invalid bin "+bin.str());
    if (mvaReaderMap_.find(bin)==mvaReaderMap_.end()) {
      mvaReaderMap_.emplace(bin, new TMVA::Reader("!Color:!Silent"));
      for (const auto& v : varIntMap_) mvaReaderMap_.at(bin)->AddVariable(v.first.c_str(), v.second);
      for (const auto& v : varFloatMap_) mvaReaderMap_.at(bin)->AddVariable(v.first.c_str(), v.second);
      for (const auto& v : specIntMap_) mvaReaderMap_.at(bin)->AddSpectator(v.first.c_str(), v.second);
      for (const auto& v : specFloatMap_) mvaReaderMap_.at(bin)->AddSpectator(v.first.c_str(), v.second);
    }
    const auto& m = *dynamic_cast<TMVA::MethodBase*>(mvaReaderMap_.at(bin)->BookMVA(method.c_str(), file.c_str()));
    for (size_t i=0; i<m.DataInfo().GetNClasses(); i++) {
      const std::string& par = m.DataInfo().GetClassInfo(i)->GetName();
      if (cutMap_.at(bin).find(par)==cutMap_.at(bin).end()) throw std::logic_error("[ERROR] setMVA: invalid particle "+par);
      if (cutMap_.at(bin).at(par).find("MVA")==cutMap_.at(bin).at(par).end()) throw std::logic_error("[ERROR] setMVA: MVA threshold not set");
      mvaClassMap_[bin][method][par] = i;
    } 
  };

  void setFunction(const std::string& funcMean, const std::string& funcError,
		   const BIN& bin, const std::string& par, const std::string& varY,
		   const std::string& varX, const std::pair<float, float>& rngX)
  {
    if (cutMap_.find(bin)==cutMap_.end()) throw std::logic_error("[ERROR] setFunction: invalid bin "+bin.str());
    if (cutMap_.at(bin).find(par)==cutMap_.at(bin).end()) throw std::logic_error("[ERROR] setFunction: invalid particle "+par);
    if (cutMap_.at(bin).at(par).find("FNC")==cutMap_.at(bin).at(par).end()) throw std::logic_error("[ERROR] setFunction: FNC threshold not set");
    const auto& cut = cutMap_.at(bin).at(par).at("FNC");
    const auto& name = bin.label()+"_"+par+"_"+varY+"_"+varX;
    funcMap_[bin][par][varX][varY].emplace("Mean" , TF1((name+"_Mean" ).c_str(), funcMean.c_str() , rngX.first, rngX.second));
    funcMap_[bin][par][varX][varY].emplace("Error", TF1((name+"_Error").c_str(), funcError.c_str(), rngX.first, rngX.second));
    funcMap_[bin][par][varX][varY].emplace("Lower", TF1((name+"_Lower").c_str(), Form("(%s) - %.1f*(%s)", funcMean.c_str(), cut, funcError.c_str()), rngX.first, rngX.second));
    funcMap_[bin][par][varX][varY].emplace("Upper", TF1((name+"_Upper").c_str(), Form("(%s) + %.1f*(%s)", funcMean.c_str(), cut, funcError.c_str()), rngX.first, rngX.second));
  };

  // Getters
  BIN findBin() const
  {
    BIN bin;
    for (const auto& c : cutMap_) {
      if (c.first.contain(varFloatMap_)) {
	bin = c.first;
	break;
      }
    }
    if (bin.empty()) {
      std::string s;
      for (const auto& d : varFloatMap_) s += Form("%s: %g, ", d.first.c_str(), *d.second);
      if (s.rfind(", ")!=std::string::npos) s = s.substr(0, s.rfind(", "));
      throw std::logic_error("[ERROR] findBin: bin not found for "+s);
    }
    return bin;
  };
  
  std::map<std::string, float>
  getValue(const BIN& bin, const std::string& par, const int& dir=0, const std::string& method="BDTG")
  {
    std::map<std::string, float> valM;
    if (!funcMap_[bin][par].empty()) {
      for (auto& x : funcMap_[bin][par]) {
	const auto& xVar = *varFloatMap_[x.first];
	valM[x.first] = 0.0;
	for (auto& v : x.second) {
	  const auto& yVar = *varFloatMap_[v.first];
	  auto chi = (yVar - v.second["Mean"].Eval(xVar))/v.second["Error"].Eval(xVar);
	  if (dir*chi > 0) chi = 0;
	  valM[x.first] += chi*chi;
	}
      }
    }
    else if (!mvaClassMap_[bin].empty()) {
      if (mvaClassMap_.at(bin).find(method)==mvaClassMap_.at(bin).end()) throw std::logic_error("[ERROR] passSelector: invalid MVA method "+method);
      if (mvaClassMap_.at(bin).at(method).find(par)==mvaClassMap_.at(bin).at(method).end()) throw std::logic_error("[ERROR] passSelector: invalid MVA particle "+par);
      valM["MVA"] = mvaReaderMap_.at(bin)->EvaluateMulticlass(mvaClassMap_.at(bin).at(method).at(par), method.c_str());
    }
    else throw std::logic_error("[ERROR] getValue: invalid bin "+bin.str());
    return valM;
  };
  
  bool passSelector(const BIN& bin, const std::string& par, const int& dir=0, const std::string& method="BDTG")
  {
    auto& c = cutMap_[bin][par];
    const auto& val =  getValue(bin, par, dir, method);
    for (const auto& v : val) {
      const auto& cut = (v.first=="MVA" ? c["MVA"] : c["FNC"]);
      if (v.second >= cut) return false;
    }
    return true;
  };

 private:
  std::map<std::string, int*> varIntMap_, specIntMap_;
  std::map<std::string, float*> varFloatMap_, specFloatMap_;
  std::map<BIN, std::map<std::string, std::map<std::string, float> > > cutMap_;
  std::map<BIN, TMVA::Reader*> mvaReaderMap_;
  std::map<BIN, std::map<std::string, std::map<std::string, size_t> > > mvaClassMap_;
  std::map<BIN, std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, TF1> > > > > funcMap_;
};


#endif // ifndef PIDSelector_h
