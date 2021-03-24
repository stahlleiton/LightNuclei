#ifndef bin_h
#define bin_h

#include "TGraphAsymmErrors.h"
#include "TString.h"
#include <iostream>
#include <map>
#include <string>
#include <initializer_list>


typedef std::string BINF_key_t;
typedef std::pair<float, float> BINF_value_t;


double sumError(const double& a, const double& b)
{
  return std::sqrt(a*a + b*b);
};


class BINF : public std::map<BINF_key_t, BINF_value_t>
{
 public:
  typedef std::map<BINF_key_t, BINF_value_t> DATA;
  BINF() : DATA() {};
  BINF(const std::initializer_list<std::pair<const BINF_key_t, BINF_value_t> >& l) : DATA(l) {};
  BINF(const std::vector<BINF>& bS) { this->add(bS); };
  BINF(const std::string& n, const float& a, const float& b) { this->add(n, a, b); };
  BINF(const TGraphAsymmErrors& g, const int& i) { this->add(g, i); };
  
  void add(const std::string& n, const float& a, const float& b, const bool& replace=true)
  {
    if (!replace && this->find(n)!=this->end()) return;
    (*this)[n] = {a, b};
  };
  
  void add(const TGraphAsymmErrors& g, const int& i, const bool& replace=true)
  {
    this->add(g.GetName(), g.GetX()[i]-g.GetEXlow()[i], g.GetX()[i]+g.GetEXhigh()[i], replace);
  };
  
  void add(const BINF& b, const bool& replace=true)
  {
    for (const auto& d : b) {
      if (d.first=="") continue;
      if (!replace && this->find(d.first)!=this->end()) continue;
      (*this)[d.first] = d.second;
    }
  };

  void add(const std::vector<BINF>& bS, const bool& replace=true)
  {
    for (const auto& b : bS) this->add(b, replace);
  };

  const std::pair<float, float>& at(const std::string& s) const
  {
    const auto& p = this->find(s);
    if (p==this->end()) throw std::runtime_error("[ERROR] Invalid bin variable "+s+"!");
    return p->second;
  };

  bool has(const std::string& s) const
  {
    return this->find(s)!=this->end();
  };

  bool contain(const std::map<std::string, float*>& vM) const
  {
    for (const auto& d : *this) {
      const auto& p = vM.find(d.first);
      if (p==vM.end() || !p->second) return false;
      if (d.second.first==d.second.second ? *p->second!=d.second.first : (*p->second < d.second.first || *p->second >= d.second.second)) return false;
    }
    return true;
  };

  bool contain(const BINF& b) const
  {
    for (const auto& d : *this) {
      const auto& p = b.find(d.first);
      if (p==b.end()) return false;
      if (p->second.first < d.second.first || p->second.second > d.second.second) return false;
    }
    return true;
  };

  std::string str(const int& t=0) const
  {
    const auto c = (t==1 ? 100. : 1.);
    std::pair<std::string, std::string> f = std::make_pair("%s: [%g, %g]", ", ");
    if (t==1) f = std::make_pair("%s_%.0f_%.0f", "_");
    std::string s;
    for (const auto& d : *this) s += Form((f.first+f.second).c_str(), d.first.c_str(), d.second.first*c, d.second.second*c);
    if (s.rfind(f.second)!=std::string::npos) s = s.substr(0, s.rfind(f.second));
    return s;
  };
  
  std::string dsCut() const
  {
    const std::map<std::string, std::string> VARM = { {"absRap", "abs(rap)"}, {"absEta", "abs(eta)"}};
    std::string s;
    for (const auto& d : *this) {
      const auto v = VARM.find(d.first);
      const auto var = (v!=VARM.end() ? v->second : d.first);
      s += Form("(%g <= %s && %s < %g) && ", d.second.first, var.c_str(), var.c_str(), d.second.second);
    }
    if (s.rfind(" && ")!=std::string::npos) s = s.substr(0, s.rfind(" && "));
    return s;
  };

  void print() const
  {
    std::cout << this->str() << std::endl;
  };
};


typedef std::map< BINF , std::map< BINF , std::map< BINF , std::array<double, 5> > > > CONT;
typedef const std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, CONT> > > > DATA_CONT;
typedef std::map< BINF , std::map< BINF , std::map< std::string , std::vector<float> > > > CONTBIN;
typedef const std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, CONTBIN> > > > DATA_CONTBIN;


#endif // #ifndef bin_h
