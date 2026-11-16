#include "gdatastream.h"
#include "node.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

constexpr double kTol = 1e-9;

void check_cond(bool cond, const std::string &msg) {
  if (!cond) {
    throw std::runtime_error("Check failed: " + msg);
  }
}

void check_near(double actual, double expected, double tol,
                const std::string &name) {
  if (!std::isfinite(actual) || std::fabs(actual - expected) > tol) {
    std::ostringstream os;
    os << name << " mismatch: actual=" << actual << ", expected=" << expected
       << ", tol=" << tol;
    throw std::runtime_error(os.str());
  }
}

void write_text_file(const std::filesystem::path &p, const std::string &content) {
  std::ofstream out(p, std::ios::out | std::ios::binary);
  if (!out) {
    throw std::runtime_error("Failed to create file: " + p.string());
  }
  out << content;
}

class TempDir {
public:
  explicit TempDir(const std::string &prefix) {
    const auto stamp =
        std::chrono::steady_clock::now().time_since_epoch().count();
    path_ = std::filesystem::temp_directory_path() /
            (prefix + "_" + std::to_string(stamp));
    std::filesystem::create_directories(path_);
  }
  ~TempDir() {
    std::error_code ec;
    std::filesystem::remove_all(path_, ec);
  }
  const std::filesystem::path &path() const { return path_; }

private:
  std::filesystem::path path_;
};

const char *kGibbsiteLst =
    R"(-j "Gibbsite-dch.json" "Gibbsite-ipm.json" "Gibbsite-dbr-0-0000.json")";

const char *kGibbsiteDchJson = R"JSON([{"_key":"Gibbsite-dch","set":"Gibbsite","dch":{"nIC":7,"nDC":25,"nPH":4,"nPS":3,"nDCs":24,"nICb":7,"nDCb":25,"nPHb":4,"nPSb":3,"nTp":1,"nPp":1,"iGrd":0,"fAalp":1,"mLook":0,"xic":[0,1,2,3,4,5,6],"xdc":[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],"xph":[0,1,2,3],"ICNL":["Al","H","N","Na","Nit","O","Zz"],"ccIC":["e","h","e","e","a","o","z"],"ICmm":[0.0269815406799316,0.00100794994831085,0.0140066995620728,0.0229897994995117,0.0140066995620728,0.015999400138855,0],"DCNL":["Al+3","AlO+","AlO2-","AlO2H@","AlOH+2","Na+","NaOH@","H2@","NO3-","NH3@","NH4+","N2@","O2@","OH-","H+","H2O@","H2","N2","O2","-OH-0.5","-OH-0.5Na+","-OH2+0.5","-OH2+0.5NO3-","Gbs","Gbs"],"ccDC":["S","S","S","S","S","S","S","S","S","S","S","S","S","S","T","W","G","G","G","X","Y","X","Y","Q","O"],"DCmm":[0.0269815406799316,0.0429809408187866,0.0589803409576416,0.0599882909059525,0.0439888907670975,0.0229897994995117,0.0399971495866775,0.0020158998966217,0.0620048999786377,0.0170305494070053,0.0180384993553162,0.0280133991241455,0.03199880027771,0.0170073500871658,0.00100794994831085,0.0180153000354767,0.0020158998966217,0.0280133991241455,0.03199880027771,0.0130075000524521,0.0359972995519638,0.0140154500007629,0.0760203499794006,0.0780035909414291,0.0780035909414291],"PHNL":["aq_gen","gas_gen","Gibbsite_an","Gibbsite"],"ccPH":["a","g","x","s"],"nDCinPH":[16,3,5,1],"A":[1,0,0,0,0,0,3,1,0,0,0,0,1,1,1,0,0,0,0,2,-1,1,1,0,0,0,2,0,1,1,0,0,0,1,2,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,2,0,0,0,0,0,0,0,1,0,0,3,-1,0,3,1,0,0,0,0,0,4,1,0,0,0,1,0,0,0,0,2,0,0,0,0,0,0,0,2,0,0,1,0,0,0,1,-1,0,1,0,0,0,0,1,0,2,0,0,0,1,0,0,2,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,0,0,1,0,0,0,0.75,-0.5,0,1,0,1,0,0.75,0.5,0,2,0,0,0,0.75,0.5,0,2,1,0,0,3.75,-0.5,1,3,0,0,0,3,0,1,3,0,0,0,3,0],"Ttol":1,"TKval":[298.15],"Psat":[3168.73733474928],"Ptol":50000,"Pval":[100000],"denW":[997.061364306278,-0.25866418191325,-0.00946964566451448,0.0450858902637571,0],"denWg":[0.767104546309221,-0,0,0,0],"epsW":[78.2451456494471,-0.355010888388417,0.00135076186526265,0.00406162824448407,0],"epsWg":[0,0,0,0,0],"V0":[-4.52428720605177e-05,3.07336406943463e-07,9.46735262367419e-06,1.30094459034099e-05,-2.72751155835722e-06,-1.20796344838775e-06,3.50925587536806e-06,2.52643576113497e-05,2.86648088757608e-05,2.44457848749428e-05,1.80844669842319e-05,3.34069131302658e-05,3.05008890803284e-05,-4.70784300725686e-06,0,1.80683971876364e-05,0.024789711565,0.024789711565,0.024789711565,1.35512971878052e-05,1.23433339595795e-05,1.35512971878052e-05,4.22161054611206e-05,3.1956000328064e-05,3.1956000328064e-05],"G0":[-483708.00144829,-660420.00050345,-827479.00091632,-864276.99998422,-692595.00090732,-261881.00017392,-418123.99998422,17729.00010995,-110905.00057747,-26669.999973696,-79395.000079017,18194.000182443,16446.000207431,-157270.00090727,0,-237181.38483664,0,0,0,-161017,-429469,-218097,-335573,-1150986,-1150986]}}])JSON";

const char *kGibbsiteIpmJson = R"JSON([{"_key":"Gibbsite-ipm","set":"Gibbsite","ipm":{"ID_key":"AluSiOMSA G  Gibbsite    0    0       1       25      0","pa_PE":1,"PV":1,"PSOL":0,"PAalp":"+","PSigm":"+","Lads":5,"FIa":1,"FIat":4,"pa_DB":1e-17,"pa_DHB":1e-13,"pa_EPS":1e-10,"pa_DK":5e-06,"pa_DS":1e-20,"pa_DF":0.01,"pa_DFM":0.01,"pa_DP":130,"pa_IIM":7000,"pa_PD":2,"pa_PRD":-4,"pa_AG":1,"pa_DGC":-0.001,"pa_PSM":1,"pa_GAR":1,"pa_GAH":1000,"pa_XwMin":1e-13,"pa_ScMin":1e-13,"pa_DcMin":1e-33,"pa_PhMin":1e-20,"pa_ICmin":1e-05,"pa_PC":2,"pa_DFYw":1e-05,"pa_DFYaq":1e-05,"pa_DFYid":1e-05,"pa_DFYr":1e-05,"pa_DFYh":1e-05,"pa_DFYc":1e-05,"pa_DFYs":1e-06,"pa_DW":1,"pa_DT":0,"pa_GAS":0.001,"pa_DG":1000,"pa_DNS":12.05,"pa_IEPS":0.005,"pKin":1,"pa_DKIN":1e-10,"pa_PLLG":30000,"pa_PSTALL":1,"tMin":0,"sMod":["HNNSNNNN","INNINNNN","ANNSNNNN"],"LsMod":[2,0,4,0,0,0,0,0,0],"PMc":[0.0642049834132195,3.72000002861023,1,1,1,0,0,0],"LsMdc":[0,0,0,0,0,0,0,0,0],"LsMdc2":[0,0,0,0,0,0,0,0,0],"LsPhl":[0,0,0,0,0,0,0,0],"SorMc":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"LsISmo":[0,0,0,0,0,0,0,0,0,0,0,8],"XSMD":[1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"LsESmo":[0,0,0,0,0,0,0,0,0,0,4,0],"mCDc":[-0.5,0,0,-0.5,-0.5,1,0,0.5,0.5,0,0,0.5,0.5,-1,0,-0.5,0,0,0,0],"kMod":["NNNNNN","NNNNNN","NNNNNN","NNNNNN"],"LsKin":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"PfFact":[1,1,1,1],"LsUpt":[0,0,0,0,0,0],"B":[0.384597684577799,111.5378443281,0.1001,0.1000001,4,57.0460687272512,0],"Pparc":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"fDQF":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"lnGmf":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6.90775527898214,6.90775527898214,6.90775527898214,6.90775527898214,0,0],"RLC":["B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B"],"RSC":["M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M","M"],"DLL":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"DUL":[1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,0],"Aalp":[0,0,20.2999992370605,0],"Sigw":[0,0,0,0],"Sigg":[0,0,0,0],"YOF":[0,0,0,0],"Nfsp":[8.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"MASDT":[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"C1":[1,0,0,0],"C2":[4,0,0,0],"C3":[0,0,0,0],"CZ":[0.5,0,0,0],"CD":[0.736736747557724,0,0,0],"XfIEC":[0,0,0,0,0,0,0,0],"SATC":[0,-1,0,0,0,-1,0,0,0,-1,0,0,0,-1,0,0,-1,-2,8,-1],"MASDJ":[8.5,0,0,0,8.5,0,0,0,8.5,0,0,0,8.5,0,0,0,0,0,0,0],"SCMC":["3CNININI"],"SACT":["C","C","C","C","N"],"DCads":["X","Y","X","Y","Q"],"dcMod":["HKF","HKF","HKF","HKF","HKF","HKF","HKF","HKF","HKF","HKF","HKF","HKF","HKF","HKF","HKF","HWF","CSN","CSN","CSN","K3C","K3C","K3C","K3C","CSC","CSC"]}}]
)JSON";

const char *kGibbsiteDbrJson = R"JSON([{"_key":"Gibbsite-dbr-0-0000","set":"Gibbsite","dbr":{"NodeHandle":0,"NodeTypeHY":0,"NodeTypeMT":0,"NodeStatusFMT":1,"NodeStatusCH":1,"IterDone":0,"TK":298.15,"P":100000,"Vs":0.0555324116408653,"Vi":0,"Ms":1.09523233353636,"Mi":0,"Hs":0,"Hi":0,"Gs":0,"IS":0.1,"pH":7,"pe":10,"Eh":0.6,"Tm":0,"dt":0,"Dif":0,"Vt":0,"vp":0,"eps":0,"Km":0,"Kf":0,"S":0,"Tr":0,"h":0,"rho":0,"al":0,"at":0,"av":0,"hDl":0,"hDt":0,"hDv":0,"nto":0,"bIC":[0.384597684577799,111.5378443281,0.1001,0.1000001,4,57.0460687272512,0],"rMB":[0,0,0,0,0,0,0],"uIC":[0,0,0,0,0,0,0],"bSP":[0,0,0,0,0,0,0],"xDC":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"gam":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"dll":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"dul":[1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,0],"aPH":[0,0,20299.9992370605,0],"xPH":[0,0,0,0],"vPS":[0,0,0],"mPS":[0,0,0],"xPA":[0,0,0],"amru":[1000000,1000000,1000000],"amrl":[0,0,0],"omPH":[0,0,0,0],"bPS":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]}}])JSON";

void test_table1_scalars_and_indices(const TNode &node) {
  check_cond(node.Get_Lads() == 5, "Get_Lads() == 5");
  check_cond(node.Get_FIa() == 1, "Get_FIa() == 1");
  check_cond(node.Get_FIat() == 4, "Get_FIat() == 4");
  check_cond(node.Get_Lads_ka(0) == 5, "Get_Lads_ka(0) == 5");
  check_cond(node.Get_Lads_ka(1) == -1, "Get_Lads_ka(1) == -1");

  long int lads = -1, fia = -1, fiat = -1, lads_ka = -1;
  node.Get_Lads(lads);
  node.Get_FIa(fia);
  node.Get_FIat(fiat);
  node.Get_Lads_ka(lads_ka, 0);
  check_cond(lads == 5 && fia == 1 && fiat == 4 && lads_ka == 5,
             "Table 1 scalar ref overloads");

  check_cond(node.ka2k(0) == 2, "ka2k(0) == 2");
  check_cond(node.ka2k(1) == -1, "ka2k(1) == -1");
  check_cond(node.k2ka(2) == 0, "k2ka(2) == 0");
  check_cond(node.k2ka(0) == -1, "k2ka(0) == -1");
  check_cond(node.k2ka("Gibbsite_an") == 0, "k2ka(\"Gibbsite_an\") == 0");
  check_cond(node.k2ka("aq_gen") == -1, "k2ka(\"aq_gen\") == -1");
  check_cond(node.k2ka("NonExistent") == -1, "k2ka(\"NonExistent\") == -1");

  long int k_ref = -1, ka_ref = -1;
  node.ka2k(k_ref, 0);
  node.k2ka(ka_ref, 2);
  check_cond(k_ref == 2 && ka_ref == 0, "ka2k/k2ka ref overloads");
  ka_ref = -1;
  node.k2ka(ka_ref, "Gibbsite_an");
  check_cond(ka_ref == 0, "k2ka(name) ref overload");

  check_cond(node.j_ka2ja(19, 0) == 0, "j_ka2ja(19, 0) == 0");
  check_cond(node.j_ka2ja(23, 0) == 4, "j_ka2ja(23, 0) == 4");
  check_cond(node.j_ka2ja(0, 0) == -1, "j_ka2ja(0, 0) == -1");
  check_cond(node.j_ka2ja("-OH-0.5Na+", 0) == 1,
             "j_ka2ja(\"-OH-0.5Na+\", 0) == 1");
  check_cond(node.j_ka2ja("NoSuchSpecies", 0) == -1,
             "j_ka2ja(\"NoSuchSpecies\", 0) == -1");

  long int ja_ref = -1, j_ref = -1;
  node.j_ka2ja(ja_ref, 20, 0);
  check_cond(ja_ref == 1, "j_ka2ja(j, ka) ref overload");
  ja_ref = -1;
  node.j_ka2ja(ja_ref, "-OH2+0.5", 0);
  check_cond(ja_ref == 2, "j_ka2ja(name, ka) ref overload");

  check_cond(node.ja2j(0) == 19, "ja2j(0) == 19");
  check_cond(node.ja2j(4) == 23, "ja2j(4) == 23");
  check_cond(node.ja2j(5) == -1, "ja2j(5) == -1");
  node.ja2j(j_ref, 3);
  check_cond(j_ref == 22, "ja2j ref overload");

  check_cond(node.Get_ist(0) == 0, "Get_ist(0) == 0");
  check_cond(node.Get_nst(0) == 1, "Get_nst(0) == 1");
  check_cond(node.Get_iss(0) == 0, "Get_iss(0) == 0");
  check_cond(node.Get_nst(4) == 0, "Get_nst(4) == 0 for carrier Q");
  check_cond(node.ja2ist(0) == 0, "ja2ist(0) == 0");
  check_cond(node.Get_ist(-1) == -1 && node.Get_nst(-1) == -1 &&
                 node.Get_iss(-1) == -1 && node.ja2ist(-1) == -1,
             "Invalid ja index returns -1");

  long int ist_r = -1, nst_r = -1, iss_r = -1, ja2ist_r = -1;
  node.Get_ist(ist_r, 0);
  node.Get_nst(nst_r, 0);
  node.Get_iss(iss_r, 0);
  node.ja2ist(ja2ist_r, 0);
  check_cond(ist_r == 0 && nst_r == 1 && iss_r == 0 && ja2ist_r == 0,
             "ist/nst/iss/ja2ist ref overloads");
}

void test_table1_getters_and_setters(TNode &node) {
  check_near(node.Get_SSA(2), 20.2999992370605, 1e-6, "Get_SSA(2)");
  check_near(node.Get_SSA_ka(0), 20.2999992370605, 1e-6, "Get_SSA_ka(0)");
  double ssa_r = 0.0, ssa_ka_r = 0.0;
  node.Get_SSA(ssa_r, 2);
  node.Get_SSA_ka(ssa_ka_r, 0);
  check_near(ssa_r, 20.2999992370605, 1e-6, "Get_SSA ref");
  check_near(ssa_ka_r, 20.2999992370605, 1e-6, "Get_SSA_ka ref");

  check_near(node.Get_zjA(0), -0.5, kTol, "Get_zjA(0)");
  check_near(node.Get_zjB(0), 0.0, kTol, "Get_zjB(0)");
  check_near(node.Get_zjC(0), 0.0, kTol, "Get_zjC(0)");
  check_near(node.Get_bzj(0), -0.5, kTol, "Get_bzj(0)");
  check_near(node.Get_zjA(1), -0.5, kTol, "Get_zjA(1)");
  check_near(node.Get_zjB(1), 1.0, kTol, "Get_zjB(1)");
  check_near(node.Get_bzj(1), 0.5, kTol, "Get_bzj(1)");
  check_near(node.Get_zjA(3), 0.5, kTol, "Get_zjA(3)");
  check_near(node.Get_zjB(3), -1.0, kTol, "Get_zjB(3)");
  check_near(node.Get_bzj(3), -0.5, kTol, "Get_bzj(3)");

  check_near(node.Get_GamCs(0, 0), 8.5, kTol, "Get_GamCs(0, 0)");
  check_near(node.Get_CapCs(0, 0), 0.0, kTol, "Get_CapCs(0, 0)");
  check_cond(node.Get_iss2ist(0, 0) == 0, "Get_iss2ist(0, 0) == 0");
  check_near(node.Get_fAsurt(0, 0), 1.0, kTol, "Get_fAsurt(0, 0)");
  check_near(node.Get_A1Tour(0, 0), 0.0, kTol, "Get_A1Tour(0, 0)");
  check_near(node.Get_A2Tour(0, 0), 0.0, kTol, "Get_A2Tour(0, 0)");
  check_near(node.Get_A3Tour(0, 0), 0.0, kTol, "Get_A3Tour(0, 0)");

  check_near(node.Get_C1(0, 0), 1.0, kTol, "Get_C1(0, 0)");
  check_near(node.Get_C2(0, 0), 4.0, kTol, "Get_C2(0, 0)");
  check_near(node.Get_C3(0, 0), 0.0, kTol, "Get_C3(0, 0)");
  check_near(node.Get_fzpSpd(0, 0), 0.5, kTol, "Get_fzpSpd(0, 0)");
  check_near(node.Get_pEtaD(0, 0), 0.0, kTol, "Get_pEtaD(0, 0)");
  check_near(node.Get_pEtaQ(0, 0), 0.0, kTol, "Get_pEtaQ(0, 0)");

  check_near(node.Get_GamC(0), 8.5, kTol, "Get_GamC(0)");
  check_near(node.Get_isoP1(0), 0.0, kTol, "Get_isoP1(0)");
  check_near(node.Get_isoP2(0), 0.0, kTol, "Get_isoP2(0)");
  check_near(node.Get_isoP3(0), 0.0, kTol, "Get_isoP3(0)");
  check_cond(node.Get_dent(0, 0) == 1, "Get_dent(0, 0) == 1");
  check_cond(node.Get_EIL_code(0, 0) == '3', "Get_EIL_code(0, 0) == '3'");
  check_cond(node.Get_SACT_code(0, 0) == 'C', "Get_SACT_code(0, 0) == 'C'");
  check_cond(node.Get_SATT(0) == 'C', "Get_SATT(0) == 'C'");

  // Exercise Set_* and Get_* round-trips, then restore original values
  const double orig_ssa = node.Get_SSA_ka(0);
  node.Set_SSA_ka(42.5, 0);
  check_near(node.Get_SSA(2), 42.5, kTol, "Set_SSA_ka -> Get_SSA");
  check_near(node.Get_SSA_ka(0), 42.5, kTol, "Set_SSA_ka -> Get_SSA_ka");
  node.Set_SSA_ka(orig_ssa, 0);

  node.Set_zjA(-0.25, 0);
  node.Set_zjB(0.75, 0);
  node.Set_zjC(-0.125, 0);
  double za = 0, zb = 0, zc = 0, bzj = 0;
  node.Get_zjA(za, 0);
  node.Get_zjB(zb, 0);
  node.Get_zjC(zc, 0);
  node.Get_bzj(bzj, 0);
  check_near(za, -0.25, kTol, "Set_zjA");
  check_near(zb, 0.75, kTol, "Set_zjB");
  check_near(zc, -0.125, kTol, "Set_zjC");
  check_near(bzj, 0.375, kTol, "Get_bzj updated");
  node.Set_zjA(-0.5, 0);
  node.Set_zjB(0.0, 0);
  node.Set_zjC(0.0, 0);
  check_near(node.Get_bzj(0), -0.5, kTol, "Restored Get_bzj(0)");

  double tmp_d = 0.0;
  long int tmp_i = -1;
  char tmp_c = '\0';

  node.Set_GamCs(11.2, 0, 1);
  node.Get_GamCs(tmp_d, 0, 1);
  check_near(tmp_d, 11.2, kTol, "Set/Get_GamCs");
  node.Set_GamCs(0.0, 0, 1);

  node.Set_CapCs(0.15, 0, 1);
  node.Get_CapCs(tmp_d, 0, 1);
  check_near(tmp_d, 0.15, kTol, "Set/Get_CapCs");
  node.Set_CapCs(0.0, 0, 1);

  node.Set_iss2ist(1, 0, 1);
  node.Get_iss2ist(tmp_i, 0, 1);
  check_cond(tmp_i == 1, "Set/Get_iss2ist");

  node.Set_fAsurt(0.3, 0, 1);
  node.Get_fAsurt(tmp_d, 0, 1);
  check_near(tmp_d, 0.3, kTol, "Set/Get_fAsurt");
  node.Set_fAsurt(0.0, 0, 1);

  node.Set_A1Tour(1.1, 0, 1);
  node.Set_A2Tour(2.2, 0, 1);
  node.Set_A3Tour(3.3, 0, 1);
  node.Get_A1Tour(tmp_d, 0, 1);
  check_near(tmp_d, 1.1, kTol, "Set/Get_A1Tour");
  node.Get_A2Tour(tmp_d, 0, 1);
  check_near(tmp_d, 2.2, kTol, "Set/Get_A2Tour");
  node.Get_A3Tour(tmp_d, 0, 1);
  check_near(tmp_d, 3.3, kTol, "Set/Get_A3Tour");
  node.Set_A1Tour(0.0, 0, 1);
  node.Set_A2Tour(0.0, 0, 1);
  node.Set_A3Tour(0.0, 0, 1);

  node.Set_C1(2.1, 0, 1);
  node.Set_C2(3.2, 0, 1);
  node.Set_C3(0.4, 0, 1);
  node.Set_fzpSpd(0.8, 0, 1);
  node.Set_pEtaD(-0.15, 0, 1);
  node.Set_pEtaQ(-0.25, 0, 1);
  node.Get_C1(tmp_d, 0, 1);
  check_near(tmp_d, 2.1, kTol, "Set/Get_C1");
  node.Get_C2(tmp_d, 0, 1);
  check_near(tmp_d, 3.2, kTol, "Set/Get_C2");
  node.Get_C3(tmp_d, 0, 1);
  check_near(tmp_d, 0.4, kTol, "Set/Get_C3");
  node.Get_fzpSpd(tmp_d, 0, 1);
  check_near(tmp_d, 0.8, kTol, "Set/Get_fzpSpd");
  node.Get_pEtaD(tmp_d, 0, 1);
  check_near(tmp_d, -0.15, kTol, "Set/Get_pEtaD");
  node.Get_pEtaQ(tmp_d, 0, 1);
  check_near(tmp_d, -0.25, kTol, "Set/Get_pEtaQ");
  node.Set_C1(0.0, 0, 1);
  node.Set_C2(0.0, 0, 1);
  node.Set_C3(0.0, 0, 1);
  node.Set_fzpSpd(0.0, 0, 1);
  node.Set_pEtaD(0.0, 0, 1);
  node.Set_pEtaQ(0.0, 0, 1);

  node.Set_GamC(9.9, 0);
  node.Set_isoP1(0.11, 0);
  node.Set_isoP2(0.22, 0);
  node.Set_isoP3(0.33, 0);
  node.Get_GamC(tmp_d, 0);
  check_near(tmp_d, 9.9, kTol, "Set/Get_GamC");
  node.Get_isoP1(tmp_d, 0);
  check_near(tmp_d, 0.11, kTol, "Set/Get_isoP1");
  node.Get_isoP2(tmp_d, 0);
  check_near(tmp_d, 0.22, kTol, "Set/Get_isoP2");
  node.Get_isoP3(tmp_d, 0);
  check_near(tmp_d, 0.33, kTol, "Set/Get_isoP3");
  node.Set_GamC(8.5, 0);
  node.Set_isoP1(0.0, 0);
  node.Set_isoP2(0.0, 0);
  node.Set_isoP3(0.0, 0);

  // Denticity: add second bond to site 1 (mapped to surface type 1) -> nst=2,
  // ja2ist=-1; then restore
  node.Set_dent(2L, 0, 1);
  node.Get_dent(tmp_i, 0, 1);
  check_cond(tmp_i == 2, "Set/Get_dent long overload");
  node.Get_dent(tmp_d, 0, 1);
  check_near(tmp_d, 2.0, kTol, "Get_dent double overload");
  check_cond(node.Get_nst(0) == 2, "Get_nst(0) == 2 when bound to 2 sites");
  check_cond(node.ja2ist(0) == -1,
             "ja2ist(0) == -1 across different surface types");
  node.Set_dent(0L, 0, 1);
  node.Set_iss2ist(0, 0, 1);
  check_cond(node.Get_nst(0) == 1 && node.ja2ist(0) == 0,
             "Restored Get_nst(0) == 1 & ja2ist(0) == 0");

  node.Set_EIL_code('T', 0, 1);
  node.Set_SACT_code('L', 0, 1);
  node.Set_SATT('L', 0);
  node.Get_EIL_code(tmp_c, 0, 1);
  check_cond(tmp_c == 'T', "Set/Get_EIL_code");
  node.Get_SACT_code(tmp_c, 0, 1);
  check_cond(tmp_c == 'L', "Set/Get_SACT_code");
  node.Get_SATT(tmp_c, 0);
  check_cond(tmp_c == 'L', "Set/Get_SATT");
  node.Set_EIL_code('N', 0, 1);
  node.Set_SACT_code('N', 0, 1);
  node.Set_SATT('C', 0);
}

void test_table2_outputs(const TNode &node) {
  check_near(node.Get_Zeff(), 1.0, 1e-4, "Get_Zeff()");
  double zeff_r = 0.0;
  node.Get_Zeff(zeff_r);
  check_near(zeff_r, node.Get_Zeff(), kTol, "Get_Zeff ref");

  check_cond(node.Get_C_dl(0, 0) > 0.0, "Get_C_dl(0, 0) > 0");
  check_cond(std::isfinite(node.Get_h_A(0, 0)) && node.Get_h_A(0, 0) >= 0.0,
             "Get_h_A(0, 0) valid");
  check_cond(std::isfinite(node.Get_h_B(0, 0)) && node.Get_h_B(0, 0) >= 0.0,
             "Get_h_B(0, 0) valid");
  check_cond(std::isfinite(node.Get_h_C(0, 0)) && node.Get_h_C(0, 0) >= 0.0,
             "Get_h_C(0, 0) valid");
  check_cond(node.Get_h_dl(0, 0) > 0.0, "Get_h_dl(0, 0) > 0");
  check_cond(node.Get_h_spz(0, 0) > 0.0, "Get_h_spz(0, 0) > 0");

  double v = 0.0;
  node.Get_C_dl(v, 0, 0);
  check_near(v, node.Get_C_dl(0, 0), kTol, "Get_C_dl ref");
  node.Get_h_A(v, 0, 0);
  check_near(v, node.Get_h_A(0, 0), kTol, "Get_h_A ref");
  node.Get_h_B(v, 0, 0);
  check_near(v, node.Get_h_B(0, 0), kTol, "Get_h_B ref");
  node.Get_h_C(v, 0, 0);
  check_near(v, node.Get_h_C(0, 0), kTol, "Get_h_C ref");
  node.Get_h_dl(v, 0, 0);
  check_near(v, node.Get_h_dl(0, 0), kTol, "Get_h_dl ref");
  node.Get_h_spz(v, 0, 0);
  check_near(v, node.Get_h_spz(0, 0), kTol, "Get_h_spz ref");

  check_cond(node.Get_Mass_s(2) > 0.0, "Get_Mass_s(2) > 0");
  check_cond(node.Get_Mols_s(2) > 0.0, "Get_Mols_s(2) > 0");
  node.Get_Mass_s(v, 2);
  check_near(v, node.Get_Mass_s(2), kTol, "Get_Mass_s ref");
  node.Get_Mols_s(v, 2);
  check_near(v, node.Get_Mols_s(2), kTol, "Get_Mols_s ref");

  const double etaA = node.Get_eta_A(0, 0);
  const double etaB = node.Get_eta_B(0, 0);
  const double etaC = node.Get_eta_C(0, 0);
  const double etaD = node.Get_eta_D(0, 0);
  check_near(etaA + etaB + etaC + etaD, 0.0, 1e-12,
             "EIL plane electroneutrality");

  node.Get_eta_A(v, 0, 0);
  check_near(v, etaA, kTol, "Get_eta_A ref");
  node.Get_eta_B(v, 0, 0);
  check_near(v, etaB, kTol, "Get_eta_B ref");
  node.Get_eta_C(v, 0, 0);
  check_near(v, etaC, kTol, "Get_eta_C ref");
  node.Get_eta_D(v, 0, 0);
  check_near(v, etaD, kTol, "Get_eta_D ref");

  const double psiA = node.Get_psi_A(0, 0);
  const double psiB = node.Get_psi_B(0, 0);
  const double psiC = node.Get_psi_C(0, 0);
  const double psiD = node.Get_psi_D(0, 0);
  const double psiZ = node.Get_psi_Z(0, 0);
  check_cond(std::isfinite(psiA) && std::fabs(psiA) > 1e-6,
             "psi_A finite and non-zero");
  check_cond(std::isfinite(psiB) && std::fabs(psiB) > 1e-6,
             "psi_B finite and non-zero");
  check_cond(std::isfinite(psiC), "psi_C finite");
  check_cond(std::isfinite(psiD) && std::fabs(psiD) > 1e-6,
             "psi_D finite and non-zero");
  check_cond(std::isfinite(psiZ) && std::fabs(psiZ) > 1e-6,
             "psi_Z finite and non-zero");
  check_cond(std::fabs(psiZ) <= std::fabs(psiD) + 1e-12 &&
                 (psiZ * psiD >= 0.0),
             "psi_Z lies between 0 and psi_D");

  node.Get_psi_A(v, 0, 0);
  check_near(v, psiA, kTol, "Get_psi_A ref");
  node.Get_psi_B(v, 0, 0);
  check_near(v, psiB, kTol, "Get_psi_B ref");
  node.Get_psi_C(v, 0, 0);
  check_near(v, psiC, kTol, "Get_psi_C ref");
  node.Get_psi_D(v, 0, 0);
  check_near(v, psiD, kTol, "Get_psi_D ref");
  node.Get_psi_Z(v, 0, 0);
  check_near(v, psiZ, kTol, "Get_psi_Z ref");

  // Activity terms for surface species ja = 0 (j = 19)
  const double lnGam_ja = node.Get_lnGam_ja(0);
  check_near(node.Get_lnGam(0), lnGam_ja, kTol, "Get_lnGam alias");
  check_near(node.Get_lnSACT_ja(0), node.Get_lnSACT(19), kTol,
             "Get_lnSACT_ja vs Get_lnSACT");
  check_near(node.Get_lnExpF_ja(0), node.Get_lnExpF(19), kTol,
             "Get_lnExpF_ja vs Get_lnExpF");
  check_near(node.Get_lnCoul_ja(0), node.Get_lnCoul(19), kTol,
             "Get_lnCoul_ja vs Get_lnCoul");
  check_cond(std::fabs(node.Get_lnCoul(19)) > 1e-6,
             "Coulomb term for charged surface species -OH-0.5 is non-zero");

  node.Get_lnGam_ja(v, 0);
  check_near(v, lnGam_ja, kTol, "Get_lnGam_ja ref");
  node.Get_lnGam(v, 0);
  check_near(v, lnGam_ja, kTol, "Get_lnGam ref");
  node.Get_lnSACT_ja(v, 0);
  check_near(v, node.Get_lnSACT(19), kTol, "Get_lnSACT_ja ref");
  node.Get_lnExpF_ja(v, 0);
  check_near(v, node.Get_lnExpF(19), kTol, "Get_lnExpF_ja ref");
  node.Get_lnCoul_ja(v, 0);
  check_near(v, node.Get_lnCoul(19), kTol, "Get_lnCoul_ja ref");
  node.Get_lnSACT(v, 19);
  check_near(v, node.Get_lnSACT(19), kTol, "Get_lnSACT ref");
  node.Get_lnExpF(v, 19);
  check_near(v, node.Get_lnExpF(19), kTol, "Get_lnExpF ref");
  node.Get_lnCoul(v, 19);
  check_near(v, node.Get_lnCoul(19), kTol, "Get_lnCoul ref");

  check_cond(node.Get_Xfts(0, 0) > 0.0, "Get_Xfts(0, 0) > 0");
  check_cond(node.Get_Dfts(0, 0) > 0.0, "Get_Dfts(0, 0) > 0");
  node.Get_Xfts(v, 0, 0);
  check_near(v, node.Get_Xfts(0, 0), kTol, "Get_Xfts ref");
  node.Get_Dfts(v, 0, 0);
  check_near(v, node.Get_Dfts(0, 0), kTol, "Get_Dfts ref");
}

class TestNode : public TNode {
public:
  using TNode::TNode;
  explicit TestNode(const TNode &other) : TNode(other) {}
  using TNode::multi_ptr;
};

void test_copy_and_binary_roundtrip(TestNode &node,
                                    const std::filesystem::path &tmp_dir) {
  TestNode copyNode(static_cast<const TNode &>(node));
  check_cond(copyNode.Get_Lads() == node.Get_Lads(), "copy Lads");
  check_cond(copyNode.Get_FIa() == node.Get_FIa(), "copy FIa");
  check_cond(copyNode.Get_FIat() == node.Get_FIat(), "copy FIat");
  check_near(copyNode.Get_h_spz(0, 0), node.Get_h_spz(0, 0), kTol,
             "copy Get_h_spz");
  check_near(copyNode.Get_zjB(1), node.Get_zjB(1), kTol, "copy Get_zjB");
  check_cond(copyNode.Get_dent(1, 0) == node.Get_dent(1, 0), "copy Get_dent");
  check_near(copyNode.Get_Xfts(0, 0), node.Get_Xfts(0, 0), kTol,
             "copy Get_Xfts");
  check_near(copyNode.Get_psi_A(0, 0), node.Get_psi_A(0, 0), kTol,
             "copy Get_psi_A");
  check_near(copyNode.Get_psi_Z(0, 0), node.Get_psi_Z(0, 0), kTol,
             "copy Get_psi_Z");

  // Binary MULTI serialization round-trip (verifies XdlZ, XSMD, mCDc, XFTS)
  const auto bin_path = tmp_dir / "multi_roundtrip.bin";
  {
    GemDataStream out_stream(bin_path.string(),
                             std::ios::out | std::ios::binary);
    node.multi_ptr()->to_file(out_stream);
    out_stream.close();
  }
  TestNode binNode;
  {
    GemDataStream in_stream(bin_path.string(), std::ios::in | std::ios::binary);
    binNode.multi_ptr()->from_file(in_stream);
    in_stream.close();
  }
  check_near(binNode.Get_h_spz(0, 0), node.Get_h_spz(0, 0), kTol,
             "binary round-trip Get_h_spz");
  check_near(binNode.Get_zjB(1), node.Get_zjB(1), kTol,
             "binary round-trip Get_zjB");
  check_cond(binNode.Get_dent(1, 0) == node.Get_dent(1, 0),
             "binary round-trip Get_dent");
  check_near(binNode.Get_Xfts(0, 0), node.Get_Xfts(0, 0), kTol,
             "binary round-trip Get_Xfts");
  check_near(binNode.Get_psi_A(0, 0), node.Get_psi_A(0, 0), kTol,
             "binary round-trip Get_psi_A");

  // Also ensure copyNode and binNode destroy cleanly when leaving scope
}

} // namespace

int main() {
  try {
    TempDir tmp("gems3k_sorption_api_test");
    const auto lst_path = tmp.path() / "Gibbsite-dat.lst";
    write_text_file(lst_path, kGibbsiteLst);
    write_text_file(tmp.path() / "Gibbsite-dch.json", kGibbsiteDchJson);
    write_text_file(tmp.path() / "Gibbsite-ipm.json", kGibbsiteIpmJson);
    write_text_file(tmp.path() / "Gibbsite-dbr-0-0000.json", kGibbsiteDbrJson);

    TestNode node;
    const long int init_rc = node.GEM_init(lst_path.string().c_str());
    check_cond(init_rc == 0, "GEM_init returned 0");

    test_table1_scalars_and_indices(node);
    test_table1_getters_and_setters(node);

    const long int run_rc = node.GEM_run(false);
    check_cond(run_rc == OK_GEM_AIA || run_rc == OK_GEM_SIA,
               "First GEM_run(false) converged");
    test_table2_outputs(node);

    const double psiA_before = node.Get_psi_A(0, 0);
    const double psiZ_before = node.Get_psi_Z(0, 0);
    const double hspz_before = node.Get_h_spz(0, 0);

    test_copy_and_binary_roundtrip(node, tmp.path());

    // Perturb C1 and fzpSpd, re-run GEM_run(false), and check physical response
    node.Set_C1(2.0, 0, 0);
    node.Set_fzpSpd(0.25, 0, 0);
    node.pCNode()->NodeStatusCH = NEED_GEM_AIA;
    const long int run_rc2 = node.GEM_run(false);
    check_cond(run_rc2 == OK_GEM_AIA || run_rc2 == OK_GEM_SIA,
               "Second GEM_run(false) after Set_C1/Set_fzpSpd converged");

    check_near(node.Get_C1(0, 0), 2.0, kTol, "Persisted Set_C1");
    check_near(node.Get_fzpSpd(0, 0), 0.25, kTol, "Persisted Set_fzpSpd");
    check_near(node.Get_h_spz(0, 0), 0.5 * hspz_before, 1e-4,
               "Get_h_spz halved when fzpSpd halved at ~same ionic strength");
    check_cond(std::fabs(node.Get_psi_A(0, 0) - psiA_before) > 1e-4,
               "Get_psi_A changed after Set_C1");
    check_cond(std::fabs(node.Get_psi_Z(0, 0) - psiZ_before) > 1e-6 &&
                   std::fabs(node.Get_psi_Z(0, 0)) >
                       std::fabs(psiZ_before) * 1.1,
               "Get_psi_Z magnitude increased after halving fzpSpd");

    std::cout << "sorption_api_test passed all checks! (psiA: " << psiA_before
              << " -> " << node.Get_psi_A(0, 0) << ", psiZ: " << psiZ_before
              << " -> " << node.Get_psi_Z(0, 0) << ")\n";
  } catch (const std::exception &ex) {
    std::cerr << "sorption_api_test FAILED: " << ex.what() << '\n';
    return 1;
  }
  return 0;
}
