#include "gdatastream.h"
#include "node.h"

#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
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

void test_gibbsite_an_indices(const TNode &node, long int &ka_out,
                              long int &k_out) {
  check_cond(node.Get_Lads() == 5, "Get_Lads() == 5");
  check_cond(node.Get_FIa() == 1, "Get_FIa() == 1");
  check_cond(node.Get_FIat() == 4, "Get_FIat() == 4");

  long int lads = -1, fia = -1, fiat = -1;
  node.Get_Lads(lads);
  node.Get_FIa(fia);
  node.Get_FIat(fiat);
  check_cond(lads == 5 && fia == 1 && fiat == 4,
             "Table 1 scalar reference overloads");

  const long int ka = node.k2ka("Gibbsite_an");
  check_cond(ka == 0, "k2ka(\"Gibbsite_an\") == 0");
  const long int k = node.ka2k(ka);
  check_cond(k == 2, "ka2k(0) == 2");
  check_cond(node.k2ka(k) == ka, "k2ka(2) == 0");
  check_cond(node.k2ka(0) == -1, "k2ka(aq_gen) == -1");
  check_cond(node.k2ka("aq_gen") == -1, "k2ka(\"aq_gen\") == -1");
  check_cond(node.k2ka("NonExistent") == -1, "k2ka(\"NonExistent\") == -1");
  check_cond(node.ka2k(1) == -1, "ka2k(1) == -1");

  long int k_ref = -1, ka_ref = -1;
  node.ka2k(k_ref, ka);
  node.k2ka(ka_ref, k);
  check_cond(k_ref == k && ka_ref == ka, "ka2k/k2ka ref overloads");
  ka_ref = -1;
  node.k2ka(ka_ref, "Gibbsite_an");
  check_cond(ka_ref == ka, "k2ka(\"Gibbsite_an\") ref overload");

  check_cond(node.Get_Lads_ka(ka) == 5, "Get_Lads_ka(Gibbsite_an) == 5");
  check_cond(node.Get_Lads_ka(1) == -1, "Get_Lads_ka(1) == -1");
  long int lads_ka = -1;
  node.Get_Lads_ka(lads_ka, ka);
  check_cond(lads_ka == 5, "Get_Lads_ka ref overload");

  const std::array<const char *, 5> expected_species = {
      "-OH-0.5", "-OH-0.5Na+", "-OH2+0.5", "-OH2+0.5NO3-", "Gbs"};
  for (long int ja = 0; ja < 5; ++ja) {
    const long int j = 19 + ja;
    check_cond(node.j_ka2ja(j, ka) == ja,
               "j_ka2ja(" + std::to_string(j) + ", ka)");
    check_cond(node.j_ka2ja(expected_species[ja], ka) == ja,
               std::string("j_ka2ja(\"") + expected_species[ja] + "\", ka)");
    check_cond(node.ja2j(ja) == j, "ja2j(" + std::to_string(ja) + ")");
    check_cond(node.Get_ist(ja) == 0, "Get_ist(" + std::to_string(ja) + ")");
    check_cond(node.Get_iss(ja) == 0, "Get_iss(" + std::to_string(ja) + ")");
    check_cond(node.ja2ist(ja) == 0, "ja2ist(" + std::to_string(ja) + ")");
    check_cond(node.Get_nst(ja) == (ja < 4 ? 1 : 0),
               "Get_nst(" + std::to_string(ja) + ")");

    long int ja_r = -1, j_r = -1, ist_r = -1, nst_r = -1, iss_r = -1,
             ja2ist_r = -1;
    node.j_ka2ja(ja_r, j, ka);
    check_cond(ja_r == ja, "j_ka2ja(j, ka) ref overload");
    ja_r = -1;
    node.j_ka2ja(ja_r, expected_species[ja], ka);
    check_cond(ja_r == ja, "j_ka2ja(name, ka) ref overload");
    node.ja2j(j_r, ja);
    node.Get_ist(ist_r, ja);
    node.Get_nst(nst_r, ja);
    node.Get_iss(iss_r, ja);
    node.ja2ist(ja2ist_r, ja);
    check_cond(j_r == j && ist_r == 0 && iss_r == 0 && ja2ist_r == 0 &&
                   nst_r == (ja < 4 ? 1 : 0),
               "species index ref overloads");
  }

  check_cond(node.j_ka2ja(0, ka) == -1, "aq species j=0 not in Gibbsite_an");
  check_cond(node.j_ka2ja("NoSuchSpecies", ka) == -1,
             "j_ka2ja(\"NoSuchSpecies\", ka) == -1");
  check_cond(node.ja2j(5) == -1, "ja2j(5) == -1");
  check_cond(node.Get_ist(-1) == -1 && node.Get_nst(-1) == -1 &&
                 node.Get_iss(-1) == -1 && node.ja2ist(-1) == -1,
             "Invalid ja index returns -1");

  ka_out = ka;
  k_out = k;
}

void test_gibbsite_an_inputs(TNode &node, long int ka, long int k) {
  check_near(node.Get_SSA(k), 20.2999992370605, 1e-6, "Get_SSA(Gibbsite_an)");
  check_near(node.Get_SSA_ka(ka), 20.2999992370605, 1e-6,
             "Get_SSA_ka(Gibbsite_an)");
  double ssa_r = 0.0, ssa_ka_r = 0.0;
  node.Get_SSA(ssa_r, k);
  node.Get_SSA_ka(ssa_ka_r, ka);
  check_near(ssa_r, 20.2999992370605, 1e-6, "Get_SSA ref");
  check_near(ssa_ka_r, 20.2999992370605, 1e-6, "Get_SSA_ka ref");

  // Plane charge distributions mCDc[ja][0..3] across all 5 Gibbsite_an species
  const std::array<std::array<double, 4>, 5> expected_z = {{
      {{-0.5, 0.0, 0.0, -0.5}}, // -OH-0.5
      {{-0.5, 1.0, 0.0, 0.5}},  // -OH-0.5Na+
      {{0.5, 0.0, 0.0, 0.5}},   // -OH2+0.5
      {{0.5, -1.0, 0.0, -0.5}}, // -OH2+0.5NO3-
      {{0.0, 0.0, 0.0, 0.0}},   // Gbs
  }};
  for (long int ja = 0; ja < 5; ++ja) {
    check_near(node.Get_zjA(ja), expected_z[ja][0], kTol,
               "Get_zjA(" + std::to_string(ja) + ")");
    check_near(node.Get_zjB(ja), expected_z[ja][1], kTol,
               "Get_zjB(" + std::to_string(ja) + ")");
    check_near(node.Get_zjC(ja), expected_z[ja][2], kTol,
               "Get_zjC(" + std::to_string(ja) + ")");
    check_near(node.Get_bzj(ja), expected_z[ja][3], kTol,
               "Get_bzj(" + std::to_string(ja) + ")");
    check_near(node.Get_GamC(ja), ja < 4 ? 8.5 : 0.0, kTol,
               "Get_GamC(" + std::to_string(ja) + ")");
    check_near(node.Get_isoP1(ja), 0.0, kTol,
               "Get_isoP1(" + std::to_string(ja) + ")");
    check_near(node.Get_isoP2(ja), 0.0, kTol,
               "Get_isoP2(" + std::to_string(ja) + ")");
    check_near(node.Get_isoP3(ja), 0.0, kTol,
               "Get_isoP3(" + std::to_string(ja) + ")");
    check_cond(node.Get_dent(ja, 0) == (ja < 4 ? 1 : 0),
               "Get_dent(" + std::to_string(ja) + ", 0)");
    check_cond(node.Get_SATT(ja) == (ja < 4 ? 'C' : 'N'),
               "Get_SATT(" + std::to_string(ja) + ")");
  }

  // Site-type and surface-type parameters for Gibbsite_an (ist = 0, iss = 0)
  check_near(node.Get_GamCs(ka, 0), 8.5, kTol, "Get_GamCs(0, 0)");
  check_near(node.Get_GamCs(ka, 0, 0), 8.5, kTol, "Get_GamCs(0, 0, 0)");
  check_near(node.Get_CapCs(ka, 0), 0.0, kTol, "Get_CapCs(0, 0)");
  check_cond(node.Get_iss2ist(ka, 0) == 0, "Get_iss2ist(0, 0) == 0");
  check_near(node.Get_fAsurt(ka, 0), 1.0, kTol, "Get_fAsurt(0, 0)");
  check_near(node.Get_A1Tour(ka, 0), 0.0, kTol, "Get_A1Tour(0, 0)");
  check_near(node.Get_A2Tour(ka, 0), 0.0, kTol, "Get_A2Tour(0, 0)");
  check_near(node.Get_A3Tour(ka, 0), 0.0, kTol, "Get_A3Tour(0, 0)");

  check_near(node.Get_C1(ka, 0), 1.0, kTol, "Get_C1(0, 0)");
  check_near(node.Get_C2(ka, 0), 0.2, 1e-6, "Get_C2(0, 0)");
  check_near(node.Get_C3(ka, 0), 0.0, kTol, "Get_C3(0, 0)");
  check_near(node.Get_fzpSpd(ka, 0), 0.5, kTol, "Get_fzpSpd(0, 0)");
  check_near(node.Get_pEtaD(ka, 0), 0.0, kTol, "Get_pEtaD(0, 0)");
  check_near(node.Get_pEtaQ(ka, 0), 0.0, kTol, "Get_pEtaQ(0, 0)");
  check_cond(node.Get_EIL_code(ka, 0) == '3', "Get_EIL_code(0, 0) == '3'");
  check_cond(node.Get_SACT_code(ka, 0) == 'C', "Get_SACT_code(0, 0) == 'C'");

  // Exercise Set_* and Get_* round-trips, restoring original values cleanly
  const double orig_ssa = node.Get_SSA_ka(ka);
  node.Set_SSA_ka(42.5, ka);
  check_near(node.Get_SSA(k), 42.5, kTol, "Set_SSA_ka -> Get_SSA");
  check_near(node.Get_SSA_ka(ka), 42.5, kTol, "Set_SSA_ka -> Get_SSA_ka");
  check_near(node.pCNode()->aPH[k], 42500.0, 1e-6,
             "Set_SSA_ka synchronized CNode aPH (m2/kg)");
  node.Set_SSA_ka(orig_ssa, ka);

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

  // Use unused site/surface slot (1) for non-perturbing setter round-trips,
  // preserving original exported values (e.g. DBL_MIN padding from GEMSGUI).
  const double orig_gamcs1 = node.Get_GamCs(ka, 1);
  const double orig_capcs1 = node.Get_CapCs(ka, 1);
  const long int orig_iss2ist1 = node.Get_iss2ist(ka, 1);
  node.Set_GamCs(11.2, ka, 1);
  node.Get_GamCs(tmp_d, ka, 1);
  check_near(tmp_d, 11.2, kTol, "Set/Get_GamCs");
  node.Set_GamCs(orig_gamcs1, ka, 1);

  node.Set_CapCs(0.15, ka, 1);
  node.Get_CapCs(tmp_d, ka, 1);
  check_near(tmp_d, 0.15, kTol, "Set/Get_CapCs");
  node.Set_CapCs(orig_capcs1, ka, 1);

  node.Set_iss2ist(1, ka, 1);
  node.Get_iss2ist(tmp_i, ka, 1);
  check_cond(tmp_i == 1, "Set/Get_iss2ist");

  const double orig_fasurt1 = node.Get_fAsurt(ka, 1);
  const double orig_a1_1 = node.Get_A1Tour(ka, 1);
  const double orig_a2_1 = node.Get_A2Tour(ka, 1);
  const double orig_a3_1 = node.Get_A3Tour(ka, 1);
  node.Set_fAsurt(0.3, ka, 1);
  node.Get_fAsurt(tmp_d, ka, 1);
  check_near(tmp_d, 0.3, kTol, "Set/Get_fAsurt");
  node.Set_fAsurt(orig_fasurt1, ka, 1);

  node.Set_A1Tour(1.1, ka, 1);
  node.Set_A2Tour(2.2, ka, 1);
  node.Set_A3Tour(3.3, ka, 1);
  node.Get_A1Tour(tmp_d, ka, 1);
  check_near(tmp_d, 1.1, kTol, "Set/Get_A1Tour");
  node.Get_A2Tour(tmp_d, ka, 1);
  check_near(tmp_d, 2.2, kTol, "Set/Get_A2Tour");
  node.Get_A3Tour(tmp_d, ka, 1);
  check_near(tmp_d, 3.3, kTol, "Set/Get_A3Tour");
  node.Set_A1Tour(orig_a1_1, ka, 1);
  node.Set_A2Tour(orig_a2_1, ka, 1);
  node.Set_A3Tour(orig_a3_1, ka, 1);

  const double orig_c1_1 = node.Get_C1(ka, 1);
  const double orig_c2_1 = node.Get_C2(ka, 1);
  const double orig_c3_1 = node.Get_C3(ka, 1);
  const double orig_cz_1 = node.Get_fzpSpd(ka, 1);
  const double orig_etad_1 = node.Get_pEtaD(ka, 1);
  const double orig_etaq_1 = node.Get_pEtaQ(ka, 1);
  node.Set_C1(2.1, ka, 1);
  node.Set_C2(3.2, ka, 1);
  node.Set_C3(0.4, ka, 1);
  node.Set_fzpSpd(0.8, ka, 1);
  node.Set_pEtaD(-0.15, ka, 1);
  node.Set_pEtaQ(-0.25, ka, 1);
  node.Get_C1(tmp_d, ka, 1);
  check_near(tmp_d, 2.1, kTol, "Set/Get_C1");
  node.Get_C2(tmp_d, ka, 1);
  check_near(tmp_d, 3.2, kTol, "Set/Get_C2");
  node.Get_C3(tmp_d, ka, 1);
  check_near(tmp_d, 0.4, kTol, "Set/Get_C3");
  node.Get_fzpSpd(tmp_d, ka, 1);
  check_near(tmp_d, 0.8, kTol, "Set/Get_fzpSpd");
  node.Get_pEtaD(tmp_d, ka, 1);
  check_near(tmp_d, -0.15, kTol, "Set/Get_pEtaD");
  node.Get_pEtaQ(tmp_d, ka, 1);
  check_near(tmp_d, -0.25, kTol, "Set/Get_pEtaQ");
  node.Set_C1(orig_c1_1, ka, 1);
  node.Set_C2(orig_c2_1, ka, 1);
  node.Set_C3(orig_c3_1, ka, 1);
  node.Set_fzpSpd(orig_cz_1, ka, 1);
  node.Set_pEtaD(orig_etad_1, ka, 1);
  node.Set_pEtaQ(orig_etaq_1, ka, 1);

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
  node.Set_iss2ist(orig_iss2ist1, ka, 1);
  check_cond(node.Get_nst(0) == 1 && node.ja2ist(0) == 0,
             "Restored Get_nst(0) == 1 & ja2ist(0) == 0");

  const char orig_eil1 = node.Get_EIL_code(ka, 1);
  const char orig_sact1 = node.Get_SACT_code(ka, 1);
  node.Set_EIL_code('T', ka, 1);
  node.Set_SACT_code('L', ka, 1);
  node.Set_SATT('L', 0);
  node.Get_EIL_code(tmp_c, ka, 1);
  check_cond(tmp_c == 'T', "Set/Get_EIL_code");
  node.Get_SACT_code(tmp_c, ka, 1);
  check_cond(tmp_c == 'L', "Set/Get_SACT_code");
  node.Get_SATT(tmp_c, 0);
  check_cond(tmp_c == 'L', "Set/Get_SATT");
  node.Set_EIL_code(orig_eil1, ka, 1);
  node.Set_SACT_code(orig_sact1, ka, 1);
  node.Set_SATT('C', 0);
}

void test_gibbsite_an_outputs(const TNode &node, long int ka, long int k) {
  check_near(node.Get_Zeff(), 1.0, 1e-4, "Get_Zeff() for NaNO3");
  double zeff_r = 0.0;
  node.Get_Zeff(zeff_r);
  check_near(zeff_r, node.Get_Zeff(), kTol, "Get_Zeff ref");

  const double c_dl = node.Get_C_dl(ka, 0);
  check_cond(std::isfinite(c_dl) && c_dl > 0.15 && c_dl < 0.25,
             "Get_C_dl(0, 0) in expected 3LM range");
  check_cond(std::isfinite(node.Get_h_A(ka, 0)) && node.Get_h_A(ka, 0) >= 0.0,
             "Get_h_A(0, 0) >= 0");
  check_cond(std::isfinite(node.Get_h_B(ka, 0)) && node.Get_h_B(ka, 0) >= 0.0,
             "Get_h_B(0, 0) >= 0");
  check_cond(std::isfinite(node.Get_h_C(ka, 0)) && node.Get_h_C(ka, 0) >= 0.0,
             "Get_h_C(0, 0) >= 0");
  check_cond(node.Get_h_dl(ka, 0) > 0.0, "Get_h_dl(0, 0) > 0");
  check_cond(node.Get_h_spz(ka, 0) > 0.0, "Get_h_spz(0, 0) > 0");

  double v = 0.0;
  node.Get_C_dl(v, ka, 0);
  check_near(v, c_dl, kTol, "Get_C_dl ref");
  node.Get_h_A(v, ka, 0);
  check_near(v, node.Get_h_A(ka, 0), kTol, "Get_h_A ref");
  node.Get_h_B(v, ka, 0);
  check_near(v, node.Get_h_B(ka, 0), kTol, "Get_h_B ref");
  node.Get_h_C(v, ka, 0);
  check_near(v, node.Get_h_C(ka, 0), kTol, "Get_h_C ref");
  node.Get_h_dl(v, ka, 0);
  check_near(v, node.Get_h_dl(ka, 0), kTol, "Get_h_dl ref");
  node.Get_h_spz(v, ka, 0);
  check_near(v, node.Get_h_spz(ka, 0), kTol, "Get_h_spz ref");

  const double mass_s = node.Get_Mass_s(k);
  const double mols_s = node.Get_Mols_s(k);
  check_near(mass_s, 30.0, 0.05, "Get_Mass_s(Gibbsite_an) in g");
  check_near(mols_s, 0.384597, 1e-4, "Get_Mols_s(Gibbsite_an) in mol");
  node.Get_Mass_s(v, k);
  check_near(v, mass_s, kTol, "Get_Mass_s ref");
  node.Get_Mols_s(v, k);
  check_near(v, mols_s, kTol, "Get_Mols_s ref");

  const double etaA = node.Get_eta_A(ka, 0);
  const double etaB = node.Get_eta_B(ka, 0);
  const double etaC = node.Get_eta_C(ka, 0);
  const double etaD = node.Get_eta_D(ka, 0);
  check_near(etaA + etaB + etaC + etaD, 0.0, 1e-12,
             "Gibbsite_an EIL plane electroneutrality");

  node.Get_eta_A(v, ka, 0);
  check_near(v, etaA, kTol, "Get_eta_A ref");
  node.Get_eta_B(v, ka, 0);
  check_near(v, etaB, kTol, "Get_eta_B ref");
  node.Get_eta_C(v, ka, 0);
  check_near(v, etaC, kTol, "Get_eta_C ref");
  node.Get_eta_D(v, ka, 0);
  check_near(v, etaD, kTol, "Get_eta_D ref");

  const double psiA = node.Get_psi_A(ka, 0);
  const double psiB = node.Get_psi_B(ka, 0);
  const double psiC = node.Get_psi_C(ka, 0);
  const double psiD = node.Get_psi_D(ka, 0);
  const double psiZ = node.Get_psi_Z(ka, 0);
  check_cond(std::isfinite(psiA) && std::fabs(psiA) > 1e-5,
             "psi_A finite and non-zero");
  check_cond(std::isfinite(psiB) && std::fabs(psiB) > 1e-5,
             "psi_B finite and non-zero");
  check_cond(std::isfinite(psiC), "psi_C finite");
  check_cond(std::isfinite(psiD) && std::fabs(psiD) > 1e-5,
             "psi_D finite and non-zero");
  check_cond(std::isfinite(psiZ) && std::fabs(psiZ) > 1e-5,
             "psi_Z finite and non-zero");
  check_cond(std::fabs(psiZ) <= std::fabs(psiD) + 1e-12 &&
                 (psiZ * psiD >= 0.0),
             "psi_Z lies between 0 and psi_D");

  node.Get_psi_A(v, ka, 0);
  check_near(v, psiA, kTol, "Get_psi_A ref");
  node.Get_psi_B(v, ka, 0);
  check_near(v, psiB, kTol, "Get_psi_B ref");
  node.Get_psi_C(v, ka, 0);
  check_near(v, psiC, kTol, "Get_psi_C ref");
  node.Get_psi_D(v, ka, 0);
  check_near(v, psiD, kTol, "Get_psi_D ref");
  node.Get_psi_Z(v, ka, 0);
  check_near(v, psiZ, kTol, "Get_psi_Z ref");

  // Activity coefficient terms across all 5 Gibbsite_an species (ja = 0..4)
  double surface_species_sum = 0.0;
  for (long int ja = 0; ja < 5; ++ja) {
    const long int j = node.ja2j(ja);
    const double lnGam_ja = node.Get_lnGam_ja(ja);
    check_near(node.Get_lnGam(ja), lnGam_ja, kTol,
               "Get_lnGam alias ja=" + std::to_string(ja));
    check_near(node.Get_lnSACT_ja(ja), node.Get_lnSACT(j), kTol,
               "Get_lnSACT_ja vs Get_lnSACT ja=" + std::to_string(ja));
    check_near(node.Get_lnExpF_ja(ja), node.Get_lnExpF(j), kTol,
               "Get_lnExpF_ja vs Get_lnExpF ja=" + std::to_string(ja));
    check_near(node.Get_lnCoul_ja(ja), node.Get_lnCoul(j), kTol,
               "Get_lnCoul_ja vs Get_lnCoul ja=" + std::to_string(ja));

    node.Get_lnGam_ja(v, ja);
    check_near(v, lnGam_ja, kTol, "Get_lnGam_ja ref");
    node.Get_lnGam(v, ja);
    check_near(v, lnGam_ja, kTol, "Get_lnGam ref");
    node.Get_lnSACT_ja(v, ja);
    check_near(v, node.Get_lnSACT(j), kTol, "Get_lnSACT_ja ref");
    node.Get_lnExpF_ja(v, ja);
    check_near(v, node.Get_lnExpF(j), kTol, "Get_lnExpF_ja ref");
    node.Get_lnCoul_ja(v, ja);
    check_near(v, node.Get_lnCoul(j), kTol, "Get_lnCoul_ja ref");
    node.Get_lnSACT(v, j);
    check_near(v, node.Get_lnSACT(j), kTol, "Get_lnSACT ref");
    node.Get_lnExpF(v, j);
    check_near(v, node.Get_lnExpF(j), kTol, "Get_lnExpF ref");
    node.Get_lnCoul(v, j);
    check_near(v, node.Get_lnCoul(j), kTol, "Get_lnCoul ref");

    if (ja < 4) {
      check_cond(std::fabs(node.Get_lnCoul_ja(ja)) > 1e-4,
                 "Coulomb term non-zero for charged species ja=" +
                     std::to_string(ja));
      check_cond(std::fabs(node.Get_lnSACT_ja(ja)) > 1e-4,
                 "SACT term non-zero for surface species ja=" +
                     std::to_string(ja));
      check_near(node.Get_lnExpF_ja(ja), 0.0, kTol,
                 "lnExpF == 0 for Langmuir SACT ja=" + std::to_string(ja));
      surface_species_sum += node.pCNode()->xDC[j];
    } else {
      // Carrier DC 'Gbs'
      check_near(node.Get_lnCoul_ja(ja), 0.0, kTol,
                 "Coulomb term == 0 for carrier Gbs");
      check_near(node.Get_lnSACT_ja(ja), 0.0, kTol,
                 "SACT term == 0 for carrier Gbs");
    }
  }

  check_near(node.Get_Xfts(0, 0), surface_species_sum, 1e-8,
             "Get_Xfts(0, 0) equals sum of Gibbsite_an surface species");
  check_near(node.Get_Dfts(0, 0), surface_species_sum, 1e-8,
             "Get_Dfts(0, 0) equals sum of occupied monodentate sites");
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

void check_roundtrip_state(const TNode &actual, const TNode &expected,
                           long int ka, long int ist, long int ja_sample,
                           const std::string &label) {
  check_cond(actual.Get_Lads() == expected.Get_Lads(), label + " Lads");
  check_cond(actual.Get_FIa() == expected.Get_FIa(), label + " FIa");
  check_cond(actual.Get_FIat() == expected.Get_FIat(), label + " FIat");
  check_near(actual.Get_h_spz(ka, ist), expected.Get_h_spz(ka, ist), kTol,
             label + " Get_h_spz");
  check_near(actual.Get_zjB(ja_sample), expected.Get_zjB(ja_sample), kTol,
             label + " Get_zjB");
  check_cond(actual.Get_dent(ja_sample, expected.Get_iss(ja_sample)) ==
                 expected.Get_dent(ja_sample, expected.Get_iss(ja_sample)),
             label + " Get_dent");
  check_near(actual.Get_Xfts(0, ist), expected.Get_Xfts(0, ist), kTol,
             label + " Get_Xfts");
  check_near(actual.Get_Dfts(expected.Get_iss(ja_sample), ist),
             expected.Get_Dfts(expected.Get_iss(ja_sample), ist), kTol,
             label + " Get_Dfts");
  check_near(actual.Get_psi_A(ka, ist), expected.Get_psi_A(ka, ist), kTol,
             label + " Get_psi_A");
  check_near(actual.Get_psi_Z(ka, ist), expected.Get_psi_Z(ka, ist), kTol,
             label + " Get_psi_Z");
}

void test_copy_and_binary_roundtrip(TestNode &node,
                                    const std::filesystem::path &tmp_dir,
                                    long int ka, long int ist,
                                    long int ja_sample) {
  TestNode copyNode(static_cast<const TNode &>(node));
  check_roundtrip_state(copyNode, node, ka, ist, ja_sample, "copy");

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
  check_roundtrip_state(binNode, node, ka, ist, ja_sample, "binary round-trip");
}

struct MontSpeciesMeta {
  const char *name;
  long int ist;
  long int iss;
  double zA;
  double zB;
  double GamC;
};

constexpr std::array<MontSpeciesMeta, 22> kMontSpecies = {{
    {"=SsO-", 1, 3, -1.0, 0.0, 0.086},
    {"=SsO-Na+", 1, 3, -1.0, 1.0, 0.086},
    {"=SsOH2+", 1, 3, 1.0, 0.0, 0.086},
    {"=SsOH2+ClO4-", 1, 3, 1.0, -1.0, 0.086},
    {"=SsOH@", 1, 3, 0.0, 0.0, 0.086},
    {"=SsONi(OH)2-", 1, 3, -1.0, 0.0, 0.086},
    {"=SsONi+", 1, 3, 1.0, 0.0, 0.086},
    {"=SsONiOH@", 1, 3, 0.0, 0.0, 0.086},
    {"=SvO-", 1, 2, -1.0, 0.0, 1.724},
    {"=SvO-Na+", 1, 2, -1.0, 1.0, 1.724},
    {"=SvOH2+", 1, 2, 1.0, 0.0, 1.724},
    {"=SvOH2+ClO4-", 1, 2, 1.0, -1.0, 1.724},
    {"=SvOH@", 1, 2, 0.0, 0.0, 1.724},
    {"=SvONi+", 1, 2, 1.0, 0.0, 1.724},
    {"=SwO-", 1, 1, -1.0, 0.0, 0.67},
    {"=SwO-Na+", 1, 1, -1.0, 1.0, 0.67},
    {"=SwOH2+", 1, 1, 1.0, 0.0, 0.67},
    {"=SwOH2+ClO4-", 1, 1, 1.0, -1.0, 0.67},
    {"=SwOH@", 1, 1, 0.0, 0.0, 0.67},
    {"=X-H+@", 0, 0, 0.0, 0.0, 0.66},
    {"=X-Na+@", 0, 0, 0.0, 0.0, 0.66},
    {"=X-Ni+2@", 0, 0, 0.0, 0.0, 0.66},
}};

void test_mont2tlm_indices(const TNode &node, long int &ka_out,
                           long int &k_out) {
  check_cond(node.Get_Lads() == 22, "Mont2tlm Get_Lads() == 22");
  check_cond(node.Get_FIa() == 1, "Mont2tlm Get_FIa() == 1");
  check_cond(node.Get_FIat() == 4, "Mont2tlm Get_FIat() == 4");

  long int lads = -1, fia = -1, fiat = -1;
  node.Get_Lads(lads);
  node.Get_FIa(fia);
  node.Get_FIat(fiat);
  check_cond(lads == 22 && fia == 1 && fiat == 4,
             "Mont2tlm scalar reference overloads");

  const long int ka = node.k2ka("Mont2tlm");
  check_cond(ka == 0, "k2ka(\"Mont2tlm\") == 0");
  const long int k = node.ka2k(ka);
  check_cond(k == 2, "Mont2tlm ka2k(0) == 2");
  check_cond(node.k2ka(k) == ka, "Mont2tlm k2ka(2) == 0");
  check_cond(node.Get_Lads_ka(ka) == 22, "Get_Lads_ka(Mont2tlm) == 22");
  check_cond(node.k2ka(0) == -1 && node.ka2k(1) == -1,
             "Mont2tlm invalid phase mappings");

  for (long int ja = 0; ja < 22; ++ja) {
    const auto &sp = kMontSpecies[ja];
    const long int j = 31 + ja;
    check_cond(node.j_ka2ja(j, ka) == ja,
               "Mont2tlm j_ka2ja(" + std::to_string(j) + ", ka)");
    check_cond(node.j_ka2ja(sp.name, ka) == ja,
               std::string("Mont2tlm j_ka2ja(\"") + sp.name + "\", ka)");
    check_cond(node.ja2j(ja) == j,
               "Mont2tlm ja2j(" + std::to_string(ja) + ")");
    check_cond(node.Get_ist(ja) == sp.ist,
               "Mont2tlm Get_ist(" + std::to_string(ja) + ")");
    check_cond(node.ja2ist(ja) == sp.ist,
               "Mont2tlm ja2ist(" + std::to_string(ja) + ")");
    check_cond(node.Get_iss(ja) == sp.iss,
               "Mont2tlm Get_iss(" + std::to_string(ja) + ")");
    check_cond(node.Get_nst(ja) == 1,
               "Mont2tlm Get_nst(" + std::to_string(ja) + ")");

    long int ja_r = -1, j_r = -1, ist_r = -1, nst_r = -1, iss_r = -1,
             ja2ist_r = -1;
    node.j_ka2ja(ja_r, sp.name, ka);
    node.ja2j(j_r, ja);
    node.Get_ist(ist_r, ja);
    node.Get_nst(nst_r, ja);
    node.Get_iss(iss_r, ja);
    node.ja2ist(ja2ist_r, ja);
    check_cond(ja_r == ja && j_r == j && ist_r == sp.ist && iss_r == sp.iss &&
                   ja2ist_r == sp.ist && nst_r == 1,
               "Mont2tlm species index ref overloads");
  }

  check_cond(node.j_ka2ja(0, ka) == -1 && node.ja2j(22) == -1 &&
                 node.Get_nst(22) == -1,
             "Mont2tlm out-of-range species guards");

  ka_out = ka;
  k_out = k;
}

void test_mont2tlm_inputs(TNode &node, long int ka, long int k) {
  check_near(node.Get_SSA(k), 800.0, 1e-6, "Get_SSA(Mont2tlm)");
  check_near(node.Get_SSA_ka(ka), 800.0, 1e-6, "Get_SSA_ka(Mont2tlm)");

  for (long int ja = 0; ja < 22; ++ja) {
    const auto &sp = kMontSpecies[ja];
    check_near(node.Get_zjA(ja), sp.zA, kTol,
               "Mont2tlm Get_zjA(" + std::to_string(ja) + ")");
    check_near(node.Get_zjB(ja), sp.zB, kTol,
               "Mont2tlm Get_zjB(" + std::to_string(ja) + ")");
    check_near(node.Get_zjC(ja), 0.0, kTol,
               "Mont2tlm Get_zjC(" + std::to_string(ja) + ")");
    check_near(node.Get_GamC(ja), sp.GamC, 1e-6,
               "Mont2tlm Get_GamC(" + std::to_string(ja) + ")");
    check_cond(node.Get_dent(ja, sp.iss) == 1,
               "Mont2tlm Get_dent(" + std::to_string(ja) + ")");
    check_cond(node.Get_SATT(ja) == 'C',
               "Mont2tlm Get_SATT(" + std::to_string(ja) + ")");
  }

  // Site-type parameters across all 4 active site types (iss = 0..3)
  constexpr std::array<long int, 4> expected_iss2ist = {{0, 1, 1, 1}};
  constexpr std::array<double, 4> expected_gamcs = {{0.66, 2.439, 2.439, 0.122}};
  constexpr std::array<double, 4> expected_capcs = {{0.87, 0.04, 0.04, 0.002}};
  for (long int iss = 0; iss < 4; ++iss) {
    check_cond(node.Get_iss2ist(ka, iss) == expected_iss2ist[iss],
               "Mont2tlm Get_iss2ist(0, " + std::to_string(iss) + ")");
    check_near(node.Get_GamCs(ka, iss), expected_gamcs[iss], 1e-5,
               "Mont2tlm Get_GamCs(0, " + std::to_string(iss) + ")");
    check_near(node.Get_CapCs(ka, iss), expected_capcs[iss], 1e-6,
               "Mont2tlm Get_CapCs(0, " + std::to_string(iss) + ")");
  }

  // Surface type 0 (NEM) and surface type 1 (3TL EIL)
  check_cond(node.Get_EIL_code(ka, 0) == 'N', "Mont2tlm ist=0 EIL == 'N'");
  check_cond(node.Get_SACT_code(ka, 0) == 'I', "Mont2tlm ist=0 SACT == 'I'");
  check_near(node.Get_fAsurt(ka, 0), 0.9876, 1e-6, "Mont2tlm ist=0 fAsurt");
  check_near(node.Get_C1(ka, 0), 0.0, kTol, "Mont2tlm ist=0 C1 == 0");
  check_near(node.Get_C2(ka, 0), 0.0, kTol, "Mont2tlm ist=0 C2 == 0");
  check_near(node.Get_C3(ka, 0), 0.0, kTol, "Mont2tlm ist=0 C3 == 0");
  check_near(node.Get_fzpSpd(ka, 0), 0.0, kTol, "Mont2tlm ist=0 fzpSpd == 0");
  check_near(node.Get_pEtaD(ka, 0), 0.66, 1e-6, "Mont2tlm ist=0 pEtaD");
  check_near(node.Get_pEtaQ(ka, 0), 0.87, 1e-6, "Mont2tlm ist=0 pEtaQ");

  check_cond(node.Get_EIL_code(ka, 1) == 'T', "Mont2tlm ist=1 EIL == 'T'");
  check_cond(node.Get_SACT_code(ka, 1) == 'I', "Mont2tlm ist=1 SACT == 'I'");
  check_near(node.Get_fAsurt(ka, 1), 0.0124, 1e-6, "Mont2tlm ist=1 fAsurt");
  check_near(node.Get_C1(ka, 1), 1.1, 1e-6, "Mont2tlm ist=1 C1 == 1.1");
  check_near(node.Get_C2(ka, 1), 0.2, 1e-6, "Mont2tlm ist=1 C2 == 0.2");
  check_near(node.Get_C3(ka, 1), 0.0, kTol, "Mont2tlm ist=1 C3 == 0");
  check_near(node.Get_fzpSpd(ka, 1), 0.5, 1e-6, "Mont2tlm ist=1 fzpSpd == 0.5");
  check_near(node.Get_pEtaD(ka, 1), 0.0, kTol, "Mont2tlm ist=1 pEtaD == 0");
  check_near(node.Get_pEtaQ(ka, 1), 0.0, kTol, "Mont2tlm ist=1 pEtaQ == 0");

  // Exercise Set/Get round-trips on surface type 0/1 & site/species slots and restore
  const double orig_ssa = node.Get_SSA_ka(ka);
  node.Set_SSA_ka(750.0, ka);
  check_near(node.Get_SSA(k), 750.0, kTol, "Mont2tlm Set_SSA_ka");
  check_near(node.pCNode()->aPH[k], 750000.0, 1e-4,
             "Mont2tlm Set_SSA_ka aPH sync");
  node.Set_SSA_ka(orig_ssa, ka);

  // Neutral species ja=19 (ist=0) charge round-trip (preserves bzj=0)
  node.Set_zjA(-0.25, 19);
  node.Set_zjB(0.75, 19);
  node.Set_zjC(-0.125, 19);
  check_near(node.Get_bzj(19), 0.375, kTol, "Mont2tlm Set_zj* -> Get_bzj");
  node.Set_zjA(0.0, 19);
  node.Set_zjB(0.0, 19);
  node.Set_zjC(0.0, 19);

  // Multi-surface-type denticity test: add bond of 3TL species ja=0 (iss=3, ist=1)
  // onto NEM site iss=0 (ist=0), verify nst=2 and ja2ist=-1, then restore
  node.Set_dent(2L, 0, 0);
  check_cond(node.Get_nst(0) == 2, "Mont2tlm cross-surface nst == 2");
  check_cond(node.ja2ist(0) == -1,
             "Mont2tlm ja2ist == -1 across NEM + 3TL surface types");
  node.Set_dent(0L, 0, 0);
  check_cond(node.Get_nst(0) == 1 && node.ja2ist(0) == 1,
             "Mont2tlm restored nst == 1 & ja2ist == 1");

  double tmp_d = 0.0;
  char tmp_c = '\0';
  const double orig_a1 = node.Get_A1Tour(ka, 1);
  node.Set_A1Tour(0.42, ka, 1);
  node.Get_A1Tour(tmp_d, ka, 1);
  check_near(tmp_d, 0.42, kTol, "Mont2tlm Set/Get_A1Tour");
  node.Set_A1Tour(orig_a1, ka, 1);

  node.Set_EIL_code('M', ka, 1);
  node.Get_EIL_code(tmp_c, ka, 1);
  check_cond(tmp_c == 'M', "Mont2tlm Set/Get_EIL_code");
  node.Set_EIL_code('T', ka, 1);
}

void test_mont2tlm_outputs(const TNode &node, long int ka, long int k) {
  check_near(node.Get_Zeff(), 1.0, 1e-4, "Mont2tlm Get_Zeff() for NaCl");
  check_near(node.Get_Mass_s(k), 0.98, 1e-3, "Get_Mass_s(Mont2tlm) in g");
  check_near(node.Get_Mols_s(k), 0.0, kTol,
             "Get_Mols_s(Mont2tlm) == 0 (site-balance implicit sorbent)");

  // Surface type 0 (NEM): zero EIL potentials, plane charges, and C_dl
  check_near(node.Get_C_dl(ka, 0), 0.0, kTol, "NEM C_dl == 0");
  check_near(node.Get_h_A(ka, 0), 0.0, kTol, "NEM h_A == 0");
  check_near(node.Get_h_B(ka, 0), 0.0, kTol, "NEM h_B == 0");
  check_near(node.Get_h_C(ka, 0), 0.0, kTol, "NEM h_C == 0");
  check_near(node.Get_h_spz(ka, 0), 0.0, kTol, "NEM h_spz == 0");
  check_near(node.Get_psi_A(ka, 0), 0.0, kTol, "NEM psi_A == 0");
  check_near(node.Get_psi_B(ka, 0), 0.0, kTol, "NEM psi_B == 0");
  check_near(node.Get_psi_C(ka, 0), 0.0, kTol, "NEM psi_C == 0");
  check_near(node.Get_psi_D(ka, 0), 0.0, kTol, "NEM psi_D == 0");
  check_near(node.Get_psi_Z(ka, 0), 0.0, kTol, "NEM psi_Z == 0");
  check_near(node.Get_eta_A(ka, 0), 0.0, kTol, "NEM eta_A == 0");
  check_near(node.Get_eta_B(ka, 0), 0.0, kTol, "NEM eta_B == 0");
  check_near(node.Get_eta_C(ka, 0), 0.0, kTol, "NEM eta_C == 0");
  check_near(node.Get_eta_D(ka, 0), 0.0, kTol, "NEM eta_D == 0");

  // Surface type 1 (3TL EIL): active EIL outputs and plane charge balance
  const double c_dl = node.Get_C_dl(ka, 1);
  check_cond(std::isfinite(c_dl) && c_dl > 0.5 && c_dl < 1.0,
             "3TL Get_C_dl(0, 1) in expected range");
  check_cond(node.Get_h_dl(ka, 1) > 0.5, "3TL Get_h_dl(0, 1) > 0");
  check_near(node.Get_h_spz(ka, 1), 0.5 * node.Get_h_dl(ka, 1), 1e-8,
             "3TL Get_h_spz == 0.5 * h_dl");

  const double etaA = node.Get_eta_A(ka, 1);
  const double etaB = node.Get_eta_B(ka, 1);
  const double etaC = node.Get_eta_C(ka, 1);
  const double etaD = node.Get_eta_D(ka, 1);
  check_cond(std::fabs(etaA) > 1e-12 && std::fabs(etaD) > 1e-12,
             "3TL non-zero EIL plane charges");
  check_near(etaA + etaB + etaC + etaD, 0.0, 1e-12,
             "3TL EIL plane electroneutrality");

  const double psiA = node.Get_psi_A(ka, 1);
  const double psiB = node.Get_psi_B(ka, 1);
  const double psiD = node.Get_psi_D(ka, 1);
  const double psiZ = node.Get_psi_Z(ka, 1);
  check_cond(std::fabs(psiA) > 1e-4 && std::fabs(psiB) > 1e-4 &&
                 std::fabs(psiD) > 1e-5 && std::fabs(psiZ) > 1e-5,
             "3TL non-zero potentials psi_A/B/D/Z");
  check_cond(std::fabs(psiZ) <= std::fabs(psiD) + 1e-12 &&
                 (psiZ * psiD >= 0.0),
             "3TL psi_Z lies between 0 and psi_D");

  // Activity terms and site/surface-type totals across all 22 species
  double nem_sum = 0.0, eil_sum = 0.0;
  std::array<double, 4> site_sums{};
  for (long int ja = 0; ja < 22; ++ja) {
    const auto &sp = kMontSpecies[ja];
    const long int j = node.ja2j(ja);
    const double xdc = node.pCNode()->xDC[j];
    site_sums[sp.iss] += xdc;
    if (sp.ist == 0)
      nem_sum += xdc;
    else
      eil_sum += xdc;

    // Site-balance SACT='I' models have zero lnGam / lnSACT / lnExpF
    check_near(node.Get_lnGam_ja(ja), 0.0, kTol,
               "Mont2tlm lnGam == 0 for ja=" + std::to_string(ja));
    check_near(node.Get_lnSACT_ja(ja), 0.0, kTol,
               "Mont2tlm lnSACT == 0 for ja=" + std::to_string(ja));
    check_near(node.Get_lnExpF_ja(ja), 0.0, kTol,
               "Mont2tlm lnExpF == 0 for ja=" + std::to_string(ja));
    check_near(node.Get_lnCoul_ja(ja), node.Get_lnCoul(j), kTol,
               "Mont2tlm Get_lnCoul_ja vs Get_lnCoul ja=" + std::to_string(ja));

    if (sp.ist == 0 || (sp.zA == 0.0 && sp.zB == 0.0)) {
      check_near(node.Get_lnCoul_ja(ja), 0.0, kTol,
                 "NEM/neutral surface species has lnCoul == 0 ja=" +
                     std::to_string(ja));
    } else {
      check_cond(std::fabs(node.Get_lnCoul_ja(ja)) > 1e-4,
                 "Charged 3TL species has non-zero lnCoul ja=" +
                     std::to_string(ja));
    }
  }

  check_near(node.Get_Xfts(0, 0), nem_sum, 1e-9,
             "Mont2tlm Get_Xfts(0, 0) equals NEM species total");
  check_near(node.Get_Xfts(0, 1), eil_sum, 1e-9,
             "Mont2tlm Get_Xfts(0, 1) equals 3TL species total");
  check_near(node.Get_Dfts(0, 0), site_sums[0], 1e-9,
             "Mont2tlm Get_Dfts(0, 0) equals NEM site 0 total");
  for (long int iss = 1; iss < 4; ++iss) {
    check_near(node.Get_Dfts(iss, 1), site_sums[iss], 1e-9,
               "Mont2tlm Get_Dfts(" + std::to_string(iss) +
                   ", 1) equals 3TL site total");
  }
}

void run_gibbsite_test() {
  const std::filesystem::path data_dir(GEMS3K_GIBBSITE_DATA_DIR);
  const std::filesystem::path lst_path = data_dir / "Gibbsite-dat.lst";
  check_cond(std::filesystem::exists(lst_path),
             "Missing exported GEMS3K fileset: " + lst_path.string());

  TestNode node;
  const long int init_rc = node.GEM_init(lst_path.string().c_str());
  check_cond(init_rc == 0, "Gibbsite GEM_init returned 0");

  // Capture reference equilibrium values loaded from Gibbsite-dbr-0-0000.json
  check_near(node.pCNode()->TK, 293.15, 1e-6, "Exported TK == 293.15");
  check_near(node.pCNode()->P, 100000.0, 1e-3, "Exported P == 1e5");
  const double ref_pH = node.pCNode()->pH;
  const double ref_IC = node.pCNode()->IC;
  const double ref_Gs = node.pCNode()->Gs;
  std::array<double, 5> ref_xDC{};
  for (long int ja = 0; ja < 5; ++ja) {
    ref_xDC[ja] = node.pCNode()->xDC[19 + ja];
  }

  // Recompute equilibrium state from the exported GEMS3K fileset
  node.pCNode()->NodeStatusCH = NEED_GEM_AIA;
  const long int run_rc_aia = node.GEM_run(false);
  check_cond(run_rc_aia == OK_GEM_AIA || run_rc_aia == OK_GEM_SIA,
             "Recomputed equilibrium via AIA converged");
  check_near(node.pCNode()->pH, ref_pH, 0.1,
             "Recomputed AIA pH consistent with exported SysEq");
  check_near(node.pCNode()->IC, ref_IC, 1e-4,
             "Recomputed AIA ionic strength consistent with exported SysEq");
  check_near(node.pCNode()->Gs, ref_Gs, 1e-3,
             "Recomputed AIA Gibbs energy consistent with exported SysEq");
  for (long int ja = 0; ja < 5; ++ja) {
    check_near(node.pCNode()->xDC[19 + ja], ref_xDC[ja], 5e-5,
               "Recomputed AIA Gibbsite_an xDC[" + std::to_string(19 + ja) +
                   "]");
  }

  // Also verify SIA recomputation starting from the converged MULTI state
  const double aia_pH = node.pCNode()->pH;
  node.pCNode()->NodeStatusCH = NEED_GEM_SIA;
  const long int run_rc_sia = node.GEM_run(false);
  check_cond(run_rc_sia == OK_GEM_SIA || run_rc_sia == OK_GEM_AIA,
             "Recomputed equilibrium via SIA converged");
  check_near(node.pCNode()->pH, aia_pH, 1e-3,
             "Recomputed SIA pH reproduces AIA equilibrium");

  long int ka = -1, k = -1;
  test_gibbsite_an_indices(node, ka, k);
  test_gibbsite_an_outputs(node, ka, k);
  test_gibbsite_an_inputs(node, ka, k);

  TempDir tmp("gems3k_sorption_api_test_gibbsite");
  test_copy_and_binary_roundtrip(node, tmp.path(), ka, 0, 1);

  const double psiA_before = node.Get_psi_A(ka, 0);
  const double psiZ_before = node.Get_psi_Z(ka, 0);
  const double hspz_before = node.Get_h_spz(ka, 0);

  // Perturb C1 and fzpSpd on Gibbsite_an surface type 0 and recompute
  node.Set_C1(2.0, ka, 0);
  node.Set_fzpSpd(0.25, ka, 0);
  node.pCNode()->NodeStatusCH = NEED_GEM_AIA;
  const long int run_rc2 = node.GEM_run(false);
  check_cond(run_rc2 == OK_GEM_AIA || run_rc2 == OK_GEM_SIA,
             "GEM_run(false) after Set_C1/Set_fzpSpd converged");

  check_near(node.Get_C1(ka, 0), 2.0, kTol, "Persisted Set_C1");
  check_near(node.Get_fzpSpd(ka, 0), 0.25, kTol, "Persisted Set_fzpSpd");
  check_near(node.Get_h_spz(ka, 0), 0.5 * hspz_before, 1e-4,
             "Get_h_spz halved when fzpSpd halved at ~same ionic strength");
  check_cond(std::fabs(node.Get_psi_A(ka, 0) - psiA_before) > 1e-4,
             "Get_psi_A changed after Set_C1");
  check_cond(std::fabs(node.Get_psi_Z(ka, 0) - psiZ_before) > 1e-6 &&
                 std::fabs(node.Get_psi_Z(ka, 0)) >
                     std::fabs(psiZ_before) * 1.1,
             "Get_psi_Z magnitude increased after halving fzpSpd");

  std::cout << "sorption_api_test passed all checks on "
               "AluSiOMSA:G:Gibbsite:1:0:1:20:0: "
            << "(pH=" << node.pCNode()->pH << ", IC=" << node.pCNode()->IC
            << ", psiA: " << psiA_before << " -> " << node.Get_psi_A(ka, 0)
            << ", psiZ: " << psiZ_before << " -> " << node.Get_psi_Z(ka, 0)
            << ")\n";
}

void run_claysorp_test() {
  const std::filesystem::path data_dir(GEMS3K_CLAYSORP_DATA_DIR);
  const std::filesystem::path lst_path = data_dir / "Mont2tlNi01-dat.lst";
  check_cond(std::filesystem::exists(lst_path),
             "Missing exported GEMS3K fileset: " + lst_path.string());

  TestNode node;
  const long int init_rc = node.GEM_init(lst_path.string().c_str());
  check_cond(init_rc == 0, "Mont2tlNi01 GEM_init returned 0");

  check_near(node.pCNode()->TK, 298.15, 1e-6, "Mont2tlNi01 TK == 298.15");
  check_near(node.pCNode()->P, 100000.0, 1e-3, "Mont2tlNi01 P == 1e5");
  const double ref_pH = node.pCNode()->pH;
  const double ref_IC = node.pCNode()->IC;
  const double ref_Gs = node.pCNode()->Gs;
  std::array<double, 22> ref_xDC{};
  for (long int ja = 0; ja < 22; ++ja) {
    ref_xDC[ja] = node.pCNode()->xDC[31 + ja];
  }

  // Recompute equilibrium state via AIA and SIA
  node.pCNode()->NodeStatusCH = NEED_GEM_AIA;
  const long int run_rc_aia = node.GEM_run(false);
  check_cond(run_rc_aia == OK_GEM_AIA || run_rc_aia == OK_GEM_SIA,
             "Mont2tlNi01 recomputed AIA converged");
  check_near(node.pCNode()->pH, ref_pH, 0.01,
             "Mont2tlNi01 AIA pH matches exported SysEq");
  check_near(node.pCNode()->IC, ref_IC, 1e-4,
             "Mont2tlNi01 AIA IC matches exported SysEq");
  check_near(node.pCNode()->Gs, ref_Gs, 1e-3,
             "Mont2tlNi01 AIA Gs matches exported SysEq");
  for (long int ja = 0; ja < 22; ++ja) {
    check_near(node.pCNode()->xDC[31 + ja], ref_xDC[ja], 5e-5,
               "Mont2tlNi01 AIA xDC[" + std::to_string(31 + ja) + "]");
  }

  const double aia_pH = node.pCNode()->pH;
  node.pCNode()->NodeStatusCH = NEED_GEM_SIA;
  const long int run_rc_sia = node.GEM_run(false);
  check_cond(run_rc_sia == OK_GEM_SIA || run_rc_sia == OK_GEM_AIA,
             "Mont2tlNi01 recomputed SIA converged");
  check_near(node.pCNode()->pH, aia_pH, 1e-3,
             "Mont2tlNi01 SIA pH reproduces AIA equilibrium");

  long int ka = -1, k = -1;
  test_mont2tlm_indices(node, ka, k);
  test_mont2tlm_outputs(node, ka, k);
  test_mont2tlm_inputs(node, ka, k);

  TempDir tmp("gems3k_sorption_api_test_claysorp");
  test_copy_and_binary_roundtrip(node, tmp.path(), ka, 1, 1);

  const double psiA1_before = node.Get_psi_A(ka, 1);
  const double psiZ1_before = node.Get_psi_Z(ka, 1);
  const double hspz1_before = node.Get_h_spz(ka, 1);

  // Perturb C1 and fzpSpd on Mont2tlm 3TL surface type (ist=1) and recompute
  node.Set_C1(2.2, ka, 1);
  node.Set_fzpSpd(0.25, ka, 1);
  node.pCNode()->NodeStatusCH = NEED_GEM_AIA;
  const long int run_rc2 = node.GEM_run(false);
  check_cond(run_rc2 == OK_GEM_AIA || run_rc2 == OK_GEM_SIA,
             "Mont2tlNi01 GEM_run(false) after Set_C1/Set_fzpSpd converged");

  check_near(node.Get_C1(ka, 1), 2.2, kTol, "Mont2tlm persisted Set_C1(1)");
  check_near(node.Get_fzpSpd(ka, 1), 0.25, kTol,
             "Mont2tlm persisted Set_fzpSpd(1)");
  check_near(node.Get_h_spz(ka, 1), 0.5 * hspz1_before, 1e-4,
             "Mont2tlm 3TL h_spz halved when fzpSpd halved");
  check_cond(std::fabs(node.Get_psi_A(ka, 1) - psiA1_before) > 1e-5,
             "Mont2tlm 3TL psi_A changed after Set_C1");
  check_cond(std::fabs(node.Get_psi_Z(ka, 1) - psiZ1_before) > 1e-6 &&
                 std::fabs(node.Get_psi_Z(ka, 1)) >
                     std::fabs(psiZ1_before) * 1.1,
             "Mont2tlm 3TL psi_Z magnitude increased after halving fzpSpd");
  // Surface type 0 (NEM) remains non-electrostatic after 3TL perturbation
  check_near(node.Get_psi_A(ka, 0), 0.0, kTol, "Mont2tlm NEM psi_A == 0");
  check_near(node.Get_psi_Z(ka, 0), 0.0, kTol, "Mont2tlm NEM psi_Z == 0");
  check_near(node.Get_h_spz(ka, 0), 0.0, kTol, "Mont2tlm NEM h_spz == 0");

  std::cout << "sorption_api_test passed all checks on "
               "ClaySorp:G:Mont2tlNi01:0:0:1:25:0: "
            << "(pH=" << node.pCNode()->pH << ", IC=" << node.pCNode()->IC
            << ", 3TL psiA: " << psiA1_before << " -> "
            << node.Get_psi_A(ka, 1) << ", psiZ: " << psiZ1_before << " -> "
            << node.Get_psi_Z(ka, 1) << ")\n";
}

} // namespace

int main() {
  try {
    run_gibbsite_test();
    run_claysorp_test();
  } catch (const std::exception &ex) {
    std::cerr << "sorption_api_test FAILED: " << ex.what() << '\n';
    return 1;
  }
  return 0;
}
