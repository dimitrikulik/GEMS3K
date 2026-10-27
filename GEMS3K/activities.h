//-------------------------------------------------------------------
// $Id: activities.h 879 2013-10-10 14:47:33Z kulik $
/// \file activities.h
/// Contains definition of the ACTIVITY structure - module for computing
/// activities for GEM IPM 3 and new speciation algorithms.
//
/// \struct ACTIVITY activity.h
/// Defines the structure of node-dependent data for
/// exchange of activity terms between the node instance and the GEM algorithms.
/// Used in TNode class.
//
// Copyright (c) 2014- by D.Kulik, A.Leal
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public TSolModLicense for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------
//
#ifndef Activity_H_

#define Activity_H_

#include "databr.h"
#include "datach.h"
#include "s_solmod.h"
#include "s_sorpmod.h"

class TNode;

typedef struct /// ACTIVITY - data structure for computing DC activities in
               /// phases
{              /// DATACH indexation throughout, no I/O file exchange

  //      Usage of this variable node GEM-in GEM-out node
  double
      // \section Chemical scalar variables
      TK, ///< Node temperature T (Kelvin)                     	         + + - -
      P,  ///< Node Pressure P (Pa)                         	             +
          ///< +      -     -
      RT, // RT product
      ln5551, // ln(H2O_mol_to_kg) = log(55.50837344)
              //
      IC,   ///< Effective aqueous ionic strength (molal)                    - -
            ///< +     +
      Zeff, ///< Effective ion valence of aqueous electrolyte
      pH,   ///< pH of aqueous solution in the activity scale (-log10 molal) - -
            ///< +     +
      pe,   ///< pe of aqueous solution in the activity scale (-log10 molal) - -
            ///< +     +
      Eh;   ///< Eh of aqueous solution (V)                                  - -
            ///< +     +

  double Pc, Tc; // T in K, P in bar
  // Here TSolMod, TSorpMod, TKinMet parameters & coefficients moved from TMULTI
  // MULTI
  //   DK, AL September 2014
  long int
      // TSolMod stuff
      *LsMod, ///< Number of interaction parameters. Max parameter order (cols
              ///< in IPx), and number of coefficients per parameter in PMc
              ///< table [3*nPS]
      *LsMdc, ///<  for multi-site models: [3*nPS] - number of nonid. params per
              ///<  component;
              /// number of sublattices nS; number of moieties nM
      *LsMdc2, ///<  new: [3*nPS] - number of DQF coeffs; reciprocal coeffs per
               ///<  end member;
               /// reserved
      *IPx,    ///< Collected indexation table for interaction parameters of
               ///< non-ideal solutions
               ///< ->LsMod[k,0] x LsMod[k,1]   over nPS

      *LsPhl,      ///< new: Number of phase links; number of link parameters;
                   ///< [nPH][2]
      (*PhLin)[2]; ///< new: indexes of linked phases and link type codes (sum
                   ///< 2*LsPhl[k][0] over nPH)

  /* TSolMod !! arrays and counters to be added (for mixed-solvent electrolyte
  phase) TW

  ncsolv, /// TW new: number of solvent parameter coefficients (columns in solvc
  array) nsolv,  /// TW new: number of solvent interaction parameters (rows in
  solvc array) *ixsolv, /// new: array of indexes of solvent interaction
  parameters [nsolv*2] *solvc, /// TW new: array of solvent interaction
  parameters [ncsolv*nsolv]

  ncdiel, /// TW new: number of dielectric constant coefficients (colums in
  dielc array) ndiel,  /// TW new: number of dielectric constant parameters
  (rows in dielc array) *ixdiel /// new: array of indexes of dielectric
  interaction parameters [ndiel*2] *dielc, /// TW new: array of dielectric
  constant parameters [ncdiel*ndiel]

  ndh,    /// TW new: number of generic DH coefficients (rows in dhc array)
  *dhc,   /// TW new: array of generic DH parameters [ndh]
  */
  double denW[5], ///< Density of water, first T, second T, first P, second P
                  ///< derivative for Tc,Pc
      denWg[5],   ///< Density of steam for Tc,Pc
      epsW[5],    ///< Diel. constant of H2O(l)for Tc,Pc
      epsWg[5];   ///< Diel. constant of steam for Tc,Pc
  double
      *PMc, ///< Collected interaction parameter coefficients for the (built-in)
            ///< non-ideal mixing models -> LsMod[k,0] x LsMod[k,2]
      *DMc, ///< Non-ideality coefficients f(TPX) for DC -> L1[k] x LsMdc[k][0]
      *MoiSN, ///< End member moiety- site multiplicity number tables ->  L1[k]
              ///< x LsMdc[k][1] x LsMdc[k][2]
      *SitFr; ///< Tables of sublattice site fractions for moieties ->
              ///< LsMdc[k][1] x LsMdc[k][2]
  double *lPhc, ///< new: Collected array of phase link parameters
                ///< (sum(LsPhl[k][1] over Fi)
      *DQFc;    ///< new: Collected array of DQF parameters for DCs in phases ->
                ///< L1[k] x LsMdc2[k][0]
  //  *rcpc,  ///< new: Collected array of reciprocal parameters for DCs in
  //  phases -> L1[k] x LsMdc2[k][1]

  // TSorpMod stuff
  long int *LsESmo, ///< new: number of EIL model layers; EIL params per layer;
                    ///< CD coefs per DC; reserved  [nPS][4]
      *LsISmo;      ///< new: number of surface sites; isotherm coeffs per site;
                    ///< isotherm coeffs per DC; max.denticity of DC [nPS][4]
  long int (
      *XSMD)[NSURFSITES]; ///< SCM Modernization: denticity of surface species
                          ///< per surface site [Lads][NSURFSITES=8]
  long int (
      *SATC)[MNCAS]; ///< Setup of surface sites and species (will be applied
                     ///< separately within each sorption phase) [Lads]
#define xSMd XSMD
#define SATX SATC
  double *Aalp; ///< Full vector of specific surface areas of phases (m2/g)
                ///< [0:FI-1]
  double (*XfIEC)[MNST]
                 [NCIEC]; ///< SCM Modernization: Fixed permanent charge
                          ///< density in 1/nm2 and eq/kg [FIs][MNST=6][NCIEC=2]
                          /// link indexes to surface type [XL_ST]; sorbent em
                          /// [XL_EM]; surf.site [XL-SI] and EDL plane [XL_SP]
  double *EImc,           ///< new: Collected EIL model coefficients k -> +=
                          ///< LsESmo[k][0]*LsESmo[k][1]
      (*mCDc)[DMCDC],     ///< SCM Modernization: Array of charges of surface
                          ///< species [Lads][DMCDC=4]
      (*MASDT)[MNST][NMSDT], ///< SCM Modernization: Area fraction and patch
                             ///< isotherm parameters [FIs][MNST=6][NMSDT=4]
      (*MASDJ)[DNFCN],       ///< Parameters of surface species in surface
                             ///< complexation models [Lads][DNFCN]
      *IsoPc, ///< new: Collected isotherm coefficients per DC k -> +=
              ///< L1[k]*LsISmo[k][2];
      *IsoSc; ///< new: Collected isotherm coeffs per site k -> +=
              ///< LsISmo[k][0]*LsISmo[k][1];
  double (
      *XcapA)[MNST]; ///< Capacitance density of 0 EDL layer, F/m2 [FIs][FIat]
  double (
      *XcapB)[MNST]; ///< Capacitance density of B EDL layer, F/m2 [FIs][FIat]
  double (
      *XcapC)[MNST]; ///< Capacitance density of Ba EDL layer F/m2 [FIs][FIat]
  double (*XcapZ)[MNST]; ///< Location of slip plane for zeta potential
                         ///< (fraction of Debye length) [FIs][FIat]
  double (
      *XcapD)[MNST]; ///< Eff. cap. density of diffuse layer, F/m2 [FIs][FIat]
  // TSorpMod & TKinMet stuff
  double *
      SorMc; ///< new: Phase-related kinetics and sorption model parameters:
             ///< [nPS][16] in the same order as from Asur until fRes2 in TPhase

  //      Usage of this variable node GEM-in GEM-out node
  double *H0; ///< DC molar enthalpies, reserved [L]
  //     double *A0;     ///< DC molar Helmholtz energies, reserved [L]
  double *S0;  ///< DC molar entropies, reserved [L]
  double *Cp0; ///< DC molar heat capacity, reserved [L]

  double *tpp_G; ///< input Gibbs energy of species in J/mol [nDC]
  double *G0;    ///< Input normalized g0_j(T,P) for DC at unified standard
                 ///< scale[nDC]               +      +      -     -
  double *G; ///< Normalized DC energy function c(j), incl. primal mu mole/mole
             ///< [nDC]            -      -      +     +

  double *lnAct; ///< ln of DC activities (mu - mu0) [nDC] -      -      +     +
  double *lnGam; ///< ln of DC activity coefficients in unified (mole-fraction)
                 ///< scale [nDC]        -      -      +     +
  double *lnGmo; ///< Copy of lnGam from previous IPM iteration (reserved)
  // TSolMod stuff (detailed output on partial energies of mixing)   lnGam
  // components
  double *lnDQFt; ///< new: DQF terms adding to overall activity coefficients
                  ///< [nDCs]               -      -      +     +
  double *lnRcpt; ///< new: reciprocal terms adding to overall activity
                  ///< coefficients [nDCs]        -      -      +     +
  double *lnExet; ///< new: excess energy terms adding to overall activity
                  ///< coefficients [nDCs]     -      -      +     +
  double *lnCnft; ///< new: configurational terms adding to overall activity
                  ///< [nDCs]                -      -      +     +
  // Takeover from MULTI
  double *Gamma; ///< DC activity coefficients in molal or other phase-specific
                 ///< scale [0:L-1]
  double *lnGmf; ///< ln of initial DC activity coefficients for correcting G0
                 ///< [0:L-1]
  double *lnGmM; ///< ln of DC pure gas fugacity (or metastability) coeffs or
                 ///< DDF correction [0:L-1]

  // TSorpMod stuff (TBD)
  double *FWGTS;  ///< new: Masses of sorbent (real or virtual) for SCMS [FIs] -
                  ///< -      +     +
  double *lnSACT; ///< new: ln isotherm-specific SACT for surface species [nDCs]
                  ///< -      -      +     +
  double *lnGammF; ///< new: Frumkin or BET non-electrostatic activity
                   ///< coefficients [nDCs]        -      -      +     +
  double *CTerms;  ///< new: Coulombic terms (electrostatic activity
                   ///< coefficients) [nDCs]         -      -      +     +

  double (*VPh)[MIXPHPROPS], ///< Volume properties for mixed phases [FIs] - -
                             ///< +     +
      (*GPh)[MIXPHPROPS], ///< Gibbs energy properties for mixed phases [FIs] -
                          ///< -      +     +
      (*HPh)[MIXPHPROPS], ///< Enthalpy properties for mixed phases [FIs] - - +
                          ///< +
      (*SPh)[MIXPHPROPS], ///< Entropy properties for mixed phases [FIs] - - + +
      (*CPh)[MIXPHPROPS], ///< Heat capacity Cp properties for mixed phases
                          ///< [FIs]           -      -      +     +
      (*APh)[MIXPHPROPS], ///< Helmholtz energy properties for mixed phases
                          ///< [FIs]           -      -      +     +
      (*UPh)[MIXPHPROPS]; ///< Internal energy properties for mixed phases [FIs]
                          ///< -      -      +     +
                          // MIXPHPROPS: see m_const_base.h

  // Coding - takeover from MULTI
  char (*sMod)[8]; ///< new: Codes for built-in mixing models of multicomponent
                   ///< phases [FIs]
  char (
      *dcMod)[6]; ///< Codes for PT corrections for dependent component data [L]
  char (*SM)[MAXDCNAME];           ///< List of DC names in the system [L]
  char (*SF)[MAXPHNAME + MAXSYMB]; ///< List of phase names in the system [FI]
  char *ICC; ///< Classifier of IC { e o h a z v i <int> } [N]
  char *
      DCC; ///< Classifier of DC { TESKWL GVCHNI JMFD QPR <0-9>  AB  XYZ O } [L]
  char *PHC;  ///< Classifier of phases { a g f p m l x d h } [FI]
  char *DCCW; ///< internal DC class codes for calculation of primal chemical
              ///< potentials [L]
  char (*SCMC)[MNST][NSCMC]; ///< SCM Modernization: Classifier of EIL models
                             ///< and SACT codes applied to surface types
                             ///< [FIs][MNST=6][NSCMC=2]
  char *SATT;   ///< Classifier of applied SACT equations (isotherm corrections)
                ///< [Lads]
  double *MM;   //< DC molar masses, g/mol [L]
  double *FWGT; //< phase (carrier) masses, g                [0:FI-1]
  double *YOF; //< Surface free energy parameter for phases (J/g) (to accomodate
               // for variable phase composition) [FI]

  // Work part (to compute concentrations and chem.potentials), takeover from
  // MULTI
  long int N, ///< Number of ICs
      NR,     ///< NR - dimensions of R matrix
      L,      ///< Number of DCs
      Ls,     ///< Total number of DCs in phases-solutions
      Lads,   ///< Total number of DCs in sorption phases
      LO,     ///< LO -   index of water-solvent in DC list
      FI,     ///< Number of phases
      FIs,    ///< Number of phases-solutions
      FIa,    ///< Number of sorption phases (with EIL and SACT)
      E,      ///< PE - flag of electroneutrality constraint { 0 1 }
      K2,     // not needed?
      *L1;    ///< l_a vector - number of DCs included into each phase [nPH];
              // copy of nDCinPH

  double *A, ///< DC stoichiometry matrix A composed of a_ji [0:N-1][0:L-1]
      *Awt,  ///< IC atomic (molar) mass, g/mole [0:N-1]
      *XF,   ///< Output total number of moles of phases Xa[0:FI-1]
      *YF,   ///< Approximation of X_a in the next IPM iteration [0:FI-1]
      *XFA,  ///< Quantity of carrier in asymmetric phases Xwa, moles [FIs]
      *YFA,  ///< Approximation of XFA in the next IPM iteration [0:FIs-1]
      *X,    ///< DC quantities at eqstate x_j, moles - primal IPM solution [L]
      *Y,    ///< Copy of x_j from previous IPM iteration [0:L-1]
      *Fx,   ///< Dual DC chemical potentials defined via u_i and a_ji [L]
      *Wx,   ///< Mole fractions Wx of DC in multi-component phases [L]
      *F,    ///< Primal DC chemical potentials defined via g0_j, Wx_j and
             ///< lnGam_j[L]
      *F0,   ///< Excess Gibbs energies for (metastable) DC, mole/mole [L]
      *U,    ///< IC chemical potentials u_i (mole/mole) - dual IPM solution [N]
      *Falp; ///< Phase stability index (PC==2) [FI]

  double *IC_m, ///< Total IC molalities in aqueous phase (excl.solvent) [0:N-1]
      *IC_lm,   ///< log total IC molalities in aqueous phase [0:N-1]
      *IC_wm;   ///< Total dissolved IC concentrations in g/kg_soln [0:N-1]

  // old sorption models - EIL models (data for EIL activity coefficients)
  double (
      *XetaA)[MNST]; ///< Total EIL charge on A (0) EIL plane, mol [FIs][FIat]
  double (*XetaB)[MNST]; ///< Total charge on B (1) EIL plane, mol [FIs][FIat]
  double (*XetaC)[MNST]; ///< New Total charge on C (2, gamma) EIL plane,
                         ///< mol [FIs][FIat]
  double (*XetaD)[MNST]; ///< Total charge on D (3) EDL plane, mol [FIs][FIat]
  double (
      *XpsiA)[MNST]; ///< Relative potential at A (0) EDL plane,V [FIs][FIat]
  double (*XpsiB)[MNST]; ///< Relative potential at B (1, beta) EDL plane,V
                         ///< [FIs][FIat]
  double (*XpsiC)[MNST]; ///< New Relative potential at C (2, gamma) EDL plane,V
                         ///< [FIs][FIat]
  double (*XpsiD)[MNST]; ///< Relative potential at D (2) plane,V [FIs][FIat]
  double (*XpsiZ)[MNST]; ///< New Zeta potential at slip plane,V [FIs][FIat]
  double (*XdlA)[MNST];  ///< Effective thickness of 0-1 A EDL layer, nm
                         ///< [FIs][FIat], tbd
  double (*XdlB)[MNST];  ///< Effective thickness of 0-2 B EDL layer, nm
                         ///< [FIs][FIat], tbd
  double (*XdlC)[MNST];  ///< Effective thickness of 0-d C EDL layer, nm
                         ///< [FIs][FIat], tbd
  double (*XdlD)[MNST];  ///< Effective thickness of diffuse layer, nm
                         ///< [FIs][FIat], tbd
  double (*XFTS)[MNST]; ///< Total number of moles of surface DC at surface type
                        ///< [FIs][FIat]
  // Old sorption models
  double (*D)[MNST];  ///< Reserved; new work array for calc. surface act.coeff.
  double (*lnSAC)[4]; ///< former lnSAT ln surface activity coeff and Coulomb's
                      ///< term  [Lads][4]

  double *GamFs; ///< Copy of activity coefficients Gamma [L]
  double *fDQF,  ///< Increments to molar G0 values of DCs from pure gas
                 ///< fugacities or DQF terms, normalized [L]
      *MU,       ///< mu_j values of differences between dual and primal DC
                 ///< chem.potentials [L]
      *EMU,      ///< Exponents of DC increment to F_a criterion for phase [L]
      *NMU,      ///< DC increments to F_a criterion for phase [L]
      *Y_la,     ///< log activity of DC in phases (mju-mji0) [0:L-1]
      *Y_m,      ///< Molalities of aqueous species and sorbates [0:Ls-1]
      *Y_w,      ///< Mass concentrations of DC in multi-component
                 ///< phases,%(ppm)[0:Ls-1]
      *Pparc,    ///< Partial pressures or fugacities of pure DC, bar (Pc by
                 ///< default) [0:L-1]
      *EZ,       ///< Formula charge of DC in multi-component phases [0:Ls-1]
      *Vol,      ///< DC molar volumes, cm3/mol [L]
      *FVOL, ///< phase volumes, cm3 comment corrected DK 04.08.2009  [0:FI-1]
      *BF,  ///< Output bulk compositions of multicomponent phases bf_ai[FIs][N]
      *BFC; ///< Total output bulk composition of all solid phases [1][N]

  double      ///< Numerical tolerances (scaled), taken over from MULTI
      XwMinM, ///< Cutoff mole amount for elimination of water-solvent { 1e-13 }
      ScMinM, ///< Cutoff mole amount for elimination of solid sorbent { 1e-13 }
      DcMinM, ///< Cutoff mole amount for elimination of solution- or surface
              ///< species { 1e-30 }
      PhMinM, ///< Cutoff mole amount for elimination of non-electrolyte
              ///< condensed phase { 1e-23 } insertion values (re-scaled to
              ///< system size)
      //    DFYwM,  ///< Insertion mole amount for water-solvent { 1e-6 }
      //    DFYaqM, ///< Insertion mole amount for aqueous and surface species {
      //    1e-6 } DFYidM, ///< Insertion mole amount for ideal solution
      //    components { 1e-6 } DFYrM,  ///< Insertion mole amount for major
      //    solution components (incl. sorbent) { 1e-6 } DFYhM,  ///< Insertion
      //    mole amount for minor solution components { 1e-6 } DFYcM,  ///<
      //    Insertion mole amount for single-component phase { 1e-6 } DFYsM,
      //    ///< Insertion mole amount used in PhaseSelect() for a condensed
      //    phase component  { 1e-7 }
      SizeFactor, ///< factor for re-scaling the cutoffs/insertions to the
                  ///< system size
      //    TMols,      ///< Input total moles in b vector before rescaling
      //    SMols,      ///< Standart total moles (upscaled) {1000}
      //    MBX,       ///< Total mass of the system, kg
      //    FX,    	    ///< Current Gibbs potential of the system in IPM,
      // IC,         ///< Effective molal ionic strength of aqueous electrolyte
      // pH,         ///< pH of aqueous solution
      // pe,         ///<  pe of aqueous solution Eh,         ///< Eh of aqueous
      // solution, V
      DHBM,    ///< balance (relative) precision criterion
      DSM,     ///< min amount of phase DS
      GWAT,    ///< used in ipm_gamma()
      YMET,    ///< reserved
      IEPS;    // SACT control  from pa_p->IEPS
  long int PC; //  Control from pa_p->PC
  double W1,   //   Control from pa_p->W1
      DNS,     // Control from pa_p->DNS
      //    PCI,        ///< Current value of Dikin criterion of IPM convergence
      //    DK>=DX DXM,        ///< IPM convergence criterion threshold DX
      //    (1e-5) lnP,        ///< log Ptotal RT,         ///< RT: 8.31451*T
      //    (J/mole/K)
      FRT, ///< F/RT, F - Faraday constant = 96485.309
      Yw,  ///< Current number of moles of solvent in aqueous phase
      // ln5551,     ///< ln(55.50837344)
      ICmin,     ///< tolerance for minimum ionic strength to calculate aqueous
                 ///< activity models
      aqsTail,   ///< v_j asymmetry correction factor for aqueous species
      lowPosNum, ///< Minimum mole amount considered in GEM calculations
                 ///< (MinPhysAmount = 1.66e-24)
      logXw,     ///< work variable
      logYFk,    ///< work variable
      YFk;       ///< Current number of moles in a multicomponent phase
} ACTIVITY;

// class TNode;

#ifndef EIL_MODEL_RESULT_DEFINED
#define EIL_MODEL_RESULT_DEFINED
// Added for SCM modernization Oct 18 2026  DK ClaySorEl project
struct EILModelResult {
  char modelcode; // Model code '4', '3', 'T', 'B', 'C', 'D', 'N'
  std::string name;
  double sigma0, sigma1, sigma2, sigmad,
      sigmanet; // Plane charges (C/m²): '0', '1', '2', 'd', 'net'
  double psi0, psi1, psi2, psid; // Plane Potentials (V): '0', '1', '2', 'd'
  double zeta;                   // Electrokinetic Zeta potential (V)
  double cap_d;                  // Capacitance of diffuse layer (F/m2)?
  double sld01, sld02, sld0d,
      dlth; // distances from 0 plane to 1, 2, d, slip plane, and
            // diffuse layer thickness counted from 0 plane
};
#endif

class TActivity {
  ACTIVITY act;

  DATACH *csd;
  //    MULTI* mp;
  DATABR *cnd;
  // Current pointer to TNode class instance
  TNode *cno;

  size_t sizeFIs;     ///< current size of phSolMod
  TSolMod **phSolMod; ///< size current FIs - number of multicomponent phases

  // new - allocation of TsorpMod
  size_t sizeFIa;       ///< current size of phSorpMod
  TSorpMod **phSorpMod; ///< size current FIa - number of adsorption phases

  bool load;
  long int IT,
      IIM; // current number of iterations of GEM algorithm (used in smoothing)
  double AG, DGC, // smoothing
      FitVar[5];  ///< Internal. FitVar[0] is total mass (g) of solids in the
                  ///< system (sum over the BFC array)
                  ///<      FitVar[1], [2] reserved
                  ///<       FitVar[4] is the AG smoothing parameter;
                  ///<       FitVar[3] is the actual smoothing coefficient

protected:
  // Internal arrays for the performance optimization  (since version 2.0.0) <-
  // MULTI
  long int sizeN; /*, sizeL, sizeAN;*/
  double *AA;
  double *BB;
  long int *arrL;
  long int *arrAN;

  void set_def();
  void Alloc_TSolMod(size_t newFIs);
  void Free_TSolMod();
  void Alloc_TSorpMod(size_t newFIs);
  void Free_TSorpMod();

  void Alloc_A_B(long int newN);
  void Free_A_B();
  void Build_compressed_xAN();
  void Free_compressed_xAN();
  void Free_internal();

public:
  void Alloc_internal();

  // Explicit constructor
  TActivity(DATACH *csd, DATABR *sbc, TNode *ptn);
  ~TActivity();

  // Returns pointer to ACTIVITY work data structure
  ACTIVITY *GetActivityDataPtr(void) { return &this->act; }

  long int CheckMassBalanceResiduals(double *Y);
  double ConvertGj_toUniformStandardState(double g0, long int j, long int k);
  double PhaseSpecificGamma(long int j, long int jb, long int je, long int k,
                            long int DirFlag = 0L); // Added 26.06.08
  void SetSmoothingFactor(long int mode); // new smoothing function (3 variants)

  // ipm_chemical.cpp
  //        void XmaxSAT_IPM2();
  //    void XmaxSAT_IPM2_reset();
  double DC_DualChemicalPotential(double U[], double AL[], long int N,
                                  long int j);
  //        void Set_DC_limits( long int Mode );
  void TotalPhasesAmounts(double X[], double XF[], double XFA[]);
  double DC_PrimalChemicalPotentialUpdate(long int j, long int k);
  double DC_PrimalChemicalPotential(double G, double logY, double logYF,
                                    double asTail, double logYw, char DCCW);
  void PrimalChemicalPotentials(double F[], double Y[], double YF[],
                                double YFA[]);
  double KarpovCriterionDC(double *dNuG, double logYF, double asTail,
                           double logYw, double Wx, char DCCW);
  void KarpovsPhaseStabilityCriteria();
  void StabilityIndexes(); // added 01.05.2010 DK
  double DC_GibbsEnergyContribution(double G, double x, double logXF,
                                    double logXw, char DCCW);
  double GX(double LM);

  void ConvertDCC();
  long int getXvolume();

  // ipm_chemical2.cpp
  void phase_bcs(long int N, long int M, long int jb, double *A, double X[],
                 double BF[]);
  void phase_bfc(long int k, long int jj);
  double bfc_mass(void);
  void CalculateConcentrationsInPhase(double X[], double XF[], double XFA[],
                                      double Factor, double MMC, double Dsur,
                                      long int jb, long int je, long int k);
  void CalculateConcentrations(double X[], double XF[], double XFA[]);
  // Added on 18.10.2026 DK  ClaySorEl project
  void IS_EtaCalc(); // Calculation of surface charge densities on surf. patches
  double virtualSorbentMass(long int jb, long int je, long int N);
  double calc_water_permittivity(double T_kelvin = 298.15,
                                 double P_bar = 1.0) const;
  double calc_debye_length(double I) const;
  double calc_cap_d_from_sigma_d(double sigma_d, double I) const;
  double calc_psi_d_from_sigma_d(double sigma_d, double I, double z_eff) const;
  double calc_gc_potential_at_distance(double psi_d, double I, double x_dist,
                                       double z_eff = 1.0) const;
  // Individual EIL Model Solvers
  EILModelResult compute_flm(char ModCode, double I, double z_eff, double C1,
                             double C2, double C3, double sig0, double sig1,
                             double sig2, double x_slip, double debye_m) const;
  EILModelResult compute_tlm(char ModCode, double I, double z_eff, double C1,
                             double C2, double sig0, double sig1, double x_slip,
                             double debye_m) const;
  EILModelResult compute_cd_music(char ModCode, double I, double z_eff,
                                  double C1, double C2, double sig0,
                                  double sig1, double x_slip,
                                  double debye_m) const;
  EILModelResult compute_bsm(char ModCode, double I, double z_eff, double C1,
                             double sig0, double x_slip, double debye_m) const;
  EILModelResult compute_ddl(char ModCode, double I, double z_eff, double sig0,
                             double x_slip) const;
  EILModelResult compute_ccm(char ModCode, double C1, double sig0) const;
  EILModelResult compute_nem(char ModCode, double sig0) const;
  EILModelResult compute_tce(char ModCode, double I, double z_eff, double C1,
                             double C2, double sig0, double sig1, double par[],
                             double x_slip, double debye_m) const;
  // Calculation of any of 6 EIL models and NEM
  long int GouyChapman(long int jb, long int je, long int k);
  //  Surface activity coefficient terms
  long int SurfaceActivityCoeff(long int jb, long int je, long int jpb,
                                long int jdb, long int k);
  long int denticity_total(long int j);
  long int denticity_num_site_types(long int j);
  //        void pm_GC_ods_link( long int k, long int jb, long int jpb, long int
  //        jdb, long int ipb ); double SmoothingFactor( ); void
  //        SetSmoothingFactor( long int mode ); // new smoothing function (3
  //        variants)
  // Main call for calculation of activity coefficients on IPM iterations
  long int CalculateActivityCoefficients(long int LinkMode);
  // Built-in activity coefficient models
  // Generic solution model calls
  void SolModCreate(
      long int jb, long int jmb, long int jsb, long int jpb, long int jdb,
      long int k, long int ipb, char ModCode, char MixCode,
      /* long int jphl, long int jlphc, */ long int jdqfc /*, long int jrcpc*/);
  void SolModParPT(long int k, char ModCode);
  void SolModActCoeff(long int k, char ModCode);
  void SolModExcessProp(long int k, char ModCode);
  void SolModIdealProp(/*long int jb,*/ long int k, char ModCode);
  void SolModStandProp(/*long int jb,*/ long int k /*, char ModCode*/);
  void SolModDarkenProp(/*long int jb,*/ long int k /*, char ModCode*/);

  // Generic access methods
  void setTemperature(double T);      // set temperature (in units of K)
  void setPressure(double P);         // set pressure (in units of Pa)
  void updateStandardGibbsEnergies(); // compute standard Gibbs energies (so far
                                      // only P,T interpolation)
  void updateStandardVolumes();
  void updateStandardEnthalpies();
  void updateStandardEntropies();
  void updateStandardHeatCapacities();
  void updateThermoData(bool norm = true);

  void setSpeciation(const double *n); // set speciation (in units of moles)
  void updateConcentrations();         // compute concentrations in all phases
  void updateActivityCoefficients();   // compute activity coefficients
  void updateChemicalPotentials();     // compute primal chemical potentials
  void updateActivities();             // compute primal activities
  void updateChemicalData();
};

#endif
// -----------------------------------------------------------------------------
// end of Activity_h
