//-------------------------------------------------------------------
// $Id$
//
/// \file m_const_base.h
/// Declarations of enums and constants from GEM-Selektor code that
/// are used in GEMS3K standalone code (this file is not used otherwise).
//
// Copyright (c) 1995-2012 S.Dmytriyeva, D.Kulik
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
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------
//

#ifndef M_CONST_BASE_H
#define M_CONST_BASE_H

const char S_OFF = '-', S_ON = '+', S_REM = '*', A_NUL = '?';

const long int MNST = 4, // fixed number of surface types (patches)
    DMCDC = 4, // new number of columns in mCDc (DC EIL plane charges) table
    DNFCN = 4, // number of columns in MASDJ table
    MNCAS = 4, // number of columns in SATC table
    NCAZP = 4, // New: number of capacitance & zeta parameters per surface type
    NMSDT = 4, // New: number of columns in MSDT table
    NSCMC = 2, // New: number of columns in SCMC table
    NFSIT = 4, // New: number of columns in FsiT table
    NCIEC = 2, // Number of columns in XfIEC input table
    NSURFSITES = 8; // New: maximum number of surface sites
                    // to consider in sorption phase

const unsigned long int MAXICNAME = 6, MAXSYMB = 4, TDBVERSION = 64;

const long int MAXDCNAME = 16, MAXPHNAME = 16, EQ_RKLEN = 58;

const int MPP_TOT = 0, // index of column with total mixed phase property
    MPP_STD = 1, // index of column with standard property sum for mixed phases
    MPP_RCP = 2, // index of column with reciprocal (Darken) property sum for
                 // mixed phases
    MPP_IDL = 3, // index of column with ideal mixing property for the phases
    MPP_EXS = 4, // index of column with excess mixing property for the phases
    MIXPHPROPS = 5; // Number of columns in the property table for mixed phases

typedef enum { // Limits on DC and phases
  // type of lmits
  NO_LIM = 'O',
  LOWER_LIM = 'L',
  UPPER_LIM = 'U',
  BOTH_LIM = 'B',
  // mode recalc of limits Set_DC_Limits()
  DC_LIM_INIT = 0,
  DC_LIM_CURRENT = 1
} DC_LIMITS;

// Work DC classifier codes  pm->DCCW
enum SolDCodes {

  DC_SINGLE = 'U',       // This DC is a single-component phase
  DC_SYMMETRIC = 'I',    // This DC is in symmetric solution phase
  DC_ASYM_SPECIES = 'S', // This is DC-solute(sorbate) in asymmetric phase
  DC_ASYM_CARRIER = 'W'  // This is carrier(solvent) DC in asymmetric phase
};

enum solmod_switches { // indexes of keys of phase (solution, sorption, kinetic)
                       // models
  SPHAS_TYP,
  DCOMP_DEP,
  SPHAS_DEP,
  SGM_MODE,
  DCE_LINK,
  MIX_TYP,

  // Link state of CalculateActivityCoefficients()
  LINK_UX_MODE,
  LINK_TP_MODE,
  LINK_PP_MODE, // LINK_PHP_MODE,
  LINK_IN_MODE, // Initialization mode for kinetics and other time-dependent
                // processes
  SORP_MOD,     // new, see also enum sorption_control
  KINR_MOD,     /// see also enum kinmet_controls
  // Posible modes of calculation of activity coefficients (private, public)
  SM_UNDEF = 'N',
  SM_TPDEP = 'T',
  SM_UXDEP = 'X',
  SM_PRIVATE_ = 'P',
  SM_PUBLIC = 'U',

  // Posible modes of calculation of activity coefficients (built-in or scripted
  // models)
  SM_STNGAM = 'S', // Built-in function for activity coefficients
  SM_NOSTGAM = 'N',

  // Codes to identify the mixing models used (during IPM iterations)
  SM_IDEAL = 'I',  // ideal solution or single-component phase
  SM_BERMAN = 'B', // built-in multicomponent multisite (a)symmetric
                   // solid-solution model
  SM_CEF = '$',    //     built-in multicomponent multisite solid-solution model
                   //     (CALPHAD)
  SM_MBW =
      '#', //     built-in Modified Bragg-Williams model, [Vinograd et al. 2018]
  SM_REDKIS =
      'G', // built-in binary Guggenheim (Redlich-Kister) solid-solution model
  SM_MARGB = 'M',   // built-in binary Margules solid-solutions (subregular)
  SM_MARGT = 'T',   // built-in ternary Margules solid-solution (regular)
  SM_VANLAAR = 'V', // built-in multi-component Van Laar solid-solution model
  SM_GUGGENM = 'K', // built-in multi-component Guggenheim solid-solution model
  SM_REGULAR = 'R', // built-in multi-component Regular solid-solution model
  SM_NRTLLIQ = 'L', // built-in multi-component NRTL model for liquid solutions
  SM_WILSLIQ =
      'W', // built-in multi-component Wilson model for liquid solutions
  SM_CGFLUID =
      'F', // built-in multi-component Churakov-Gottschalk (CG) fluid EoS model
  SM_PRFLUID =
      'P', // built-in Peng-Robinson-Stryjek-Vera (PRSV) fluid EoS model
  SM_PCFLUID = '5', // built-in perturbed-chain statistical-association (PCSAFT)
                    // fluid EoS model (reserved)
  SM_STFLUID = '6', // built-in Sterner-Pitzer (STP) fluid EoS model
  SM_PR78FL = '7',  // built-in Peng-Robinson (PR78) fluid EoS model
  SM_CORKFL = '8',  // built-in compensated Redlich-Kwong (CORK) fluid EoS model
  SM_REFLUID = '9', // built-in reference EoS fluid model (reserved)
  SM_SRFLUID = 'E', // built-in Soave-Redlich-Kwong (SRK) fluid EoS model
  SM_AQDAV = 'D',   // built-in Davies model (with 0.3) for aqueous electrolytes
  SM_AQDH1 =
      '1', // built-in Debye-Hueckel limiting law for aqueous electrolytes
  SM_AQDH2 =
      '2', // built-in 2-term Debye-Hueckel model for aqueous electrolytes
  SM_AQDH3 = '3',  // built-in 3-term Debye-Hueckel model for aqueous
                   // electrolytes (Karpov version)
  SM_AQDHH = 'H',  // built-in 3-term Debye-Hueckel model for aqueous
                   // electrolytes (Helgeson version)
  SM_AQDHS = 'Y',  // built-in 3-term Debye-Hueckel model for aqueous
                   // electrolytes (Shvarov version)
  SM_AQSIT = 'S',  // built-in SIT model for aqueous electrolytes
  SM_AQEXUQ = 'Q', // built-in extended UNIQUAC model for aqueous electrolytes
  SM_AQPITZ = 'Z', // built-in Pitzer HMW model for aqueous electrolytes
  SM_AQMIX =
      'C', // built-in mixed-solvent aqueous Debye-Hueckel model (reserved)
  SM_AQELVIS = 'J', // built-in modified extended UNIQUAC model (ELVIS) for
                    // aqueous electrolyte
  SM_DONNAN = 'X',  // ion exchange (Donnan volume model) (reserved)
  SM_SURCOM = 'A',  // models of surface complexation at solid-aqueous interface
  SM_USERDEF = 'U', // user-defined mixing model (scripts in Phase record)
  SM_OTHER = 'O', // other built-in phase-specific models of non-ideal solutions
                  // (selected by phase name)

  // Codes to identify specific mixing rules and temperature functions in EoS
  // and activity models - set in ph->sol_t[MIX_TYP] (ph->sol_t[6])
  MR_UNDEF = 'N', // Default mixing rule or form of interaction parameter
                  // coefficients; NEM for adsorption 'A'
  MR_WAAL = 'W',  // Basic Van der Waals mixing rules in cubic EoS models
  MR_CONST =
      'C', // Constant one-term interaction parameter kij; CCM for sorption 'A'
  MR_TEMP = 'T', // Temperature-dependent one-term interaction parameter kij
                 // (Jaubert et al. 2005); TLM for adsorption 'A'
  MR_LJ = 'J',   // Lemmon-Jacobsen mixing rule (Lemmon and Jacobsen, 1999)
  MR_KW1 = 'K',  // Kunz-Wagner mixing rule (Kunz and Wagner, 2007)
  MR_PITZ5 =
      '5', // 5-term Pitzer model temperature dependence (TOUGHREACT variant)
  MR_PITZ6 =
      '6', // 6-term Pitzer model temperature dependence (FREZCHEM variant)
  MR_PITZ8 = '8',  // 8-term Pitzer model temperature dependence
  MR_B_RCPT = 'R', // Use CEF reciprocal non-ideality terms in Berman multi-site
                   // ss model
  MR_SORP_SMB = 'B', // new for SM_SURCOM  ='A': this sorption model is based on
                     // site mole balances for all surface types and species,
                     // disables SACTs but can use EIL models; has implicit
                     // sorbent described by its current mass and SSA
  MR_SORP_ELEM = 'E' // this sorption model has only surface species in
  // elemental stoichiometry on all surface types, with sorbent endmembers and
  // SSA given for the whole sorbent, can use both SACTs and EIL models

  /* new the type of EIL model is now shifted to ph->SCMC[i][0] per surface type
  MR_A_DLM = 'D',  // Diffuse-layer electrostatic model (DDLM) for adsorption
                   // SM_SURCOM = 'A'
  MR_A_BSM = 'B',  // Basic Stern electrostatic model (BSM) for adsorption 'A'
  MR_A_CDLM =
      'M', // CD-MUSIC (3-layer) electrostatic model (3LM) for adsorption 'A'
  MR_A_ETLM = 'E' // Extended TLM electrostatic model (ETLM) for adsorption 'A'
  */
};

typedef enum { // classes of independent components IC, used in ccIC/ICC code
               // list
  IC_ELEMENT = 'e',   // chemical element (except oxygen and hydrogen)
  IC_OXYGEN = 'o',    // oxygen
  IC_HYDROGEN = 'h',  // hydrogen (natural mixture of isotopes) H
  IC_PROTIUM = 'p',   // protium (reserved) Hp
  IC_DEYTERIUM = 'd', // deuterium (reserved) D
  IC_TRITIUM = 't',   // tritium (reserved) T
  IC_FORMULA = 'f',   // formula unit (eg. for Sio - a symbol of SiO2)
  IC_METALION = 'm',  // metal ion (cation), reserved
  IC_LIGAND = 'l',    // ligand (anion), reserved
  IC_ADDIT = 'a',     // IC with unknown stoichiometry (eg: Hum, humic acid)
  IC_ISOTOPE = 'i',   // isotope of chemical element from 1 to 250
  IC_OXYGEN16 = 'q',  // q: oxygen 16O (reserved)
  IC_OXYGEN18 = 'r',  // r: oxygen 18O (reserved)
  IC_CHARGE = 'z',    // z: electrical charge
  IC_VOLUME = 'v',    // volume
  IC_SSITE = 's'      // sorption site for site balance constraint (ClaySor)
} IC_CLASSES;

typedef enum { // Classifications of DC
  // Type of input data for
  SRC_DCOMP = 'd',  // the key points to existing PDB record in DCOMP chain
  SRC_REACDC = 'r', // the key points to existing PDB record in REACDC chain
  SRC_NEWDC = 'n',  // the key new reaction-defined component
  SRC_NEWISO = 'i', // the same as n, but this component is an isotopic form
  SRC_FICT = 'f',   // fictive species
  // Aqueous electrolyte phase
  DC_AQ_PROTON = 'T',   // hydrogen ion H+
  DC_AQ_ELECTRON = 'E', // electron (as a DC)
  DC_AQ_SPECIES = 'S',  // other aqueous species (ions, complexes and ion pairs)
  DC_AQ_SURCOMP = 'K',  // Surface complex represented as aqueous species
  DC_AQ_SOLVENT = 'W',  // water H2O (major solvent)
  DC_AQ_SOLVCOM = 'L',  // other components of a solvent (eg. alcohol)
  // Gas phase ( G code can be used for all gases; V,C,H,N codes are reserved)
  DC_GAS_COMP = 'G', // other gases
  DC_GAS_H2O = 'V',  // H2O steam
  DC_GAS_CO2 = 'C',  // CO2 (carbon dioxide)
  DC_GAS_H2 = 'H',   // H2 hydrogen
  DC_GAS_N2 = 'N',   // N2 nitrogen
  // Solid/liquid non-electrolyte multicomponent phases
  DC_SOL_IDEAL = 'I', // end-member of ideal solution
  DC_SOL_MINOR =
      'J', // junior independent end member (for initial approximation)
  DC_SOL_MAJOR =
      'M', // major independent end member (for initial approximation)
  DC_SOL_MINDEP =
      'F', // junior dependent end member (for initial approximation)
  DC_SOL_MAJDEP = 'D', // major dependent end member (for initial approximation)
  // Sorption phases and poly(oligo)electrolytes
  DC_SUR_CARRIER = 'Q', // Principal end-member of solid carrier
  DC_SUR_MINAL = 'P',   // Minor end-member of solid carrier
  DC_PEL_CARRIER = 'R', // Carrier of poly(oligo)electrolyte
  //  DC_SSC_A0 = '0',      // Strong surface complex on site type 0 (A plane)
  //  DC_SSC_A1 = '2',      // Strong surface complex on site type 2 (A plane)
  //  DC_SSC_A2 = '4',      // Strong surface complex on site type 4 (A plane)
  //  DC_SSC_A3 = '6',      // Strong surface complex on site type 6 (A plane)
  //  DC_SSC_A4 = '8',      // Strong surface complex on site type 8 (A plane)
  //  DC_WSC_A0 = '1',      // Weak surface complex on site type 1 (B plane)
  //  DC_WSC_A1 = '3',      // Weak surface complex on site type 3 (B plane)
  //  DC_WSC_A2 = '5',      // Weak surface complex on site type 5 (B plane)
  //  DC_WSC_A3 = '7',      // Weak surface complex on site type 7 (B plane)
  //  DC_WSC_A4 = '9',      // Weak surface complex on site type 9 (B plane)
  DC_IESC_A = 'A', // Strong exchange ion const-charge A plane
  DC_IEWC_B = 'B', // Weak exchange ion const-charge B or C plane
  // Aliases for multi-site model   Revision 17.10.2026  DK
  DC_SUR_GROUP = '@',   // Surface amphoteric group on A plane -> '0'  was 'X'
  DC_SUR_COMPLEX = 'X', // Strong surface complex on only A plane -> '0' was 'Y'
  DC_SUR_IPAIR = 'Y',   // Weak sur complex on A, B or C plane -> '1' was 'Z'
  DC_SUR_SPECIES = 'Z', // Generic CD-MUSIC surface species (surf.complex)

  // Single-component phases
  DC_SCP_CONDEN = 'O', // DC forming a single-component phase

  // New SCMs with site mole balances (added 16.11.2017 by DK)
  DC_SCM_SPECIES = 'U'

} DC_CLASSES;

#define SORPTION_DC "0246813579ABXYZPQ"

// This code defines standard state and reference scale of concentrations
// for components of this phase. It is used by many subroutines
// during calculations of equilibrium states
enum PH_CLASSES {  // Possible values
  PH_AQUEL = 'a',  // aqueous electrolyte
  PH_GASMIX = 'g', // mixture of gases
  PH_FLUID = 'f',  // fluid phase
  PH_PLASMA = 'p', // plasma
  PH_SINCOND =
      's', // condenced solid phase, also multicomponent (solid solution)
  PH_SINDIS = 'd',   // dispersed solid phase, also multicomponent
  PH_HCARBL = 'h',   // mixture of condensed hydrocarbons
  PH_LIQUID = 'l',   // non-electrolyte liquid (melt)
  PH_SIMELT = 'm',   // silicate (magmatic) melt or non-aqueous electrolyte
  PH_IONEX = 'i',    // ion exchange on permanent charge with site balances e.g.
                     // B&B Clay, no SACT and no EIL terms
  PH_ADSORPT = 'q',  // surface complexation with surface site balances e.g. B&B
                     // ClaySor, formerly 'z', no SACT and no EIL terms
  PH_SORPTION = 'x', // Sorption phase (sorbent+sorbates) without site balances
                     // and with SACT and EIL terms
  PH_POLYEL = 'y',   // colloidal poly- (oligo)electrolyte, membrane etc. TBD
  PH_ADSSBEIL = 'z'  // new: surface complexation with site balances & EIL
                     // models e.g. ClaySorEl, no SACT terms
};

enum sorption_control {
  // EDL interface models - separate for site(patch) types in v. 3.1
  // to be set in ph->SCMC list for each surface type
  SC_NEM = 'N', // Non-electrostatic model, charge densities on zero plane
                // calculated, but potentials and Coulomb terms are neglected
  // Main EIL models (revision 2026, cf. [Hiemstra & Lutzenkirchen 2025])
  SC_4LM = '4', // Four-layer model
  SC_3LM = '3', // Three-layer EIL model for CD-MUSIC
  SC_TLM = 'T', // Triple-layer model
  SC_BSM = 'S', // Basic Stern model
  SC_DDL = 'D', // Double-layer model
  SC_CCM = 'C', // Constant-capacitance model
  SC_TCE = 'E', // Tournassat EIL for clay edge zero plane potential
  // Auxiliary models (experimental), now disabled since 18.10.2026
  SC_MTL = 'M',
  SC_MXC = 'E',
  SC_NNE = 'X',
  SC_IEV = 'I',
  SC_ISCCM = 'P',
  SC_UNUSED = 'U',   // same as NEM but if set on all surface types, the EIL
                     // structures are not used at all
  SC_NOT_USED = 'N', // Compatibility alias for SC_NEM

  // New - column index in the EIL capacitance and zeta parameter table
  EIL_CAPDA = 0, // C1 (CA) capacitance density parameter (Stern layer models)
                 // default 1.0
  EIL_CAPDB = 1, // C2 (CB) parameter (4LM, 3LM, TLM) default 0.2
  EIL_CAPDC = 2, // C3 (CC) parameter (4LM) default 0.2
  EIL_SLPLF = 3, // fraction of Debye length to locate the slip plane for zeta
                 // potential, default 0.5 for hydrophilic and 0.17 for
                 // hydrophobic surfaces, range -0.2 to 2.0

  // Methods of Surface Activity Terms calculation
  SAT_COMP = 'C',
  SAT_NCOMP = 'N',
  SAT_SOLV = 'S',
  SAT_INDEF = 'I',

  // New methods for surface activity coefficient terms (2004)
  SACT_INDEF = 'I',
  SACT_LAN_NCOMP = 'L',
  SACT_LAN_COMP = 'C',
  SACT_QCA_NCOMP = 'Q',
  SACT_QCA1_NCOMP = '1',
  SACT_QCA2_NCOMP = '2',
  SACT_QCA3_NCOMP = '3',
  SACT_QCA4_NCOMP = '4',
  SACT_FFG_NCOMP = 'F',
  SACT_FFG_COMP = 'R',
  SACT_PIV_NCOMP = 'P',
  SACT_FREU_NCOMP = 'f',
  SACT_BET_NCOMP = 'B',
  SACT_VIR_NCOMP = 'W',
  // TBD: codes for surface type isotherms with species bound to different site
  // types

  // Assignment of surtype to carrier (end-member)
  CCA_VOL = 'V',
  CCA_0 = '0',
  CCA_1,
  CCA_2,
  CCA_3,
  CCA_4,
  CCA_5,
  CCA_6,
  CCA_7,
  CCA_8,
  CCA_9,
  SPL_0 = '0',
  SPL_1,
  SPL_2,
  SPL_3,
  SPL_B = 'b',
  SPL_D = 'd',
  SPL_C = 'c',
  SDU_N = 'n',
  SDU_m = 'm',
  SDU_M = 'M',
  SDU_g = 'g',
  CST_0 = '0',
  CST_1,
  CST_2,
  CST_3,
  CST_4,
  CST_5, // surface type index
  CSI_0 = '0',
  CSI_1,
  CSI_2,
  CSI_3,
  CSI_4,
  CSI_5, // surface site index

  // New: Column index in SCMC array (rows for surface types/patches)
  SCMC_EILM_CODE = 0,
  SCMC_SACT_CODE = 1,
  // New: Column index in MASDT array (rows for surface types/patches)
  MSDT_FSUTSA = 0, // Fraction of this surface type rel. to sorbent surface area
  MSDT_TCE_PA1 = 1, // Tournassat parameter A1 for clay edge
  MSDT_TCE_PA2 = 2, // Tournassat parameter A2 for clay edge
  MSDT_TCE_PA3 = 3, // Tournassat parameter A3 for clay edge
  // NEW: Column index in CapT array (capacitances F/m2 of EIL layers)
  CAPT_A_CAPD = 0, // Capacitance density of A layer (0-1 plane)
  CAPT_B_CAPD = 1, // Capacitance density of B layer (1-2 or 1-d plane)
  CAPT_C_CAPD = 2, // Capacitance density of C layer (2-d plane, 4LM only)
  CAPT_F_ZETA = 3, // fraction of Debye length locating slip plane for zeta
  // New: column index in FsiT array (rows for site types)
  FSIT_S_DENSITY = 0,  // Density of sites of this type (1/nm2)
  FSIT_S_CAPACITY = 1, // Site capacity in mol/(kg sorbent) e.g. in ClaySor

  // Number of parameters per surface species in the MaSdj array
  // MCAS = 4 = DNFCN ; position index    added by KD 25.10.2004,
  // mod. 11.10.2026
  COV_DENS = 0, // Max. full-site-coverage density of species 1/nm2
  ISOTH_PAR_1,  // SACT isotherm parameter 1 (FFG, BET)
  ISOTH_PAR_2,  // SACT isotherm parameter 2 (BET, Pivovarov)
  ISOTH_PAR_3,  // SACT isotherm parameter 3 (BET etc.)

  // Column indices of surface species allocation table SATC [MCAS]
  SURSP_ALLOC_SURFTYPE = 0,
  SURSP_INDIV_SACT_CODE,
  SURSP_ALLOC_SORBENT_EM,
  SURSP_COV_DENS_UNITS,

  // SCM Modernization: Column indices of mCDc table of charge distribution
  CHARGE_ON_A = 0,
  CHARGE_ON_B,
  CHARGE_ON_C,
  CHARGE_BAL,

  // SCM Modernization: Column indices of MaSdj table of adsorption parameters
  // [DNFCN=4]
  PI_DEN = 0, // Complete coverage density (1/nm2)
  PI_P1 = 1,  // SACT isotherm parameter 1 (Frumkin interaction alpha)
  PI_P2 = 2,  // SACT isotherm parameter 2 (q in QCA/BET/Pivovarov)
  PI_P3 = 3   // SACT isotherm parameter 3 (b in BET or S in Pivovarov)
};

enum ph_kinmet_controls { /// TKinMet: codes to control kinetic rate models

  KM_UNDEF = 'N', /// not defined, no account for
  KinProCode = 2,
  KM_PRO_MWR = 'M',  /// Kinetics of generic dissolution/precipitation (no
                     /// uptake, ionex, adsorption)
  KM_PRO_UPT = 'U',  /// Kinetics of uptake/entrapment (of minor/trace element)
                     /// into solid solution
  KM_PRO_IEX = 'X',  /// Kinetics of ion exchange (clays, C-S-H, zeolites, ...)
  KM_PRO_ADS = 'A',  /// Kinetics of adsorption (on MWI), redox
  KM_PRO_NUPR = 'P', /// Kinetics of nucleation and precipitation
  KinModCode = 3,
  KM_MOD_TST = 'T', /// Generic TST dissolution/precipitation model following
                    /// Shott ea 2012
  KM_MOD_PAL =
      'P', /// Dissolution/precipitation model of the form (Palandri 2004)
  KM_MOD_WOL = 'W',  /// Carbonate growth model following (Wolthers 2012)
  KM_MOD_NUGR = 'U', /// Mineral nucleation and growth model with
                     /// nuclei/particle size distr. (TBD)
  KinSorpCode = 4,
  KM_UPT_ENTRAP =
      'E',           ///	Unified entrapment model (Thien,Kulik,Curti 2012)
  KM_UPT_UPDP = 'M', ///	DePaolo (2011) uptake kinetics model
  KM_UPT_SEMO = 'G', ///  Growth (surface) entrapment model (Watson 2004)
  KM_IEX_FAST = 'F', ///  Fast ion exchange kinetics (e.g. montmorillonite, CSH)
  KM_IEX_SLOW = 'L', ///  Slow ion exchange kinetics (e.g. illite, zeolites)
  KM_ADS_INHIB = 'I', ///  Adsorption inhibition
  KM_NUCL_SSMP = 'P', ///  Solid solution nucleation model (Prieto 2013)
  KinLinkCode = 5,
  KM_LNK_SURF = 'S', ///   Link to (fraction of) solid substrate surface
  KM_LNK_PVOL = 'P', ///    Link to (fraction of) solid substrate (pore) volume
  KM_LNK_MASS = 'M', ///	Link to (fraction of) solid substrate mass
  KinSizedCode = 6,
  KM_SIZED_ETM = 'T', ///  Empirical f(time) cubic polynomial f = a + bt +ct^2 +
                      ///  dt^3 (default)
  KM_SIZED_ESI =
      'S', ///  Empirical f(lgSI) cubic polynomial f = a + bt +ct^2 + dt^3
  KM_SIZED_ESA = 'A', ///  Empirical f(sarea-change) cubic polynomial f = a + bt
                      ///  +ct^2 + dt^3
  KM_SIZED_EVOL = 'V', ///  Empirical f(volume-change) cubic polynomial f = a +
                       ///  bt +ct^2 + dt^3
  KM_SIZED_MASS = 'M', ///  Empirical f(mass-change) cubic polynomial f = a + bt
                       ///  +ct^2 + dt^3
  KM_SIZED_MOL = 'X',  ///  Empirical f(amount-change) cubic polynomial f = a +
                       ///  bt +ct^2 + dt^3
  KM_SIZED_UNI = 'U',  /// 	Uniform particle/pore size distribution
  KM_SIZED_BIN = 'B',  /// 	Binodal particle/pore size distribution
  KM_SIZED_FUN = 'F',  ///    Empirical distribution function
  KinResCode = 7,
  KM_RES_SURF_N = 'A', /// surface-scaled rate constant (k in mol/m2/s), default
  KM_RES_SURF_M = 'M', /// surface-scaled rate constant (k in kg/m2/s)
  KM_RES_PVS_N = 'V',  /// pore-volume-scaled rate constant (k in mol/m3/s)
  KM_RES_PVS_M = 'W',  /// pore-volume-scaled rate constant (k in kg/m3/s)
  KM_RES_ABS_N = 'F',  /// absolute (unscaled) rate constant (k in mol/s)
  KM_RES_ABS_M = 'G',  /// absolute (unscaled) rate constant (k in kg/s)
  KM_LIN_RATE = 'L'    /// linear growth/dissolution rate constant (v in m/s)

};

typedef enum { // Units of measurement of quantities and concentrations
  // number of components and phases
  QUAN_MKMOL = 'Y',
  QUAN_MMOL = 'h',
  QUAN_MOL = 'M', // NUMBER OF MOLES
  QUAN_MGRAM = 'y',
  QUAN_GRAM = 'g',
  QUAN_KILO = 'G', // MASS
  // concentrations of components and phases
  CON_MOLFR = 'n',
  CON_MOLPROC = 'N',
  CON_pMOLFR = 'f', // MOLE FRACTION
  CON_VOLFR = 'v',
  CON_VOLPROC = 'V',
  CON_pVOLFR = 'u', // VOLUME FRACTION
  CON_MOLAL = 'm',
  CON_MMOLAL = 'i',
  CON_pMOLAL = 'p', // MOLALITY
  CON_MOLAR = 'L',
  CON_MMOLAR = 'j',
  CON_pMOLAR = 'q', // MOLARITY
  CON_WTFR = 'w',
  CON_WTPROC = '%',
  CON_PPM = 'P', // MASS FRACTION
  CON_AQWFR = 'C',
  CON_AQWPROC = 'c',
  CON_AQPPM = 'a', // CONCENTRATION
  // aqueous species
  CON_AQGPL = 'd',
  CON_AQMGPL = 'e',
  CON_AQMKGPL = 'b', // VOLUME CONCENTRATION

  // Units of measurement of pressure Pr, P  { b B p P A }'
  PVT_BAR = 'b',   // bar - default, 1 bar = 1e5 Pa
  PVT_KBAR = 'B',  // kbar, 1 kbar = 1000 bar
  PVT_PASC = 'p',  // Pascal (Pa)
  PVT_KPASC = 'P', // MPa, 1 MPa = 10 bar = 1e6 Pa
  PVT_ATM = 'A',   // atm, 1 atm = 1.013 bar
  // Attention: Only b code can be used in this version!

  // Units of measurement of molar volume  { c j a L m }'
  PVT_CM3 = 'c',  // cm3, cm3/mole
  PVT_LITR = 'L', // liters (L) - volume of the system only, 1 L = 1000 cm3
  PVT_JBAR = 'j', // J/bar, 10 cm3/mole = 1 J/bar
  PVT_CBAR = 'a', // (cal/bar), 41.84 cm3/mole = 1 cal/bar
  // m  - reserved.
  // Attention: only j code can be used in this version!

  // Units of measurement of reference temperature Tr { C K F }'
  PVT_CELS = 'C',   // degrees Celsius (C)
  PVT_KELVIN = 'K', // Kelvins (K), 0 C = 273.15 K
  PVT_FAREN = 'F',  // degrees Fahrenheit (F)
  // Attention: Only C code can be used in this version of GEM IPM algorithm.

  // Units of measurement of energy values { j c J C n N }
  TDAT_JOUL = 'j',  // Joules (J/mole)
  TDAT_KJOUL = 'J', // kilojoules (kJ/mole)
  TDAT_CAL = 'c',   // calories (cal/mole); 1 cal = 4.184 J;
  TDAT_KCAL = 'C',  // kilocalories (kcal/mole)
  TDAT_NORM = 'N'   // normalized (mole/mole, J/mole/RT, cal/mole/RT)
                    // Attention: Only j code can be used in this version!
} SPPA_UNITS;

#endif // M_CONST_BASE_H
