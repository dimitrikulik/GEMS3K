//-------------------------------------------------------------------
// $Id: s_sorpmod.cpp 725 2026-10-08 15:43:37Z kulik $
//
/// \file s_sorpmod.cpp  Implementation of TSorpMod base class
/// for surface complexation models on a surface type (patch)
///  (TSCM_NEM, TSCM_CCM, TSCM_DLM, TSCM_BSM, TSCM_TLM, TSCM_CDM, TSCM_ETLM,
///  TSCM_BET) and T_NEM_calc derived class for NEM SCMs
//
// Copyright (C) 2017, 2026  D.Kulik
// <GEMS Development Team, mailto:dkulik@congineer.com>
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

#include "s_sorpmod.h"
// #include "s_solmod.h"

//=============================================================================================
// TSorpMod base class for multicomponent solid, liquid, fluid and aqueous
// mixtures (c) DK/TW 2006 (version 1) (c) DK/TW November 2010 (version 2)
//=============================================================================================

/// generic constructor (new)
TSorptionPhaseMod::TSorptionPhaseMod(SorptionPhaseData *sdph)
    : ModCode(sd->Mod_Code), MixCode(sd->Mix_Code), PhaseName(sd->phaseName),
      NComp(sd->NSpecies), NPar(sd->NParams), NPcoef(sd->NPcoefs),
      MaxOrd(sd->MaxOrder), NP_DC(sd->NPperDC), /*NPTP_DC(NPTPperDC),*/
      NSub(sd->NSublat), NMoi(sd->NMoiet), R_CONST(8.31451), Tk(sd->T_k),
      Pbar(sd->P_bar) {
  DC_Names = sd->arSM;
  //   NlPh = sd->NlPhs;
  //    NlPc = sd->NlPhC;
  NDQFpc = sd->NDQFpDC;
  //    NrcPpc = sd->NrcPpDC;
  //   lPhcf = sd->lPhc;
  DQFcf = sd->arDQFc; // read-only
                      //    rcpcf = sd->rcpc;  // read-only
                      //    if(rcpcf == NULL) NrcPpc = 0;
  // PhLin = sd->arPhLin;
  //    PhLin = new long int[NlPh][2];
  //    for (long int i=0; i<NlPh; i++)
  //    {
  //        PhLin[i][0] = sd->arPhLin[2*i];
  //        PhLin[i][1] = sd->arPhLin[2*i+1];
  //    }
  // pointer assignments
  aIPx = sd->arIPx; // Direct access to index list and parameter coeff arrays!
  aIPc = sd->arIPc;
  aIP = new double[NPar];
  for (int ii = 0; ii < NPar; ii++)
    aIP[ii] = 0.0;

  aGEX = sd->arGEX; // Reciprocal energies, Darken terms, pure fugacities of DC
                    // (corr. to TP)
  aPparc = sd->arPparc; // Partial pressures NComp
  aMoiSN = sd->arMoiSN; // End member moiety- site multiplicity number tables
                        // NComp x NSub x NMoi
  aSitFR = sd->arSitFr; // Tables of sublattice site fractions for moieties NSub
                        // x NMoi
  aDCc = sd->arDCc;
  DC_Codes = sd->DC_Codes;
  x = sd->arWx;
  phVOL = sd->aphVOL;
  aVol = sd->arVol;
  lnGamma = sd->arlnGam;
  lnGamConf = sd->arlnCnft; // new double[NComp];
  // lnGamCorr = sd->arlnCnft;  // new double[NComp];
  lnGamRecip = sd->arlnRcpt; // new double[NComp];
  lnGamEx = sd->arlnExet;    // new double[NComp];
  lnGamDQF = sd->arlnDQFt;   // new double[NComp];
  CTerm = sd->arCTermt;      // Coulombic terms (new)
  // Arrays for lnGamma components - should we zero off?
  for (long int i = 0; i < NComp; i++) {
    lnGamConf[i] = 0.0;
    lnGamRecip[i] = 0.0;
    lnGamEx[i] = 0.0;
    // lnGamCorr[i] = 0.0;
    lnGamDQF[i] = 0.0;
    CTerm[i] = 0.0;
  }
  // initialize phase properties
  Gex = 0.0;
  Hex = 0.0;
  Sex = 0.0;
  CPex = 0.0;
  Vex = 0.0;
  Aex = 0.0;
  Uex = 0.0;
  Gid = 0.0;
  Hid = 0.0;
  Sid = 0.0;
  CPid = 0.0;
  Vid = 0.0;
  Aid = 0.0;
  Uid = 0.0;
  Gdq = 0.0;
  Hdq = 0.0;
  Sdq = 0.0;
  CPdq = 0.0;
  Vdq = 0.0;
  Adq = 0.0;
  Udq = 0.0;
  Grs = 0.0;
  Hrs = 0.0;
  Srs = 0.0;
  CPrs = 0.0;
  Vrs = 0.0;
  Ars = 0.0;
  Urs = 0.0;
  alloc_multisite();
  init_multisite();
}

//=============================================================================================
// Implementation of NEM SCM
// References: Kulik et al., 2018 (NAB 18-27)
// (c) DKulik November 2026
//=============================================================================================

// Generic constructor for the TSCM_NEM class
T_NEM_calc::T_NEM_calc(long int NCmp, double Pp, double Tkp)
    : TSurfaceTypeMod(sds) {}

T_NEM_calc::~T_NEM_calc() {}

// Updates pure surface species properties
long int T_NEM_calc::SurfaceSpecies() { return 0; }

// Calculates T,P corrected SCM or EIL parameters
long int T_NEM_calc::PTSurfSpecParam() { return 0; }

// Calculates SACT (isotherm) corrections (no site balance models only)
long int T_NEM_calc::IsothermCorrMod() { return 0; }
// For site balance based SCMs, comp.Langmuir is ensured for each site type

// Calculates SCM surf.plane charges and potentials with Coulomb term
// corrections
long int T_NEM_calc::EILCorrMod() {
  // NEM ignores Coulomb corrections, but surface charge densities on planes
  // still can be computed

  long int j, retCode = 0;

  if (MixCode != MR_UNDEF_)
    return 1; // This should be 'N' code for NEM

  for (j = 0; j < NComp; j++)
    CTerm[j] = 0.0; // For NEM, Coulombic terms equal 0

  return retCode;

  return 0;
}
// All these functions apply to one surface type (patch) as the whole instance
// of TSorpMod base class does

/// Calculates ideal configurational terms in case of multi-site mixing
/// to preserve values computed in Phase scripts.
/// Only increments lnGamma[j] - may need to be cleaned before running MixMod
long int T_NEM_calc::SurfTypeMixMod() {}

// double TSCM_NEM::SCM_conf_entropy() {
//   long int j;
//   double si = 0.0;
//
//   for (j = 0; j < NComp; j++) {
//     if (x[j] > 1.0e-32)
//       si += x[j] * log(x[j]);
//   }
//   auto Sid1 = (-1.) * R_CONST * si;
//   return Sid1;
//  }

//--------------------- End of s_sorption.cpp ---------------------------
