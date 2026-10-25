//-------------------------------------------------------------------
// $Id$
//
/// \file ipm_chemical2.cpp
/// Implementation of chemistry-specific functions (concentrations,
/// activity coefficients, adsorption models etc.)
/// for the IPM convex programming Gibbs energy minimization algorithm
//
// Copyright (c) 1992-2012  D.Kulik, S.Dmitrieva, K.Chudnenko
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

#include "m_const_base.h"
#include "ms_multi.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculating bulk stoichiometry of a multicomponent phase
//
void TMultiBase::phase_bcs(long int N, long int M, long int jb, double *A,
                           double X[], double BF[]) {
  long int ii, i, j;
  double Xx;

  if (!A || !X || !BF)
    return;
  for (i = 0; i < N; i++)
    BF[i] = 0.;

  for (j = 0; j < M; j++) {
    Xx = X[j];
    if (fabs(Xx) < 1e-16) // was 1e-12
      continue;
    for (ii = arrL[j + jb]; ii < arrL[j + jb + 1]; ii++) {
      i = arrAN[ii];
      BF[i] += A[i + j * N] * Xx;
    }
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Adds phase to total bulk stoichiometry of all solid phases in the system
// Done on request by TW in November 2006
//
void TMultiBase::phase_bfc(long int k, long int jj) {
  long int ii, i, j;
  double Xx;

  if (pm.PHC[k] == PH_AQUEL || pm.PHC[k] == PH_GASMIX ||
      pm.PHC[k] == PH_FLUID || pm.PHC[k] == PH_PLASMA ||
      pm.PHC[k] == PH_SIMELT || pm.PHC[k] == PH_LIQUID)
    return;
  for (j = 0; j < pm.L1[k]; j++) {
    Xx = pm.X[j + jj];
    if (fabs(Xx) < 1e-12)
      continue;
    for (ii = arrL[j + jj]; ii < arrL[j + jj + 1]; ii++) {
      i = arrAN[ii];
      pm.BFC[i] += pm.A[i + (jj + j) * pm.N] * Xx;
    }
  }
}

/// Returns mass of all solid phases in grams (from the BFC vector)
double TMultiBase::bfc_mass(void) {
  double TotalMass = 0.;
  for (long int i = 0; i < pm.N; i++)
    TotalMass += pm.BFC[i] * pm.Awt[i];
  return TotalMass;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#define a(j, i) ((*(pm.A + (i) + (j) * pm.N)))
/// Calculation of dual chemical potentials, activities, and primal
/// concentrations for DCs (indexed jb to je) in a k-th phase.
// Input arrays X, XF, XFA,  input factors: Factor, MMC
//
//  Do we need this all in GEMIPM ?
//
void TMultiBase::CalculateConcentrationsInPhase(double X[], double XF[],
                                                double XFA[], double Factor,
                                                double MMC, double /*Dsur*/,
                                                long int jb, long int je,
                                                long int k) {
  long int j, ii, i;
  double Muj, /* DsurT=0.0,*/ SPmol, lnFmol = 4.016535;

  if (pm.PHC[0] == PH_AQUEL) { // mole fraction to molality conversion
    if (!k)
      lnFmol = log(1000. / MMC); // aq species
    else
      lnFmol = log(H2O_mol_to_kg); // 4.016535; 	   // other species
  }

  for (j = jb; j < je;
       j++) { // loop over DC - with important bugfixes from 02.04.2003
    Muj = DC_DualChemicalPotential(pm.U, pm.A + j * pm.N, pm.NR, j);
    pm.Fx[j] = Muj; // *pmp->RT; Fixed DK 12.03.2012

    //        if( X[j] <= pm.lowPosNum )
    if (X[j] <= pm.DcMinM) { // zeroing off
      pm.Wx[j] = 0.0;
      pm.VL[j] = log(pm.lowPosNum);
      pm.Y_w[j] = 0.0;
      pm.lnGam[j] = 0.0;
      if (pm.PHC[0] == PH_AQUEL)
        pm.Y_m[j] = 0.0;
      switch (pm.DCC[j]) // choice of expressions
      { // since 10.03.2008, changed the concept of DualTh activity
      case DC_SCP_CONDEN:
        pm.Y_la[j] = ln_to_lg * (Muj - pm.G0[j]);
        break;
      case DC_AQ_ELECTRON:
      case DC_AQ_PROTON:
      case DC_AQ_SPECIES:
      case DC_AQ_SURCOMP:
        pm.Y_la[j] = ln_to_lg * (Muj - pm.G0[j] + lnFmol);
        break;
      case DC_AQ_SOLVENT:
      case DC_AQ_SOLVCOM:
        pm.Y_la[j] = ln_to_lg * (Muj - pm.G0[j]);
        break;
      case DC_GAS_COMP:
      case DC_GAS_H2O:
      case DC_GAS_CO2: // gases
      case DC_GAS_H2:
      case DC_GAS_N2:
        pm.Y_la[j] = ln_to_lg * (Muj - pm.G0[j]);
        if (pm.P > 1e-29)
          pm.Y_la[j] += log10(pm.P);
        break;
      case DC_SOL_IDEAL:
      case DC_SOL_MINOR:
      case DC_SOL_MAJOR:
      case DC_SOL_MINDEP:
      case DC_SOL_MAJDEP:
      case DC_SCM_SPECIES:
        pm.Y_la[j] = ln_to_lg * (Muj - pm.G0[j]);
        break;
      case DC_SUR_GROUP:
        pm.Y_la[j] =
            ln_to_lg * (Muj - pm.G0[j]); // + lnFmol ); corr. 06.10.10 DK
        break;
      case DC_SUR_COMPLEX:
      case DC_SUR_IPAIR:
      case DC_SUR_SPECIES:
      case DC_IESC_A:
      case DC_IEWC_B:
        pm.Y_la[j] =
            ln_to_lg * (Muj - pm.G0[j]); // + lnFmol ); corr. 06.10.10 DK
        break;
      case DC_PEL_CARRIER:
      case DC_SUR_MINAL:
      case DC_SUR_CARRIER: // sorbent
        pm.Y_la[j] = ln_to_lg * (Muj - pm.G0[j]);
        break;
      default:
        break; // error in DC class code
      }
      continue;
    }
    // calculation of the mole fraction
    pm.Wx[j] = X[j] / XF[k];
    if (X[j] > std::min(pm.lowPosNum, pm.DcMinM))
      pm.VL[j] =
          log(pm.Wx[j]); // this is used nowhere except in some scripts. Remove?
    else
      pm.VL[j] = log(pm.lowPosNum); // debugging 29.11.05 KD
    pm.Y_la[j] = 0.0;
    switch (pm.DCC[j]) // choice of expressions
    {
    case DC_SCP_CONDEN:
      pm.Wx[j] = 1;
      pm.VL[j] = 0.0;
      if (pm.LO) {                 //  bugfix DK 08.02.10
        pm.Y_m[j] = X[j] * Factor; // molality
      }
      pm.Y_w[j] = // mass % of the system
          1e2 * X[j] * pm.MM[j] / pm.MBX;
      pm.Y_la[j] = ln_to_lg * (Muj - pm.G0[j]);
      pm.FVOL[k] += pm.Vol[j] * X[j];
      break;
    case DC_AQ_ELECTRON:
      pm.Y_m[j] = 0.0;
      pm.Y_la[j] = 0.0 - pm.pe;
      pm.Y_w[j] = 0.0;
      break;
    case DC_AQ_PROTON: // in molal scale!
      pm.pH = -ln_to_lg * (Muj - pm.G0[j] + lnFmol);
      [[fallthrough]];
    case DC_AQ_SPECIES:
      pm.IC += 0.5 * X[j] * Factor *
               (pm.EZ[j] * pm.EZ[j]); // increment to effective molal ionic
                                      // strength Ieff = 0.5*sum(c_i * z_i^2)
      pm.Zeff += X[j] * Factor * (pm.EZ[j] * pm.EZ[j] * pm.EZ[j] * pm.EZ[j]);
      //      Increment to sum(c_i * z_i^4) of effective mean-field valency
      //      z_eff = sqrt(sum(c_i * z_i^4) / sum(c_i * z_i^2))
      [[fallthrough]];
    case DC_AQ_SURCOMP:
      SPmol = X[j] * Factor; // molality
      //           pm.IC += 0.5* SPmol *(pm.EZ[j]*pm.EZ[j]); // Bugfix
      //           DK 21.10.2011 'K' species don't count here!
      pm.FVOL[k] += pm.Vol[j] * X[j]; // fixed 04.02.03 KD
      pm.Y_m[j] = SPmol;
      pm.Y_la[j] = ln_to_lg * (Muj - pm.G0[j] + lnFmol);
      pm.Y_w[j] = 1e6 * X[j] * pm.MM[j] / pm.FWGT[k];
      if (pm.DCC[j] != DC_AQ_SURCOMP) { // Bugfix 21.10.2011  DK - excluding 'K'
                                        // species from total IC molality
        //  Optimized for performance - calculation inline
        for (i = arrL[j]; i < arrL[j + 1]; i++) {
          ii = arrAN[i];
          if (ii >= pm.NR)
            continue;
          pm.IC_m[ii] += SPmol * a(j, ii); // total aqueous molality
          pm.IC_wm[ii] += X[j] * a(j, ii); // total aqueous mass concentration
        }
      }
      break;
    case DC_AQ_SOLVENT: // mole fractions of solvent
    case DC_AQ_SOLVCOM:
      //            pm.Y_m[j] = X[j]/XFA[k];  Replaced by DK 12.03.2012
      pm.Y_m[j] = H2O_mol_to_kg;
      pm.Y_w[j] = 1e3 * X[j] * pm.MM[j] / pm.FWGT[k];
      pm.FVOL[k] += pm.Vol[j] * X[j];
      pm.Y_la[j] = ln_to_lg * (Muj - pm.G0[j]);
      break;
    case DC_GAS_COMP:
    case DC_GAS_H2O:
    case DC_GAS_CO2: // gases
    case DC_GAS_H2:
    case DC_GAS_N2:
      pm.FVOL[k] += pm.Vol[j] * X[j];
      pm.Y_la[j] = ln_to_lg * (Muj - pm.G0[j]);
      if (pm.P > 1e-9)
        pm.Y_la[j] += log10(pm.P);
      break;
    case DC_SOL_IDEAL:
    case DC_SOL_MINOR: // solution end member
    case DC_SOL_MAJOR:
    case DC_SOL_MINDEP:
    case DC_SOL_MAJDEP:
    case DC_SCM_SPECIES:
      pm.FVOL[k] += pm.Vol[j] * X[j];
      pm.Y_la[j] = ln_to_lg * (Muj - pm.G0[j]);
      if (pm.LO) {                 //  bugfix DK 16.04.2012
        pm.Y_m[j] = X[j] * Factor; // molality
      }
      break;
    case DC_SUR_GROUP:           // adsorption:
      pm.Y_m[j] = X[j] * Factor; // molality
      pm.Y_w[j] =                // mg/g sorbent
          1e3 * X[j] * pm.MM[j] / (MMC * XFA[k]);
      pm.Y_la[j] = ln_to_lg * (Muj - pm.G0[j]); // + lnFmol ); corr. 06.10.10
      pm.FVOL[k] += pm.Vol[j] * X[j];           // fixed 11.03.2008 KD
      break;
      // case DC_SUR_GROUP:  // Check if commenting away of this is correct
      // !!!
    case DC_SUR_COMPLEX:
    case DC_SUR_IPAIR:
    case DC_SUR_SPECIES:
    case DC_IESC_A:
    case DC_IEWC_B:
      pm.Y_m[j] = X[j] * Factor; // molality
      pm.Y_w[j] =                // mg/g sorbent
          1e3 * X[j] * pm.MM[j] / (MMC * XFA[k]);
      //           DsurT = MMC * pm.Aalp[k] * pa->p.DNS*1.66054e-6;
      pm.Y_la[j] =
          ln_to_lg *
          (Muj - pm.G0[j]); // + lnFmol - pm.GEX[j] + Dsur + DsurT/( 1.0+DsurT )
      pm.FVOL[k] += pm.Vol[j] * X[j]; // fixed 11.03.2008   KD
      break;
    case DC_PEL_CARRIER:
    case DC_SUR_MINAL:
    case DC_SUR_CARRIER:         // sorbent
      pm.Y_m[j] = X[j] * Factor; // molality
      pm.Y_w[j] = 0.0;
      if (pm.YF[0] >= pm.DSM)
        pm.Y_w[j] = // mg of sorbent per kg aq solution
            1e6 * X[j] * pm.MM[j] / pm.FWGT[0];
      pm.Y_la[j] =
          ln_to_lg * (Muj - pm.G0[j]); // - pm.GEX[j] + Dsur - 1. + 1./(1.+Dsur)
                                       // - DsurT + DsurT/(1+DsurT)
      pm.FVOL[k] += pm.Vol[j] * X[j];
      break;
    default:
      break; // error in DC class code
    };
  } // j
  pm.Zeff = sqrt(pm.Zeff / (2.0 * pm.IC));
  // Effective valency z_eff = sqrt(sum(c_i * z_i^4) / sum(c_i * z_i^2))
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of derived values (concentrations etc.) on IPM iteration
///  from X,XF, and XFA vectors. Also calculates pH, pe, Eh
// This function has to be rewritten using new set of built-in
// chemical functions.
//
void TMultiBase::CalculateConcentrations(double X[], double XF[],
                                         double XFA[]) {
  long int k, ii, i, j, ist, jj, jja;
  double Factor = 0.0, Dsur = 0.0, MMC = 0.0;
  const BASE_PARAM *pa_p = base_param();

  //    if( pm.Ls < 2 || !pm.FIs )  Temporary disabled  09.03.2010 DK
  //        return;

  for (i = 0; i < pm.N; i++)
    pm.BFC[i] = 0.; // cleanup for phase_bfc() calculation

  for (j = 0; j < pm.Ls; j++) {
    pm.Wx[j] = 0.;
    pm.VL[j] = 0.;
  }
  j = 0;
  pm.VXc = 0.0;
  for (k = 0; k < pm.FI; k++) { // cycle by phases
    i = j + pm.L1[k];
    pm.FWGT[k] = 0.0;
    pm.FVOL[k] = 0.0;
    //   Dsur = 0.0;

    if (XF[k] > pm.DSM && !(pm.PHC[k] == PH_SORPTION && XFA[k] <= pa_p->ScMin))
      phase_bfc(k, j);

    if (k >= pm.FIs || pm.L1[k] == 1) { // this is a single- component phase
      pm.Wx[j] = 1.0;                   // SD 04/05/2010
      if (XF[k] < pm.DSM) {             // This phase to be zeroed off
        if (pm.LO)
          pm.Y_m[j] = 0.0;
        pm.Wx[j] = 0.0;
        pm.Y_w[j] = 0.0;
        pm.Fx[j] = DC_DualChemicalPotential(pm.U, pm.A + j * pm.N, pm.NR, j);
        pm.Y_la[j] = ln_to_lg * (pm.Fx[j] - pm.G0[j]); // -pm.GEX[j]
        //                pm.Fx[j] *= pm.RT;     // el-chem potential
        goto NEXT_PHASE;
      }
      pm.VL[j] = 0.0;
      if (pm.LO && XFA[0] > 0)
        pm.Y_m[j] = X[j] * 1000. / 18.01528 / XFA[0]; // molality
      pm.Y_w[j] =                                     // mass % in the system
          1e2 * X[j] * pm.MM[j] / pm.MBX;
      pm.Fx[j] = DC_DualChemicalPotential(pm.U, pm.A + j * pm.N, pm.NR, j);
      pm.Y_la[j] = ln_to_lg * (pm.Fx[j] - pm.G0[j]); // - pm.GEX[j]
      //            pm.Fx[j] *= pm.RT;     // el-chem potential
      pm.FWGT[k] += X[j] * pm.MM[j];
      pm.FVOL[k] += X[j] * pm.Vol[j];
      goto NEXT_PHASE;
    }
    if (pm.LO && !k) { // aqueous phase present
      pm.IC = 0.;      // fix 20.06.13 DK
      pm.Zeff = 0.;    // added 17.10.2026 DK
      for (ii = 0; ii < pm.N; ii++) {
        pm.IC_m[ii] = 0.0;
        pm.IC_lm[ii] = 0.0;
        pm.IC_wm[ii] = 0.0;
      }
    }
    if (XF[k] <= pm.DSM ||
        (pm.PHC[k] == PH_AQUEL && (XFA[k] <= pm.XwMinM || XF[k] <= pm.DSM)) ||
        (pm.PHC[k] == PH_SORPTION && XFA[k] <= pm.ScMinM)) {
      for (jj = 0; jj < pm.N; jj++)
        pm.BF[k * pm.N + jj] = 0.;

      for (jj = j; jj < i; jj++) // Loop added 10.03.01  KD (GTDEMO)
      {
        pm.Wx[j] = 0.0;
        if (pm.LO)
          pm.Y_m[jj] = 0.0;
        pm.Y_w[jj] = 0.0;
        pm.Fx[jj] = DC_DualChemicalPotential(pm.U, pm.A + jj * pm.N, pm.NR, jj);
        pm.Y_la[jj] = ln_to_lg * (pm.Fx[jj] - pm.G0[jj]);
        if (pm.PHC[k] ==
            PH_AQUEL) // || pm.PHC[k] == PH_SORPTION ) corr. 06.10.10 DK
          pm.Y_la[jj] += 1.74438;
        if (pm.PHC[k] == PH_GASMIX || pm.PHC[k] == PH_PLASMA)
          pm.Y_la[jj] += log10(pm.P);
        //                pm.Fx[jj] *= pm.RT;     // el-chem potential
        pm.lnGam[jj] = 0.0;
      }
      goto NEXT_PHASE;
    }
    // calculate bulk stoichiometry of a multicomponent phase
    phase_bcs(pm.N, pm.L1[k], j, pm.A + j * pm.N, pm.X + j, pm.BF + k * pm.N);

    switch (pm.PHC[k]) {
    case PH_AQUEL:
      MMC = 0.0; // molar mass of carrier
      //            Dsur = XFA[k]/XF[k] - 1.0; // Asymm.corr. - aq only!
      //            if( XFA[k] > pm.lowPosNum )
      if (XFA[k] > pm.XwMinM) {
        for (jj = j; jj < i; jj++)
          if (pm.DCC[jj] == DC_AQ_SOLVENT || pm.DCC[jj] == DC_AQ_SOLVCOM)
            MMC += pm.MM[jj] * X[jj] / XFA[k];
      } else
        MMC = 18.01528; // Assuming water-solvent
      //            if( (XFA[k] > pm.lowPosNum) && (MMC > pm.lowPosNum) )
      if ((XFA[k] > pm.XwMinM) && (MMC > pm.lowPosNum))
        Factor = 1000. / MMC / XFA[k]; // molality
      else
        Factor = 0.0;
      //            pm.IC=0.;
      pm.pe = ln_to_lg * pm.U[pm.N - 1];
      pm.Eh = 0.000086 * pm.U[pm.N - 1] * pm.T;
      [[fallthrough]];
    case PH_GASMIX:
    case PH_FLUID:
    case PH_PLASMA:
    case PH_SIMELT:
    case PH_HCARBL:
    case PH_IONEX:
    case PH_ADSORPT:
    case PH_SINCOND:
    case PH_SINDIS:
    case PH_LIQUID:
      pm.YFk = XF[k];
      for (jj = j; jj < i; jj++) {
        if (X[jj] > pm.DcMinM) // fixed 30.08.2009 DK
          pm.FWGT[k] += X[jj] * pm.MM[jj];
      }
      break;
    case PH_POLYEL:
    case PH_SORPTION:  // only sorbent end-members!
      pm.YFk = XFA[k]; // not XF[k] as any other phase?
      MMC = 0.0;

      for (ist = 0; ist < pm.FIat; ist++)
        pm.XFTS[k][ist] = 0.0;
      //           if( XFA[k] < pm.lowPosNum ) XFA[k] = pm.lowPosNum;
      if (XFA[k] < pm.ScMinM)
        XFA[k] = pm.ScMinM;
      for (jj = j; jj < i; jj++) {
        jja = jj - (pm.Ls - pm.Lads);
        if (pm.DCC[jj] == DC_SUR_CARRIER || pm.DCC[jj] == DC_SUR_MINAL ||
            pm.DCC[jj] == DC_PEL_CARRIER) {
          MMC += pm.MM[jj] * X[jj] / XFA[k];
          // Only sorbent mass
          pm.FWGT[k] += X[jj] * pm.MM[jj];
        } else { // Doubtful, do we really need it? Total amount of species on
                 // surf site type or on surface type (patch)?
          /*!!!!!*/ ist = pm.SATC[jja][XL_ST]; // Check carfefully!
          pm.XFTS[k][ist] += X[jj];
        }
      }
      pm.FWGTS[k] = pm.FWGT[k]; // Use for surfaces and charges only pm.FWGTS[k]
      pm.logYFk = log(pm.YFk);
      //            Dsur = XFA[k]/XF[k] - 1.0;  // Also for sorption phases
      //            if( Dsur <= -1.0 ) Dsur = -0.999999;
      break;
    case PH_ADSSBEIL: // New: sorption phase with site balances and without
                      // carrier DC
      if (XF[k] < pm.ScMinM)
        // this is in error because of site mole balance, so the amount cannot
        // be zero
        break;
      pm.YFk = XF[k];
      if (pm.FWGTS[k] < pm.ScMinM) {
        // this is probably the first calculation for this chemical system, so
        // we get the virtual sorbent mass from site IC amounts
        pm.FWGTS[k] = virtualSorbentMass(j, i, pm.N);
      }
      pm.FWGT[k] = pm.FWGTS[k];
      pm.logYFk = log(pm.YFk);
      break;
    default:
      return; // Phase class code error!
    }
    // calculation of species concentrations in k-th phase
    CalculateConcentrationsInPhase(X, XF, XFA, Factor, MMC, Dsur, j, i, k);

  NEXT_PHASE:
    pm.VXc += pm.FVOL[k];
    if (pm.PHC[k] == PH_AQUEL && XF[k] > pm.DSM && XFA[k] > pm.XwMinM)
      for (ii = 0; ii < pm.NR; ii++) {
        if (pm.LO) {
          if (pm.IC_m[ii] >= pa_p->DB)
            pm.IC_lm[ii] = ln_to_lg * log(pm.IC_m[ii]);
          else
            pm.IC_lm[ii] = 0;
          if (pm.FWGT[k] >= pa_p->DB)
            pm.IC_wm[ii] *= pm.Awt[ii] * 1000. / pm.FWGT[k];
          else
            pm.IC_wm[ii] = 0;
        }
      }
    j = i;
  } // k
}

// Calculation of virtual sorbent mass (g) from site IC amounts in the system
// First, we detect from DC stoichiometry which site ICs are used in this phase,
// then calculate the sorbent mass from those IC amounts and their mol masses
double TMultiBase::virtualSorbentMass(long int jb, long int je, long int N) {
  const BASE_PARAM *pa_p = base_param();
  long int j, i, nDConSite[N];
  double SMASS = 0.0, Btot[N];
  for (i = 0; i < N; i++) { // initial cleaning
    nDConSite[i] = 0;
    Btot[i] = 0.0;
  }
  for (j = jb; j < je; j++) {
    // screening which site ICs are used in DCs in this phase
    for (i = 0; i < N; i++) {
      if (a(j, i) <= pa_p->DB)
        continue;
      if (pm.ICC[i] == IC_SSITE ||
          pm.ICC[i] ==
              IC_ADDIT) { // IC_ADDIT is provisorial and should be removed!
                          // Now kept here because of error E37FPun when
                          // opening a project caused by IC_SSITE
        nDConSite[i] += 1;
        if (Btot[i] < pa_p->DB)
          Btot[i] = pm.B[i];
      }
    }
  }
  // Now calculating the mass of sorbent in the system for the sorption phase
  for (i = 0; i < N; i++) {
    if (Btot[i] < pa_p->DB || nDConSite[i] == 0)
      continue;
    SMASS += Btot[i] * pm.Awt[i];
  }
  return SMASS;
  //  #define a(j, i) ((*(pm.A + (i) + (j) * pm.N)))
}

//--------------------------------------------------------------------------------
/// Calculation of surface charge densities on multi-surface sorption phase
// Revised by DK on 17.10.2026 using updated input data propagation
void TMultiBase::IS_EtaCalc() {
  long int k, i, ist, j = 0, ja;
  double XetaS = 0., XetaW = 0., etaA, Ez, CDA, CDB, CDC;

  for (k = 0; k < pm.FIs; k++) { // loop over all solution phases
    i = j + pm.L1[k];
    if (pm.FIat > 0 && (pm.PHC[k] == PH_POLYEL || pm.PHC[k] == PH_SORPTION ||
                        pm.PHC[k] == PH_ADSSBEIL))
      for (ist = 0; ist < MNST; ist++) { // Replaced pm.FIat with MNST
        pm.XetaA[k][ist] = 0.0;          // plane charges in eq (moles)
        pm.XetaB[k][ist] = 0.0;
        pm.XetaC[k][ist] = 0.0;
        pm.XetaD[k][ist] = 0.0; // added 12.09.05  KD
      }

    if (pm.XF[k] <= pm.DSM ||
        (pm.PHC[k] == PH_AQUEL && (pm.X[pm.LO] <= pm.XwMinM //  pa->p.XwMin
                                   || pm.XF[k] <= pm.DHBM)) ||
        (pm.PHC[k] == PH_SORPTION &&
         pm.XF[k] <= pm.ScMinM) // no site balances, sorbent DC
        || (pm.PHC[k] == PH_ADSSBEIL &&
            pm.XF[k] <= pm.ScMinM) // site balances, no sorbent DC
        )                          //  pa->p.ScMin) )
      goto NEXT_PHASE;

    switch (pm.PHC[k]) { // initialization according to the phase class
    case PH_AQUEL:       // aqueous solution
      pm.Yw = pm.XFA[k];
      XetaW = 0.0;
      break;
    case PH_PLASMA:
    case PH_SIMELT:
      XetaS = 0.0;
      break;
    case PH_POLYEL:
    case PH_SORPTION:
    case PH_ADSSBEIL:
      for (ist = 0; ist < MNST; ist++) {
        etaA = 0.0; // Units of etaA, etaB, ... in moles (eq)
        if (fabs(pm.XfIEC[k][ist][FSIT_S_CAPACITY]) > pm.DHBM) {
          etaA += pm.XfIEC[k][ist][FSIT_S_CAPACITY] // eq/kg
                  * pm.FWGT[k] / 1e3; // must be converted to eq (moles)
        } else {
          // capacity in eq/kg not given, trying to get the amount of perm
          // charge from density in 1/nm2, mass and SSA
          if (fabs(pm.XfIEC[k][ist][FSIT_S_DENSITY]) > pm.DHBM)
            etaA += pm.XfIEC[k][ist][FSIT_S_DENSITY] * pm.FWGT[k] * pm.Aalp[k] *
                    1e-18 * pm.MASDT[k][ist][MSDT_FSUTSA]; // convert nm2 to m2
        }
        pm.XetaA[k][ist] = etaA; // Moles of permanent charge
      }
      break;
    default:
      break;
    }
    for (; j < i; j++) { // loop over DC for calculating total phase charge
      if (pm.X[j] <= pm.lowPosNum * 100.)
        continue; // Skipping too low concentrations
      ja = j - (pm.Ls - pm.Lads);

      switch (pm.DCC[j]) // select expressions for species classes
      {
      case DC_AQ_ELECTRON:
      case DC_AQ_PROTON:
      case DC_AQ_SPECIES:
      case DC_AQ_SURCOMP:
        XetaW += pm.X[j] * pm.EZ[j];
      case DC_AQ_SOLVCOM:
      case DC_AQ_SOLVENT:
        break;
      case DC_PEL_CARRIER:
      case DC_SUR_MINAL:
      case DC_SUR_CARRIER: // charge of carrier: ???????
                           // pm.XetaA[k] += pm.X[j]*pm.EZ[j];
        break;
      // surface species
      case DC_SUR_GROUP:
      case DC_SUR_COMPLEX:
      case DC_SUR_IPAIR:
      case DC_SUR_SPECIES:
      case DC_IESC_A:
      case DC_IEWC_B:
        // SCM Modernization: charges from mCDc[ja][0..1]
        CDA = pm.mCDc ? pm.mCDc[ja][CHARGE_ON_A]
                      : 0.0; // species charge that goes into 0 plane
        CDB = pm.mCDc ? pm.mCDc[ja][CHARGE_ON_B]
                      : 0.0; // species charge that goes into b (1) plane
        CDC = pm.mCDc ? pm.mCDc[ja][CHARGE_ON_C]
                      : 0.0; // species charge that goes into c (2) plane
        //
        Ez = pm.EZ[j]; // take formula charge if no plane charges are provided
        if (fabs(CDA) < pm.DHBM && fabs(CDB) < pm.DHBM && fabs(CDC) < pm.DHBM)
          CDA = Ez; // Only as 0-plane charge
        // FInd out to which surface type this surface species belongs
        ist = pm.SATC[ja][SURSP_ALLOC_SURFTYPE];
        pm.XetaA[k][ist] += pm.X[j] * CDA;
        pm.XetaB[k][ist] += pm.X[j] * CDB;
        pm.XetaC[k][ist] += pm.X[j] * CDC;
        break;
      default: // Any other phase: charge balance over the phase
        XetaS += pm.X[j] * pm.EZ[j];
        break;
      }
    } // j
    // Charge balance for d-plane on all surface types
    if (pm.FIat > 0 && (pm.PHC[k] == PH_POLYEL || pm.PHC[k] == PH_SORPTION ||
                        pm.PHC[k] == PH_ADSSBEIL)) {
      for (ist = 0; ist < MNST; ist++) {
        pm.XetaD[k][ist] =
            0.0 - pm.XetaA[k][ist] - pm.XetaB[k][ist] - pm.XetaC[k][ist];
      }
    }
  NEXT_PHASE:
    j = i;
    if (pm.LO && !k && pm.FIat > 0) {
      pm.XetaA[k][0] = XetaW;
      pm.XetaB[k][0] = XetaW;
      pm.XetaC[k][0] = XetaW;
      pm.XetaD[k][0] = XetaW;
    }
    if ((pm.PHC[k] == PH_PLASMA || pm.PHC[k] == PH_SIMELT) && pm.FIat)
      pm.XetaA[k][0] = XetaS; // Doubtful if k runs only over sorption phases?
    // it seems that in this function, it runs over all solution phases...
  } // k
    // compare pm.Xetaf[k]+pm.XetaA[k]+pm.XetaB[k] and XetaS
    // Test XetaW
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of the surface potential pm[q].XpsiD[k] on diffuse.
/// layer plane on k-th sorption phase. From total charge pm.Xeta[k]
/// ( in moles ) using Gouy-Chapman equation.
/// Strictly valid at PSI < 30 mV. Modified by DK 5 Jan 2000
/// to add a Basic Stern EIL model, then revised on 18.10.2026 to process
/// EIL models plus NEM using the same data structure
///    Added 13.03.2008 by DK: returns int value showing (if not 0)
///    that some extreme values were reached for charge densities or
///    electric potentials (for detecting bad PIA), 0 otherwise
//
long int TMultiBase::GouyChapman(long int, long int, long int k) {
  long int ist, status = 0;
  double sigA = 0., sigB = 0., sigC = 0., sigD = 0., sigmaD = 0., debye_m,
         x_slip;
  double SigmA[MNST], SigmB[MNST], SigmC[MNST], SigmD[MNST];
  double psiA, psiB, psiC, psiD, psiZ, convFactor, ConvFactor[MNST];
  double PsiA[MNST], PsiB[MNST], PsiC[MNST], PsiD[MNST], PsiZ[MNST];
  double f1, f3, A, Sig, F2RT, I, Zeff, Cap, rhoW, epsW, SMASSg = 0.0;
  EILModelResult emr;

  if (pm.XF[k] < pm.ScMinM) // TProfil::pm->pa.p.ScMin )
    return status;          // no sorbent or sorbed component amounts

  // Determining the sorbent mass in grams
  if (pm.PHC[k] == PH_POLYEL || pm.PHC[k] == PH_SORPTION ||
      pm.PHC[k] == PH_ADSSBEIL)
    SMASSg = pm.FWGTS[k];
  // newly implemented for the sorbent mass, both with or without site balances
  // Now calculating charge densities in C/m2
  for (ist = 0; ist < MNST;
       ist++) // loop over surface types // ist<pm.FIat=MNST
  {
    sigA = 0.0, sigB = 0.0, sigC = 0.0, sigD = 0.0;
    convFactor = 1.;
    psiA = 0.0, psiB = 0.0, psiC = 0.0, psiD = 0.0, psiZ = 0.0; // Cleanup

    // SCM Modernization: SCMC[k][ist][0] and MASDT [k][ist][MSDT_FSUTSA]
    if (pm.SCMC[k][ist][0] == SC_NOT_USED ||
        pm.MASDT[k][ist][MSDT_FSUTSA] < 1e-9)
      continue; // potentially error
    convFactor = F_CONSTANT / (pm.FWGTS[k] * pm.Aalp[k] *
                               pm.MASDT[k][ist][MSDT_FSUTSA]); // C/mol/(m2)
    // Calculation of charge densities (now limited to total charges > max.
    // balance residual)
    if (fabs(pm.XetaA[k][ist]) > pm.DHBM)   // moles (eq)
      sigA = pm.XetaA[k][ist] * convFactor; // result in C/m2
    if (fabs(pm.XetaB[k][ist]) > pm.DHBM)   // moles (eq)
      sigB = pm.XetaB[k][ist] * convFactor; // result in C/m2
    if (fabs(pm.XetaC[k][ist]) > pm.DHBM)   // moles
      sigC = pm.XetaC[k][ist] * convFactor; // result in C/m2
    if (fabs(pm.XetaD[k][ist]) > pm.DHBM)   // moles
      sigD = pm.XetaD[k][ist] * convFactor; // C/m2

    // Limit maximum charge densities to prevent GEM divergence
    if (fabs(sigA) > 1.4) {
      ipm_logger->debug(
          "EIL charge density plane A {} truncated to +- 1.4 C/m2  "
          "IT={} k= {} ist= {}",
          sigA, pm.IT, k, ist);
      SigmA[ist] = sigA < 0.0 ? -1.4 : 1.4;
      status = 60;
    } else
      SigmA[ist] = sigA;
    if (fabs(sigB) > 1.0) {
      ipm_logger->debug(
          "EIL charge density plane B {} truncated to +- 1.0 C/m2  "
          "IT={} k= {} ist= {}",
          sigB, pm.IT, k, ist);
      SigmB[ist] = sigB < 0.0 ? -2.0 : 2.0;
      status = 61;
    } else
      SigmB[ist] = sigB;
    if (fabs(sigC) > 1.0) {
      ipm_logger->debug(
          "EIL charge density plane C {} truncated to +- 1.0 C/m2  "
          "IT={} k= {} ist= {}",
          sigC, pm.IT, k, ist);
      SigmC[ist] = sigC < 0.0 ? -2.0 : 2.0;
      status = 62;
    } else
      SigmC[ist] = sigC;
    if (fabs(sigD) > 0.7) {
      ipm_logger->debug(
          "EIL charge density plane D {} truncated to +- 0.7 C/m2  "
          "IT={} k= {} ist= {}",
          sigD, pm.IT, k, ist);
      SigmD[ist] = sigD < 0.0 ? -1.4 : 1.4;
      status = 63;
    } else
      SigmD[ist] = sigD;

    F2RT = pm.FRT / 2.;
    I = pm.IC;
    Zeff = pm.Zeff;
    rhoW = pm.denW[0];
    epsW = pm.epsW[0];
    // Calculation of Debye length in metres
    debye_m = calc_debye_length(I);
    x_slip = pm.XcapZ[k][ist] * debye_m; // spli plane shift factor f * debye_m

    // Calculating EIL models for this surface type (patch)
    switch (pm.SCMC[k][ist][SCMC_EILM_CODE]) {
    case SC_NEM: // Non-electrostatic, charge densities on zero plane
      // calculated, but potentials and Coulomb terms are zeros
      emr = compute_nem(SC_NEM, SigmA[ist]);
      break;
    case SC_4LM: // Four-layer model
      emr = compute_flm(SC_4LM, I, Zeff, pm.XcapA[k][ist], pm.XcapB[k][ist],
                        pm.XcapC[k][ist], SigmA[ist], SigmB[ist], SigmC[ist],
                        x_slip, debye_m);
      break;
    case SC_3LM: // Three-layer EIL model for CD-MUSIC
      emr =
          compute_cd_music(SC_3LM, I, Zeff, pm.XcapA[k][ist], pm.XcapB[k][ist],
                           SigmA[ist], SigmB[ist], x_slip, debye_m);
      break;
    case SC_TLM: // Triple-layer model
      emr = compute_tlm(SC_TLM, I, Zeff, pm.XcapA[k][ist], pm.XcapB[k][ist],
                        SigmA[ist], SigmB[ist], x_slip, debye_m);
      break;
    case SC_BSM: // Basic Stern model
      emr = compute_bsm(SC_BSM, I, Zeff, pm.XcapA[k][ist], SigmA[ist], x_slip,
                        debye_m);
      break;
    case SC_CCM: // Constant-capacitance model
      emr = compute_ccm(SC_CCM, pm.XcapA[k][ist], SigmA[ist]);
      break;
    case SC_DDL: // Double-layer model
      emr = compute_ddl(SC_DDL, I, Zeff, SigmA[ist], x_slip);
      break;
    case SC_TCE: // Tournassat 2013 clay edge potential model with charge
                 // spillover TBD
      emr = compute_tce(SC_TCE, I, Zeff, pm.XcapA[k][ist], pm.XcapB[k][ist],
                        SigmA[ist], SigmB[ist], pm.MASDT[k][ist], x_slip,
                        debye_m);
      break;
    default:
      break;
    }
    // Unpack results from EILModelResult emr fields to MULTI data for this
    // phase and surface type
    pm.XpsiA[k][ist] = emr.psi0;
    pm.XpsiB[k][ist] = emr.psi1;
    pm.XpsiC[k][ist] = emr.psi2;
    pm.XpsiD[k][ist] = emr.psid;
    pm.XpsiZ[k][ist] = emr.zeta;  // zeta potential
    pm.XcapD[k][ist] = emr.cap_d; // Eff. capacitance density of diffuse layer
    pm.XdlA[k][ist] =
        emr.sld01; // Thickness of inner Stern layer (0 to 1 plane) nm
    pm.XdlB[k][ist] =
        emr.sld02; // Thickness of outer Stern layer (1 to 2 plane) nm
    pm.XdlC[k][ist] =
        emr.sld0d; // Eff. thickness of Stern layer (0 to 3-plane or d-plane)
    pm.XdlD[k][ist] = emr.dlth; // Eff. thickness of diffuse layer (nm)
                                // or 2*Debye_length ?
  } // ist   end of loop over surface types

  return status;
}
/*
  SigA = 0.;
  SigB = 0.;
  SigD = 0.;
  SigDDL = 0.; // Cleanup 07.12.2009 DK
  pm.XcapD[k][ist] = 0.0;
  pm.XpsiD[k][ist] = 0.0;
  if (fabs(XetaA[ist]) < pm.DHBM && // pm.lowPosNum*1e6 &&
      fabs(XetaB[ist]) < pm.DHBM && // pm.lowPosNum*1e6 &&
      fabs(XetaD[ist]) < pm.DHBM)   // pm.lowPosNum*1e6 )
    goto GEMU_CALC;                 // skipping at near-zero charge
  // calculating charge density at diffuse layer
  switch (pm.SCMC[k][ist][0]) {
  case SC_CCM:   // Constant-Capacitance Model Schindler, extension Nilsson
  case SC_ISCCM: // 2026 SC_ISCCM uses dependcy of capacity with IS C =
                 // sqrt(IS)/alpha, alpha is given instead of capacity
    SigA = pm.XfIEC[k][ist][0] + XetaA[ist];
    SigDDL = -SigA - XetaB[ist];
    SigB = XetaB[ist];
    break;
  case SC_DDLM: // Generalized Double Layer Model [Dzombak and Morel, 1990]
    SigA = pm.XfIEC[k][ist][0] + XetaA[ist];
    SigDDL = -SigA;
    SigB = 0.0;
    break;
  case SC_TLM: // Triple-Layer Model [Hayes and Leckie, 1987]
    SigA = pm.XfIEC[k][ist][0] + XetaA[ist];
    SigB = XetaB[ist];
    SigDDL = -SigA - XetaB[ist];
    break;
  case SC_MTL: // Modified Triple-Layer Model [Robertson, 1997]
    SigA = pm.XfIEC[k][ist][0] + XetaA[ist];
    SigB = XetaB[ist];
    SigDDL = -SigA - XetaB[ist];
    break;
  case SC_BSM: // Basic Stern model: [Christl and Kretzschmar, 1999]
    SigA = pm.XfIEC[k][ist][0] + XetaA[ist];
    SigB = XetaB[ist];
    SigDDL = -SigA - XetaB[ist];
    break;
  case SC_3LM: // Three-layer model: Hiemstra ea 1996; Tadanier & Eick 2002
    SigA = pm.XfIEC[k][ist][0] + XetaA[ist];
    SigB = XetaB[ist];
    SigD = XetaD[ist];
    SigDDL = -SigA - SigB - SigD;
    break;
  case SC_MXC: // BSM for Ion-Exchange on permanent-charge surface
    SigA = pm.XfIEC[k][ist][0] + XetaA[ist];
    SigB = XetaB[ist];
    SigDDL = -SigA - XetaB[ist];
    break;
  case SC_NNE: // Non-Electrostatic
    break;
  default:
    continue;
  }
  //        if( fabs( SigD ) > 1 )
  //        {
  //           ipm_logger->debug("EDL charge density D {} truncated to +- 1
  //           C/m2  IT={} k= {} ist= {}", SigD, pm.IT, k, ist); SigD = SigD <
  //           0.0 ? -1.0: 1.0;
  //        }
  // Gouy-Chapman equation
  // parameters of diffuse layer using [Damaskin, 1987,p.192-195]
  A = 1e-9;
  F2RT = pm.FRT / 2.;
  Sig = SigDDL; //  - XetaW[ist] ;
  I = pm.IC;
  Zeff = pm.Zeff;
  rhoW = pm.denW[0];
  epsW = pm.epsW[0];
  if (I > 1e-7)
  // Aq solution density Ro included acc. to [Machesky ea., 1999]
  // Only for the basic Stern model for now
  {
    double Ro = 1.0;
    if (pm.SCMC[k][ist][0] == SC_BSM && pm.FVOL[0] > 1e-16)
      Ro = pm.FWGT[0] / pm.FVOL[0];
    A = sqrt(2000. * 8.854e-12 * pm.epsW[0] * pm.RT * I * Ro);
  }
  Cap = F2RT * sqrt(4. * A * A + Sig * Sig);

  // SD: workaround because of problems with log argument
  f3 = sqrt(1. + Sig * Sig / (4. * A * A)) - Sig / (2. * A);
  if (f3 < 1) {
    f1 = exp(-3. * F2RT);
    if (f3 < f1)
      f3 = f1;
  } else {
    f1 = exp(3. * F2RT);
    if (f3 > f1)
      f3 = f1;
  }
  ipm_logger->trace(" f1={} f2={}", f1, f3);
  PSIo = log(f3) / F2RT;
  //          PSIo = log( sqrt( 1.+Sig*Sig/(4.*A*A) ) - Sig/(2.*A) ) / F2RT;
  //          Cap0 = fabs(Sig/PSIo);
  //          Del = A*1e9/(2.*I*F)/cosh(PSIo*F2RT);
  //          pm.XdlD[k] = Del;
  pm.XcapD[k][ist] = Cap;
  pm.XpsiD[k][ist] = PSIo;
  PsiD = PSIo;
  // Truncating diffuse plane potential to avoid divergence
  if (fabs(PsiD) > 0.4) {
    ipm_logger->debug(
        "All EIL models: PsiD ={} truncated to +- 0.4 V IT={} k= {} ist= {}",
        PsiD, pm.IT, k, ist);
    PsiD = PsiD < 0 ? -0.4 : 0.4;
    status = 63;
  }
GEMU_CALC: // calculating potentials at EIL planes
  switch (pm.SCMC[k][ist][SCMC_EILM_CODE]) {
  case SC_NEM: // Non-electrostatic, charge densities on zero plane
               // calculated, but potentials and Coulomb terms are zeros
    break;
  case SC_4LM: // Four-layer model
    break;
  case SC_3LM: // Three-layer EIL model for CD-MUSIC
  case SC_TLM: // Triple-layer model
    break;
  case SC_BSM: // Basic Stern model
    break;
  case SC_CCM: // Constant-capacitance model
    break;
  case SC_DDL: // Double-layer model
    break;
  default:
    break;
  }
  // calculating potentials at EIL planes
  switch (pm.SCMC[k][ist][0]) {
  case SC_DDLM: // Generalized Double Layer Model [Dzombak & Morel 1990]
    pm.XpsiA[k][ist] = PsiD;
    pm.XpsiB[k][ist] = PsiD;
    break;
  case SC_CCM:   // Constant-Capacitance Model Schindler, ext. Nilsson
  case SC_ISCCM: // 2026 SC_ISCCM uses dependcy of capacity with IS C =
                 // sqrt(IS)/alpha, alpha is given instead of capacity
    if (pm.XcapA[k][ist] >
        0.0) { // Classic CCM Schindler with inner-sphere species only
      if (pm.SCMC[k][ist][0] == SC_ISCCM)
        PsiA = SigA / (sqrt(I) / pm.XcapA[k][ist]);
      else
        PsiA = SigA / pm.XcapA[k][ist];
      if (fabs(PsiA) > 0.7) // truncated 0-plane potential
      {
        PsiA = PsiA < 0 ? -0.7 : 0.7;
        status = 64;
      }
      pm.XpsiA[k][ist] = PsiA;
    }
    if (pm.XcapB[k][ist] <
        0.001) { // Extended CCM model [Nilsson ea 1996] as TLM with PsiD = 0
      if (approximatelyZero(pm.XcapB[k][ist])) {
        PsiB = 0.0; //-SigB<0? -0.3: 0.3;
        ipm_logger->debug("Extended CCM model XcapB = {} IT={} k= {} ist= {}",
                          pm.XcapB[k][ist], pm.IT, k, ist);
      } else {
        if (pm.SCMC[k][ist][0] == SC_ISCCM)
          PsiB = -SigB / (sqrt(I) / pm.XcapB[k][ist]);
        else
          PsiB = -SigB / pm.XcapB[k][ist];
      }
      if (fabs(PsiB) > 0.3) // truncated B-plane potential
      {
        PsiB = PsiB < 0 ? -0.3 : 0.3;
        status = 65;
      }
      if (approximatelyZero(pm.XcapA[k][ist])) {
        PsiA = 0.0; // SigA<0? -0.7: 0.7;
        ipm_logger->debug("Extended CCM model XcapA = {} IT={} k= {} ist= {}",
                          pm.XcapA[k][ist], pm.IT, k, ist);
      } else {
        if (pm.SCMC[k][ist][0] == SC_ISCCM)
          PsiA = PsiB + SigA / (sqrt(I) / pm.XcapA[k][ist]);
        else
          PsiA = PsiB + SigA / pm.XcapA[k][ist];
      }
      if (fabs(PsiA) > 0.7) {
        PsiA = PsiA < 0 ? -0.7 : 0.7;
        status = 66;
      }
      pm.XpsiA[k][ist] = PsiA;
      pm.XpsiB[k][ist] = PsiB;
    }
    break;
  case SC_MTL: // Modified Triple Layer Model for X- Robertson | Kulik
               // PsiD = 0.0; // test
    if (approximatelyZero(pm.XcapB[k][ist])) {
      PsiB = 0.0; // SigB<0 ? -0.6: 0.6;
      ipm_logger->debug(
          "Modified Triple Layer Model XcapB = {} IT={} k= {} ist= {}",
          pm.XcapB[k][ist], pm.IT, k, ist);
    } else {
      PsiB = PsiD - SigDDL / pm.XcapB[k][ist];
    }
    if (fabs(PsiB) > 0.6) // truncated B-plane potential
    {
      ipm_logger->debug(
          "EIL (MTL) PsiB = {} truncated to +- 0.6 V IT={} k= {} ist= {}",
PsiB, pm.IT, k, ist); PsiB = PsiB < 0 ? -0.6 : 0.6; status = 67;
    }
    if (approximatelyZero(pm.XcapA[k][ist])) {
      PsiA = 0.0; // SigA<0? -1.1: 1.1;
      ipm_logger->debug(
          "Modified Triple Layer Model XcapA = {} IT={} k= {} ist= {}",
          pm.XcapA[k][ist], pm.IT, k, ist);
    } else {
      PsiA = PsiB + SigA / pm.XcapA[k][ist];
    }
    if (fabs(PsiA) > 1.1) // truncated 0 plane potential
    {
      ipm_logger->debug(
          "EIL (MTL) PsiA ={} truncated to +- 1.1 V IT={} k= {} ist= {}",
PsiA, pm.IT, k, ist); PsiA = PsiA < 0 ? -1.1 : 1.1; status = 68;
    }
    pm.XpsiA[k][ist] = PsiA;
    pm.XpsiB[k][ist] = PsiB;
    break;
  case SC_TLM: // Triple-Layer Model   [Hayes 1987]
    if (approximatelyZero(pm.XcapB[k][ist])) {
      PsiB = 0.0; // SigB<0 ? -0.6: 0.6;
      ipm_logger->debug("Triple-Layer Model XcapB = {} IT={} k= {} ist= {}",
                        pm.XcapB[k][ist], pm.IT, k, ist);
    } else {
      PsiB = PsiD - SigDDL / pm.XcapB[k][ist];
    }
    if (fabs(PsiB) > 0.6) // // truncated B-plane potential
    {
      ipm_logger->debug(
          "EIL (TLM) PsiB ={} truncated to +- 0.6 V IT={} k= {} ist= {}",
PsiB, pm.IT, k, ist); PsiB = PsiB < 0 ? -0.6 : 0.6; status = 69;
    }
    if (approximatelyZero(pm.XcapA[k][ist])) {
      PsiA = 0.0; // SigA<0? -1.1: 1.1;
      ipm_logger->debug("Triple-Layer Model XcapA = {} IT={} k= {} ist= {}",
                        pm.XcapA[k][ist], pm.IT, k, ist);
    } else {
      PsiA = PsiB + SigA / pm.XcapA[k][ist];
    }
    if (fabs(PsiA) > 1.1) // truncated 0-plane potential
    {
      ipm_logger->debug(
          "EIL (TLM) PsiA ={} truncated to +- 1.1 V IT={} k= {} ist= {}",
PsiA, pm.IT, k, ist); PsiA = PsiA < 0 ? -1.1 : 1.1; status = 70;
    }
    pm.XpsiA[k][ist] = PsiA;
    pm.XpsiB[k][ist] = PsiB;
    break;
  case SC_3LM: // Three-Layer Model [Hiemstra & van Riemsdijk 1996]
               //            PsiB = PsiD + SigD / pm.XcapB[k][ist];
    ipm_logger->debug("EIL (3LM) PsiB(D) ={} IT={} k= {} ist= {}", PsiB,
pm.IT, k, ist); if (approximatelyZero(pm.XcapB[k][ist])) { PsiB = 0.0; //
SigA<0 ? -0.6: 0.6; ipm_logger->debug("Three-Layer Model XcapB = {} IT={} k=
{} ist= {}", pm.XcapB[k][ist], pm.IT, k, ist); } else { PsiB = PsiD + (SigA +
SigB) / pm.XcapB[k][ist]; // Compare!
    }
    ipm_logger->debug("EIL (3LM) PsiB(AB) ={} IT={} k= {} ist= {}", PsiB,
pm.IT, k, ist); if (fabs(PsiB) > 0.6) // truncated B-plane potential
    {
      ipm_logger->debug(
          "EIL (3LM) PsiB ={} truncated to +- 0.6 V IT={} k= {} ist= {}",
PsiB, pm.IT, k, ist); PsiB = PsiB < 0 ? -0.6 : 0.6; status = 71;
    }
    if (approximatelyZero(pm.XcapA[k][ist])) {
      PsiA = 0.0; // SigA<0? -1.1: 1.1;
      ipm_logger->debug("Three-Layer Model XcapA = {} IT={} k= {} ist= {}",
                        pm.XcapA[k][ist], pm.IT, k, ist);
    } else {
      PsiA = PsiB + SigA / pm.XcapA[k][ist];
    }
    if (fabs(PsiA) > 1.1) // truncated 0-plane potential
    {
      ipm_logger->debug(
          "EIL (3LM) PsiA ={} truncated to +- 1.1 V IT={} k= {} ist= {}",
PsiA, pm.IT, k, ist); PsiA = PsiA < 0 ? -1.1 : 1.1; status = 72;
    }
    pm.XpsiA[k][ist] = PsiA;
    pm.XpsiB[k][ist] = PsiB;
    break;
  case SC_BSM: // Basic Stern model, Christl & Kretzschmar, 1999
    PsiB = PsiD;
    if (fabs(PsiB) > 0.6) // truncated B-plane potential
    {
      ipm_logger->debug(
          "EIL (BSM) PsiB ={} truncated to +- 0.6 V IT={} k= {} ist= {}",
PsiB, pm.IT, k, ist); PsiB = PsiB < 0 ? -0.6 : 0.6; status = 73;
    }
    if (approximatelyZero(pm.XcapA[k][ist])) {
      PsiA = 0.0; // SigA<0? -1.1: 1.1;
      ipm_logger->debug("Basic Stern model XcapA = {} IT={} k= {} ist= {}",
                        pm.XcapA[k][ist], pm.IT, k, ist);
    } else {
      PsiA = PsiB + SigA / pm.XcapA[k][ist];
    }
    if (fabs(PsiA) > 1.1) // truncated 0-plane potential
    {
      ipm_logger->debug(
          "EIL (BSM) PsiA ={} truncated to +- 1.1 V IT={} k= {} ist= {}",
PsiA, pm.IT, k, ist); PsiA = PsiA < 0 ? -1.1 : 1.1; status = 74;
    }
    pm.XpsiA[k][ist] = PsiA;
    pm.XpsiB[k][ist] = PsiB;
    break;
  case SC_MXC: // BSM for permanent charge surfaces
    PsiB = PsiD;
    if (fabs(PsiB) > 0.6) // truncated B-plane potential
    {
      ipm_logger->debug(
          "EIL (MXC) PsiB ={} truncated to +- 0.6 V IT={} k= {} ist= {}",
PsiB, pm.IT, k, ist); PsiB = PsiB < 0 ? -0.6 : 0.6; status = 75;
    }
    if (approximatelyZero(pm.XcapA[k][ist])) {
      PsiA = 0.0; // SigA<0? -1.1: 1.1;
      ipm_logger->debug(
          "BSM for permanent charge surfaces XcapA = {} IT={} k= {} ist= {}",
          pm.XcapA[k][ist], pm.IT, k, ist);
    } else {
      PsiA = PsiB + SigA / pm.XcapA[k][ist];
    }
    if (fabs(PsiA) > 1.1) // truncated 0-plane potential
    {
      ipm_logger->debug(
          "EIL (MXC) PsiA ={} truncated to +- 1.1 V IT={} k= {} ist= {}",
PsiA, pm.IT, k, ist); PsiA = PsiA < 0 ? -1.1 : 1.1; status = 76;
    }
    pm.XpsiA[k][ist] = PsiA;
    pm.XpsiB[k][ist] = PsiB;
    break;
  case SC_NNE: // Non-Electrostatic
    pm.XpsiA[k][ist] = 0.0;
    pm.XpsiB[k][ist] = 0.0;
    break;
  default:
    break;
  }
*/
// } // ist   end of loop over surface types
// return status;
// }

// Added on 18.10.2026  by DK in ClaySorEl project
// Dynamic relative permittivity of water via Malmberg & Maryott (1956)
// polynomial.
double TMultiBase::calc_water_permittivity(double T_kelvin,
                                           double P_bar) const {
  double T_C = T_kelvin - 273.15;
  double epsW1, epsW2, epsW;
  epsW1 = 87.740 - 0.40008 * T_C + 9.398e-4 * (T_C * T_C) -
          1.410e-6 * (T_C * T_C * T_C);
  // for P > 1 bar, get water permittivity from SUPCRT
  if (P_bar > 1.0) {
    epsW2 = pm.epsW[0];
    epsW = epsW2;
  } else
    epsW = epsW1;
  return epsW;
}

// Calculates Debye screening length (kappa^-1) in meters from ionic strength I
// (mol/L).
double TMultiBase::calc_debye_length(double I) const {
  double EPS_0 = 8.8541878128e-12;
  double FARADAY = 96485.3321;
  double RT = pm.RT;
  double I_m3 = I * 1000.0;
  //  double Rho = pm.FWGT[0] / pm.FVOL[0]; // density of aq phase
  double rhoW = pm.denW[0]; // density of water at T,P
  double epsW =
      calc_water_permittivity(pm.Tc, pm.Pc); // permittivity of water at T,P
  double debye_m;                            // Debye length in m
  debye_m =
      std::sqrt((EPS_0 * epsW * RT) / (2.0 * FARADAY * FARADAY * I_m3 * rhoW));
  // debye_m = std::sqrt((EPS_0 * epsW * RT) /
  //                 (2.0 * FARADAY * FARADAY * I_m3 * Rho));
  return debye_m;
}

// Calculates capacitance density of the diffuse layer from ionic strength I
//
double TMultiBase::calc_cap_d_from_sigma_d(double sigma_d, double I) const {
  double I_m3 = I * 1000.0;
  double EPS_0 = 8.8541878128e-12;
  double rhoW = pm.denW[0];
  double epsW =
      calc_water_permittivity(pm.Tc, pm.Pc); // permittivity of water at T,P
  double RT = pm.RT;
  double Sig = sigma_d;
  double A = 1.0;
  if (I > 1e-7) {
    double Rho = 1.0;
    if (pm.FVOL[0] > 1e-16)
      Rho = pm.FWGT[0] / pm.FVOL[0];
    A = sqrt(2.0 * EPS_0 * epsW * pm.RT * I_m3 * Rho); // or Rho = rhoW ?
    // Aq solution density Rho included acc. to [Machesky ea., 1999]
  }
  double Cap_d = pm.FRT / 2.0 * sqrt(4. * A * A + Sig * Sig);
  return Cap_d;
}

// Inverted Gouy-Chapman relation generalized for mixed aqueous electrolytes
// with effective valency z_eff.
double TMultiBase::calc_psi_d_from_sigma_d(double sigma_d, double I,
                                           double z_eff) const {
  double I_m3 = I * 1000.0;
  double EPS_0 = 8.8541878128e-12;
  double RT = pm.RT;
  double F_RT = pm.FRT; // FARADAY / (R_GAS * T);
  double rhoW = pm.denW[0];
  double eps_r = calc_water_permittivity(pm.Tc, pm.Pc);
  double gc_coeff = std::sqrt(8.0 * EPS_0 * eps_r * RT * I_m3 * rhoW);
  double arg = sigma_d / (-gc_coeff);
  double psiD = (2.0 / (z_eff * F_RT)) * std::asinh(arg);
  // Truncating diffuse plane potential to avoid divergence of GEM
  if (fabs(psiD) > 0.5) {
    //   ipm_logger->debug(
    //       "All EIL models: PsiD ={} truncated to +- 0.5 V IT={} k= {} ist=
    //       {}", psi_d, pm.IT, k, ist);
    psiD = psiD < 0.0 ? -0.5 : 0.5;
    //    status = 64;
  }
  return psiD;
}

// Non-linear Poisson-Boltzmann potential decay through the diffuse layer at a
// given distance x_dist.
double TMultiBase::calc_gc_potential_at_distance(double psi_d, double I,
                                                 double x_dist,
                                                 double z_eff) const {
  if (std::abs(psi_d) < 1e-7)
    return 0.0;
  // double F_RT = FARADAY / (R_GAS * T);
  double F_RT = pm.FRT;
  double gamma_0 = std::tanh(z_eff * F_RT * psi_d / 4.0);
  double debye_len = calc_debye_length(I);
  double ratio = gamma_0 * std::exp(-x_dist / debye_len);
  ratio = std::clamp(ratio, -0.9999999, 0.9999999);
  double psi_x = (4.0 / (z_eff * F_RT)) * std::atanh(ratio);
  return psi_x;
}

EILModelResult TMultiBase::compute_flm(char ModCode, double I, double z_eff,
                                       double C1, double C2, double C3,
                                       double sig0, double sig1, double sig2,
                                       double x_slip, double debye_m) const {
  double sigma0 = sig0;
  double sigma1 = sig1;
  double sigma2 = sig2;
  double sigmad = -(sigma0 + sigma1 + sigma2);
  double sigmanet = sigma0 + sigma1 + sigma2 + sigmad;

  double cap_d = calc_cap_d_from_sigma_d(sigmad, I);
  double psi_d = calc_psi_d_from_sigma_d(sigmad, I, z_eff);
  double psi_2 = psi_d + (-(sigmad) / C3);
  double psi_1 = psi_2 + ((sigma0 + sigma1) / C2);
  double psi_0 = psi_1 + (sigma0 / C1);
  double zeta = calc_gc_potential_at_distance(psi_d, I, x_slip, z_eff);

  EILModelResult res;
  res.modelcode = ModCode;
  res.name = "Four-Layer Model (FLM, 4LM)";
  res.sigma0 = sigma0;
  res.sigma1 = sigma1;
  res.sigma2 = sigma2;
  res.sigmad = sigmad;
  res.sigmanet = sigmanet;
  res.psi0 = psi_0;
  res.psi1 = psi_1;
  res.psi2 = psi_2;
  res.psid = psi_d;
  res.zeta = zeta;
  res.cap_d = cap_d;
  // res.sld01 = dist01; ....
  return res;
}

EILModelResult TMultiBase::compute_tlm(char ModCode, double I, double z_eff,
                                       double C1, double C2, double sig0,
                                       double sig1, double x_slip,
                                       double debye_m) const {
  double sigma0 = sig0;
  double sigma_beta = sig1;
  double sigmad = -(sigma0 + sigma_beta);
  double sigmanet = sigma0 + sigma_beta + sigmad;

  double cap_d = calc_cap_d_from_sigma_d(sigmad, I);
  double psi_d = calc_psi_d_from_sigma_d(sigmad, I, z_eff);
  double psi_beta = psi_d - (sigmad / C2);
  double psi_0 = psi_beta + (sigma0 / C1);
  double zeta = calc_gc_potential_at_distance(psi_d, I, x_slip, z_eff);

  EILModelResult res;
  res.modelcode = ModCode;
  res.name = "Triple-Layer Model (TLM)";
  res.sigma0 = sigma0;
  res.sigma1 = sigma_beta;
  res.sigma2 = 0.0;
  res.sigmad = sigmad;
  res.sigmanet = sigmanet;
  res.psi0 = psi_0;
  res.psi1 = psi_beta;
  res.psi2 = 0.0;
  res.psid = psi_d;
  res.zeta = zeta;
  res.cap_d = cap_d;
  // res.sld01 = dist01; ....
  return res;
}

EILModelResult TMultiBase::compute_cd_music(char ModCode, double I,
                                            double z_eff, double C1, double C2,
                                            double sig0, double sig1,
                                            double x_slip,
                                            double debye_m) const {
  double sigma0 = sig0;
  double sigma1 = sig1;
  double sigmad = -(sigma0 + sigma1);
  double sigmanet = sigma0 + sigma1 + sigmad;

  double cap_d = calc_cap_d_from_sigma_d(sigmad, I);
  double psi_2 = calc_psi_d_from_sigma_d(sigmad, I, z_eff);
  double psi_1 = psi_2 + ((sigma0 + sigma1) / C2);
  double psi_0 = psi_1 + (sigma0 / C1);
  double zeta = calc_gc_potential_at_distance(psi_2, I, x_slip, z_eff);

  EILModelResult res;
  res.modelcode = ModCode;
  res.name = "CD-MUSIC Model (3LM)";
  res.sigma0 = sigma0;
  res.sigma1 = sigma1;
  res.sigma2 = 0.0;
  res.sigmad = sigmad;
  res.sigmanet = sigmanet;
  res.psi0 = psi_0;
  res.psi1 = psi_1;
  res.psi2 = 0.0;
  res.psid = psi_2;
  res.zeta = zeta;
  res.cap_d = cap_d;
  // res.sld01 = dist01; ....
  return res;
}

EILModelResult TMultiBase::compute_bsm(char ModCode, double I, double z_eff,
                                       double C1, double sig0, double x_slip,
                                       double debye_m) const {
  double sigma0 = sig0;
  double sigmad = -sigma0;
  double sigmanet = sigma0 + sigmad;

  double cap_d = calc_cap_d_from_sigma_d(sigmad, I);
  double psi_d = calc_psi_d_from_sigma_d(sigmad, I, z_eff);
  double psi_0 = psi_d + (sigma0 / C1);
  double zeta = calc_gc_potential_at_distance(psi_d, I, x_slip, z_eff);

  EILModelResult res;
  res.modelcode = ModCode;
  res.name = "Basic Stern Model (BSM)";
  res.sigma0 = sigma0;
  res.sigma1 = 0.0;
  res.sigma2 = 0.0;
  res.sigmad = sigmad;
  res.sigmanet = sigmanet;
  res.psi0 = psi_0;
  res.psi1 = 0.0;
  res.psi2 = 0.0;
  res.psid = psi_d;
  res.zeta = zeta;
  res.cap_d = cap_d;
  // res.sld01 = dist01; ....
  return res;
}

EILModelResult TMultiBase::compute_ddl(char ModCode, double I, double z_eff,
                                       double sig0, double x_slip) const {
  double sigma0 = sig0;
  double sigmad = -sigma0;
  double sigmanet = sigma0 + sigmad;

  double cap_d = calc_cap_d_from_sigma_d(sigmad, I);
  double psi_d = calc_psi_d_from_sigma_d(sigmad, I, z_eff);
  double psi_0 = psi_d;
  double zeta = calc_gc_potential_at_distance(psi_d, I, x_slip, z_eff);

  EILModelResult res;
  res.modelcode = ModCode;
  res.name = "Diffuse Double Layer (DDL) model";
  res.sigma0 = sigma0;
  res.sigma1 = 0.0;
  res.sigma2 = 0.0;
  res.sigmad = sigmad;
  res.sigmanet = sigmanet;
  res.psi0 = psi_0;
  res.psi1 = 0.0;
  res.psi2 = 0.0;
  res.psid = psi_d;
  res.zeta = zeta; // perhaps not reliable in case of DDL model
  res.cap_d = cap_d;
  // res.sld01 = dist01; ....
  return res;
}

EILModelResult TMultiBase::compute_ccm(char ModCode, double C1,
                                       double sig0) const {
  double sigma0 = sig0;
  double sigmad = -sigma0;
  double sigmanet = sigma0 + sigmad;
  double psi_0 = sigma0 / C1;
  double cap_d = 0.0;
  double psi_d = 0.0;
  double zeta = 0.0;

  EILModelResult res;
  res.modelcode = ModCode;
  res.name = "Constant Capacitance model (CCM)";
  res.sigma0 = sigma0;
  res.sigma1 = 0.0;
  res.sigma2 = 0.0;
  res.sigmad = sigmad;
  res.sigmanet = sigmanet;
  res.psi0 = psi_0;
  res.psi1 = 0.0;
  res.psi2 = 0.0;
  res.psid = psi_d;
  res.zeta = zeta;
  res.cap_d = cap_d;
  // res.sld01 = dist01; ....
  return res;
}

EILModelResult TMultiBase::compute_nem(char ModCode, double sig0) const {
  double sigma0 = sig0;
  double sigmanet = sigma0;

  EILModelResult res;
  res.modelcode = ModCode;
  res.name = "Non-electrostatic model (NEM)";
  res.sigma0 = sigma0;
  res.sigma1 = 0.0;
  res.sigma2 = 0.0;
  res.sigmad = 0.0;
  res.sigmanet = sigmanet;
  res.psi0 = 0.0;
  res.psi1 = 0.0;
  res.psi2 = 0.0;
  res.psid = 0.0;
  res.zeta = 0.0;
  res.cap_d = 0.0;
  // res.sld01 = 0.0; ....
  return res;
}

EILModelResult TMultiBase::compute_tce(char ModCode, double I, double z_eff,
                                       double C1, double C2, double sig0,
                                       double sig1, double par[], double x_slip,
                                       double debye_m) const {

  double A1 = par[1];
  double A2 = par[2];
  double A3 = par[3];
  double sigma0 = sig0;
  double sigma_beta = sig1;
  double sigmad = -(sigma0 + sigma_beta);
  double sigmanet = sigma0 + sigma_beta + sigmad;

  double cap_d = calc_cap_d_from_sigma_d(sigmad, I);
  double psi_d = calc_psi_d_from_sigma_d(sigmad, I, z_eff);
  double psi_beta = psi_d - (sigmad / C2);
  double psi_0 = psi_beta + (sigma0 / C1);
  double zeta = calc_gc_potential_at_distance(psi_d, I, x_slip, z_eff);

  EILModelResult res;
  res.modelcode = ModCode;
  res.name = "TLM for clay edge with charge spillover";
  res.sigma0 = sigma0;
  res.sigma1 = sigma_beta;
  res.sigma2 = 0.0;
  res.sigmad = sigmad;
  res.sigmanet = sigmanet;
  res.psi0 = psi_0;
  res.psi1 = psi_beta;
  res.psi2 = 0.0;
  res.psid = psi_d;
  res.zeta = zeta;
  res.cap_d = cap_d;
  // res.sld01 = dist01; ....
  return res;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
///  Calculation of new surface activity coefficient terms SACT (Kulik, 2004)
//
///  Revised by KD in April 2004 (PSI) to introduce new activity
///  coefficient terms SACT rigorously derived from Langmuir and QCA
///  isotherms (Kulik 2006, Radiochimica Acta).
///  Revised again by DK in October 2026 in ClaySorEl prolect for compatibility
///  with new data structures in Phase and MULTI.
//
///  SACT are placed into pm.lnGam[j], as other activity coefficients except
///  relative surface potentials (Coulombic terms) kept separately.
//
///  Old (obsolete) SAT calculations (Kulik 2000, 2002) are retained.
//
///  pm.lnSAC[*][3] vector is now used to keep original DUL[j] to restore
///  them after IPM-2 refinements for surface complexes.
///    Added 13.03.2008 by DK: returns int value showing (if true)
///    that some extreme values were obtained for some SACTs,
///    0 otherwise (for detecting bad PIA)
//
long int TMultiBase::SurfaceActivityCoeff(long int jb, long int je, long int,
                                          long int, long int k) {
  long int status = 0;
  long int i, ii, j, ja, ist = 0, iss, dentt, dentsn, Cj,
                         iSite[MNST]; // SCM Modernization: MNST
  double XS0, xj0, XVk, XSk, XSkC, xj, Mm, rIEPS, ISAT, XSs, SATst, xjn, q1, q2,
      GamC, alpF, qBet, bBet, cN, eF, lnGamjo, lnDiff, lnFactor;
  const BASE_PARAM *pa_p = base_param();

  if (pm.PHC[k] == PH_SORPTION || pm.PHC[k] == PH_POLYEL) {
    if (pm.XF[k] <= pm.DSM) // No sorption phase retained by the IPM
      return status;
    if (pm.XFA[k] <= pa_p->ScMin) // no amount of sorbent remains
      return status;
  } else {
    // SACTS do not apply for sorption phases with site balances [yet]
    if (pm.PHC[k] == PH_ADSSBEIL || pm.PHC[k] == PH_ADSORPT)
      return status;
  }
  // This needs a deeper reconsideration!
  // Initialization of work arrays: D[NSURFSITES][MNST]
  for (i = 0; i < MNST; i++) // SCM Modernization: loop over MNST surf patches
  {
    iSite[i] = -1;
    for (ii = 0; ii < NSURFSITES; ii++) // SCM Modernization: site types
      pm.D[ii][i] = 0.0; // cleaning the matrix for sites totals in surface
  }

  // Extraction of site indices for the neutral >OH group
  for (j = jb; j < je; j++) {
    ja = j - (pm.Ls - pm.Lads);
    if (pm.SATT[ja] == SAT_SOLV) {
      ist = pm.SATC[ja][XL_ST]; // / MSPN;
      iSite[ist] = j;
    }
    // Counting current sites totals (these species do not compete)
    if (pm.DCC[j] == DC_SUR_CARRIER || pm.DCC[j] == DC_SUR_MINAL ||
        pm.DCC[j] == DC_PEL_CARRIER || pm.SATT[ja] == SAT_SOLV)
      continue;
    // Get ist (index of surface type)
    ist = pm.SATC[ja][XL_ST];   // / MSPN;
    if (ist < 0 || ist >= MNST) // SCM Modernization: MNST
      ist = 0;                  // default: zero surface type
    // Get iss (index of site on surface type)
    iss = pm.SATC[ja][XL_SI];
    if (iss < 0 || iss >= MNST) // SCM Modernization: MNST
      iss = 0; // default: zero site is the weekest and the most abundant one?
    pm.D[iss][ist] += pm.X[j]; // adding to total amount on a site
  }

  for (j = jb; j < je; j++) { // Main loop for DCs - surface complexes
    lnGamjo = pm.lnGmo[j];    // bugfix 16.03.2008 DK
    if (pm.X[j] < std::min(pm.DcMinM, pm.lowPosNum))
      continue; // This surface DC has been killed by IPM
                //        OSAT = pm.lnGmo[j]; // added 6.07.01 by KDA
    ja = j - (pm.Ls - pm.Lads); // Assuming that sorption phase(s) are the last
                                // ones in the list of phases-solutions
    rIEPS = pa_p->IEPS;         // default 1e-3 (for old SAT - 1e-9)
                                //        dent = 1;  // default - monodentate
    switch (pm.DCC[j])          // code of species class
    {
    default: // pm.lnGam[j] = 0.0;  not a surface species
      continue;
    case DC_IEWC_B:
    case DC_IESC_A:
    case DC_SUR_GROUP:
    case DC_SUR_COMPLEX:
    case DC_SUR_IPAIR:
    case DC_SUR_SPECIES:
      // Get ist (index of surface type)
      ist = pm.SATC[ja][XL_ST]; // / MSPN;
      // get iss - index of site on surf.type
      iss = pm.SATC[ja][XL_SI];
      if (iss < 0 || iss >= MNST) // SCM Modernization: MNST
        iss = 0; // zero site is the weekest and the most abundant one
      // Cj - index of carrier DC
      Cj = pm.SATC[ja][XL_EM];
      if (Cj < 0) { // Assigned to the whole sorbent
        XVk = pm.XFA[k];
        Mm = pm.FWGT[k] / XVk;
      } else { // Assigned to one of the sorbent end-members
        XVk = pm.X[Cj];
        if (XVk < pm.ScMinM) // pm.DSM*0.1 )
          continue;          // This end-member is zeroed off by IPM
        Mm = pm.MM[Cj] * XVk / pm.XFA[k]; // mol.mass
      }
      XSk = pm.XFTS[k][ist]; // Total moles of sorbates on ist surface type
      XSs = pm.D[iss][ist];  // Total moles of SC on iss site type
      xj = pm.X[j];          // Current moles of this surface species

      // Extracting SACT isotherm parameters
      if (pm.MASDJ) {
        cN = pm.MASDJ[ja][PI_P2];    // Frumkin/Pivovarov water coord. number
        dentt = denticity_total(ja); // total denticity of this surf species
        dentsn = denticity_num_site_types(ja); // Number of site types to bind
        GamC = pm.MASDJ[ja][PI_DEN]; // Complete coverage density (individual,
                                     // or as set for the whole surf type)
        alpF = pm.MASDJ[ja][PI_P1];  // Frumkin lateral interaction energy alpha
        qBet = pm.MASDJ[ja][PI_P2];  // q parameter (QCA/BET) or cN in Pivovarov
        bBet = pm.MASDJ[ja][PI_P3];  // Parameter b in BET or S in Pivovarov
      } else {
        cN = 0.0;
        alpF = 1.0;
        qBet = 0.0;
        bBet = 0.0;
        dentt = 1; // bet = 1.0;
        dentsn = 1;
      }
      switch (pm.SATT[ja]) // selection of the SACT model
      {
      case SACT_LAN_NCOMP: // Non-competitive individual monodentate Langmuir
                           // TBD
        break;
      case SACT_LAN_COMP: // Competitive monodentate Langmuir on a surface type
                          // and site type
        XSkC = XSs / XVk / Mm * 1e6 / pm.MASDT[k][ist][MSDT_FSUTSA] /
               pm.Aalp[k] / 1.66054; // per nm2
        XS0 = (fabs(pm.MASDJ[ja][PI_DEN]) / pm.Aalp[k] / 1.66054);
        // max. density per nm2
        if (pa_p->PC <= 2)
          rIEPS = pa_p->IEPS * XS0; // relative IEPS
        if (XSkC < 0.0)
          XSkC = 0.0;
        if (XSkC >= XS0) // Setting limits
          XSkC = XS0 - 2.0 * rIEPS;
        q1 = XS0 - XSkC;
        if ((pa_p->PC == 3 && !pm.W1) || pa_p->PC != 3) {
          q2 = rIEPS * XS0;
          if (q1 > q2)
            q2 = q1;
        } else {
          q2 = q1;
          if (q2 <= 1e-22)
            q2 = 1e-22;
        }
        ISAT = log(XS0) - log(q2);
        pm.lnGam[j] = ISAT;
        pm.lnSACT[j] = ISAT;
        pm.lnSAC[ja][0] = ISAT;
        break;
      case SACT_QCA_NCOMP: // dent++;   bidentate is default for QCA
        xj0 = (fabs(pm.MASDJ[ja][PI_DEN]) / pm.Aalp[k] / 1.66054);
        // Max site density per nm2
        xj = XSs / XVk / Mm / pm.MASDT[k][ist][MSDT_FSUTSA] * 1e6 // xj
             / pm.Aalp[k] / 1.66054; // Density per nm2 on site type iss
        if (pa_p->PC <= 2)
          rIEPS = pa_p->IEPS * xj0; // relative IEPS
        if (xj >= xj0 / static_cast<double>(dentt))
          xj = xj0 / static_cast<double>(dentt) - rIEPS; // upper limit
        //                ISAT = 0.0;
        q2 = xj0 - xj * dentt; // Computing differences in QCA gamma
        q1 = xj0 - xj;
        if (q1 < 1e-22)
          q1 = 1e-22;
        if (q2 < 1e-22)
          q2 = 1e-22;
        ISAT = log(xj0) + log(q1) * (dentt - 1) - log(q2) * dentt;
        pm.lnGam[j] = ISAT;
        pm.lnSACT[j] = ISAT;
        pm.lnSAC[ja][0] = ISAT;
        break;
        // (Non)competitive QCA-L for 1 to 4 dentate species
      case SACT_QCA1_NCOMP: // dentt = 4;     // code '4'
      case SACT_QCA2_NCOMP: // dentt = 3;     // code '3'
      case SACT_QCA3_NCOMP: // dentt = 2;     // code '2'
      case SACT_QCA4_NCOMP: // dentt = 1;     // code '1'
        break;
      case SACT_FFG_COMP: // Competitive Frumkin (FFG) isotherm
        dentt = 1;        // monodentate for now
        xj0 = (fabs(pm.MASDJ[ja][PI_DEN]) / pm.Aalp[k] / 1.66054);
        // Max site density per nm2
        xj = XSs / XVk / Mm / pm.MASDT[k][ist][MSDT_FSUTSA] * 1e6 //  xj
             / pm.Aalp[k] / 1.66054; // Current species density per nm2
        if (pa_p->PC <= 2)
          rIEPS = pa_p->IEPS * xj0; // relative IEPS
        if (xj >= xj0 / dentt)
          xj = xj0 / dentt - rIEPS; // upper limit
        q2 = xj0 - xj * dentt;      // Computing differences in QCA gamma
        q1 = xj0 - xj;
        if (q1 < 1e-22)
          q1 = 1e-22;
        if (q2 < 1e-22)
          q2 = 1e-22;
        ISAT = log(xj0) + log(q1) * (dentt - 1) - log(q2) * dentt;
        // Calculation of the Frumkin exponential term
        if (fabs(alpF) < 1e-9 || fabs(cN) < 1e-9)
          eF = 0.0;
        else { // changed from cN to -cN on 16.09.2005 by KD
          eF = -cN * alpF * xj * dentt / xj0; // Fi = Fi'/(kT) Bockris p.938
        }
        pm.lnGam[j] = ISAT + eF;
        pm.lnSACT[j] = ISAT;
        pm.lnGammF[j] = eF;
        pm.lnSAC[ja][0] = ISAT;
        pm.lnSAC[ja][1] = eF;
        break;
      case SACT_FFG_NCOMP: // Non-Competitive Frumkin isotherm
                           // for permanent charge surfaces
        XSkC = xj / XVk / Mm / pm.MASDT[k][ist][MSDT_FSUTSA] * 1e6 /
               pm.Aalp[k] / 1.66054; // per nm2
        XS0 = (pm.MASDJ[ja][PI_DEN] / pm.Aalp[k] / 1.66054);
        // max.dens.per nm2
        if (pa_p->PC <= 2)
          rIEPS = pa_p->IEPS * XS0; // relative IEPS
        if (XSkC < 0.0)
          XSkC = 0.0;
        if (XSkC >= XS0) // Limits
          XSkC = XS0 - 2.0 * rIEPS;
        q1 = XS0 - XSkC;
        if ((pa_p->PC == 3 && !pm.W1) || pa_p->PC != 3) {
          q2 = rIEPS * XS0;
          if (q1 > q2)
            q2 = q1;
        } else {
          q2 = q1;
          if (q2 <= 1e-22)
            q2 = 1e-22;
        }
        ISAT = log(XS0) - log(q2);
        // Calculation of the Frumkin exponential term (competitive)
        if (fabs(alpF) < 1e-9 || fabs(cN) < 1e-9)
          eF = 0.0;
        else { // changed from cN to -cN on 16.09.2005 by KD
          eF = -cN * alpF * XSkC / XS0; // Fi = Fi'/(kT) Bockris p.938
        }
        pm.lnGam[j] = ISAT + eF;
        pm.lnSACT[j] = ISAT;
        pm.lnGammF[j] = eF;
        pm.lnSAC[ja][0] = ISAT;
        pm.lnSAC[ja][1] = eF;
        break;
      case SACT_PIV_NCOMP: // (Non)competitive Pivovarov isotherm
        dentt = 1;         // monodentate for now
        xj0 = (fabs(pm.MASDJ[ja][PI_DEN]) / pm.Aalp[k] / 1.66054);
        // Max site density per nm2
        xj = XSs / XVk / Mm / pm.MASDT[k][ist][MSDT_FSUTSA] * 1e6 / pm.Aalp[k] /
             1.66054; // Current density per nm2
        if (pa_p->PC <= 2)
          rIEPS = pa_p->IEPS * xj0; // relative IEPS
        if (xj >= xj0 / dentt)
          xj = xj0 / dentt - rIEPS; // upper limit
        q2 = xj0 - xj * dentt;      // Computing differences in gamma
        q1 = xj0 - xj;
        if (q1 < 1e-22)
          q1 = 1e-22;
        if (q2 < 1e-22)
          q2 = 1e-22;
        ISAT = log(xj0) + log(q1) * (dentt - 1) - log(q2) * dentt;
        // Calculation of the Frumkin exponential term
        if (fabs(alpF) < 1e-9 || fabs(cN) < 1e-9)
          eF = 0.0;
        else {
          double pivovar;
          eF = cN * alpF * xj / xj0; // Fi = Fi'/(kT) Bockris p.938
          // Calculation of the Pivovarov 98 exponential correction term
          pivovar = xj0 / (xj0 + xj * (cN - 1));
          eF *= pivovar;
        }
        pm.lnGam[j] = ISAT + eF;
        pm.lnSACT[j] = ISAT;
        pm.lnGammF[j] = eF;
        pm.lnSAC[ja][0] = ISAT;
        pm.lnSAC[ja][1] = eF;
        break;
      case SACT_VIR_NCOMP: // Non-Competitive virial isotherm
        dentt = 1;         // monodentate for now
        xj0 = (fabs(pm.MASDJ[ja][PI_DEN]) / pm.Aalp[k] / 1.66054);
        // Max site density per nm2
        xj = XSs / XVk / Mm / pm.MASDT[k][ist][MSDT_FSUTSA] * 1e6 / pm.Aalp[k] /
             1.66054; // Current density per nm2
        if (pa_p->PC <= 2)
          rIEPS = pa_p->IEPS * xj0; // relative IEPS
        if (xj >= xj0 / dentt)
          xj = xj0 / dentt - rIEPS; // upper limit
        ISAT = 0.0;
        if (fabs(alpF) < 1e-9 || fabs(cN) < 1e-9)
          eF = 0.0;
        else { // changed from cN to -cN on 16.09.2005 by KD
          eF = -cN * alpF * xj / xj0; // Fi = Fi'/(kT) Bockris p.938
        }
        pm.lnGam[j] = ISAT + eF;
        pm.lnSACT[j] = ISAT;
        pm.lnGammF[j] = eF;
        pm.lnSAC[ja][0] = ISAT;
        pm.lnSAC[ja][1] = eF;
        break;
      case SACT_FREU_NCOMP: // Non-competitive Freundlich isotherm TBD
        break;
      case SACT_BET_NCOMP: // Non-competitive BET for surface precipitation
        ISAT = 0.0;
        eF = 0.0;
        // To be completed
        //
        //
        pm.lnGam[j] = ISAT + eF;
        pm.lnSACT[j] = ISAT;
        pm.lnGammF[j] = eF;
        pm.lnSAC[ja][0] = ISAT;
        break;
        // TBD: cases for surface type isotherms with species bound to different
        // site
        // types
      case SACT_INDEF: // No SACT calculation whatsoever
        pm.lnGam[j] = 0.0;
        pm.lnSACT[j] = 0.0;
        pm.lnGammF[j] = 0.0;
        pm.lnSAC[ja][0] = 0;
        break;
      default: // pm.lnGam[j] = 0.0;
        break;
      }
      // Make optional ....
      switch (pm.SATT[ja]) // selection of the old SAT model
      {
      //  Obsolete old SAT calculations (retained from previous versions)
      case SAT_COMP: // Competitive SAT (obsolete) on a surface type
        if (iSite[ist] < 0)
          xjn = 0.0;
        else
          xjn = pm.X[iSite[ist]]; // neutral site does not compete!
        // SCM Modernization: MASDT patch fraction and density
        XS0 = fabs(pm.MASDJ[ja][PI_DEN]) * XVk * Mm / 1e6 *
              pm.MASDT[k][ist][MSDT_FSUTSA]; // expected total in moles
        if (pa_p->PC <= 2)
          rIEPS = pa_p->IEPS * XS0; // relative IEPS
        XSkC = XSk - xjn - xj;      // occupied by the competing species;
                                    // this sorbate cannot compete to itself
        if (XSkC < 0.0)
          XSkC = 0.0;
        if (XSkC >= XS0) // Limits
          XSkC = XS0 - 2.0 * rIEPS;
        xj0 = XS0 - XSkC; // expected moles of this sorbate
        if (xj >= xj0)
          xj = xj0 - rIEPS; // Limits
        if (xj * 2 <= xj0)
          ISAT = 0.0; // ideal case
        else {
          q1 = xj0 - xj;
          q2 = rIEPS * XS0;
          if ((pa_p->PC == 3 && !pm.W1) || pa_p->PC != 3) {
            if (q1 > q2)
              q2 = q1;
          } else {
            q2 = q1;
            if (q2 <= 1e-33)
              q2 = 1e-33;
          }
          ISAT = log(xj) - log(q2);
        }
        pm.lnGam[j] = ISAT;
        pm.lnSAC[ja][0] = ISAT;
        break;
      case SAT_NCOMP: // Non-competitive truncated Langmuir SAT (obsolete)
        // rIEPS = pa->p.IEPS * 2;
        xj0 = fabs(pm.MASDJ[ja][PI_DEN]) * XVk * Mm / 1e6 *
              pm.MASDT[k][ist][MSDT_FSUTSA]; // in moles
        if (pa_p->PC <= 2)
          rIEPS = pa_p->IEPS * xj0; // relative IEPS
        if (xj >= xj0)
          xj = xj0 - rIEPS;  // upper limit
        if (xj * 2.0 <= xj0) // Linear adsorption - to improve !
          ISAT = 0.0;
        else {
          q1 = xj0 - xj; // limits: rIEPS to 0.5*xj0
          q2 = xj0 * rIEPS;
          if (pa_p->PC == 3 && pm.W1)
            ISAT = log(xj) - log(q1);
          else {
            if (q1 > q2)
              ISAT = log(xj) - log(q1);
            else
              ISAT = log(xj) - log(q2);
          }
        }
        pm.lnGam[j] = ISAT;
        pm.lnSAC[ja][0] = ISAT;
        break;
      case SAT_SOLV: // Neutral surface site (e.g. >O0.5H@ group)
                     // applies to the whole surface type!
        XSs = 0.0;   // calc total moles on all sites on surface type
        for (i = 0; i < MNST; i++)
          XSs += pm.D[i][ist];
        XSkC = XSs / XVk / Mm * 1e6 // total non-solvent surf.species
               / pm.MASDT[k][ist][MSDT_FSUTSA] / pm.Aalp[k] /
               1.66054; // per nm2
        XS0 = fabs(pm.MASDJ[ja][PI_DEN]);
        SATst = pa_p->DNS * 1.66054 * pm.Aalp[k] / XS0;
        XS0 = XS0 / pm.Aalp[k] / 1.66054;
        if (pa_p->PC <= 2)
          rIEPS = pa_p->IEPS * XS0; // relative IEPS
        if (XSkC < 0.0)
          XSkC = 0.0;
        if (XSkC >= XS0) // Limits
          XSkC = XS0 - 2.0 * rIEPS;
        q1 = XS0 - XSkC;
        if ((pa_p->PC == 3 && !pm.W1) || pa_p->PC != 3) {
          q2 = rIEPS * XS0;
          if (q1 > q2)
            q2 = q1;
        } else {
          q2 = q1;
          if (q2 <= 1e-22)
            q2 = 1e-22;
        }
        ISAT = log(XS0) - log(q2);
        pm.lnGam[j] = ISAT + log(SATst);
        pm.lnSAC[ja][0] = ISAT;
        pm.lnSAC[ja][1] = log(SATst);
        break;
      }
    }

    if (lnGamjo > pm.DHBM) { // Workaround DK 07.12.2009
      lnFactor = 0.2;
      lnDiff = pm.lnGam[j] - lnGamjo;
      if (fabs(lnDiff) > fabs(lnFactor)) // e times
      {
        if (fabs(lnDiff) > 6.907755) // 1000 times
          status = 101;              // the SACT has changed too much;
                        // the threshold pa->p.IEPS needs adjustment!
        // Smoothing (provisional)
        pm.lnGam[j] = lnDiff > 0 ? lnGamjo + lnDiff - lnFactor
                                 : lnGamjo + lnDiff + lnFactor;
      }
    }
  } // j
  return status;
}

// Total denticity of surface species (binding to one or more site types)
long int TMultiBase::denticity_total(long int j) {
  long int i, dent_total = 0; //, dent_types = 0;
  for (i = 0; i < NSURFSITES; i++) {
    //   if (pm.XSMD[j][i] > 0)
    //     dentypes++;
    dent_total += pm.XSMD[j][i];
  }
  return dent_total;
}

// Number of different site types binding to this surface species
long int TMultiBase::denticity_num_site_types(long int j) {
  long int i, dent_types = 0; // , dent_total = 0;
  for (i = 0; i < NSURFSITES; i++) {
    if (pm.XSMD[j][i] > 0)
      dent_types++;
    //    dent_total += pm.XSMD[j][i];
  }
  return dent_types;
}

//--------------------- End of ipm_chemical2.cpp ---------------------------
