//-------------------------------------------------------------------
// $Id: s_activity2.cpp 845 2013-06-20 15:58:57Z kulik $
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

#include "node.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculating bulk stoichiometry of a multicomponent phase
//
void TActivity::phase_bcs(long int N, long int M, long int jb, double *A,
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
void TActivity::phase_bfc(long int k, long int jj) {
  long int ii, i, j;
  double Xx;

  if (act.PHC[k] == PH_AQUEL || act.PHC[k] == PH_GASMIX ||
      act.PHC[k] == PH_FLUID || act.PHC[k] == PH_PLASMA ||
      act.PHC[k] == PH_SIMELT || act.PHC[k] == PH_LIQUID)
    return;
  for (j = 0; j < act.L1[k]; j++) {
    Xx = act.X[j + jj];
    if (fabs(Xx) < 1e-12)
      continue;
    for (ii = arrL[j + jj]; ii < arrL[j + jj + 1]; ii++) {
      i = arrAN[ii];
      act.BFC[i] += act.A[i + (jj + j) * act.N] * Xx;
    }
  }
}

/// Returns mass of all solid phases in grams (from the BFC vector)
double TActivity::bfc_mass(void) {
  double TotalMass = 0.;
  for (long int i = 0; i < act.N; i++)
    TotalMass += act.BFC[i] * act.Awt[i];
  return TotalMass;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#define a(j, i) ((*(act.A + (i) + (j) * act.N)))
/// Calculation of dual chemical potentials, activities, and primal
/// concentrations for DCs (indexed jb to je) in a k-th phase.
// Input arrays X, XF, XFA,  input factors: Factor, MMC
//
void TActivity::CalculateConcentrationsInPhase(double X[], double XF[],
                                               double /*XFA*/[], double Factor,
                                               double MMC, double /*Dsur*/,
                                               long int jb, long int je,
                                               long int k) {
  long int j, i, ii;
  double Muj, /* DsurT=0.0,*/ SPmol, lnFmol = 4.016535;

  if (act.PHC[0] == PH_AQUEL) { // mole fraction to molality conversion
    if (!k)
      lnFmol = log(1000. / MMC); // aq species
    else
      lnFmol = log(H2O_mol_to_kg); // 4.016535; 	   // other species
  }

  for (j = jb; j < je;
       j++) { // loop over DC - with important bugfixes from 02.04.2003
    Muj = DC_DualChemicalPotential(act.U, act.A + j * act.N, act.NR, j);
    act.Fx[j] = Muj; // *pmp->RT; Fixed DK 12.03.2012

    //        if( X[j] <= act.lowPosNum )
    if (X[j] <= act.DcMinM) { // zeroing off
      act.Wx[j] = 0.0;
      //            act.VL[j] = log( lowPosNum );
      //            act.Y_w[j] = 0.0;
      act.lnGam[j] = 0.0;
      if (act.PHC[0] == PH_AQUEL)
        act.Y_m[j] = 0.0;
      switch (act.DCC[j]) // choice of expressions
      { // since 10.03.2008, changed the concept of DualTh activity
      case DC_SCP_CONDEN:
        act.Y_la[j] = ln_to_lg * (Muj - act.G0[j]);
        break;
      case DC_AQ_ELECTRON:
      case DC_AQ_PROTON:
      case DC_AQ_SPECIES:
      case DC_AQ_SURCOMP:
        act.Y_la[j] = ln_to_lg * (Muj - act.G0[j] + lnFmol);
        break;
      case DC_AQ_SOLVENT:
      case DC_AQ_SOLVCOM:
        act.Y_la[j] = ln_to_lg * (Muj - act.G0[j]);
        break;
      case DC_GAS_COMP:
      case DC_GAS_H2O:
      case DC_GAS_CO2: // gases
      case DC_GAS_H2:
      case DC_GAS_N2:
        act.Y_la[j] = ln_to_lg * (Muj - act.G0[j]);
        if (act.P > 1e-29)
          act.Y_la[j] += log10(act.P);
        break;
      case DC_SOL_IDEAL:
      case DC_SOL_MINOR:
      case DC_SOL_MAJOR:
      case DC_SOL_MINDEP:
      case DC_SOL_MAJDEP:
      case DC_SCM_SPECIES:
        act.Y_la[j] = ln_to_lg * (Muj - act.G0[j]);
        break;
      case DC_SUR_GROUP:
        act.Y_la[j] =
            ln_to_lg * (Muj - act.G0[j]); // + lnFmol ); corr. 06.10.10 DK
        break;
      case DC_SUR_SPECIES:
      case DC_SUR_COMPLEX:
      case DC_SUR_IPAIR:
      case DC_IESC_A:
      case DC_IEWC_B:
        act.Y_la[j] =
            ln_to_lg * (Muj - act.G0[j]); // + lnFmol ); corr. 06.10.10 DK
        break;
      case DC_PEL_CARRIER:
      case DC_SUR_MINAL:
      case DC_SUR_CARRIER: // sorbent
        act.Y_la[j] = ln_to_lg * (Muj - act.G0[j]);
        break;
      default:
        break; // error in DC class code
      }
      continue;
    }
    // calculation of the mole fraction
    act.Wx[j] = X[j] / XF[k];
    //        if( X[j] > DcMinM )
    //            act.VL[j] = log( act.Wx[j] );     // this is used nowhere
    //            except in some scripts. Remove?
    //        else act.VL[j] = log( DcMinM );   // was lowPosNum
    act.Y_la[j] = 0.0;
    switch (act.DCC[j]) // choice of expressions
    {
    case DC_SCP_CONDEN:
      act.Wx[j] = 1;
      //            act.VL[j] = 0.0;
      if (act.LO) {                 //  bugfix DK 08.02.10
        act.Y_m[j] = X[j] * Factor; // molality
      }
      //            act.Y_w[j] = // mass % of the system
      //                          1e2 * X[j] * act.MM[j] / act.MBX;
      act.Y_la[j] = ln_to_lg * (Muj - act.G0[j]);
      act.FVOL[k] += act.Vol[j] * X[j];
      break;
    case DC_AQ_ELECTRON:
      act.Y_m[j] = 0.0;
      act.Y_la[j] = 0.0 - act.pe;
      //            act.Y_w[j] = 0.0;
      break;
    case DC_AQ_PROTON: // in molal scale!
      act.pH = -ln_to_lg * (Muj - act.G0[j] + lnFmol);
      [[fallthrough]];
    case DC_AQ_SPECIES:
      act.IC += 0.5 * X[j] * Factor *
                (act.EZ[j] *
                 act.EZ[j]); // increment to effective molal ionic strength
                             // strength Ieff = 0.5*sum(c_i * z_i^2)
      act.Zeff +=
          X[j] * Factor * (act.EZ[j] * act.EZ[j] * act.EZ[j] * act.EZ[j]);
      //      Increment to sum(c_i * z_i^4) of effective mean-field valency
      //      z_eff = sqrt(sum(c_i * z_i^4) / sum(c_i * z_i^2))
      [[fallthrough]];
    case DC_AQ_SURCOMP:
      SPmol = X[j] * Factor; // molality
      //           act.IC += 0.5* SPmol *(act.EZ[j]*act.EZ[j]); // Bugfix
      //           DK 21.10.2011 'K' species don't count here!
      act.FVOL[k] += act.Vol[j] * X[j]; // fixed 04.02.03 KD
      act.Y_m[j] = SPmol;
      act.Y_la[j] = ln_to_lg * (Muj - act.G0[j] + lnFmol);
      act.Y_w[j] = 1e6 * X[j] * act.MM[j] / act.FWGT[k];
      if (act.DCC[j] !=
          DC_AQ_SURCOMP) { // Bugfix 21.10.2011  DK - excluding 'K' species from
                           //            total IC molality
                           //  Optimized for performance - calculation inline
        for (i = arrL[j]; i < arrL[j + 1]; i++) {
          ii = arrAN[i];
          if (ii >= act.NR)
            continue;
          act.IC_m[ii] += SPmol * a(j, ii); // total aqueous molality
          act.IC_wm[ii] += X[j] * a(j, ii); // total aqueous mass concentration
        }
      }
      break;
    case DC_AQ_SOLVENT: // mole fractions of solvent
    case DC_AQ_SOLVCOM:
      //            act.Y_m[j] = X[j]/XFA[k];  Replaced by DK 12.03.2012
      act.Y_m[j] = H2O_mol_to_kg;
      //            act.Y_w[j] = 1e3*X[j]*act.MM[j]/act.FWGT[k];
      act.FVOL[k] += act.Vol[j] * X[j];
      act.Y_la[j] = ln_to_lg * (Muj - act.G0[j]);
      break;
    case DC_GAS_COMP:
    case DC_GAS_H2O:
    case DC_GAS_CO2: // gases
    case DC_GAS_H2:
    case DC_GAS_N2:
      act.FVOL[k] += act.Vol[j] * X[j];
      act.Y_la[j] = ln_to_lg * (Muj - act.G0[j]);
      if (act.P > 1e-9)
        act.Y_la[j] += log10(act.P);
      break;
    case DC_SOL_IDEAL:
    case DC_SOL_MINOR: // solution end member
    case DC_SOL_MAJOR:
    case DC_SOL_MINDEP:
    case DC_SOL_MAJDEP:
    case DC_SCM_SPECIES:
      act.FVOL[k] += act.Vol[j] * X[j];
      act.Y_la[j] = ln_to_lg * (Muj - act.G0[j]);
      if (act.LO) {                 //  bugfix DK 16.04.2012
        act.Y_m[j] = X[j] * Factor; // molality
      }
      break;
    case DC_SUR_GROUP:            // adsorption:
      act.Y_m[j] = X[j] * Factor; // molality
                                  //            act.Y_w[j] =  // mg/g sorbent
      //                1e3 * X[j] * act.MM[j] / (MMC*XFA[k]);
      act.Y_la[j] = ln_to_lg * (Muj - act.G0[j]); // + lnFmol ); corr. 06.10.10
      act.FVOL[k] += act.Vol[j] * X[j];           // fixed 11.03.2008 KD
      break;
    case DC_SUR_COMPLEX:
    case DC_SUR_IPAIR:
    case DC_SUR_SPECIES:
    case DC_IESC_A:
    case DC_IEWC_B:
      act.Y_m[j] = X[j] * Factor; // molality
      act.Y_w[j] =                // mg/g sorbent
          1e3 * X[j] * act.MM[j] / (MMC * act.XFA[k]);
      //           DsurT = MMC * act.Aalp[k] * pa->p.DNS*1.66054e-6;
      act.Y_la[j] =
          ln_to_lg *
          (Muj -
           act.G0[j]); // + lnFmol - act.GEX[j] + Dsur + DsurT/( 1.0+DsurT )
      act.FVOL[k] += act.Vol[j] * X[j]; // fixed 11.03.2008   KD
      break;
    case DC_PEL_CARRIER:
    case DC_SUR_MINAL:
    case DC_SUR_CARRIER:          // sorbent
      act.Y_m[j] = X[j] * Factor; // molality
                                  //            act.Y_w[j] = 0.0;
      if (act.YF[0] >= act.DSM)
        //              act.Y_w[j] = // mg of sorbent per kg aq solution
        //                1e6 * X[j] * act.MM[j] / act.FWGT[0];
        act.Y_la[j] =
            ln_to_lg *
            (Muj - act.G0[j]); // - act.GEX[j] + Dsur - 1. + 1./(1.+Dsur) -
                               // DsurT + DsurT/(1+DsurT)
      act.FVOL[k] += act.Vol[j] * X[j];
      break;
    default:
      break; // error in DC class code
    };
  } // j
  act.Zeff = sqrt(act.Zeff / (2.0 * act.IC));
  // Effective aq valency z_eff = sqrt(sum(c_i * z_i^4) / sum(c_i * z_i^2))
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of derived values (concentrations etc.) on IPM iteration
///  from X,XF, and XFA vectors. Also calculates pH, pe, Eh
// This function has to be rewritten using new set of built-in
// chemical functions.
//
void TActivity::CalculateConcentrations(double X[], double XF[], double XFA[]) {
  long int k, i, j, jj, jja, ist, ii;
  double Factor = 0.0, Dsur = 0.0, MMC = 0.0, VXc; //, YFk;

  //    if( act.Ls < 2 || !act.FIs )  Temporary disabled  09.03.2010 DK
  //        return;

  for (i = 0; i < act.N; i++)
    act.BFC[i] = 0.; // cleanup for phase_bfc() calculation

  for (j = 0; j < act.Ls; j++) {
    act.Wx[j] = 0.;
    //        act.VL[j] = 0.;
  }
  j = 0;
  VXc = 0.0;
  for (k = 0; k < act.FI; k++) { // cycle by phases
    i = j + act.L1[k];
    act.FWGT[k] = 0.0;
    act.FVOL[k] = 0.0;
    //   Dsur = 0.0;

    if (XF[k] > act.DSM && !(act.PHC[k] == PH_SORPTION && XFA[k] <= act.ScMinM))
      phase_bfc(k, j);

    if (k >= act.FIs || act.L1[k] == 1) { // this is a single- component phase
      act.Wx[j] = 1.0;                    // SD 04/05/2010
      if (XF[k] < act.DSM) {              // This phase to be zeroed off
        if (act.LO)
          act.Y_m[j] = 0.0;
        act.Wx[j] = 0.0;
        //                act.Y_w[j] = 0.0;
        act.Fx[j] =
            DC_DualChemicalPotential(act.U, act.A + j * act.N, act.NR, j);
        act.Y_la[j] = ln_to_lg * (act.Fx[j] - act.G0[j]); // -act.GEX[j]
        //                act.Fx[j] *= act.RT;     // el-chem potential
        goto NEXT_PHASE;
      }
      //            act.VL[j] = 0.0;
      if (act.LO && XFA[0] > 0)
        act.Y_m[j] = X[j] * 1000. / 18.01528 / XFA[0]; // molality
      //            act.Y_w[j] = // mass % in the system
      //                1e2 * X[j] * act.MM[j] / act.MBX;
      act.Fx[j] = DC_DualChemicalPotential(act.U, act.A + j * act.N, act.NR, j);
      act.Y_la[j] = ln_to_lg * (act.Fx[j] - act.G0[j]); // - act.GEX[j]
      //            act.Fx[j] *= act.RT;     // el-chem potential
      act.FWGT[k] += X[j] * act.MM[j];
      act.FVOL[k] += X[j] * act.Vol[j];
      goto NEXT_PHASE;
    }
    if (act.LO && !k) { // aqueous phase present
      act.IC = 0.;      // fix 20.06.13 DK
      act.Zeff = 0.;    // added 17.10.2026 DK
      for (ii = 0; ii < act.N; ii++) {
        act.IC_m[ii] = 0.0;
        act.IC_lm[ii] = 0.0;
        act.IC_wm[ii] = 0.0;
      }
    }
    if (XF[k] <= act.DSM ||
        (act.PHC[k] == PH_AQUEL &&
         (XFA[k] <= act.XwMinM || XF[k] <= act.DSM)) ||
        (act.PHC[k] == PH_SORPTION && XFA[k] <= act.ScMinM)) {
      for (jj = 0; jj < act.N; jj++)
        act.BF[k * act.N + jj] = 0.;

      for (jj = j; jj < i; jj++) // Loop added 10.03.01  KD (GTDEMO)
      {
        act.Wx[j] = 0.0;
        if (act.LO)
          act.Y_m[jj] = 0.0;
        //                act.Y_w[jj] = 0.0;
        act.Fx[jj] =
            DC_DualChemicalPotential(act.U, act.A + jj * act.N, act.NR, jj);
        act.Y_la[jj] = ln_to_lg * (act.Fx[jj] - act.G0[jj]);
        if (act.PHC[k] ==
            PH_AQUEL) // || act.PHC[k] == PH_SORPTION ) corr. 06.10.10 DK
          act.Y_la[jj] += 1.74438;
        if (act.PHC[k] == PH_GASMIX || act.PHC[k] == PH_PLASMA)
          act.Y_la[jj] += log10(act.P);
        //                act.Fx[jj] *= act.RT;     // el-chem potential
        act.lnGam[jj] = 0.0;
      }
      goto NEXT_PHASE;
    }
    // calculate bulk stoichiometry of a multicomponent phase
    phase_bcs(act.N, act.L1[k], j, act.A + j * act.N, act.X + j,
              act.BF + k * act.N);

    switch (act.PHC[k]) {
    case PH_AQUEL:
      MMC = 0.0; // molar mass of carrier
      //            Dsur = XFA[k]/XF[k] - 1.0; // Asymm.corr. - aq only!
      //            if( XFA[k] > act.lowPosNum )
      if (XFA[k] > act.XwMinM) {
        for (jj = j; jj < i; jj++)
          if (act.DCC[jj] == DC_AQ_SOLVENT || act.DCC[jj] == DC_AQ_SOLVCOM)
            MMC += act.MM[jj] * X[jj] / XFA[k];
      } else
        MMC = 18.01528; // Assuming water-solvent
      //            if( (XFA[k] > act.lowPosNum) && (MMC > act.lowPosNum) )
      if ((XFA[k] > act.XwMinM) && (MMC > act.XwMinM))
        Factor = 1000. / MMC / XFA[k]; // molality
      else
        Factor = 0.0;
      //            act.IC=0.;
      act.pe = ln_to_lg * act.U[act.N - 1];
      act.Eh = 0.000086 * act.U[act.N - 1] * act.Tc;
      [[fallthrough]];
    case PH_GASMIX:
    case PH_FLUID:
    case PH_PLASMA:
    case PH_SIMELT:
    case PH_HCARBL:
    case PH_SINCOND:
    case PH_SINDIS:
    case PH_LIQUID:
      //    case PH_IONEX:
      //   case PH_ADSORPT:
      // YFk = XF[k];
      for (jj = j; jj < i; jj++) {
        if (X[jj] > act.DcMinM) // fixed 30.08.2009 DK
          act.FWGT[k] += X[jj] * act.MM[jj];
      }
      break;
    case PH_POLYEL:
    case PH_SORPTION:   // only sorbent end-members!
      act.YFk = XFA[k]; // not XF[k] as any other phase?
      MMC = 0.0;

      for (ist = 0; ist < act.FIat; ist++)
        act.XFTS[k][ist] = 0.0;
      //           if( XFA[k] < act.lowPosNum ) XFA[k] = act.lowPosNum;
      if (XFA[k] < act.ScMinM)
        XFA[k] = act.ScMinM;
      for (jj = j; jj < i; jj++) {
        jja = jj - (act.Ls - act.Lads);
        if (act.DCC[jj] == DC_SUR_CARRIER || act.DCC[jj] == DC_SUR_MINAL ||
            act.DCC[jj] == DC_PEL_CARRIER) {
          MMC += act.MM[jj] * X[jj] / XFA[k];
          // Only sorbent mass
          act.FWGT[k] += X[jj] * act.MM[jj];
        } else { // Doubtful, do we really need it? Total amount of species on
                 // surf site type or on surface type (patch)?
          /*!!!!!*/ ist = act.SATC[jja][XL_ST]; // Check carfefully!
          act.XFTS[k][ist] += X[jj];
        }
      }
      act.FWGTS[k] =
          act.FWGT[k]; // Use for surfaces and charges only act.FWGTS[k]
      act.logYFk = log(act.YFk);
      //            Dsur = XFA[k]/XF[k] - 1.0;  // Also for sorption phases
      //            if( Dsur <= -1.0 ) Dsur = -0.999999;
      break;
    case PH_ADSSBEIL: // New: sorption phase with site balances and without
                      // carrier DC
      if (XF[k] < act.ScMinM)
        // this is in error because of site mole balance, so the amount cannot
        // be zero
        break;
      act.YFk = XF[k];
      if (act.FWGTS[k] < act.ScMinM) {
        // this is probably the first calculation for this chemical system, so
        // we get the virtual sorbent mass from site IC amounts
        act.FWGTS[k] = virtualSorbentMass(j, i, act.N);
      }
      act.FWGT[k] = act.FWGTS[k];
      act.logYFk = log(act.YFk);
      break;
    default:
      return; // Phase class code error!
    }
    // calculation of species concentrations in k-th phase
    CalculateConcentrationsInPhase(X, XF, XFA, Factor, MMC, Dsur, j, i, k);

  NEXT_PHASE:
    VXc += act.FVOL[k];
    if (act.PHC[k] == PH_AQUEL && XF[k] > act.DSM && XFA[k] > act.XwMinM)
      for (ii = 0; ii < act.NR; ii++) {
        if (act.LO) {
          if (act.IC_m[ii] >= act.DHBM)
            act.IC_lm[ii] = ln_to_lg * log(act.IC_m[ii]);
          else
            act.IC_lm[ii] = 0;
          if (act.FWGT[k] >= act.DHBM)
            act.IC_wm[ii] *= act.Awt[ii] * 1000. / act.FWGT[k];
          else
            act.IC_wm[ii] = 0;
        }
      }
    j = i;
  } // k
}

// Calculation of virtual sorbent mass (g) from site IC amounts in the system
// First, we detect from DC stoichiometry which site ICs are used in this phase,
// then calculate the sorbent mass from those IC amounts and their mol masses
double TActivity::virtualSorbentMass(long int jb, long int je, long int N) {
  long int j, i;
  std::vector<long int> nDConSite(N, 0);
  std::vector<double> Btot(N, 0.0);
  double SMASS = 0.0;
  double db = act.DHBM > 0.0 ? act.DHBM : 1e-12;
  double *b_vec = cnd ? cnd->bIC : nullptr;

  for (j = jb; j < je; j++) {
    // screening which site ICs are used in DCs in this phase
    for (i = 0; i < N; i++) {
      if (a(j, i) <= db)
        continue;
      if (act.ICC[i] == IC_SSITE ||
          act.ICC[i] ==
              IC_ADDIT) { // IC_ADDIT is provisorial and should be removed!
                          // Now kept here because of error E37FPun when
                          // opening a project caused by IC_SSITE
        nDConSite[i] += 1;
        if (Btot[i] < db && b_vec)
          Btot[i] = b_vec[i];
      }
    }
  }
  // Now calculating the mass of sorbent in the system for the sorption phase
  for (i = 0; i < N; i++) {
    if (Btot[i] < db || nDConSite[i] == 0)
      continue;
    SMASS += Btot[i] * act.Awt[i];
  }
  return SMASS;
}

/// Calculation of surface charge densities on multi-surface sorption phase
// Revised by DK on 17.10.2026 using updated input data propagation
void TActivity::IS_EtaCalc() {
  long int k, i, ist, j = 0, ja;
  double XetaS = 0., XetaW = 0., etaA, Ez, CDA, CDB, CDC;

  for (k = 0; k < act.FIs; k++) { // loop over all solution phases
    i = j + act.L1[k];
    if (act.FIat > 0 && (act.PHC[k] == PH_POLYEL || act.PHC[k] == PH_SORPTION ||
                         act.PHC[k] == PH_ADSSBEIL))
      for (ist = 0; ist < MNST; ist++) { // Replaced act.FIat with MNST
        act.XetaA[k][ist] = 0.0;         // plane charges in eq (moles)
        act.XetaB[k][ist] = 0.0;
        act.XetaC[k][ist] = 0.0;
        act.XetaD[k][ist] = 0.0; // added 12.09.05  KD
      }

    if (act.XF[k] <= act.DSM ||
        (act.PHC[k] == PH_AQUEL && (act.X[act.LO] <= act.XwMinM //  pa->p.XwMin
                                    || act.XF[k] <= act.DHBM)) ||
        (act.PHC[k] == PH_SORPTION &&
         act.XF[k] <= act.ScMinM) // no site balances, sorbent DC
        || (act.PHC[k] == PH_ADSSBEIL &&
            act.XF[k] <= act.ScMinM) // site balances, no sorbent DC
        )                            //  pa->p.ScMin) )
      goto NEXT_PHASE;

    switch (act.PHC[k]) { // initialization according to the phase class
    case PH_AQUEL:        // aqueous solution
      act.Yw = act.XFA[k];
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
        if (fabs(act.XfIEC[k][ist][FSIT_S_CAPACITY]) > act.DHBM) {
          etaA += act.XfIEC[k][ist][FSIT_S_CAPACITY] // eq/kg
                  * act.FWGT[k] / 1e3; // must be converted to eq (moles)
        } else {
          // capacity in eq/kg not given, trying to get the amount of perm
          // charge from density in 1/nm2, mass and SSA
          if (fabs(act.XfIEC[k][ist][FSIT_S_DENSITY]) > act.DHBM)
            etaA += act.XfIEC[k][ist][FSIT_S_DENSITY] * act.FWGT[k] *
                    act.Aalp[k] * 1e-18 *
                    act.MASDT[k][ist][MSDT_FSUTSA]; // convert nm2 to m2
        }
        act.XetaA[k][ist] = etaA; // Moles of permanent charge
      }
      break;
    default:
      break;
    }
    for (; j < i; j++) { // loop over DC for calculating total phase charge
      if (act.X[j] <= act.lowPosNum * 100.)
        continue; // Skipping too low concentrations
      ja = j - (act.Ls - act.Lads);

      switch (act.DCC[j]) // select expressions for species classes
      {
      case DC_AQ_ELECTRON:
      case DC_AQ_PROTON:
      case DC_AQ_SPECIES:
      case DC_AQ_SURCOMP:
        XetaW += act.X[j] * act.EZ[j];
      case DC_AQ_SOLVCOM:
      case DC_AQ_SOLVENT:
        break;
      case DC_PEL_CARRIER:
      case DC_SUR_MINAL:
      case DC_SUR_CARRIER: // charge of carrier: ???????
                           // act.XetaA[k] += act.X[j]*act.EZ[j];
        break;
      // surface species
      case DC_SUR_GROUP:
      case DC_SUR_COMPLEX:
      case DC_SUR_IPAIR:
      case DC_SUR_SPECIES:
      case DC_IESC_A:
      case DC_IEWC_B:
        // SCM Modernization: charges from mCDc[ja][0..1]
        CDA = act.mCDc ? act.mCDc[ja][CHARGE_ON_A]
                       : 0.0; // species charge that goes into 0 plane
        CDB = act.mCDc ? act.mCDc[ja][CHARGE_ON_B]
                       : 0.0; // species charge that goes into b (1) plane
        CDC = act.mCDc ? act.mCDc[ja][CHARGE_ON_C]
                       : 0.0; // species charge that goes into c (2) plane
        //
        Ez = act.EZ[j]; // take formula charge if no plane charges are provided
        if (fabs(CDA) < act.DHBM && fabs(CDB) < act.DHBM &&
            fabs(CDC) < act.DHBM)
          CDA = Ez; // Only as 0-plane charge
        // FInd out to which surface type this surface species belongs
        ist = act.SATC[ja][SURSP_ALLOC_SURFTYPE];
        act.XetaA[k][ist] += act.X[j] * CDA;
        act.XetaB[k][ist] += act.X[j] * CDB;
        act.XetaC[k][ist] += act.X[j] * CDC;
        break;
      default: // Any other phase: charge balance over the phase
        XetaS += act.X[j] * act.EZ[j];
        break;
      }
    } // j
    // Charge balance for d-plane on all surface types
    if (act.FIat > 0 && (act.PHC[k] == PH_POLYEL || act.PHC[k] == PH_SORPTION ||
                         act.PHC[k] == PH_ADSSBEIL)) {
      for (ist = 0; ist < MNST; ist++) {
        act.XetaD[k][ist] =
            0.0 - act.XetaA[k][ist] - act.XetaB[k][ist] - act.XetaC[k][ist];
      }
    }
  NEXT_PHASE:
    j = i;
    if (act.LO && !k && act.FIat > 0) {
      act.XetaA[k][0] = XetaW;
      act.XetaB[k][0] = XetaW;
      act.XetaC[k][0] = XetaW;
      act.XetaD[k][0] = XetaW;
    }
    if ((act.PHC[k] == PH_PLASMA || act.PHC[k] == PH_SIMELT) && act.FIat)
      act.XetaA[k][0] = XetaS; // Doubtful if k runs only over sorption phases?
    // it seems that in this function, it runs over all solution phases...
  } // k
    // compare act.Xetaf[k]+act.XetaA[k]+act.XetaB[k] and XetaS
    // Test XetaW
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of the surface potential pm[q].XpsiD[k] on diffuse.
/// layer plane on k-th sorption phase. From total charge act.Xeta[k]
/// ( in moles ) using Gouy-Chapman equation.
/// Strictly valid at PSI < 30 mV. Modified by DK 5 Jan 2000
/// to add a Basic Stern EIL model, then revised on 18.10.2026 to process
/// EIL models plus NEM using the same data structure
///    Added 13.03.2008 by DK: returns int value showing (if not 0)
///    that some extreme values were reached for charge densities or
///    electric potentials (for detecting bad PIA), 0 otherwise
//
long int TActivity::GouyChapman(long int, long int, long int k) {
  long int ist, status = 0;
  double sigA = 0., sigB = 0., sigC = 0., sigD = 0., sigmaD = 0., debye_m,
         x_slip;
  double SigmA[MNST], SigmB[MNST], SigmC[MNST], SigmD[MNST];
  double psiA, psiB, psiC, psiD, psiZ, convFactor, ConvFactor[MNST];
  double PsiA[MNST], PsiB[MNST], PsiC[MNST], PsiD[MNST], PsiZ[MNST];
  double f1, f3, A, Sig, F2RT, I, Zeff, Cap, rhoW, epsW, SMASSg = 0.0;
  EILModelResult emr;

  if (act.XF[k] < act.ScMinM) // TProfil::pm->pa.p.ScMin )
    return status;            // no sorbent or sorbed component amounts

  // Determining the sorbent mass in grams
  if (act.PHC[k] == PH_POLYEL || act.PHC[k] == PH_SORPTION ||
      act.PHC[k] == PH_ADSSBEIL)
    SMASSg = act.FWGTS[k];
  // newly implemented for the sorbent mass, both with or without site balances
  // Now calculating charge densities in C/m2
  for (ist = 0; ist < MNST;
       ist++) // loop over surface types // ist<act.FIat=MNST
  {
    sigA = 0.0, sigB = 0.0, sigC = 0.0, sigD = 0.0;
    convFactor = 1.;
    psiA = 0.0, psiB = 0.0, psiC = 0.0, psiD = 0.0, psiZ = 0.0; // Cleanup

    // SCM Modernization: SCMC[k][ist][0] and MASDT [k][ist][MSDT_FSUTSA]
    if (act.SCMC[k][ist][0] == SC_NOT_USED ||
        act.MASDT[k][ist][MSDT_FSUTSA] < 1e-9)
      continue; // potentially error
    convFactor = F_CONSTANT / (act.FWGTS[k] * act.Aalp[k] *
                               act.MASDT[k][ist][MSDT_FSUTSA]); // C/mol/(m2)
    // Calculation of charge densities (now limited to total charges > max.
    // balance residual)
    if (fabs(act.XetaA[k][ist]) > act.DHBM)  // moles (eq)
      sigA = act.XetaA[k][ist] * convFactor; // result in C/m2
    if (fabs(act.XetaB[k][ist]) > act.DHBM)  // moles (eq)
      sigB = act.XetaB[k][ist] * convFactor; // result in C/m2
    if (fabs(act.XetaC[k][ist]) > act.DHBM)  // moles
      sigC = act.XetaC[k][ist] * convFactor; // result in C/m2
    if (fabs(act.XetaD[k][ist]) > act.DHBM)  // moles
      sigD = act.XetaD[k][ist] * convFactor; // C/m2

    // Limit maximum charge densities to prevent GEM divergence
    if (fabs(sigA) > 1.4) {
      //      ipm_logger->debug(
      //          "EIL charge density plane A {} truncated to +- 1.4 C/m2  "
      //          "IT={} k= {} ist= {}",
      //          sigA, act.IT, k, ist);
      SigmA[ist] = sigA < 0.0 ? -1.4 : 1.4;
      status = 60;
    } else
      SigmA[ist] = sigA;
    if (fabs(sigB) > 1.0) {
      //      ipm_logger->debug(
      //          "EIL charge density plane B {} truncated to +- 1.0 C/m2  "
      //          "IT={} k= {} ist= {}",
      //          sigB, act.IT, k, ist);
      SigmB[ist] = sigB < 0.0 ? -2.0 : 2.0;
      status = 61;
    } else
      SigmB[ist] = sigB;
    if (fabs(sigC) > 1.0) {
      //      ipm_logger->debug(
      //          "EIL charge density plane C {} truncated to +- 1.0 C/m2  "
      //          "IT={} k= {} ist= {}",
      //          sigC, act.IT, k, ist);
      SigmC[ist] = sigC < 0.0 ? -2.0 : 2.0;
      status = 62;
    } else
      SigmC[ist] = sigC;
    if (fabs(sigD) > 0.7) {
      //     ipm_logger->debug(
      //         "EIL charge density plane D {} truncated to +- 0.7 C/m2  "
      //         "IT={} k= {} ist= {}",
      //         sigD, act.IT, k, ist);
      SigmD[ist] = sigD < 0.0 ? -1.4 : 1.4;
      status = 63;
    } else
      SigmD[ist] = sigD;

    F2RT = act.FRT / 2.;
    I = act.IC;
    Zeff = act.Zeff;
    rhoW = act.denW[0];
    epsW = act.epsW[0];
    // Calculation of Debye length in metres
    debye_m = calc_debye_length(I);
    x_slip = act.XcapZ[k][ist] * debye_m; // spli plane shift factor f * debye_m

    // Calculating EIL models for this surface type (patch)
    switch (act.SCMC[k][ist][SCMC_EILM_CODE]) {
    case SC_NEM: // Non-electrostatic, charge densities on zero plane
      // calculated, but potentials and Coulomb terms are zeros
      emr = compute_nem(SC_NEM, SigmA[ist]);
      break;
    case SC_4LM: // Four-layer model
      emr = compute_flm(SC_4LM, I, Zeff, act.XcapA[k][ist], act.XcapB[k][ist],
                        act.XcapC[k][ist], SigmA[ist], SigmB[ist], SigmC[ist],
                        x_slip, debye_m);
      break;
    case SC_3LM: // Three-layer EIL model for CD-MUSIC
      emr = compute_cd_music(SC_3LM, I, Zeff, act.XcapA[k][ist],
                             act.XcapB[k][ist], SigmA[ist], SigmB[ist], x_slip,
                             debye_m);
      break;
    case SC_TLM: // Triple-layer model
      emr = compute_tlm(SC_TLM, I, Zeff, act.XcapA[k][ist], act.XcapB[k][ist],
                        SigmA[ist], SigmB[ist], x_slip, debye_m);
      break;
    case SC_BSM: // Basic Stern model
      emr = compute_bsm(SC_BSM, I, Zeff, act.XcapA[k][ist], SigmA[ist], x_slip,
                        debye_m);
      break;
    case SC_CCM: // Constant-capacitance model
      emr = compute_ccm(SC_CCM, act.XcapA[k][ist], SigmA[ist]);
      break;
    case SC_DDL: // Double-layer model
      emr = compute_ddl(SC_DDL, I, Zeff, SigmA[ist], x_slip);
      break;
    case SC_TCE: // Tournassat 2013 clay edge potential model with charge
                 // spillover TBD
      emr = compute_tce(SC_TCE, I, Zeff, act.XcapA[k][ist], act.XcapB[k][ist],
                        SigmA[ist], SigmB[ist], act.MASDT[k][ist], x_slip,
                        debye_m);
      break;
    default:
      break;
    }
    // Unpack results from EILModelResult emr fields to MULTI data for this
    // phase and surface type
    act.XpsiA[k][ist] = emr.psi0;
    act.XpsiB[k][ist] = emr.psi1;
    act.XpsiC[k][ist] = emr.psi2;
    act.XpsiD[k][ist] = emr.psid;
    act.XpsiZ[k][ist] = emr.zeta;  // zeta potential
    act.XcapD[k][ist] = emr.cap_d; // Eff. capacitance density of diffuse layer
    act.XdlA[k][ist] =
        emr.sld01; // Thickness of inner Stern layer (0 to 1 plane) nm
    act.XdlB[k][ist] =
        emr.sld02; // Thickness of outer Stern layer (1 to 2 plane) nm
    act.XdlC[k][ist] =
        emr.sld0d; // Eff. thickness of Stern layer (0 to 3-plane or d-plane)
    act.XdlD[k][ist] = emr.dlth; // Eff. thickness of diffuse layer (nm)
                                 // or 2*Debye_length ?
  } // ist   end of loop over surface types

  return status;
}
// Added on 18.10.2026  by DK in ClaySorEl project
// Dynamic relative permittivity of water via Malmberg & Maryott (1956)
// polynomial.
double TActivity::calc_water_permittivity(double T_kelvin, double P_bar) const {
  double T_C = T_kelvin - 273.15;
  double epsW1, epsW2, epsW;
  epsW1 = 87.740 - 0.40008 * T_C + 9.398e-4 * (T_C * T_C) -
          1.410e-6 * (T_C * T_C * T_C);
  // for P > 1 bar, get water permittivity from SUPCRT
  if (P_bar > 1.0) {
    epsW2 = act.epsW[0];
    epsW = epsW2;
  } else
    epsW = epsW1;
  return epsW;
}

// Calculates Debye screening length (kappa^-1) in meters from ionic strength I
// (mol/L).
double TActivity::calc_debye_length(double I) const {
  double EPS_0 = 8.8541878128e-12;
  double FARADAY = 96485.3321;
  double RT = act.RT;
  double I_m3 = I * 1000.0;
  //  double Rho = act.FWGT[0] / act.FVOL[0]; // density of aq phase
  double rhoW = act.denW[0]; // density of water at T,P
  double epsW =
      calc_water_permittivity(act.Tc, act.Pc); // permittivity of water at T,P
  double debye_m;                              // Debye length in m
  debye_m =
      std::sqrt((EPS_0 * epsW * RT) / (2.0 * FARADAY * FARADAY * I_m3 * rhoW));
  // debye_m = std::sqrt((EPS_0 * epsW * RT) /
  //                 (2.0 * FARADAY * FARADAY * I_m3 * Rho));
  return debye_m;
}

// Calculates capacitance density of the diffuse layer from ionic strength I
//
double TActivity::calc_cap_d_from_sigma_d(double sigma_d, double I) const {
  double I_m3 = I * 1000.0;
  double EPS_0 = 8.8541878128e-12;
  double rhoW = act.denW[0];
  double epsW =
      calc_water_permittivity(act.Tc, act.Pc); // permittivity of water at T,P
  double RT = act.RT;
  double Sig = sigma_d;
  double A = 1.0;
  if (I > 1e-7) {
    double Rho = 1.0;
    if (act.FVOL[0] > 1e-16)
      Rho = act.FWGT[0] / act.FVOL[0];
    A = sqrt(2.0 * EPS_0 * epsW * act.RT * I_m3 * Rho); // or Rho = rhoW ?
    // Aq solution density Rho included acc. to [Machesky ea., 1999]
  }
  double Cap_d = act.FRT / 2.0 * sqrt(4. * A * A + Sig * Sig);
  return Cap_d;
}

// Inverted Gouy-Chapman relation generalized for mixed aqueous electrolytes
// with effective valency z_eff.
double TActivity::calc_psi_d_from_sigma_d(double sigma_d, double I,
                                          double z_eff) const {
  double I_m3 = I * 1000.0;
  double EPS_0 = 8.8541878128e-12;
  double RT = act.RT;
  double F_RT = act.FRT; // FARADAY / (R_GAS * T);
  double rhoW = act.denW[0];
  double eps_r = calc_water_permittivity(act.Tc, act.Pc);
  double gc_coeff = std::sqrt(8.0 * EPS_0 * eps_r * RT * I_m3 * rhoW);
  double arg = sigma_d / (-gc_coeff);
  double psiD = (2.0 / (z_eff * F_RT)) * std::asinh(arg);
  // Truncating diffuse plane potential to avoid divergence of GEM
  if (fabs(psiD) > 0.5) {
    //   ipm_logger->debug(
    //       "All EIL models: PsiD ={} truncated to +- 0.5 V IT={} k= {} ist=
    //       {}", psi_d, act.IT, k, ist);
    psiD = psiD < 0.0 ? -0.5 : 0.5;
    //    status = 64;
  }
  return psiD;
}

// Non-linear Poisson-Boltzmann potential decay through the diffuse layer at a
// given distance x_dist.
double TActivity::calc_gc_potential_at_distance(double psi_d, double I,
                                                double x_dist,
                                                double z_eff) const {
  if (std::abs(psi_d) < 1e-7)
    return 0.0;
  // double F_RT = FARADAY / (R_GAS * T);
  double F_RT = act.FRT;
  double gamma_0 = std::tanh(z_eff * F_RT * psi_d / 4.0);
  double debye_len = calc_debye_length(I);
  double ratio = gamma_0 * std::exp(-x_dist / debye_len);
  ratio = std::clamp(ratio, -0.9999999, 0.9999999);
  double psi_x = (4.0 / (z_eff * F_RT)) * std::atanh(ratio);
  return psi_x;
}

EILModelResult TActivity::compute_flm(char ModCode, double I, double z_eff,
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

EILModelResult TActivity::compute_tlm(char ModCode, double I, double z_eff,
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

EILModelResult TActivity::compute_cd_music(char ModCode, double I, double z_eff,
                                           double C1, double C2, double sig0,
                                           double sig1, double x_slip,
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

EILModelResult TActivity::compute_bsm(char ModCode, double I, double z_eff,
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

EILModelResult TActivity::compute_ddl(char ModCode, double I, double z_eff,
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

EILModelResult TActivity::compute_ccm(char ModCode, double C1,
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

EILModelResult TActivity::compute_nem(char ModCode, double sig0) const {
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

EILModelResult TActivity::compute_tce(char ModCode, double I, double z_eff,
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
///  SACT are placed into act.lnGam[j], as other activity coefficients except
///  relative surface potentials (Coulombic terms) kept separately.
//
///  Old (obsolete) SAT calculations (Kulik 2000, 2002) are retained.
//
///  act.lnSAC[*][3] vector is now used to keep original DUL[j] to restore
///  them after IPM-2 refinements for surface complexes.
///    Added 13.03.2008 by DK: returns int value showing (if true)
///    that some extreme values were obtained for some SACTs,
///    0 otherwise (for detecting bad PIA)
//
long int TActivity::SurfaceActivityCoeff(long int jb, long int je, long int,
                                         long int, long int k) {
  long int status = 0;
  long int i, ii, j, ja, ist = 0, iss, dentt, dentsn, Cj,
                         iSite[MNST]; // SCM Modernization: MNST
  double XS0, xj0, XVk, XSk, XSkC, xj, Mm, rIEPS, ISAT, XSs, SATst, xjn, q1, q2,
      GamC, alpF, qBet, bBet, cN, eF, lnGamjo, lnDiff, lnFactor;
  // const BASE_PARAM *pa_p = base_param();

  if (act.PHC[k] == PH_SORPTION || act.PHC[k] == PH_POLYEL) {
    if (act.XF[k] <= act.DSM) // No sorption phase retained by the IPM
      return status;
    if (act.XFA[k] <= act.ScMinM) // no amount of sorbent remains
      return status;
  } else {
    // SACTS do not apply for sorption phases with site balances [yet]
    if (act.PHC[k] == PH_ADSSBEIL || act.PHC[k] == PH_ADSORPT)
      return status;
  }
  // This needs a deeper reconsideration!
  // Initialization of work arrays: D[NSURFSITES][MNST]
  for (i = 0; i < MNST; i++) // SCM Modernization: loop over MNST surf patches
  {
    iSite[i] = -1;
    for (ii = 0; ii < NSURFSITES; ii++) // SCM Modernization: site types
      act.D[ii][i] = 0.0; // cleaning the matrix for sites totals in surface
  }

  // Extraction of site indices for the neutral >OH group
  for (j = jb; j < je; j++) {
    ja = j - (act.Ls - act.Lads);
    if (act.SATT[ja] == SAT_SOLV) {
      ist = act.SATC[ja][XL_ST]; // / MSPN;
      iSite[ist] = j;
    }
    // Counting current sites totals (these species do not compete)
    if (act.DCC[j] == DC_SUR_CARRIER || act.DCC[j] == DC_SUR_MINAL ||
        act.DCC[j] == DC_PEL_CARRIER || act.SATT[ja] == SAT_SOLV)
      continue;
    // Get ist (index of surface type)
    ist = act.SATC[ja][XL_ST];  // / MSPN;
    if (ist < 0 || ist >= MNST) // SCM Modernization: MNST
      ist = 0;                  // default: zero surface type
    // Get iss (index of site on surface type)
    iss = act.SATC[ja][XL_SI];
    if (iss < 0 || iss >= MNST) // SCM Modernization: MNST
      iss = 0; // default: zero site is the weekest and the most abundant one?
    act.D[iss][ist] += act.X[j]; // adding to total amount on a site
  }

  for (j = jb; j < je; j++) { // Main loop for DCs - surface complexes
    lnGamjo = act.lnGmo[j];   // bugfix 16.03.2008 DK
    if (act.X[j] < std::min(act.DcMinM, act.lowPosNum))
      continue; // This surface DC has been killed by IPM
                //        OSAT = act.lnGmo[j]; // added 6.07.01 by KDA
    ja = j - (act.Ls - act.Lads); // Assuming that sorption phase(s) are the
                                  // last ones in the list of phases-solutions
    rIEPS = act.IEPS;   // pa_p->IEPS;      // default 1e-3 (for old SAT - 1e-9)
                        //         dent = 1;  // default - monodentate
    switch (act.DCC[j]) // code of species class
    {
    default: // act.lnGam[j] = 0.0;  not a surface species
      continue;
    case DC_IEWC_B:
    case DC_IESC_A:
    case DC_SUR_GROUP:
    case DC_SUR_COMPLEX:
    case DC_SUR_IPAIR:
    case DC_SUR_SPECIES:
      // Get ist (index of surface type)
      ist = act.SATC[ja][XL_ST]; // / MSPN;
      // get iss - index of site on surf.type
      iss = act.SATC[ja][XL_SI];
      if (iss < 0 || iss >= MNST) // SCM Modernization: MNST
        iss = 0; // zero site is the weekest and the most abundant one
      // Cj - index of carrier DC
      Cj = act.SATC[ja][XL_EM];
      if (Cj < 0) { // Assigned to the whole sorbent
        XVk = act.XFA[k];
        Mm = act.FWGT[k] / XVk;
      } else { // Assigned to one of the sorbent end-members
        XVk = act.X[Cj];
        if (XVk < act.ScMinM) // act.DSM*0.1 )
          continue;           // This end-member is zeroed off by IPM
        Mm = act.MM[Cj] * XVk / act.XFA[k]; // mol.mass
      }
      XSk = act.XFTS[k][ist]; // Total moles of sorbates on ist surface type
      XSs = act.D[iss][ist];  // Total moles of SC on iss site type
      xj = act.X[j];          // Current moles of this surface species

      // Extracting SACT isotherm parameters
      if (act.MASDJ) {
        cN = act.MASDJ[ja][PI_P2];   // Frumkin/Pivovarov water coord. number
        dentt = denticity_total(ja); // total denticity of this surf species
        dentsn = denticity_num_site_types(ja); // Number of site types to bind
        GamC = act.MASDJ[ja][PI_DEN]; // Complete coverage density (individual,
                                      // or as set for the whole surf type)
        alpF = act.MASDJ[ja][PI_P1]; // Frumkin lateral interaction energy alpha
        qBet = act.MASDJ[ja][PI_P2]; // q parameter (QCA/BET) or cN in Pivovarov
        bBet = act.MASDJ[ja][PI_P3]; // Parameter b in BET or S in Pivovarov
      } else {
        cN = 0.0;
        alpF = 1.0;
        qBet = 0.0;
        bBet = 0.0;
        dentt = 1; // bet = 1.0;
        dentsn = 1;
      }
      switch (act.SATT[ja]) // selection of the SACT model
      {
      case SACT_LAN_NCOMP: // Non-competitive individual monodentate Langmuir
                           // TBD
        break;
      case SACT_LAN_COMP: // Competitive monodentate Langmuir on a surface type
                          // and site type
        XSkC = XSs / XVk / Mm * 1e6 / act.MASDT[k][ist][MSDT_FSUTSA] /
               act.Aalp[k] / 1.66054; // per nm2
        XS0 = (fabs(act.MASDJ[ja][PI_DEN]) / act.Aalp[k] / 1.66054);
        // max. density per nm2
        if (act.PC <= 2)
          rIEPS = act.IEPS * XS0; // relative IEPS
        if (XSkC < 0.0)
          XSkC = 0.0;
        if (XSkC >= XS0) // Setting limits
          XSkC = XS0 - 2.0 * rIEPS;
        q1 = XS0 - XSkC;
        if ((act.PC == 3 && !act.W1) || act.PC != 3) {
          q2 = rIEPS * XS0;
          if (q1 > q2)
            q2 = q1;
        } else {
          q2 = q1;
          if (q2 <= 1e-22)
            q2 = 1e-22;
        }
        ISAT = log(XS0) - log(q2);
        act.lnGam[j] = ISAT;
        act.lnSACT[j] = ISAT;
        act.lnSAC[ja][0] = ISAT;
        break;
      case SACT_QCA_NCOMP: // dent++;   bidentate is default for QCA
        xj0 = (fabs(act.MASDJ[ja][PI_DEN]) / act.Aalp[k] / 1.66054);
        // Max site density per nm2
        xj = XSs / XVk / Mm / act.MASDT[k][ist][MSDT_FSUTSA] * 1e6 // xj
             / act.Aalp[k] / 1.66054; // Density per nm2 on site type iss
        if (act.PC <= 2)
          rIEPS = act.IEPS * xj0; // relative IEPS
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
        act.lnGam[j] = ISAT;
        act.lnSACT[j] = ISAT;
        act.lnSAC[ja][0] = ISAT;
        break;
        // (Non)competitive QCA-L for 1 to 4 dentate species
      case SACT_QCA1_NCOMP: // dentt = 4;     // code '4'
      case SACT_QCA2_NCOMP: // dentt = 3;     // code '3'
      case SACT_QCA3_NCOMP: // dentt = 2;     // code '2'
      case SACT_QCA4_NCOMP: // dentt = 1;     // code '1'
        break;
      case SACT_FFG_COMP: // Competitive Frumkin (FFG) isotherm
        dentt = 1;        // monodentate for now
        xj0 = (fabs(act.MASDJ[ja][PI_DEN]) / act.Aalp[k] / 1.66054);
        // Max site density per nm2
        xj = XSs / XVk / Mm / act.MASDT[k][ist][MSDT_FSUTSA] * 1e6 //  xj
             / act.Aalp[k] / 1.66054; // Current species density per nm2
        if (act.PC <= 2)
          rIEPS = act.IEPS * xj0; // relative IEPS
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
        act.lnGam[j] = ISAT + eF;
        act.lnSACT[j] = ISAT;
        act.lnGammF[j] = eF;
        act.lnSAC[ja][0] = ISAT;
        act.lnSAC[ja][1] = eF;
        break;
      case SACT_FFG_NCOMP: // Non-Competitive Frumkin isotherm
                           // for permanent charge surfaces
        XSkC = xj / XVk / Mm / act.MASDT[k][ist][MSDT_FSUTSA] * 1e6 /
               act.Aalp[k] / 1.66054; // per nm2
        XS0 = (act.MASDJ[ja][PI_DEN] / act.Aalp[k] / 1.66054);
        // max.dens.per nm2
        if (act.PC <= 2)
          rIEPS = act.IEPS * XS0; // relative IEPS
        if (XSkC < 0.0)
          XSkC = 0.0;
        if (XSkC >= XS0) // Limits
          XSkC = XS0 - 2.0 * rIEPS;
        q1 = XS0 - XSkC;
        if ((act.PC == 3 && !act.W1) || act.PC != 3) {
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
        act.lnGam[j] = ISAT + eF;
        act.lnSACT[j] = ISAT;
        act.lnGammF[j] = eF;
        act.lnSAC[ja][0] = ISAT;
        act.lnSAC[ja][1] = eF;
        break;
      case SACT_PIV_NCOMP: // (Non)competitive Pivovarov isotherm
        dentt = 1;         // monodentate for now
        xj0 = (fabs(act.MASDJ[ja][PI_DEN]) / act.Aalp[k] / 1.66054);
        // Max site density per nm2
        xj = XSs / XVk / Mm / act.MASDT[k][ist][MSDT_FSUTSA] * 1e6 /
             act.Aalp[k] / 1.66054; // Current density per nm2
        if (act.PC <= 2)
          rIEPS = act.IEPS * xj0; // relative IEPS
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
        act.lnGam[j] = ISAT + eF;
        act.lnSACT[j] = ISAT;
        act.lnGammF[j] = eF;
        act.lnSAC[ja][0] = ISAT;
        act.lnSAC[ja][1] = eF;
        break;
      case SACT_VIR_NCOMP: // Non-Competitive virial isotherm
        dentt = 1;         // monodentate for now
        xj0 = (fabs(act.MASDJ[ja][PI_DEN]) / act.Aalp[k] / 1.66054);
        // Max site density per nm2
        xj = XSs / XVk / Mm / act.MASDT[k][ist][MSDT_FSUTSA] * 1e6 /
             act.Aalp[k] / 1.66054; // Current density per nm2
        if (act.PC <= 2)
          rIEPS = act.IEPS * xj0; // relative IEPS
        if (xj >= xj0 / dentt)
          xj = xj0 / dentt - rIEPS; // upper limit
        ISAT = 0.0;
        if (fabs(alpF) < 1e-9 || fabs(cN) < 1e-9)
          eF = 0.0;
        else { // changed from cN to -cN on 16.09.2005 by KD
          eF = -cN * alpF * xj / xj0; // Fi = Fi'/(kT) Bockris p.938
        }
        act.lnGam[j] = ISAT + eF;
        act.lnSACT[j] = ISAT;
        act.lnGammF[j] = eF;
        act.lnSAC[ja][0] = ISAT;
        act.lnSAC[ja][1] = eF;
        break;
      case SACT_FREU_NCOMP: // Non-competitive Freundlich isotherm TBD
        break;
      case SACT_BET_NCOMP: // Non-competitive BET for surface precipitation
        ISAT = 0.0;
        eF = 0.0;
        // To be completed
        //
        //
        act.lnGam[j] = ISAT + eF;
        act.lnSACT[j] = ISAT;
        act.lnGammF[j] = eF;
        act.lnSAC[ja][0] = ISAT;
        break;
        // TBD: cases for surface type isotherms with species bound to different
        // site
        // types
      case SACT_INDEF: // No SACT calculation whatsoever
        act.lnGam[j] = 0.0;
        act.lnSACT[j] = 0.0;
        act.lnGammF[j] = 0.0;
        act.lnSAC[ja][0] = 0;
        break;
      default: // act.lnGam[j] = 0.0;
        break;
      }
      // Make optional ....
      switch (act.SATT[ja]) // selection of the old SAT model
      {
      //  Obsolete old SAT calculations (retained from previous versions)
      case SAT_COMP: // Competitive SAT (obsolete) on a surface type
        if (iSite[ist] < 0)
          xjn = 0.0;
        else
          xjn = act.X[iSite[ist]]; // neutral site does not compete!
        // SCM Modernization: MASDT patch fraction and density
        XS0 = fabs(act.MASDJ[ja][PI_DEN]) * XVk * Mm / 1e6 *
              act.MASDT[k][ist][MSDT_FSUTSA]; // expected total in moles
        if (act.PC <= 2)
          rIEPS = act.IEPS * XS0; // relative IEPS
        XSkC = XSk - xjn - xj;    // occupied by the competing species;
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
          if ((act.PC == 3 && !act.W1) || act.PC != 3) {
            if (q1 > q2)
              q2 = q1;
          } else {
            q2 = q1;
            if (q2 <= 1e-33)
              q2 = 1e-33;
          }
          ISAT = log(xj) - log(q2);
        }
        act.lnGam[j] = ISAT;
        act.lnSAC[ja][0] = ISAT;
        break;
      case SAT_NCOMP: // Non-competitive truncated Langmuir SAT (obsolete)
        // rIEPS = pa->p.IEPS * 2;
        xj0 = fabs(act.MASDJ[ja][PI_DEN]) * XVk * Mm / 1e6 *
              act.MASDT[k][ist][MSDT_FSUTSA]; // in moles
        if (act.PC <= 2)
          rIEPS = act.IEPS * xj0; // relative IEPS
        if (xj >= xj0)
          xj = xj0 - rIEPS;  // upper limit
        if (xj * 2.0 <= xj0) // Linear adsorption - to improve !
          ISAT = 0.0;
        else {
          q1 = xj0 - xj; // limits: rIEPS to 0.5*xj0
          q2 = xj0 * rIEPS;
          if (act.PC == 3 && act.W1)
            ISAT = log(xj) - log(q1);
          else {
            if (q1 > q2)
              ISAT = log(xj) - log(q1);
            else
              ISAT = log(xj) - log(q2);
          }
        }
        act.lnGam[j] = ISAT;
        act.lnSAC[ja][0] = ISAT;
        break;
      case SAT_SOLV: // Neutral surface site (e.g. >O0.5H@ group)
                     // applies to the whole surface type!
        XSs = 0.0;   // calc total moles on all sites on surface type
        for (i = 0; i < MNST; i++)
          XSs += act.D[i][ist];
        XSkC = XSs / XVk / Mm * 1e6 // total non-solvent surf.species
               / act.MASDT[k][ist][MSDT_FSUTSA] / act.Aalp[k] /
               1.66054; // per nm2
        XS0 = fabs(act.MASDJ[ja][PI_DEN]);
        SATst = act.DNS * 1.66054 * act.Aalp[k] / XS0;
        XS0 = XS0 / act.Aalp[k] / 1.66054;
        if (act.PC <= 2)
          rIEPS = act.IEPS * XS0; // relative IEPS
        if (XSkC < 0.0)
          XSkC = 0.0;
        if (XSkC >= XS0) // Limits
          XSkC = XS0 - 2.0 * rIEPS;
        q1 = XS0 - XSkC;
        if ((act.PC == 3 && !act.W1) || act.PC != 3) {
          q2 = rIEPS * XS0;
          if (q1 > q2)
            q2 = q1;
        } else {
          q2 = q1;
          if (q2 <= 1e-22)
            q2 = 1e-22;
        }
        ISAT = log(XS0) - log(q2);
        act.lnGam[j] = ISAT + log(SATst);
        act.lnSAC[ja][0] = ISAT;
        act.lnSAC[ja][1] = log(SATst);
        break;
      }
    }

    if (lnGamjo > act.DHBM) { // Workaround DK 07.12.2009
      lnFactor = 0.2;
      lnDiff = act.lnGam[j] - lnGamjo;
      if (fabs(lnDiff) > fabs(lnFactor)) // e times
      {
        if (fabs(lnDiff) > 6.907755) // 1000 times
          status = 101;              // the SACT has changed too much;
                        // the threshold pa->p.IEPS needs adjustment!
        // Smoothing (provisional)
        act.lnGam[j] = lnDiff > 0 ? lnGamjo + lnDiff - lnFactor
                                  : lnGamjo + lnDiff + lnFactor;
      }
    }
  } // j
  return status;
}

// Total denticity of surface species (binding to one or more site types)
long int TActivity::denticity_total(long int j) {
  long int i, dent_total = 0; //, dent_types = 0;
  for (i = 0; i < NSURFSITES; i++) {
    //   if (act.XSMD[j][i] > 0)
    //     dentypes++;
    dent_total += act.XSMD[j][i];
  }
  return dent_total;
}

// Number of different site types binding to this surface species
long int TActivity::denticity_num_site_types(long int j) {
  long int i, dent_types = 0; // , dent_total = 0;
  for (i = 0; i < NSURFSITES; i++) {
    if (act.XSMD[j][i] > 0)
      dent_types++;
    //    dent_total += act.XSMD[j][i];
  }
  return dent_types;
}
//--------------------- End of s_activity2.cpp ---------------------------
