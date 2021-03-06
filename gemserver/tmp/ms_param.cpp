//-------------------------------------------------------------------
// $Id$
/// \file ms_param.cpp
/// Implementation  of default settings for the Interior Points Method
/// (IPM) module for convex programming Gibbs energy minimization
/// as well as I/O functions for IPM files
//
// Copyright (c) 1992,2012 S.Dmytriyeva, D.Kulik, K.Chudnenko
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



#include  <sstream>
#include <math.h>
//#include "m_param.h"
#include "num_methods.h"
#include "gdatastream.h"
#include "node.h"

#ifdef IPMGEMPLUGIN

//TProfil* TProfil::pm;

const double R_CONSTANT = 8.31451,
              NA_CONSTANT = 6.0221367e23,
                F_CONSTANT = 96485.309,
                  e_CONSTANT = 1.60217733e-19,
                    k_CONSTANT = 1.380658e-23,
// Conversion factors
                      cal_to_J = 4.184,
                        C_to_K = 273.15,
                          lg_to_ln = 2.302585093,
                            ln_to_lg = 0.434294481,
                             H2O_mol_to_kg = 55.50837344,
                               Min_phys_amount = 1.66e-24;

//enum volume_code {  // Codes of volume parameter
//    VOL_UNDEF, VOL_CALC, VOL_CONSTR
//};

SPP_SETTING pa_ = {
    " GEMS-GUI v.3.1 r.2150 (rc) " " GEMS3K v.3.1 r.690 (rc) ",
    {    // Typical default set (03.04.2012) new PSSC( logSI ) & uDD()
         2,  /* PC */  2,     /* PD */   -5,   /* PRD */
         1,  /* PSM  */ 130,  /* DP */   1,   /* DW */
         0, /* DT */     30000,   /* PLLG */   1,  /* PE */  7000, /* IIM */
         1000., /* DG */   1e-13,  /* DHB */  1e-20,  /* DS */
         1e-6,  /* DK */  0.01,  /* DF */  0.01,  /* DFM */
         1e-5,  /* DFYw */  1e-5,  /* DFYaq */    1e-5,  /* DFYid */
         1e-5,  /* DFYr,*/  1e-5,  /* DFYh,*/   1e-5,  /* DFYc,*/
         1e-6, /* DFYs, */  1e-17,  /* DB */   1.,   /* AG */
         0.,   /* DGC */   1.0,   /* GAR */  1000., /* GAH */
         1e-3, /* GAS */   12.05,  /* DNS */   1e-13,  /* XwMin, */
         1e-13,  /* ScMin, */  1e-33, /* DcMin, */   1e-20, /* PhMin, */
         1e-5,  /* ICmin */   1e-10,  /* EPS */   1e-3,  /* IEPS */
         1e-10,  /* DKIN  */ 0,  /* tprn */
    },
}; // SPP_SETTING



void BASE_PARAM::write(fstream& oss)
{
  short arr[10];

  arr[0] = PC;
  arr[1] = PD;
  arr[2] = PRD;
  arr[3] = PSM;
  arr[4] = DP;
  arr[5] = DW;
  arr[6] = DT;
  arr[7] = PLLG;
  arr[8] = PE;
  arr[9] = IIM;

  oss.write( (char*)arr, 10*sizeof(short) );
  oss.write( (char*)&DG, 28*sizeof(double) );
  oss.write( (char*)&tprn, sizeof(char*) );
}

void SPP_SETTING::write(fstream& oss)
{
    oss.write( ver, TDBVERSION );
    p.write( oss );
}

TProfil::TProfil( TMultiBase* amulti )
{
    pa= pa_;
    multi = amulti;
    pmp = multi->GetPM();
}


/// GEM IPM calculation of equilibrium state in MULTI
double TProfil::ComputeEquilibriumState( /*long int& RefinLoops_,*/ long int& NumIterFIA_, long int& NumIterIPM_ )
{
   //long int RefineLoops = RefinLoops_;
   //RefineLoops = 0L; // Provisional
   return multi->CalculateEquilibriumState( /*RefineLoops,*/ NumIterFIA_, NumIterIPM_ );
}

///// Writing structure MULTI (GEM IPM work structure) to binary file
//void TProfil::outMulti( GemDataStream& ff, std::string& /*path*/  )
//{
//	 short arr[10];

//	  arr[0] = pa.p.PC;
//	  arr[1] = pa.p.PD;
//	  arr[2] = pa.p.PRD;
//	  arr[3] = pa.p.PSM;
//	  arr[4] = pa.p.DP;
//	  arr[5] = pa.p.DW;
//	  arr[6] = pa.p.DT;
//	  arr[7] = pa.p.PLLG;
//	  arr[8] = pa.p.PE;
//	  arr[9] = pa.p.IIM;

//	ff.writeArray( arr, 10 );
//    ff.writeArray( &pa.p.DG, 28 );
//    multi->to_file( ff );
//}

///// Writing structure MULTI ( free format text file  )
//void TProfil::outMultiTxt( const char *path, bool append  )
//{
//    multi->to_text_file( path, append );
//}

///// Reading structure MULTI (GEM IPM work structure) from binary file
//void TProfil::readMulti( GemDataStream& ff, DATACH  *dCH )
//{
//     short arr[10];

//	 ff.readArray( arr, 10 );
//	  pa.p.PC = arr[0];
//	  pa.p.PD = arr[1];
//	  pa.p.PRD = arr[2];
//	  pa.p.PSM = arr[3];
//	  pa.p.DP = arr[4];
//	  pa.p.DW = arr[5];
//	  pa.p.DT = arr[6];
//	  pa.p.PLLG = arr[7];
//	  pa.p.PE = arr[8];
//	  pa.p.IIM = arr[9];

//      ff.readArray( &pa.p.DG, 28 );
//      multi->from_file( ff );

//      // copy intervals for minimizatiom
//      if(  dCH->nPp > 1  )
//      {
//         pmp->Pai[0] = dCH->Pval[0];
//         pmp->Pai[1] = dCH->Pval[dCH->nPp-1];
//         pmp->Pai[2] = (pmp->Pai[1]-pmp->Pai[0])/(double)dCH->nPp;
//      }
//      pmp->Pai[3] = dCH->Ptol;
//      if(  dCH->nTp > 1  )
//      {
//         pmp->Tai[0] = dCH->TKval[0];
//         pmp->Tai[1] = dCH->TKval[dCH->nTp-1];
//         pmp->Tai[2] = (pmp->Tai[1]-pmp->Tai[0])/(double)dCH->nTp;
//      }
//      pmp->Tai[3] = dCH->Ttol;

//  }




//-------------------------------------------------------------------------
// internal functions

void strip(string& str)
{
  string::size_type pos1 = str.find_first_not_of(' ');
  string::size_type pos2 = str.find_last_not_of(' ');
  str = str.substr(pos1 == string::npos ? 0 : pos1,
    pos2 == string::npos ? str.length() - 1 : pos2 - pos1 + 1);
}

void replaceall(string& str, char ch1, char ch2)
{
  for(size_t ii=0; ii<str.length(); ii++ )
   if( str[ii] == ch1 )
            str[ii] = ch2;
}

/// read string as: "<characters>"
istream& f_getline(istream& is, std::string& str, char delim)
{
    char ch;
    is.get(ch);
    str="";

    while( is.good() && ( ch==' ' || ch=='\n' || ch== '\t') )
        is.get(ch);
    if(ch == '\"')
        is.get(ch);
    while( is.good() &&  ch!=delim && ch!= '\"' )
    {
        str += ch;
        is.get(ch);
    }
    while( is.good() &&  ch!=delim )
            is.get(ch);

   return is;
}

std::string u_makepath(const std::string& dir,
           const std::string& name, const std::string& ext)
{
    std::string Path(dir);
    if( dir != "")
      Path += "/";
    Path += name;
    Path += ".";
    Path += ext;

    return Path;
}

void u_splitpath(const std::string& aPath, std::string& dir,
            std::string& name, std::string& ext)
{
    std::string Path = aPath;
    replaceall( Path, '\\', '/');
    size_t pos = Path.rfind("/");
    if( pos != std::string::npos )
        dir = Path.substr(0, pos), pos++;
    else
        dir = "",    pos = 0;

    size_t pose = Path.rfind(".");
    if( pose != std::string::npos )
    {
        ext = Path.substr( pose+1, std::string::npos );
        name = Path.substr(pos, pose-pos);
    }
    else
    {
        ext = "";
        name = Path.substr(pos, std::string::npos);
    }
}

const long int bGRAN = 200;

/// Get Path of file and Reading list of file names from it, return number of files
char  (* f_getfiles(const char *f_name, char *Path,
		long int& nElem, char delim ))[fileNameLength]
{
  long int ii, bSize = bGRAN;
  char  (*filesList)[fileNameLength];
  char  (*filesListNew)[fileNameLength];
  filesList = new char[bSize][fileNameLength];
  std::string name;

// Get path
   std::string path_;
   std::string flst_name = f_name;
    replaceall( flst_name, '\\', '/');
   unsigned long int pos = flst_name.rfind("/");
   path_ = "";
   if( pos < std::string::npos )
      path_ = flst_name.substr(0, pos+1);
   strncpy( Path, path_.c_str(), 256-fileNameLength);
   Path[255] = '\0';

//  open file stream for the file names list file
   fstream f_lst( f_name/*flst_name.c_str()*/, ios::in );
   ErrorIf( !f_lst.good(), f_name, "Fileopen error");

// Reading list of names from file
  nElem = 0;
  while( !f_lst.eof() )
  {
	f_getline( f_lst, name, delim);
    if( nElem >= bSize )
    {    bSize = bSize+bGRAN;
         filesListNew = new char[bSize][fileNameLength];
         for( ii=0; ii<nElem-1; ii++ )
		   strncpy( filesListNew[ii], filesList[ii], fileNameLength);
	     delete[] filesList;
		 filesList =  filesListNew;
	}
    strncpy( filesList[nElem], name.c_str(), fileNameLength);
    filesList[nElem][fileNameLength-1] = '\0';
    nElem++;
  }

  // Realloc memory for reading size
  if( nElem != bSize )
  {
    filesListNew = new char[nElem][fileNameLength];
    for(  ii=0; ii<nElem; ii++ )
	  strncpy( filesListNew[ii], filesList[ii], fileNameLength);
	delete[] filesList;
	filesList =  filesListNew;
  }

  return filesList;
}
#endif

///// Writes Multi to a json/key-value string
///// \param brief_mode - Do not write data items that contain only default values
///// \param with_comments - Write files with comments for all data entries
//std::string TProfil::gemipm_to_string( bool addMui, bool with_comments, bool brief_mode )
//{
//    std::stringstream ss;
//    multi->to_text_file_gemipm( ss, addMui, with_comments, brief_mode );
//    return ss.str();
//}

///// Reads Multi structure from a json/key-value string
//bool TProfil::gemipm_from_string( const std::string& data,  DATACH  *dCH )
//{
//    std::stringstream ss;
//    ss.str(data);
//    multi->from_text_file_gemipm( ss, dCH );
//    return true;
//}

///// Reading structure MULTI (GEM IPM work structure) from text file
//void TProfil::readMulti( const char* path, DATACH  *dCH )
//{
//    // reads sizes and constants from txt file
//    fstream ff( path, ios::in );
//    ErrorIf( !ff.good() , path, "Fileopen error");
//    multi->from_text_file_gemipm( ff, dCH);
//}

// ------------------ End of ms_param.cpp -----------------------




