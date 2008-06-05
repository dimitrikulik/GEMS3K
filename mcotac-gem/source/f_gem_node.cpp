#include <iomanip>

#include "node.h"

// (1) Initialization of TNode data exchange interface
// This function reads DATACH and MULTI text input files
// returns 0 if success or 1 if error
//  Must be called once before the beginning of the coupled RMT calculation
//
#ifdef __unix
#ifdef __PGI
    extern "C" int f_gem_init_( char* string_,
 #else
   extern "C" int f_gem_init( char* string_,
#endif
#else
   extern "C" int __stdcall F_GEM_INIT( char* string_,
#endif
                  unsigned int length_ )
{

  gstring string_cto_i1c( string_, 0, length_);
  string_cto_i1c.strip();

	// Creating TNode structure accessible trough node pointer
  TNode* node = new TNode();

   // Here we read the files needed as input for initializing GEMIPM2K
   // The easiest way to prepare them is to use GEMS-PSI code (GEM2MT module)
   if( node->GEM_init( string_cto_i1c.c_str() ) ){
       printf("\nerror during gems init\n");
       return 1;  // error occured during reading the files
     }
   return 0;
}


// (2) This function gets some elements of the DATACH data structure
// that can be used in the part of the program external to GEMIPM2K
// returns 0 if success or 1 if error
// Parameter list may be extended in future with other DCH elements
//
#ifdef __unix
#ifdef __PGI
    extern "C" int  f_gem_get_dch_( int& p_nICb, int& p_nDCb, int& p_nPHb, int& p_nPSb,int& p_nPH)
//    extern "C" int  f_gem_get_dch_( int& p_nICb, int& p_nDCb, int& p_nPHb, float* p_A )
#else
   extern "C" int  f_gem_get_dch( int& p_nICb, int& p_nDCb, int& p_nPHb, int& p_nPSb,int& p_nPH)
//   extern "C" int  f_gem_get_dch( int& p_nICb, int& p_nDCb, int& p_nPHb, float* p_A )
#endif
#else
  extern "C" int  __stdcall  F_GEM_GET_DCH(  // All parameters are return values
   int& p_nICb,   // Number of Independent Components (ICs) in chemical system
   int& p_nDCb,   // Number of Dependent Components (DCs) in chemical system
   int& p_nPHb,   // Number of Phases (PHs) in chemical system
   int& p_nPSb,   // number of Phases-solutions (<= nPS) used in the data bridge
   int& p_nPH,   // Number of Phases (PHs) in chemical system

//   float* p_A     // Stoichiometry matrix A with nDCb rows and nICb columns
                  // before calling  F_GEM_GET_DCH(), a memory block of the size
                  // greater than sizeof(float)*nDCb*nICb must be allocated and
                  // a pointer to is passed in p_A
)
#endif
{

   p_nICb = (int) TNode::na->pCSD()->nICb;
   p_nDCb = (int) TNode::na->pCSD()->nDCb;
   p_nPHb = (int) TNode::na->pCSD()->nPHb;
   p_nPSb = (int) TNode::na->pCSD()->nPSb;
   p_nPH  = (int) TNode::na->pCSD()->nPH;

//    for( int i=0; i<p_nICb; i++ )
//       for( int j=0; j<p_nDCb; j++ )
//		   p_A[j*p_nICb+i] = nodeCH_A(j, i);
   return 0;
}

// (2) This function gets some elements of the DATACH data structure
// that can be used in the part of the program external to GEMIPM2K
// returns 0 if success or 1 if error
// Parameter list may be extended in future with other DCH elements
//
#ifdef __unix
#ifdef __PGI
    extern "C" int  f_gem_get_pa_( float* p_A )
//    extern "C" int  f_gem_get_dch_( int& p_nICb, int& p_nDCb, int& p_nPHb, float* p_A )
#else
   extern "C" int  f_gem_get_pa( float* p_A )
//   extern "C" int  f_gem_get_dch( int& p_nICb, int& p_nDCb, int& p_nPHb, float* p_A )
#endif
#else
  extern "C" int  __stdcall  F_GEM_GET_PA(  // All parameters are return values

   float* p_A;     // Stoichiometry matrix A with nDCb rows and nICb columns
                  // before calling  F_GEM_GET_DCH(), a memory block of the size
                  // greater than sizeof(float)*nDCb*nICb must be allocated and
                  // a pointer to is passed in p_A
)
#endif
{
   int p2_nICb;   // Number of Independent Components (ICs) in chemical system
   int p2_nDCb;   // Number of Dependent Components (DCs) in chemical system
 
   p2_nICb = TNode::na->pCSD()->nICb;
   p2_nDCb = TNode::na->pCSD()->nDCb;
   
    for( int i=0; i<p2_nICb; i++ )
       for( int j=0; j<p2_nDCb; j++ )
		   p_A[j*p2_nICb+i] = nodeCH_A(j, i);
   return 0;
}

// (3) Reads DATABR input file that describes a single node composition,
//   runs the GEM calculation and provides via parameters the GEM input and
//   output data
//
//  returns 0 if success or 1 if error
//  Usually is called once at the beginning of coupled modeling for initializing
//  all nodes
//
#ifdef __unix
#ifdef __PGI
   extern "C" int   f_gem_read_node_( char* string_,  
#else
  extern "C" int   f_gem_read_node( char* string_,
#endif
#else
  extern "C" int __stdcall   F_GEM_READ_NODE(
  char* string_,        // path (file name) of the DATABR file
  unsigned int length_, // length of the block string_
#endif
   int& p_NodeHandle,    // Node identification handle
   int& p_NodeTypeHY,    // Node type (hydraulic); see typedef NODETYPE
   int& p_NodeTypeMT,    // Node type (mass transport); see typedef NODETYPE
   int& p_NodeStatusFMT, // Node status code FMT; see typedef NODECODEFMT
   int& p_NodeStatusCH,  // Node status code CH;  see typedef NODECODECH
   int& p_IterDone,      // Number of iterations performed by IPM (if not need GEM)
// Chemical scalar variables
   double &p_T,     // Temperature T, K                        +      +      -     -
   double &p_P,     // Pressure P, bar                         +      +      -     -
   double &p_Vs,    // Volume V of reactive subsystem, cm3     -      -      +     +
   double &p_Vi,    // Volume of inert subsystem  ?            +      -      -     +
   double &p_Ms,    // Mass of reactive subsystem, kg          +      +      -     -
   double &p_Mi,    // Mass of inert part, kg    ?             +      -      -     +
   double &p_Gs,    // Gibbs energy of reactive subsystem (J)  -      -      +     +
   double &p_Hs,    // Enthalpy of reactive subsystem (J)      -      -      +     +
   double &p_Hi,    // Enthalpy of inert subsystem (J) ?       +      -      -     +
   double &p_IC,    // Effective molal aq ionic strength           -      -      +     +
   double &p_pH,    // pH of aqueous solution                      -      -      +     +
   double &p_pe,    // pe of aqueous solution                      -      -      +     +
   double &p_Eh,    // Eh of aqueous solution, V                   -      -      +     +
// Dynamic data - dimensions see in DATACH.H and DATAMT.H structures
// exchange of values occurs through lists of indices, e.g. xDC, xPH
   double  *p_bIC,  // bulk mole amounts of IC[nICb]                +      +      -     -
   double  *p_rMB,  // MB Residuals from GEM IPM [nICb]             -      -      +     +
   double  *p_uIC,  // IC chemical potentials (mol/mol)[nICb]       -      -      +     +
   double  *p_xDC,    // DC mole amounts at equilibrium [nDCb]      -      -      +     +
   double  *p_gam,    // activity coeffs of DC [nDCb]               -      -      +     +
   double  *p_dul,  // upper kinetic restrictions [nDCb]            +      +      -     -
   double  *p_dll,  // lower kinetic restrictions [nDCb]            +      +      -     -
   double  *p_aPH,  // Specific surface areas of phases (m2/g)      +      +      -     -
   double  *p_xPH,  // total mole amounts of phases [nPHb]          -      -      +     +
   double  *p_vPS,  // phase volume, cm3/mol        [nPSb]          -      -      +     +
   double  *p_mPS,  // phase (carrier) mass, g      [nPSb]          -      -      +     +
   double  *p_bPS,  // bulk compositions of phases  [nPSb][nICb]    -      -      +     +
   double  *p_xPA  // amount of carrier in phases  [nPSb] ??       -      -      +     +
  // What else?
//   double  *p_dRes1
#ifdef __unix
   , unsigned int length_
#endif
)
{
   short NodeHandle, NodeTypeHY,  NodeTypeMT, NodeStatusFMT, NodeStatusCH, IterDone;

   NodeHandle = (short) p_NodeHandle;
   NodeTypeHY = (short) p_NodeTypeHY;
   NodeTypeMT = (short) p_NodeTypeMT;
   NodeStatusFMT = (short) p_NodeStatusFMT;
   NodeStatusCH = (short) p_NodeStatusCH;
   IterDone = (short) p_IterDone;

   gstring dbr_file_name( string_, 0, length_);
   dbr_file_name.strip();

 // Read DATABR structure from file
    TNode::na->GEM_read_dbr( dbr_file_name.c_str() );

 // re-calculating equilibrium by calling GEMIPM
 //	 short NodeStatusCH = TNode::na->GEM_run();
 //     if( !( NodeStatusCH == OK_GEM_AIA || NodeStatusCH == OK_GEM_PIA ) )
 //        return 1;

 // Extracting chemical data into FMT part
     TNode::na->GEM_restore_MT( NodeHandle, NodeStatusCH, p_T,
                                p_P, p_Vs, p_Ms, p_bIC, p_dul, p_dll, p_aPH );

  // Extracting GEMIPM output data to FMT part
     TNode::na->GEM_to_MT( NodeHandle, NodeStatusCH, IterDone,
                           p_Vs, p_Ms, p_Gs, p_Hs, p_IC, p_pH, p_pe, p_Eh,
                           p_rMB, p_uIC, p_xDC, p_gam, p_xPH, p_vPS, p_mPS, p_bPS, p_xPA );

 p_NodeHandle = (int) NodeHandle;
 p_NodeTypeHY = (int) NodeTypeHY;
 p_NodeTypeMT = (int) NodeTypeMT;
 p_NodeStatusFMT = (int) NodeStatusFMT;
 p_NodeStatusCH = (int) NodeStatusCH;
 p_IterDone = (int) IterDone;

 return 0;
}

// (4) Provides GEM input data to GEMIPM2K kernel (by copying them from
//   parameters to DATABR work structure);
//   runs the GEM-IPM calculation;
//   and provides via parameters
//   the GEM output data back to the calling program
//
//  returns  0 if success or 1 if error
//  Is called on each external iteration for each node
//
#ifdef __unix
#ifdef __PGI
    extern "C" int  f_gem_calc_node_(
 #else
   extern "C" int  f_gem_calc_node(
#endif
#else
   extern "C" int  __stdcall   F_GEM_CALC_NODE(
#endif
   int& p_NodeHandle,    // Node identification handle
   int& p_NodeTypeHY,    // Node type (hydraulic); see typedef NODETYPE
   int& p_NodeTypeMT,    // Node type (mass transport); see typedef NODETYPE
   int& p_NodeStatusFMT, // Node status code FMT; see typedef NODECODEFMT
   int& p_NodeStatusCH,  // Node status code CH;  see typedef NODECODECH
   int& p_IterDone,      // Number of iterations performed by IPM (if not need GEM)
// Chemical scalar variables
   double &p_T,     // Temperature T, K                        +      +      -     -
   double &p_P,     // Pressure P, bar                         +      +      -     -
   double &p_Vs,    // Volume V of reactive subsystem, cm3     -      -      +     +
   double &p_Vi,    // Volume of inert subsystem  ?            +      -      -     +
   double &p_Ms,    // Mass of reactive subsystem, kg          +      +      -     -
   double &p_Mi,    // Mass of inert part, kg    ?             +      -      -     +
   double &p_Gs,    // Gibbs energy of reactive subsystem (J)  -      -      +     +
   double &p_Hs,    // Enthalpy of reactive subsystem (J)      -      -      +     +
   double &p_Hi,    // Enthalpy of inert subsystem (J) ?       +      -      -     +
   double &p_IC,    // Effective molal aq ionic strength           -      -      +     +
   double &p_pH,    // pH of aqueous solution                      -      -      +     +
   double &p_pe,    // pe of aqueous solution                      -      -      +     +
   double &p_Eh,    // Eh of aqueous solution, V                   -      -      +     +
// Dynamic data - dimensions see in DATACH.H and DATAMT.H structures
// exchange of values occurs through lists of indices, e.g. xDC, xPH
   double  *p_bIC,  // bulk mole amounts of IC[nICb]                +      +      -     -
   double  *p_rMB,  // MB Residuals from GEM IPM [nICb]             -      -      +     +
   double  *p_uIC,  // IC chemical potentials (mol/mol)[nICb]       -      -      +     +
   double  *p_xDC,    // DC mole amounts at equilibrium [nDCb]      -      -      +     +
   double  *p_gam,    // activity coeffs of DC [nDCb]               -      -      +     +
   double  *p_dul,  // upper kinetic restrictions [nDCb]           +      +      -     -
   double  *p_dll,  // lower kinetic restrictions [nDCb]           +      +      -     -
   double  *p_aPH,  // Specific surface areas of phases (m2/g)      +      +      -     -
   double  *p_xPH,  // total mole amounts of phases [nPHb]          -      -      +     +
   double  *p_vPS,  // phase volume, cm3/mol        [nPSb]          -      -      +     +
   double  *p_mPS,  // phase (carrier) mass, g      [nPSb]          -      -      +     +
   double  *p_bPS,  // bulk compositions of phases  [nPSb][nICb]    -      -      +     +
   double  *p_xPA  // amount of carrier in phases  [nPSb] ??       -      -      +     +
  // What else?
  // double  *p_dRes1
 )
{
  int iRet = 0;
  int ii;
  bool uPrimalSol = true;
 

  short NodeHandle, NodeStatusCH, IterDone;

  NodeHandle =(short) p_NodeHandle;
  NodeStatusCH =(short) p_NodeStatusCH;
  IterDone =(short) p_IterDone;
  // (2) ----------------------------------------------
  // Work loop for the coupled FMT-GEM modelling


  if (uPrimalSol) {
   // sends only speciation changed by mcotac 
//   TNode::na->GEM_from_MT( NodeHandle, NodeStatusCH,
//             p_T, p_P, p_Vs, p_Ms, p_bIC, p_dul, p_dll,  p_aPH, p_xDC,p_gam);
	for ( ii=0 ; ii < TNode::na->pCSD()->nICb ; ii++ ) p_bIC[ii]=0;
		
   TNode::na->GEM_from_MT( NodeHandle, NodeStatusCH,
             p_T, p_P, p_Vs, p_Ms, p_bIC, p_dul, p_dll,  p_aPH, p_xDC);

//	for ( ii=0 ; ii < TNode::na->pCSD()->nICb ; ii++ ) {
//	      printf("bIC: %i %g \n",ii, TNode::na->pCNode()->bIC[ii]);
//		}
             TNode::na->pCNode()->bIC[(TNode::na->pCSD()->nICb)-1]=0.0;

   }
   else
   {
  // Setting input data for GEMIPM
 
   // sends only speciation changed by mcotac 
	for ( ii=0 ; ii < TNode::na->pCSD()->nICb ; ii++ ) 
             p_bIC[ii]=0;
   TNode::na->GEM_from_MT( NodeHandle, NodeStatusCH,
             p_T, p_P, p_Vs, p_Ms, p_bIC, p_dul, p_dll,  p_aPH, p_xDC);
          TNode::na->pCNode()->bIC[(TNode::na->pCSD()->nICb)-1]=0.0;

 
//  TNode::na->GEM_from_MT( NodeHandle, NodeStatusCH,
//             p_T, p_P, p_Vs, p_Ms, p_bIC, p_dul, p_dll,  p_aPH, p_xDC);
//   TNode::na->GEM_from_MT( NodeHandle, NodeStatusCH,
//             p_T, p_P, p_Vs, p_Ms, p_bIC, p_dul, p_dll,  p_aPH );
   }
  TNode::na->GEM_write_dbr( "Test-clay-cement-dbr.dat");
 // Calling GEMIPM calculation
   iRet = TNode::na->GEM_run(uPrimalSol);
   if( !( iRet == OK_GEM_AIA || iRet == OK_GEM_SIA ) )
   {
	  return 1;
   }

  // Extracting GEMIPM output data to FMT part
   TNode::na->GEM_to_MT( NodeHandle, NodeStatusCH, IterDone,
       p_Vs, p_Ms, p_Gs, p_Hs, p_IC, p_pH, p_pe, p_Eh, p_rMB, p_uIC,
       p_xDC, p_gam, p_xPH, p_vPS, p_mPS, p_bPS, p_xPA  );

 p_NodeHandle =(int) NodeHandle;
 p_NodeStatusCH =(int)  NodeStatusCH;
 p_IterDone = (int) IterDone;

/**************************************************************
  int ii;
  int nIC, nDC, nPH;
 // Extracting data bridge dimensionalities
   nIC = TProfil::pm->multi->CSD->nICb;
   nDC = TProfil::pm->multi->CSD->nDCb;
   nPH = TProfil::pm->multi->CSD->nPHb;

  fstream f_log("MAIF_CALC.txt", ios::out|ios::app );

  f_log << "Node = " <<  p_NodeHandle << "( " << readF << ", " <<
         indN << ", " << indM << ", " <<  indK << " )";
  f_log << "  NodeStatus = " <<  p_NodeStatusCH;
  f_log << "  ReturnStatus = " <<  iRet;
  f_log << "  IterDone = " <<  p_IterDone << endl;

  f_log << "p_xDC" << endl;
  for(ii=0; ii<nDC; ii++ )
    f_log << setprecision(15) << scientific << p_xDC[ii]<< ", ";
  f_log <<  endl;

  f_log << "p_gam" << endl;
  for(ii=0; ii<nDC; ii++ )
    f_log << setprecision(15) << scientific << p_gam[ii]<< ", ";
  f_log <<  endl;

  f_log << "p_xPH" << endl;
  for(ii=0; ii<nPH; ii++ )
    f_log << setprecision(15) << scientific << p_xPH[ii]<< ", ";
  f_log <<  endl;

  f_log << "p_uIC" << endl;
  for(ii=0; ii<nIC; ii++ )
    f_log << setprecision(15) << scientific << p_uIC[ii]<< ", ";
  f_log <<  endl;

  f_log << "p_rMB" << endl;
  for(ii=0; ii<nIC; ii++ )
    f_log << setprecision(15) << scientific << p_rMB[ii]<< ", ";
  f_log <<  endl;

  f_log <<  endl;
  f_log <<  endl;

*************************************************************/

  return 0;
}

// (5) Writes a DATABR text file (with file path name provided in string_)
//    from the currently available content of the DATABR work structure.
// can be used for interruption of coupled modeling or for debugging purposes
//
#ifdef __unix
 extern "C" void f_gem_write_node( char* string_, unsigned int length_ )
#else
 extern "C"  void  __stdcall   F_GEM_WRITE_NODE(
   char* string_,        // path (file name) of the DATABR file
   unsigned int length_ // length of the block fname_
)
#endif
{
  gstring fname( string_, 0, length_);
  fname.strip();
  TNode::na->GEM_write_dbr( fname.c_str()  );
}

// (6) For detailed examination of GEM work data structure:
// writes GEMIPM internal MULTI data structure into text file
// path name fname in debugging format (different from MULTI input format).
// This file cannot be read back with F_GEM_INIT()!
//
#ifdef __unix
 extern "C" void f_gem_print_ipm( char* string_, unsigned int length_ )
#else
 extern "C"  void  __stdcall   F_GEM_PRINT_IPM(
   char* string_,        // path (file name) of the DATABR file
   unsigned int length_ // length of the block fname_
)
#endif
{
  gstring fname( string_, 0, length_);
  fname.strip();
  TNode::na->GEM_print_ipm( fname.c_str()  );
}


// end of f_gem_node.cpp
// -------------------------------------------------------------------------------------
