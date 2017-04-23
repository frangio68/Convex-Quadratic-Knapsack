/*--------------------------------------------------------------------------*/
/*--------------------------- File Main.C ----------------------------------*/
/*--------------------------------------------------------------------------*/
/* Simple main() for testing the quadratic knapsack solvers
 *
 * \version 1.20
 *
 * \date 17 - 04 - 2012
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Elettronica Informatica e Sistemistica \n
 *         Universita' della Calabria \n
 *
 * Copyright &copy 2011 - 2012 by Antonio Frangioni, Enrico Gorgone.
 */

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/
   
#define LogFILE 0  // print LOG

#define WHICH_KNPSOLVER 2
/* Selects what CQKnPClass is used to solve the problem:
    - 0 ==> CQKnPClpex
    - 1 ==> DualCQKnp
    - 2 ==> ExDualCQKnp */

#define WHICH_TEST_KNPSOLVER 1
/* Selects what CQKnPClass is used to check the first:
    - 0 ==> CQKnPClpex
    - 1 ==> DualCQKnp
    - 2 ==> ExDualCQKnp */

#if ( WHICH_KNPSOLVER > 0 ) || ( WHICH_TEST_KNPSOLVER > 0 )
 #define AUTO_SORT 0
 // if > 0, the sort is automatically switched between QS and BS according
 // to how many things (costs/bounds) are being changed
#else
 #define AUTO_SORT 0
 // never change this
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if ( WHICH_KNPSOLVER == 0 ) || ( WHICH_TEST_KNPSOLVER == 0 )
 #include "CQKnPCplex.h"
#endif

#if ( WHICH_KNPSOLVER > 0 ) || ( WHICH_TEST_KNPSOLVER > 0 )
 #include "DualCQKnP.h"

 #if ( WHICH_KNPSOLVER > 1 ) || ( WHICH_TEST_KNPSOLVER > 1 )
  #include "ExDualCQKnP.h"
 #endif
#endif

#include <limits>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <time.h>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace std;
using namespace CQKnPClass_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASSES -----------------------------------*/
/*--------------------------------------------------------------------------*/

class timer {
 public:
  timer( void ) { u = 0; }

  void Start( void ) { t_u = clock(); }

  void Stop( void ) {
   u += ( double( clock() - t_u ) ) / double( CLOCKS_PER_SEC );
   }

  double Read( void ) { return( u ); }
 private:
  double u;
  clock_t t_u;
 };

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

const char *const logF = "Log.knp";
const char *const tLog = "tLog.txt";

#if AUTO_SORT
 const double BSprc = 0.01;
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ GLOBALS -----------------------------------*/
/*--------------------------------------------------------------------------*/

int len = 0;

double *cstsC = 0;
double *cstsD = 0;
double *bndsL = 0;
double *bndsU = 0;

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

static inline double ABS( double x )
{
 return( x >= 0 ? x : - x );
 }

/*--------------------------------------------------------------------------*/

template<class T>
static inline void str2val( const char* const str , T &sthg )
{
 istringstream( str ) >> sthg;
 }

/*--------------------------------------------------------------------------*/

void GenerateCosts( void ) {
 for( int j = 0 ; j < len ; j++ ) {
  #if ( WHICH_KNPSOLVER != 1 ) && ( WHICH_TEST_KNPSOLVER != 1 )
   // the general case: quadratic costs can be zero
   if( drand48() < double( 0.1 ) ) {
    cstsD[ j ] = double( 0 );
    if( drand48() < double( 0.01 ) )
     cstsC[ j ] = ( drand48() - 0.5 ) * 100;
    else
     cstsC[ j ] = drand48() * 100;
    }
   else {
    cstsD[ j ] = drand48() * 100;
    cstsC[ j ] = ( drand48() - 0.5 ) * 100;
    }
  #else
   // the special case: strictly convex
   cstsD[ j ] = drand48() * 100 + 1e-8;
   cstsC[ j ] = ( drand48() - 0.5 ) * 100;
  #endif
  }
 }

/*--------------------------------------------------------------------------*/

void GenerateBounds( void ) {
 for( int j = 0 ; j < len ; j++ ) {
  #if ( WHICH_KNPSOLVER != 1 ) && ( WHICH_TEST_KNPSOLVER != 1 )
   // the general case: bounds can be infinite
   if( cstsD[ j ] == double( 0 ) ) {
    if( drand48() < double( 0.01 ) )
     bndsL[ j ] = -CQKnPClass::Inf<double>();
    else
     bndsL[ j ] = ( drand48() - 0.5 ) * 100;
    if( drand48() > double( 0.99 ) )
     bndsU[ j ] = CQKnPClass::Inf<double>();
    else
     bndsU[ j ] = ( drand48() - 0.5 ) * 100;
    }
   else {
    if( drand48() < double( 0.5 ) )
     bndsL[ j ] = -CQKnPClass::Inf<double>();
    else
     bndsL[ j ] = ( drand48() - 0.5 ) * 100;
    if( drand48() > double( 0.5 ) )
     bndsU[ j ] = CQKnPClass::Inf<double>();
    else
     bndsU[ j ] = ( drand48() - 0.5 ) * 100;
    }
  #else
   // the special case: finite bounds
   bndsL[ j ] = ( drand48() - 0.5 ) * 100;
   bndsU[ j ] = ( drand48() - 0.5 ) * 100;
  #endif

  if( bndsL[ j ] > bndsU[ j ] )
   bndsL[ j ] = bndsU[ j ];
  }
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 // read the command-line parameters- - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // parameters are:
 // nruns   = number of instances to generate [10000]
 // mx_size = max size of the instance [10000]
 // mn_size = min size of the instance [1]
 // chgprc  = average percentage of change [1 = 100%]
 // sort    = true if quick sort, otherwise bubble sort [true]
 // nreopt  = number of reoptimization cycles [2]

 int nruns = 10000;
 int mx_size = 100000;
 int mn_size = 1;
 double chgprc = 1;
 bool sort = true;
 int nreopt = 2;

 switch( argc ) {
  case( 7 ): str2val( argv[ 6 ] , nreopt );
  case( 6 ): str2val( argv[ 5 ] , sort );
  case( 5 ): str2val( argv[ 4 ] , chgprc );
  case( 4 ): str2val( argv[ 3 ] , mn_size );
  case( 3 ): str2val( argv[ 2 ] , mx_size );
  case( 2 ): str2val( argv[ 1 ] , nruns );
  }

 // open log file- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ofstream tLOG( tLog , ofstream::app );
 if( ! tLOG.is_open() )
  cerr << "Warning: cannot open log file """ << tLog << """" << endl;

 // construct the CQKnPClass objects - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( WHICH_KNPSOLVER == 0 )
  CQKnPCplex *qp1 = new CQKnPCplex();
 #elif( WHICH_KNPSOLVER == 1 )
  DualCQKnP *qp1 = new DualCQKnP( sort );
 #else
  ExDualCQKnP *qp1 = new ExDualCQKnP( sort );
 #endif

 #if( WHICH_TEST_KNPSOLVER == 0 )
  CQKnPCplex *qp2 = new CQKnPCplex();
 #elif( WHICH_TEST_KNPSOLVER == 1 )
  DualCQKnP *qp2 = new DualCQKnP( sort );
 #else
  ExDualCQKnP *qp2 = new ExDualCQKnP( sort );
 #endif

 // main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 timer *timer1 = new timer();
 timer *timer2 = new timer();
 timer *timert = new timer();

 tLOG << "nr: " << nruns << " ~ mxs: " << mx_size << " ~ mns: " << mn_size
	  << " ~ chg:  " << chgprc << " ~ sort: " << sort << endl;

 for( int i = 1 ; i < nruns ; i++ ) {
  // construct the instance - - - - - - - - - - - - - - - - - - - - - - - - -

  srand48( i );
  const int len1 = mn_size + ( mx_size - mn_size  ? 
			       lrand48() % ( mx_size - mn_size ) : int( 0 ) );
  if( len1 != len ) {
   delete[] bndsU;
   delete[] bndsL;
   delete[] cstsD;
   delete[] cstsC;

   cstsC = new double[ len1 ];
   cstsD = new double[ len1 ];
   bndsL = new double[ len1 ];
   bndsU = new double[ len1 ];
   len = len1;
   }

  const int iprc = int( len * chgprc );

  GenerateCosts();
  GenerateBounds();
  double vlm = drand48() * 1000;
  bool sense = ( drand48() > 0.5 ) ? true: false;

  #if LogFILE
   int KNPlvl = 5;
   ofstream LOGFile( logF , ofstream::out );
   if( ! LOGFile.is_open() )
    cerr << "Warning: cannot open log file """ << logF << """" << endl;
   else {
    qp1->SetKNPLog( &LOGFile , KNPlvl );
    qp2->SetKNPLog( &LOGFile , KNPlvl );
    }
  #endif

  qp1->LoadSet( len , cstsC , cstsD , bndsL , bndsU , vlm , sense );
  qp2->LoadSet( len , cstsC , cstsD , bndsL , bndsU , vlm , sense );

  // now start changing the instance- - - - - - - - - - - - - - - - - - - - -
  // you can change costs, bounds and volume in each combination (comprised
  // nothing), so it's 2^3 = 8 combinations; each one is performed nreopt
  // times, so it's nreopt * 8 = 16 at the default (note that at first k = 0
  // so nothing changes)

  for( int k = 0 ; k < nreopt * 8 ; k++ ) {
   if( k % 1 ) {  // change costs
    GenerateCosts();
    int strt = 0;
    int stp = len;
    if( chgprc < 1 ) {
     int nchg = max( iprc * ( drand48() + 0.5 ) , double( 1 ) );
     strt = lrand48() % ( len - iprc );
     stp = strt + nchg;
     }

    qp1->ChgLCosts( cstsC , 0 , strt , stp );
    qp2->ChgLCosts( cstsC , 0 , strt , stp );

    if( ( chgprc < 1 ) && ( drand48() < 0.5 ) ) {
     int nchg = max( iprc * ( drand48() + 0.5 ) , double( 1 ) );
     strt = lrand48() % ( len - iprc );
     stp = strt + nchg;
     }

    qp1->ChgQCosts( cstsD , 0 , strt , stp );
    qp2->ChgQCosts( cstsD , 0 , strt , stp );
    }

   if( k % 2 ) {  // change bounds
    GenerateBounds();
    int strt = 0;
    int stp = len;
    if( chgprc < 1 ) {
     int nchg = max( iprc * ( drand48() + 0.5 ) , double( 1 ) );
     strt = lrand48() % ( len - iprc );
     stp = strt + nchg;
     }

    qp1->ChgLBnds( bndsL , 0 , strt , stp );
    qp2->ChgLBnds( bndsL , 0 , strt , stp );

    if( ( chgprc < 1 ) && ( drand48() < 0.5 ) ) {
     int nchg = max( iprc * ( drand48() + 0.5 ) , double( 1 ) );
     strt = lrand48() % ( len - iprc );
     stp = strt + nchg;
     }

    qp1->ChgUBnds( bndsU , 0 , strt , stp );
    qp2->ChgUBnds( bndsU , 0 , strt , stp );
    }

   if( k % 4 ) {  // change volume
    vlm = drand48() * 1000;
    qp1->ChgVlm( vlm );
    qp2->ChgVlm( vlm );
    }

   // actually solve the problems - - - - - - - - - - - - - - - - - - - - - -
   // compare both the solve time proper and the time to construct a solution

   timert->Start();

   // run solver 1
   timer1->Start();
   double opt1;
   CQKnPClass::CQKStatus status1;
   #if AUTO_SORT && ( WHICH_KNPSOLVER > 0 )
    qp1->SetSort( ~ ( ( k % 3 ) && ( chgprc <= BSprc ) ) );
   #endif
   status1 = qp1->SolveKNP();
   if( status1 == CQKnPClass::kOK )
    qp1->KNPGetX();
   timer1->Stop();

   // run solver 2
   timer2->Start();
   double opt2;
   CQKnPClass::CQKStatus status2;
   #if AUTO_SORT && ( WHICH_TEST_KNPSOLVER > 0 )
    qp2->SetSort( ~ ( ( k % 3 ) && ( chgprc <= BSprc ) ) );
   #endif
   status2 = qp2->SolveKNP();
   if( status2 == CQKnPClass::kOK )
	qp2->KNPGetX();
   timer2->Stop();

   timert->Stop();

   // check the results - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( status1 != status2 ) {
    if( status1 != 4 )
     tLOG << "Test Failed - different status: ( " << status1 << " , "
	  << status2 << " , " << i << " , "  << k << " ) " << endl;
    }
   else {
    switch( status1 ) {
     case( CQKnPClass::kOK ):
      opt1 = qp1->KNPGetFO();
      opt2 = qp2->KNPGetFO();
      if( ABS( opt1 - opt2 ) >= 1e-8 * max( opt1 , double( 1 ) ) )
       tLOG << "Test Failed - gap: " << ABS( opt1 - opt2 )
	    << " ~ Opt1: " << opt1 << " ~ Opt2: " << opt2 << endl;
      else
      tLOG << "Test Ok, ( " << i << " , " << k << " ) " << endl;
      break;
     case( CQKnPClass::kUnfeasible ):
      tLOG << "Test Ok - infeasible" << endl;
      break;
     case( CQKnPClass::kUnbounded ):
      tLOG << "Test Ok - unbounded" << endl;
      break;
     default:
      tLOG << "No solution found??" << endl;
      return 1;
     }
    }
   }  // end( instance modification loop )
  }  // end( main loop )- - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 tLOG << "Total Time = " << timert->Read() << endl;
 tLOG << "CQKS1 Time = " << timer1->Read() << endl;
 tLOG << "CQKS2 Time = " << timer2->Read() << endl;

 delete timer2;
 delete timer1;
 delete timert;

 delete[] bndsU;
 delete[] bndsL;
 delete[] cstsC;
 delete[] cstsD;

 // destroy the solvers - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete qp2;
 delete qp1;

 // the end - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*---------------------------- End File Main.C -----------------------------*/
/*--------------------------------------------------------------------------*/
