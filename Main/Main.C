/*--------------------------------------------------------------------------*/
/*--------------------------- File Main.C ----------------------------------*/
/*--------------------------------------------------------------------------*/
/* Simple main() for testing the quadratic knapsack solvers
 *
 * \version 1.01
 *
 * \date 19 - 12 - 2012
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
   
#define FileCsv 0
/**< Use a csv file to print final report. */

#define WHICH_KNPSOLVER 1
/* Selects what CQKnPClass is used to solve the problem:
    - 0 ==> CQKnPClpex
    - 1 ==> DualCQKnp
    - 2 ==> ExDualCQKnp */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if WHICH_KNPSOLVER == 0
 #include "CQKnPCplex.h"
#elif WHICH_KNPSOLVER == 1
 #include "DualCQKnP.h"
#else
 #include "ExDualCQKnP.h"
#endif

#include <limits>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include "CQKnPClass.h"
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

#if( FileCsv )
 const char *const tLog = "tLog.csv";
#else
 const char *const tLog = "tLog.txt";
#endif

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
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 // read the command-line parameters- - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // parameters are:
 // filename
 // sort = true if quick sort, otherwise bubble sort [true]

 bool sort = true;

 switch( argc ) {
  case( 3 ): str2val( argv[ 2 ] , sort );
  case( 2 ): break;
  }

 // enter the try-block - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 try{

  // open log file - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ofstream tLOG( tLog , ofstream::app );
  if( ! tLOG.is_open() )
   cerr << "Warning: cannot open log file """ << tLog << """" << endl;

  // construct the CQKnPClass objects- - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if( WHICH_KNPSOLVER == 0 )
   CQKnPCplex *qp = new CQKnPCplex();
  #elif( WHICH_KNPSOLVER == 1 )
   DualCQKnP *qp = new DualCQKnP( sort );
  #else
   ExDualCQKnP *qp = new ExDualCQKnP( sort );
  #endif

  // read the instance - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ifstream instance( argv[ 1 ] );
  if( ! instance.is_open() )
   cerr << "Warning: cannot open instance file """ << argv[ 1 ] << """" << endl;

  qp->ReadInstance( instance );

  timer *timert = new timer();

  // actually solve the problems - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  timert->Start();

  double opt;
  CQKnPClass::CQKStatus status;

  status = qp->SolveKNP();
  if( status == CQKnPClass::kOK )
   qp->KNPGetX();

  timert->Stop();

  // check the results- - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if FileCsv == 0
   switch( status ) {
    case( CQKnPClass::kOK ):
     opt = qp->KNPGetFO();
     tLOG << "Value = " << opt << "\t" << "Ok" << endl;
     break;
    case( CQKnPClass::kUnfeasible ):
     tLOG << "Infeasible" << endl;
     break;
    case( CQKnPClass::kUnbounded ):
     tLOG << "Unbounded" << endl;
     break;
    default:
     tLOG << "No solution found??" << endl;
     return 1;
    }
   tLOG << "\t" << "Time = " << timert->Read() << endl;
  #else
   tLOG << argv[ 1 ] << "\t" << timert->Read() << "\t";
   switch( status ) {
    case( CQKnPClass::kOK ):
     opt = qp->KNPGetFO();
     tLOG << opt << "\t" << "Ok" << endl;
     break;
    case( CQKnPClass::kUnfeasible ):
     tLOG << "Infeasible" << endl;
     break;
    case( CQKnPClass::kUnbounded ):
     tLOG << "Unbounded" << endl;
     break;
    default:
     tLOG << "No solution found??" << endl;
     return 1;
    }
  #endif

  delete timert;

  // destroy the solver - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  delete qp;

  tLOG.close();
  instance.close();

  // the end - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  return( 0 );
  }
 catch( exception &e ) {
  cerr << e.what() << endl;
  return( 1 );
  }
 }  // end( main )

/*--------------------------------------------------------------------------*/
/*---------------------------- End File Main.C -----------------------------*/
/*--------------------------------------------------------------------------*/
