/*--------------------------------------------------------------------------*/
/*-------------------------- File CQKnPCplex.C -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Continuous Quadratic Knapsack Problems (CQKnP) solver based on calls to
 * the CPLEX callable library.
 *
 * \version 1.03
 *
 * \date 22 - 12 - 2012
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
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CQKnPCplex.h"
#include "math.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace CQKnPClass_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( CQKnPClass_LOG )
 #define KLOG( l , x ) if( KNPLLvl > l ) *KNPLog << x
#else
 #define KLOG( l , x )
 #define Log1()
#endif

/*--------------------------------------------------------------------------*/
/*-------------------- IMPLEMENTATION OF CQKnPCplex ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

CQKnPCplex::CQKnPCplex( const int alg , const double eps ) : CQKnPClass()
{
 int cplexStatus;

 // cplex:initializing env (if necessary); consume only a license
 if( ! (users++) )
  if( cplexEnv == 0 ) {
   cplexEnv = CPXopenCPLEX( &cplexStatus );
   if( cplexStatus )
    throw(
     CQKException( "CQKnPCplex(): Error initializing cplex environment" ) );
   }
  else
   throw( CQKException(
	 "CQKnPCplex(): Error: cplex environment already initialized??" ) );

 switch ( alg ) {
  case ( 1 ): SetCplexParam( CPX_PARAM_QPMETHOD , CPX_ALG_PRIMAL ); break;
  case ( 2 ): SetCplexParam( CPX_PARAM_QPMETHOD , CPX_ALG_DUAL ); break;
  default: SetCplexParam( CPX_PARAM_QPMETHOD , CPX_ALG_AUTOMATIC ); break;
  }

 SetCplexParam( CPX_PARAM_BAREPCOMP , eps );

 XSol = 0;
 lp = 0;

 }  // end( CQKnPCplex() )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void CQKnPCplex::LoadSet( const int pn ,
			  const double *pC , const double *pD ,
			  const double *pA , const double *pB ,
			  const double pV , const bool sns )
{
 if( lp ) {
  if( CPXfreeprob( cplexEnv , &lp ) )
   throw( CQKException(
                  "CQKnPCplex::LoadSet: Error in freeing cplex problem" ) );
  lp = 0;
  }

 if( XSol ) {
  delete[] XSol;
  XSol = 0;
  }

 n = pn;
 if( ! n )  // just sit down in the corner and wait
  return;
 else
  XSol = new double[ n ];

 int cplexStatus;

 // cplex:instantiating problem object  - - - - - - - - - - - - - - - - - - -

 lp = CPXcreateprob( cplexEnv , &cplexStatus , "CQKnP" );
 if( cplexStatus )
  throw( CQKException(
	        "CQKnPCplex::LoadSet: Error initializing cplex problem" ) );

 // knapsack constraint - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int *matbeg = new int[ n ];
 int *matcnt = new int[ n ];
 int *matind = new int[ n ];
 double *lb = new double[ n ];
 double *ub = new double[ n ];
 double *matval = new double[ n ];
 double *costsL = new double[ n ];

 for( int i = 0 ; i < n ; i++ ) {
  matbeg[ i ] = i;
  matcnt[ i ] = 1;
  matind[ i ] = 0;
  matval[ i ] = 1;
  
  if( ( ! pA ) || ( pA[ i ] <= -Inf<double>() ) )
   lb[ i ] = - CPX_INFBOUND;
  else
   lb[ i ] = pA[ i ];

  if( ( ! pB ) || ( pB[ i ] >= Inf<double>() ) )
   ub[ i ] = + CPX_INFBOUND;
  else
   ub[ i ] = pB[ i ];

  if( ! pC )
   costsL[ i ] = 0;
  else
   costsL[ i ] = pC[ i ];
  }

 // cplex: passing data block - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // copy the LP part of the problem data into the lp  - - - - - - - - - - - -

 const char sense = sns ? 'E' : 'L';

 cplexStatus = CPXcopylp( cplexEnv , lp , n , 1 , 1 , costsL , &pV , &sense ,
			  matbeg , matcnt , matind , matval , lb , ub , 0 );
 if( cplexStatus )
  throw( CQKException(
                "CQKnPCplex::LoadSet: error initializing cplex problem" ) );

 // copy the quadratic part of objective function  - - - - - - - - - - - - -

 if( pD ) {  // ... if any
  for( int i = 0 ; i < n ; i++ )
   costsL[ i ] = 2 * pD[ i ];

  cplexStatus = CPXcopyqpsep( cplexEnv , lp , costsL );
  if( cplexStatus )
   throw( CQKException(
                  "CQKnPCplex::LoadSet: error setting quadratic costs" ) );
  }

 delete[] costsL;
 delete[] matval;
 delete[] ub;
 delete[] lb;
 delete[] matind;
 delete[] matcnt;
 delete[] matbeg;

 }  // end( CQKnPCplex::LoadSet() )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

CQKnPClass::CQKStatus CQKnPCplex::SolveKNP( void )
{
 CPXqpopt( cplexEnv , lp );

 #if CQKnPClass_LOG
  if( CPXwriteprob( cplexEnv , lp , "Q2Knp.lp" , 0 ) )
   throw( CQKException(
                  "CQKnPCplex::LoadSet: error writing problem Q2Knp.lp" ) );
 #endif

 switch( CPXgetstat( cplexEnv , lp ) ) {
  case( CPX_STAT_OPTIMAL ):    status = kOK;
                               break;
  case( CPX_STAT_INFEASIBLE ): status = kUnfeasible;
                               break;
  case( CPX_STAT_UNBOUNDED ):  status = kUnbounded;
                               break;
  case( CPX_STAT_NUM_BEST ):   status = kStopped;
                               break;
  default:                     status = kError;
  }

 return( CQKStatus( status ) );

 } // end( CQKnPCplex::SolveKNP() )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

const double *CQKnPCplex::KNPGetX( void )
{
 if( status != kOK )
  throw( CQKException( "CQKnPCplex::KNPGetX(): problem not solved" ) );

 if( CPXgetx( cplexEnv , lp , XSol , 0 , n - 1 ) )
  throw( CQKException(
		  "CQKnPCplex::KNPGetX(): error retrieving cplex solution" ) );
 Log1();

 for( int j = 0 ; j < n ; j++ )
  if( XSol[ j ] == CPX_INFBOUND )
   XSol[ j ] =  Inf<double>();
  else
   if( XSol[ j ] == - CPX_INFBOUND )
    XSol[ j ] =  - Inf<double>();

 return( XSol );

 }  // end( CQKnPCplex::KNPGetX )

/*--------------------------------------------------------------------------*/

double CQKnPCplex::KNPGetPi( void )
{
 double pi;
 if( CPXgetpi( cplexEnv , lp , &pi , 0 , 0 ) )
  throw( CQKException(
           "CQKnPCplex::KNPGetPi(): error retrieving cplex dual solution" ) );

 return( pi );

 }  // end( CQKnPCplex::KNPGetPi() )

/*--------------------------------------------------------------------------*/

double CQKnPCplex::KNPGetFO( void )
{
 double fo;

 switch( status ) {
 case( kOK ):
 case( kStopped ):
  if( CPXgetobjval( cplexEnv , lp , &fo ) )
   throw( CQKException(
          "CQKnPCplex::KNPGetFO(): error retrieving cplex optimal value" ) );
  KLOG( 1 , std::endl << "Opt. value: " << fo << std::endl );
  break;
 case( kUnfeasible ): fo = + Inf<double>(); break;
 case( kUnbounded ):  fo = - Inf<double>(); break;
 default:
  throw( CQKException( "CQKnPCplex::KNPGetFO(): problem not solved" ) );
 }

 return( fo );

 }  // end( CQKnPCplex::KNPGetFO() )

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

void CQKnPCplex::KNPLCosts( double *csts , const int *nms , int strt , int stp )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; )
   CPXgetcoef( cplexEnv , lp , -1 , i , csts++ );
  }
 else
  CPXgetobj( cplexEnv , lp , csts , strt , stp - 1 );

 }  // end( CQKnPCplex::KNPLCosts )

/*--------------------------------------------------------------------------*/

void CQKnPCplex::KNPQCosts( double *csts , const int *nms , int strt , int stp )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; ) {
   CPXgetqpcoef( cplexEnv , lp , i , i , csts );
   *(csts++) *= double( 0.5 );
   }
  }
 else
  for( int i = strt ;  i < stp ; i++ ) {
   CPXgetqpcoef( cplexEnv , lp , i , i , csts );
   *(csts++) *= double( 0.5 );
   }
 }  // end( CQKnPCplex::KNPQCosts )

/*--------------------------------------------------------------------------*/

void CQKnPCplex::KNPLBnds( double *bnds , const int *nms , int strt , int stp )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; )
   CPXgetlb( cplexEnv , lp , bnds++,  i , i );
  }
 else
  CPXgetlb( cplexEnv , lp , bnds , strt , stp - 1 );

 }  // end( CQKnPCplex::KNPLBnds )

/*--------------------------------------------------------------------------*/

void CQKnPCplex::KNPUBnds( double *bnds , const int *nms , int strt , int stp )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; )
   CPXgetub( cplexEnv , lp , bnds++,  i , i );
  }
 else
  CPXgetub( cplexEnv , lp , bnds , strt , stp - 1 );

 }  // end( CQKnPCplex::KNPUBnds )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void CQKnPCplex::ChgLCosts( const double *csts , const int *nms ,
			    int strt , int stp )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {

  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; )
   CPXchgobj( cplexEnv , lp , 1 , &i , csts++ );
  }
 else
  for( int i = strt ; i < stp ; i++ )
   CPXchgobj( cplexEnv , lp , 1 , &i , &csts[ i ] );

 }  // end( CQKnPCplex::ChgLCosts )

/*--------------------------------------------------------------------------*/

void CQKnPCplex::ChgQCosts( const double *csts , const int *nms ,
			    int strt , int stp )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; )
   CPXchgqpcoef( cplexEnv , lp , i , i , *(csts++) * double( 2 ) );
  }
 else
  for( int i = strt ; i < stp ; i++ )
   CPXchgqpcoef( cplexEnv , lp , i , i , csts[ i ] * double( 2 ) );

 }  // end( CQKnPCplex::ChgQCosts )

/*--------------------------------------------------------------------------*/

void CQKnPCplex::ChgLBnds( const double *bnds , const int *nms ,
			   int strt , int stp  )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 const char lu = 'L';
 double lb;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; ) {
   if( *bnds <= -Inf<double>() )
    lb = -CPX_INFBOUND;
   else
    lb = *bnds;
   CPXchgbds( cplexEnv , lp , 1 , &i , &lu , &lb );
   bnds++;
   }
  }
 else
  for( int i = strt ; i < stp ; i++ ) {
   if( *bnds <= -Inf<double>() )
    lb = -CPX_INFBOUND;
   else
    lb = *bnds;
   CPXchgbds( cplexEnv , lp , 1 , &i , &lu , &lb );
   bnds++;
   }

 }  // end( CQKnPCplex::ChgLBnds )

/*--------------------------------------------------------------------------*/

void CQKnPCplex::ChgUBnds( const double *bnds , const int *nms ,
			   int strt , int stp  )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 const char lu = 'U';
 double ub;

 if( nms ) {
  while( *nms < strt )
  nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; ) {
   if( *bnds >= Inf<double>() )
    ub = CPX_INFBOUND;
   else
    ub = *bnds;
   CPXchgbds( cplexEnv , lp , 1 , &i , &lu , &ub );
   bnds++;
   }
  }
 else
  for( int i = strt ; i < stp ; i++ ) {
   if( *bnds >= Inf<double>() )
    ub = CPX_INFBOUND;
   else
    ub = *bnds;
   CPXchgbds( cplexEnv , lp , 1 , &i , &lu , &ub );
   bnds++;
   }

 } // end ( CQKnPCplex::ChgUBnds )

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

CQKnPCplex::~CQKnPCplex()
{
 if( lp ) {
  if( CPXfreeprob( cplexEnv , &lp ) )
   throw( CQKException( "~CQKnPCplex(): Error in freeing cplex problem" ) );
  lp = 0;
  }

 delete[] XSol;

 if( ! (--users) )
  if( CPXcloseCPLEX( &cplexEnv ) )
   throw( CQKException( "~CQKnPCplex(): Error in closing cplex environment" )
	  );
  else
   cplexEnv = 0;

 } // end( CQKnPCplex::~CQKnPCplex() )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

#if CQKnPClass_LOG

inline void CQKnPCplex::Log1( void  )
{
 if( KNPLLvl > 1 ) {
  KNPLog->precision( 4 );
  *KNPLog << std::endl << "Opt. sol.: ";
  for( int i = 0 ; i < n ; i++ ) {
   if( ! ( i % 4 ) )
    *KNPLog << std::endl;

   *KNPLog << "x( " <<  i  << " ) [ " << XSol[ i ] << " ] ";
   }

  *KNPLog << std::endl;
  }
 }

#endif

/*-------------------------------------------------------------------------*/
/*------------------initialize static members------------------------------*/
/*-------------------------------------------------------------------------*/

CPXENVptr CQKnPCplex::cplexEnv = 0;
int CQKnPCplex::users = 0;

/*--------------------------------------------------------------------------*/
/*---------------------- End File CQKnPCplex.C -----------------------------*/
/*--------------------------------------------------------------------------*/
