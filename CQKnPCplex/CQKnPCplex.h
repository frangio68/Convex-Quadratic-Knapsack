/*--------------------------------------------------------------------------*/
/*-------------------------- File CQKnPCplex.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Continuous Quadratic Knapsack Problems (CQKnP) solver based on calls to
 * the CPLEX callable library. Conforms to the standard interface for CQKnP
 * solver defined by the abstract base class CQKnpClass.
 *
 * \version 1.04
 *
 * \date 22 - 12 - 2012
 *
 * \author Antonio Frangioni \n
 * \author Operations Research Group \n
 * \author Dipartimento di Informatica \n
 * \author Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 * \author Dipartimento di Elettronica Informatica e Sistemistica \n
 * \author Universita' della Calabria \n
 *
 * Copyright &copy 2011 - 2012 by Antonio Frangioni, Enrico Gorgone.
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __QCQKnPCplex
 #define __QCQKnPCplex /* self-identification: #endif at the end of the file*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CQKnPClass.h"
#include "cplex.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE and USING ----------------------------*/
/*--------------------------------------------------------------------------*/

namespace CQKnPClass_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** Solver of Continuous Quadratic Knapsack Problems (CQKnP) based on calls
    to the \c Cplex callable library.

    Derives from CQKnpClass and therefore it conforms to its interface. */

class CQKnPCplex : public CQKnPClass {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

  CQKnPCplex( const int alg = 0 , const double eps = 1e-6 );

/**< The parameter alg controls which algorithm CPLEX uses to solve continuous
    quadratic knapsack problem. Possible values of this parameter are:

     0   [default] et CPLEX choose.
     1   Primal simplex;
     2   Dual simplex

     [See the documentation of CPX_PARAM_QPMETHOD in the Cplex manual for
      details.]

     The parameter Eps defines the convergence tolerance required to construct
     the solution [default value is 1e-6].

     [See the documentation of CPX_PARAM_BAREPCOMP in the Cplex manual for
      details.]

     The choices can be changed at any time with SetAlg() and SetEps(),
     respectively [see below]. */

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

   void LoadSet( const int pn = 0 ,
		 const double *pC = 0 , const double *pD = 0 ,
		 const double *pA = 0 , const double *pB = 0 ,
		 const double pV = 0 , const bool sns = true );

/*--------------------------------------------------------------------------*/

   inline void SetCplexParam( int whichparam , int value );

   inline void SetCplexParam( int whichparam , double value );

/**< The first two methods allow to set the many algorithmic parameters of
   Cplex; see the documentation of CPXsetintparam() and CPXsetdblparam() in
   the Cplex manual for details. */

/*--------------------------------------------------------------------------*/

   inline void SetAlg( const int alg = 0 );

/**< Allows to change the algorithm  to be used in the next calls to
    SolveKNP(); see the comments to the constructor for details. */

/*--------------------------------------------------------------------------*/

   inline void SetEps( double eps = 1e-6 );

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

   CQKnPClass::CQKStatus SolveKNP( void );

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   const double *KNPGetX( void );

   double KNPGetPi( void );

   double KNPGetFO( void );

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   void KNPLCosts( double *CostC , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() );

   void KNPQCosts( double *CostD , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() );

   inline double KNPLCost( const int i );

   inline double KNPQCost( const int i );

   void KNPLBnds( double *bnds , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() );

   void KNPUBnds( double *bnds , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() );

   inline double KNPLBnd( const int i );

   inline double KNPUBnd( const int i );

   inline double KNPVlm( void );

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

   void ChgLCosts( const double *csts , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() );

   void ChgQCosts( const double *csts , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() );

   void ChgLBnds( const double *bnds , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() );

   void ChgUBnds( const double *bnds, const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() );

/*--------------------------------------------------------------------------*/

   inline void ChgLCost( int item , const double cst );

   inline void ChgQCost( int item , const double cst );

   inline void ChgLBnd( int item , const double bnd );

   inline void ChgUBnd( int item , const double bnd );

   inline void ChgVlm( const double NVlm );

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   ~CQKnPCplex();

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*------------------------ PRIVATE DATA STRUCTURES  ------------------------*/
/*--------------------------------------------------------------------------*/

  #if( CQKnPClass_LOG )
   inline void Log1( void  );
  #endif

/*--------------------------------------------------------------------------*/
/*------------------------ PRIVATE DATA STRUCTURES  ------------------------*/
/*--------------------------------------------------------------------------*/

   double *XSol;               ///< primal solution

   static CPXENVptr cplexEnv;  ///< Cplex environment (static)
   static int users;           ///< number of active instances
   CPXLPptr lp;                ///< Cplex problem

 }; // end( class CQKnPCplex )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void CQKnPCplex::SetCplexParam( int whichparam , int value )
{
 CPXsetintparam( cplexEnv , whichparam , value );
 }

/*---------------------------------------------------------------------------*/

inline void CQKnPCplex::SetCplexParam( int whichparam , double value )
{
 CPXsetdblparam( cplexEnv , whichparam , value );
 }

/*---------------------------------------------------------------------------*/

inline void CQKnPCplex::SetAlg( const int alg ) {

 switch ( alg ) {
  case ( 1 ): SetCplexParam( CPX_PARAM_QPMETHOD , CPX_ALG_PRIMAL ); break;
  case ( 2 ): SetCplexParam( CPX_PARAM_QPMETHOD , CPX_ALG_DUAL ); break;
  default: SetCplexParam( CPX_PARAM_QPMETHOD , CPX_ALG_AUTOMATIC ); break;
  }
 }

/*---------------------------------------------------------------------------*/

inline void CQKnPCplex::SetEps( const double eps ) {

 SetCplexParam( CPX_PARAM_BAREPCOMP , eps );
 }

/*--------------------------------------------------------------------------*/

inline double CQKnPCplex::KNPLCost( const int i )
{
 double c;
 CPXgetcoef( cplexEnv , lp , -1 , i , &c );
 return( c );
 }

/*---------------------------------------------------------------------------*/

inline double CQKnPCplex::KNPQCost( const int i )
{
 double d;
 CPXgetqpcoef( cplexEnv , lp , i , i , &d );
 d *= double( 0.5 );
 return( d );
 }

/*--------------------------------------------------------------------------*/

inline double CQKnPCplex::KNPLBnd( const int i )
{
 double a;
 CPXgetlb( cplexEnv , lp , &a , i , i );
 return( a );
 }

/*--------------------------------------------------------------------------*/

inline double CQKnPCplex::KNPUBnd( const int i )
{
 double b;
 CPXgetub( cplexEnv , lp , &b , i , i );
 return( b );
 }

/*--------------------------------------------------------------------------*/

inline double CQKnPCplex::KNPVlm( void )
{
 double rhs;
 CPXgetrhs( cplexEnv , lp , &rhs , 0 , 0 );
 return( rhs );
 }

/*--------------------------------------------------------------------------*/

inline void CQKnPCplex::ChgLCost( int item , const double cst )
{
 CPXchgcoef( cplexEnv , lp , -1 , item , cst );
 }

/*--------------------------------------------------------------------------*/

inline void CQKnPCplex::ChgQCost( int item , const double cst )
{
 CPXchgqpcoef( cplexEnv , lp , item , item , cst * double(2) );
 }

/*--------------------------------------------------------------------------*/

inline void CQKnPCplex::ChgLBnd( int item , const double bnd )
{
 const char lu = 'L';
 double lb = ( bnd <= -Inf<double>() ? -CPX_INFBOUND : bnd );
 CPXchgbds( cplexEnv , lp , 1 , &item , &lu , &lb );
 }

/*--------------------------------------------------------------------------*/

inline void CQKnPCplex::ChgUBnd( int item , const double bnd )
{
 const char lu = 'U';
 double ub = ( bnd >= Inf<double>() ? CPX_INFBOUND : bnd );
 CPXchgbds( cplexEnv , lp , 1 , &item , &lu , &ub );
 }

/*--------------------------------------------------------------------------*/

inline void CQKnPCplex::ChgVlm( double NVlm )
{
 int index = 0;
 CPXchgrhs( cplexEnv , lp , 1 , &index , &NVlm );
 }

/*--------------------------------------------------------------------------*/

 };  // end( namespace CQKnPClass_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* CQKnPCplex.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File CQKnPCplex.h --------------------------*/
/*--------------------------------------------------------------------------*/
