/*--------------------------------------------------------------------------*/
/*----------------------------- File CQKnPClone.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of implementation of the template class CQKnPClone.
 *
 * \version 1.00
 *
 * \date 23 - 12 - 2012
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Operations Research Group \n
 *         Dipartimento di Elettronica Informatica e Sistemistica \n
 *         Universita' della Calabria \n
 *
 * Copyright(C) 1992 - 2012 Antonio Frangioni
 */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __CQKnPClone
 #define __CQKnPClone  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CQKnPClass.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS CQKnPClone ------------------------------*/
/*--------------------------------------------------------------------------*/
/** Class for cross-testing solvers of the Continuous Quadratic Knapsack
    Problem implemented as derived object of the class CQKnPClass [see
    CQKnPClass.h].

    The class is template over two objects, a "Master" one and a "Slave" one,
    that are assumed to be both derived from CQKnPClass; it derives from the
    Master class, and creates an object of the Slave class. Every call to a
    method is "reflected" on both objects. In some cases, checks are performed
    that the results agree.

    The rationale for having CQKnPClone to derive from the Master class is
    that some methods (typically, these that read back the data of the problem)
    do not need to be redefined, as being CQKnPClone derived from Master the
    Master:: implementation of the method is automatically used (of course
    this implies that no checks are done that Master and Slave agree on which
    data they have, which perhaps should be done). */

template<class Master, class Slave>
class CQKnPClone : public Master {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** Construct both the Master and the Slave object. */

   CQKnPClone( void ) : Master() {
    SlvCQK = new Slave;
    }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** Load data in both the Master and the Slave object. */

   void LoadSet( const int pn = 0 ,
		 const double *pC = 0 , const double *pD = 0 ,
		 const double *pA = 0 , const double *pB = 0 ,
		 const double pV = 0 , const bool sns = true ) {
    Master::LoadSet( pn , pC , pD , pA , pB , pV , sns );
    SlvCQK->LoadSet( pn , pC , pD , pA , pB , pV , sns );
    }

/*--------------------------------------------------------------------------*/

#if CQKnPClass_LOG
   void SetKNPLog( std::ostream *outs = 0 , const char lvl = 0 ) {
    Master::SetKNPLog( outs , lvl );
    SlvCQK->SetKNPLog( outs , lvl );
    }
#endif

/*--------------------------------------------------------------------------*/
/** Set tolerance in both the Master and the Slave object. */

   void SetEps( const double eps = 1e-6 ) {
    Master::SetEps( eps );
    SlvCQK->SetEps( eps );
    }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/
/** Solve the problem in both the Master and the Slave object, check that the
    solution status agrees. */

   CQKnPClass::CQKStatus SolveKNP( void ) {
    CQKnPClass::CQKStatus status1 = Master::SolveKNP( );
    CQKnPClass::CQKStatus status2 = SlvCQK->SolveKNP( );
    if( status1 != status2 )
     throw( CQKnPClass::CQKException( 
                                 "CQKnPClone::SolveKNP: different status" ) );
    return( status1 );
    }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** Recover optimal value in both the Master and the Slave object, check that
    they agree. */

   double KNPGetFO( void ) {
    const double F1 = Master::KNPGetFO();
    const double F2 = SlvCQK->KNPGetFO();
    const double epsilon = 1e-6 * max( ABS( F1 ) , double( 1 ) );
    if( ABS( F1 - F2 ) >= epsilon )
     throw( CQKnPClass::CQKException( 
                                "CQKnPClone::KNPGetFO: different values" ) );
    return( Master::KNPGetFO() );
    }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR CHANGING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** Change linear costs in both the Master and the Slave object. */

   void ChgLCosts( const double *csts , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() ) {
    Master::ChgLCosts( csts , nms , strt , stp );
    SlvCQK->ChgLCosts( csts , nms , strt , stp );
    }

/*--------------------------------------------------------------------------*/
/** Change quadratic costs in both the Master and the Slave object. */

   void ChgQCosts( const double *csts , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() ) {
    Master::ChgQCosts( csts , nms , strt , stp );
    SlvCQK->ChgQCosts( csts , nms , strt , stp );
    }

/*--------------------------------------------------------------------------*/
/** Change a linear cost in both the Master and the Slave object. */

   void ChgLCost( int i , const double cst ) {
    Master::ChgLCost( i , cst );
    SlvCQK->ChgLCost( i , cst );
    }

/*--------------------------------------------------------------------------*/
/** Change a quadratic cost in both the Master and the Slave object. */

   void ChgQCost( int i , const double cst ) {
    Master::ChgQCost( i , cst );
    SlvCQK->ChgQCost( i , cst );
    }

/*--------------------------------------------------------------------------*/
/** Change lower bounds in both the Master and the Slave object. */

   void ChgLBnds( const double *bnds , const int *nms = 0 ,
		  int strt = 0 , int stp = Inf<int>() ) {
    Master::ChgLBnds( bnds , nms , strt , stp );
    SlvCQK->ChgLBnds( bnds , nms , strt , stp );
    }

/*--------------------------------------------------------------------------*/
/** Change upper bounds in both the Master and the Slave object. */

   void ChgUBnds( const double *bnds , const int *nms = 0 ,
		  int strt = 0 , int stp = Inf<int>() ) {
    Master::ChgUBnds( bnds , nms , strt , stp );
    SlvCQK->ChgUBnds( bnds , nms , strt , stp );
    }

/*--------------------------------------------------------------------------*/
/** Change a lower bound in both the Master and the Slave object. */

   void ChgLBnd( int i , const double bnd )  {
    Master::ChgLBnd( i , bnd );
    SlvCQK->ChgLBnd( i , bnd );
    }

/*--------------------------------------------------------------------------*/
/** Change an upper bound in both the Master and the Slave object. */

   void ChgUBnd( int i , const double bnd ) {
    Master::ChgUBnd( i , bnd );
    SlvCQK->ChgUBnd( i , bnd );
    }

/*--------------------------------------------------------------------------*/
/** Change volume in both the Master and the Slave object. */

   void ChgVlm( const double NVlm ) {
    Master::ChgVlm( NVlm );
    SlvCQK->ChgVlm( NVlm );
    }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Destroy the Slave object together with the Master. */

   ~CQKnPClone() {
    delete SlvCQK;
    }

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

   Slave *SlvCQK;     /**< pointer to the "slave" object; it is public in
			 order to allow calling the methods of the specialized
			 interface of the slave class */

/*--------------------------------------------------------------------------*/

 };   // end( class CQKnPClone )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* CQKnPClone.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File CQKnPClone.h -----------------------------*/
/*--------------------------------------------------------------------------*/
