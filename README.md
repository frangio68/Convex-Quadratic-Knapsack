# Convex-Quadratic-Knapsack
Solution algorithms for Convex Quadratic Separable Continuous Knapsack Problems

This is the CQKnPClass project, a small collection of solvers for Continuous
(Convex, Separable) Quadratic Knapsack Problems (CQKnP). The project comprises
the abstract base class CQKnPClass and a few derived classes:

-  CQKnPClone implementing a "fake" CQKnP solver that takes two "real" ones and
   does everything on both; useful for testing the solvers (either for correctness
   or for efficiency) when used within "complex" approaches

-  CQKnPCplex implementing a CQKnP solver conforming to the CQKnPClass interface
   based on calls to the commercial (but now free for academic purposes)
   Cplex solver from IBM
   
-  DualCQKnP implementing a CQKnP solver based on the standard dual-ascent
   approach, restricted to instances where \e all items have strictly
   positive quadratic costs and finite bounds (both lower and upper), so
   in fact it does not fully implement the interface
   
-  ExDualCQKnP, derived from DualCQKnP and extending it to support all the
   features of the problem (possibly zero quadratic costs, possibly infinite
   upper and lower bounds), at the cost of being slightly less efficient, so
   DualCQKnP should be preferred for the instances that it can solve.

See doc/ for the Doxygen documentation of the code.

To install the file, go into Main and type "make". If you want to use the
CQKnPCplex class, which comes commented out by default, uncomment it from
lib/makefile-o and edit extlib/makefile-libCPX to insert the right Cplex path
libraries.

Two test Main files are provided into Main/, see the documentation for details.

More information about the implemented algorithms can be found at

  http://pages.di.unipi.it/frangio/abstracts.html#EJOR13
  




