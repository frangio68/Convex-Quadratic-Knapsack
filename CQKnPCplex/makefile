##############################################################################
#                                                                            #
#   makefile of CQKnPCplex                                                   #
#   all modules are given in input                                           #
#                                                                            #
#   Input:  $(CC)        = compiler command                                  #
#           $(SW)        = compiler options                                  #
#           $(CQKH)      = the include file(s) for CQKnPClass.h              #
#           $(CQKINC)    = the -I$(CQKnPClass directory)                     #
#           $(libCPXINC) = the -I$(include directories)	for libCPX           #
#           $(CPKDIR)    = the directory where the source is                 #
#                                                                            #
#   Output: $(CPKOBJ)    = the final object file(s)                          #
#           $(CPKLIB)    = external libreries + -L<libdirs>		     #
#           $(CPKH)      = the .h files to include                           #
#           $(CPKINC)    = the -I$(include directories)                      #
#                                                                            #
#                                VERSION 1.10	                             #
#                               23 - 04 - 2017                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                           Operations Research Group                        #
#                          Dipartimento di Informatica	                     #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(CPKDIR)*.o $(CPKDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

CPKOBJ = $(CPKDIR)CQKnPCplex.o
CPKLIB = $(libCPXLIB)
CPKH   = $(CPKDIR)CQKnPCplex.h
CPKINC = -I$(CPKDIR) $(libCPXINC)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(CPKOBJ): $(CPKDIR)CQKnPCplex.C $(CPKH) $(CQKH)
	$(CC) -c $*.C -o $@ $(libCPXINC) $(CQKINC) $(SW)

############################ End of makefile #################################
