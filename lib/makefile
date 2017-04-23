##############################################################################
################################# makefile ###################################
##############################################################################
#                                                                            #
#   makefile of libCQK.a                                                     #
#   all modules are given in input                                           #
#                                                                            #
#   Input:  $(CC)        = compiler command	                             #
#           $(SW)        = compiler options                                  #
#           $(OPTH)      = the OPT****.h files to include                    #
#           $(OPTINC)    = the -I$(OPT****.h directory)                      #
#           $(libCPXINC) = the -I$(include directories)	for libCPX           #
#           $(libCQKDIR) = the basic directory for the library               #
#                                                                            #
#   Output: $(libCQKOBJ) = the library (that must be built)                  #
#           $(libCQKLIB) = the library + external libreries + -L<libdirs>    #
#           $(libCQKH)   = all the .h files of the CQK library               #
#           $(libCQKINC) = all the -I$(include dirs) of the CQK library      #
#                                                                            #
#                                VERSION 2.00                                #
#                               17 - 05 - 2012                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                            Operations Research Group                       #
#                           Dipartimento di Informatica	                     #
#                            Universita' di Pisa                             #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(libCQKOBJ)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object/library
# *LIB external libreries + -L<libdirs>
# *H   is the list of all include files
# *INC is the -I<include directories>

# libCQK - outputs the list of OBJ
include $(libCQKDIR)lib/makefile-o

# main module (archiving phase) - - - - - - - - - - - - - - - - - - - - - - -

# this is local
lcllibCQKOBJ := $(libCQKOBJ)

# chamnge module name
libCQKOBJ := $(libCQKDIR)lib/libCQK.a

$(libCQKOBJ): $(lcllibCQKOBJ)
	ar -ruv $(libCQKOBJ) $(lcllibCQKOBJ)

############################ End of makefile #################################
