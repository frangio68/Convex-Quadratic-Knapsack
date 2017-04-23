##############################################################################
#                                                                            #
#   makefile of CQKnPSolve                                                   #
#                                                                            #
# can use two different main files: Main.C and MainRnd.C                     #
#                                                                            #
##############################################################################

# module
NAME = CQKnPSolve

# basic directory
DIR = ./

# debug switches
#SW = -g -I$(DIR)

# production switches
SW = -O3 -I$(DIR)

# libreries
LIB = -lm 

# compiler
CC = g++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

all default: $(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)*.o $(DIR)*~ $(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *LIB is the final object/library
# *H   is the list of all include files
# *INC is the -I<include directories>

# the whole libCQK, with all required external libreries
libCQKDIR = ../
include $(libCQKDIR)lib/makefile-c

# Main
#MAIN = MainRnd
MAIN = Main

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

$(NAME): $(libCQKOBJ) $(MAIN).o
	$(CC) -o $(NAME) $(MAIN).o $(libCQKOBJ) $(libCQKLIB) $(LIB) $(SW)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(DIR)$(MAIN).o: $(DIR)$(MAIN).C $(libCQKH) $(OPTH)
	$(CC) -c $*.C -o $@ $(libCQKINC) $(OPTINC) $(SW)

############################ End of makefile #################################