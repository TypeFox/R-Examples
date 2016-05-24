#------------------------------------------------------------------------------#
#                             R interface to GLPK                              #
#------------------------------------------------------------------------------#

#  glpk.R
#  R interface to GLPK.
#
#  Copyright (C) 2011-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of glpkAPI.
#
#  GlpkAPI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  GlpkAPI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with glpkAPI  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
#                       global variables (from glpk.h [4.42])                  #
#------------------------------------------------------------------------------#


# simplex integer parameters (not from glpk.h, self defined)
MSG_LEV  <- 101
METH     <- 102
PRICING  <- 103
R_TEST   <- 104
IT_LIM   <- 105
TM_LIM   <- 106
OUT_FRQ  <- 107
OUT_DLY  <- 108
PRESOLVE <- 109

# simplex double parameters (not from glpk.h, self defined)
TOL_BND  <- 201
TOL_DJ   <- 202
TOL_PIV  <- 203
OBJ_LL   <- 204
OBJ_UL   <- 205

# additional interior integer parameters (not from glpk.h, self defined)
ORD_ALG  <- 301

# basis factorization integer control parameters (not from glpk.h, self defined)
TYPE     <- 401
LU_SIZE  <- 402
PIV_LIM  <- 403
SUHL     <- 404
NFS_MAX  <- 405
NRS_MAX  <- 406
RS_SIZE  <- 407

# basis factorization double control parameters (not from glpk.h, self defined)
PIV_TOL  <- 501
EPS_TOL  <- 502
MAX_GRO  <- 503
UPD_TOL  <- 504

# MIP integer parameters (not from glpk.h, self defined)
BR_TECH  <- 601
BT_TECH  <- 602
PP_TECH  <- 603
FP_HEUR  <- 604
GMI_CUTS <- 605
MIR_CUTS <- 606
COV_CUTS <- 607
CLQ_CUTS <- 608
CB_SIZE  <- 609
BINARIZE <- 610

# MIP double parameters (not from glpk.h, self defined)
TOL_INT  <- 701
TOL_OBJ  <- 702
MIP_GAP  <- 703

# use callback routine glpkCallback in MIP
CB_FUNC  <- 651


#------------------------------------------------------------------------------#

# LP/MIP problem object
# optimization direction flag:
GLP_MIN <- 1  # minimization
GLP_MAX <- 2  # maximization

# kind of structural variable:
GLP_CV <- 1  # continuous variable
GLP_IV <- 2  # integer variable
GLP_BV <- 3  # binary variable

# type of auxiliary/structural variable:
GLP_FR <- 1  # free variable
GLP_LO <- 2  # variable with lower bound
GLP_UP <- 3  # variable with upper bound
GLP_DB <- 4  # double-bounded variable
GLP_FX <- 5  # fixed variable

# status of auxiliary/structural variable:
GLP_BS <- 1  # basic variable
GLP_NL <- 2  # non-basic variable on lower bound
GLP_NU <- 3  # non-basic variable on upper bound
GLP_NF <- 4  # non-basic free variable
GLP_NS <- 5  # non-basic fixed variable

# scaling options:
GLP_SF_GM   <- 0x01  # perform geometric mean scaling
GLP_SF_EQ   <- 0x10  # perform equilibration scaling
GLP_SF_2N   <- 0x20  # round scale factors to power of two
GLP_SF_SKIP <- 0x40  # skip if problem is well scaled
GLP_SF_AUTO <- 0x80  # choose scaling options automatically

# solution indicator:
GLP_SOL <- 1  # basic solution
GLP_IPT <- 2  # interior-point solution
GLP_MIP <- 3  # mixed integer solution

# solution status:
GLP_UNDEF  <- 1  # solution is undefined
GLP_FEAS   <- 2  # solution is feasible
GLP_INFEAS <- 3  # solution is infeasible
GLP_NOFEAS <- 4  # no feasible solution exists
GLP_OPT    <- 5  # solution is optimal
GLP_UNBND  <- 6  # solution is unbounded


#------------------------------------------------------------------------------#
# basis factorization control parameters
# type
GLP_BF_FT  <- 0x01  # LUF + Forrest-Tomlin
GLP_BF_BG  <- 0x02  # LUF + Schur compl. + Bartels-Golub
GLP_BF_GR  <- 0x03  # LUF + Schur compl. + Givens rotation
GLP_BF_LUF <- 0x00  # plain LU-factorization
GLP_BF_BTF <- 0x10  # block triangular LU-factorization

#------------------------------------------------------------------------------#
# simplex method control parameters
# msg_lev           # message level:
GLP_MSG_OFF <- 0    # no output
GLP_MSG_ERR <- 1    # warning and error messages only
GLP_MSG_ON  <- 2    # normal output
GLP_MSG_ALL <- 3    # full output
GLP_MSG_DBG <- 4    # debug output

# meth              # simplex method option:
GLP_PRIMAL <- 1     # use primal simplex
GLP_DUALP  <- 2     # use dual; if it fails, use primal
GLP_DUAL   <- 3     # use dual simplex

# pricing           # pricing technique:
GLP_PT_STD <- 0x11  # standard (Dantzig rule)
GLP_PT_PSE <- 0x22  # projected steepest edge

# r_test            # ratio test technique:
GLP_RT_STD <- 0x11  # standard (textbook)
GLP_RT_HAR <- 0x22  # two-pass Harris' ratio test


#------------------------------------------------------------------------------#
# interior-point solver control parameters
# ord_alg            # ordering algorithm:
GLP_ORD_NONE   <- 0  # natural (original) ordering
GLP_ORD_QMD    <- 1  # quotient minimum degree (QMD)
GLP_ORD_AMD    <- 2  # approx. minimum degree (AMD)
GLP_ORD_SYMAMD <- 3  # approx. minimum degree (SYMAMD)


#------------------------------------------------------------------------------#
# integer optimizer control parameters
# br_tech        # branching technique:
GLP_BR_FFV <- 1  # first fractional variable
GLP_BR_LFV <- 2  # last fractional variable
GLP_BR_MFV <- 3  # most fractional variable
GLP_BR_DTH <- 4  # heuristic by Driebeck and Tomlin
GLP_BR_PCH <- 5  # hybrid pseudocost heuristic

# bt_tech        # backtracking technique:
GLP_BT_DFS <- 1  # depth first search
GLP_BT_BFS <- 2  # breadth first search
GLP_BT_BLB <- 3  # best local bound
GLP_BT_BPH <- 4  # best projection heuristic

# pp_tech         # preprocessing technique:
GLP_PP_NONE <- 0  # disable preprocessing
GLP_PP_ROOT <- 1  # preprocessing only on root level
GLP_PP_ALL  <- 2  # preprocessing on all levels


#------------------------------------------------------------------------------#
# additional row attributes
# the row origin flag:
GLP_RF_REG  <- 0  # regular constraint
GLP_RF_LAZY <- 1  # "lazy" constraint
GLP_RF_CUT  <- 2  # cutting plane constraint

# the row class descriptor:
# klass
GLP_RF_GMI <- 1  # Gomory's mixed integer cut
GLP_RF_MIR <- 2  # mixed integer rounding cut
GLP_RF_COV <- 3  # mixed cover cut
GLP_RF_CLQ <- 4  # clique cut


#------------------------------------------------------------------------------#
# enable/disable flag:
GLP_ON  <- 1  # enable something
GLP_OFF <- 0  # disable something


#------------------------------------------------------------------------------#
# reason codes:
GLP_IROWGEN <- 0x01  # request for row generation
GLP_IBINGO  <- 0x02  # better integer solution found
GLP_IHEUR   <- 0x03  # request for heuristic solution
GLP_ICUTGEN <- 0x04  # request for cut generation
GLP_IBRANCH <- 0x05  # request for branching
GLP_ISELECT <- 0x06  # request for subproblem selection
GLP_IPREPRO <- 0x07  # request for preprocessing


#------------------------------------------------------------------------------#
# branch selection indicator:
GLP_NO_BRNCH <- 0  # select no branch
GLP_DN_BRNCH <- 1  # select down-branch
GLP_UP_BRNCH <- 2  # select up-branch


#------------------------------------------------------------------------------#
# return codes:
GLP_EBADB   <- 0x01  # invalid basis
GLP_ESING   <- 0x02  # singular matrix
GLP_ECOND   <- 0x03  # ill-conditioned matrix
GLP_EBOUND  <- 0x04  # invalid bounds
GLP_EFAIL   <- 0x05  # solver failed
GLP_EOBJLL  <- 0x06  # objective lower limit reached
GLP_EOBJUL  <- 0x07  # objective upper limit reached
GLP_EITLIM  <- 0x08  # iteration limit exceeded
GLP_ETMLIM  <- 0x09  # time limit exceeded
GLP_ENOPFS  <- 0x0A  # no primal feasible solution
GLP_ENODFS  <- 0x0B  # no dual feasible solution
GLP_EROOT   <- 0x0C  # root LP optimum not provided
GLP_ESTOP   <- 0x0D  # search terminated by application
GLP_EMIPGAP <- 0x0E  # relative mip gap tolerance reached
GLP_ENOFEAS <- 0x0F  # no primal/dual feasible solution
GLP_ENOCVG  <- 0x10  # no convergence
GLP_EINSTAB <- 0x11  # numerical instability
GLP_EDATA   <- 0x12  # invalid data
GLP_ERANGE  <- 0x13  # result out of range


#------------------------------------------------------------------------------#
# condition indicator:
GLP_KKT_PE <- 1  # primal equalities
GLP_KKT_PB <- 2  # primal bounds
GLP_KKT_DE <- 3  # dual equalities
GLP_KKT_DB <- 4  # dual bounds
GLP_KKT_CS <- 5  # complementary slackness


#------------------------------------------------------------------------------#
# MPS file format:
GLP_MPS_DECK <- 1  # fixed (ancient)
GLP_MPS_FILE <- 2  # free (modern)


#------------------------------------------------------------------------------#

# return codes of glpk optimizations
return_codeGLPK <- function(code) {
    if (code == 0)                { return( "solution process was successful" ) }
    else if (code == GLP_EBADB)   { return( "invalid basis" ) }
    else if (code == GLP_ESING)   { return( "singular matrix" ) }
    else if (code == GLP_ECOND)   { return( "ill-conditioned matrix" ) }
    else if (code == GLP_EBOUND)  { return( "invalid bounds" ) }
    else if (code == GLP_EFAIL)   { return( "solver failed" ) }
    else if (code == GLP_EOBJLL)  { return( "objective lower limit reached" ) }
    else if (code == GLP_EOBJUL)  { return( "objective upper limit reached" ) }
    else if (code == GLP_EITLIM)  { return( "iteration limit exceeded" ) }
    else if (code == GLP_ETMLIM)  { return( "time limit exceeded" ) }
    else if (code == GLP_ENOPFS)  { return( "no primal feasible solution" ) }
    else if (code == GLP_ENODFS)  { return( "no dual feasible solution" ) }
    else if (code == GLP_EROOT)   { return( "root LP optimum not provided" ) }
    else if (code == GLP_ESTOP)   { return( "search terminated by application" ) }
    else if (code == GLP_EMIPGAP) { return( "relative mip gap tolerance reached" ) }
    else if (code == GLP_ENOFEAS) { return( "no primal/dual feasible solution" ) }
    else if (code == GLP_ENOCVG)  { return( "no convergence" ) }
    else if (code == GLP_EINSTAB) { return( "numerical instability" ) }
    else if (code == GLP_EDATA)   { return( "invalid data" ) }
    else if (code == GLP_ERANGE)  { return( "result out of range" ) }
    else { return(paste("Failed to obtain solution, unknown error code:", code)) }
}


# solution status codes of glpk optimizations
status_codeGLPK <- function(code) {
    if (code == GLP_OPT)         { return( "solution is optimal" ) }
    else if (code == GLP_UNDEF)  { return( "solution is undefined" ) }
    else if (code == GLP_FEAS)   { return( "solution is feasible" ) }
    else if (code == GLP_INFEAS) { return( "solution is infeasible" ) }
    else if (code == GLP_NOFEAS) { return( "no feasible solution exists" ) }
    else if (code == GLP_UNBND)  { return( "solution is unbounded" ) }
    else { return(paste("unknown status code:", code)) }
}

