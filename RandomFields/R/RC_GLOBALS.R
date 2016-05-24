## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
## Sebastian Gross
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  




## @FILE-STARP******************************************************************
# @NAME		ZF_GLOBALS
# @DESCRIPTION	Any value that is used throughout the randomfield package
#               has to appear here
# @AUTHOR	Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
#               Martin Schlather
# @DATE		26.08.2011 (Gross), 2012 -- 2013 (Schlather)
#
# @FILE-END*********************************************************************




###############################################################################
##                        DEFINITIONS OF SYMBOLS                             ##
###############################################################################

# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_SYMBOLS_PLUS/MAL/SYMBOLS
# @DESCRIPTION	The + operator in any valid model formula
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @GLOBAL-END*******************************************************************
ZF_SYMBOLS_PLUS <- '+'
ZF_PLUS <- c("RMplus", ZF_SYMBOLS_PLUS)

ZF_SYMBOLS_MULT <- '*'
ZF_MULT <- c("RMmult", ZF_SYMBOLS_MULT)

## Special Models
DOLLAR <- c("$", "RMS")
ZF_DOLLAR <- rev(DOLLAR)

ZF_SYMBOLS_C <- "R.c"
ZF_SYMBOLS_CONST <- "R.const"

ZF_COVARIATE <- "RMcovariate"

ZF_CARTCOORD_NAMES <- c("x", "y", "z", "T")
ZF_GENERAL_COORD_NAME <- c("coords.x", "coords.T")## check general_coordinates if changed
ZF_EARTHCOORD_NAMES <- c("longitude", "latitude", "height", "time")


# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_SYMBOLS_AT
# @DESCRIPTION	The former @ operator in any valid model formula used to create fixed effects must ffs not be "*"
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @GLOBAL-END*******************************************************************
ZF_SYMBOLS_AT <- "@"

# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_SYMBOLS_L_PAR
# @DESCRIPTION	Left parenthesis
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @GLOBAL-END*******************************************************************
ZF_SYMBOLS_L_PAR <- "("

# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_SYMBOLS_R_PAR
# @DESCRIPTION	Right parenthesis
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @GLOBAL-END*******************************************************************
ZF_SYMBOLS_R_PAR <- ")"



###############################################################################
##                        DEFINITIONS OF MODELNAMES                          ##
###############################################################################

# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_FIXED and other names of special models
# @DESCRIPTION	The function name of fixed effects
# @AUTHOR       Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
#               Martin Schlather
# @DATE		29.08.2011 (Gross) and 2013--2013 (Schlather)
# @GLOBAL-END*******************************************************************
ZF_FIXED <- "RMfixed"
ZF_INTERNALMIXED <- "internalRMmixed"
ZF_TREND <- c("RMtrend", "trend")
ZF_TRENDFCT <- paste(ZF_TREND[1], "(", sep="")
ZF_DISTR <- c('RRdistr', 'Distr')
ZF_USER <- c('RMuser', 'U')
ZF_COORD <- "RMcoord"
ZF_MODEL <- "RMmodel"

ZF_MIXED <- c( "RMmixed", "mixed") 
ZF_NUGGET <- c("RMnugget", "nugget")

ZF_MODELEXT <- "RMmodelFit"


# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_MODEL_FACTORY
# @DESCRIPTION	Each covariance model is an object of this class
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		29.08.2011
# @GLOBAL-END*******************************************************************
ZF_MODEL_FACTORY <- "RMmodelgenerator"

# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_DEFAULT_STRING
# @DESCRIPTION	var, scale, Aniso, proj  get this value assigned by functions
#               of class 'RMmodelgenerator', if no such arguments are passed
#               by the user. NULL is passed to C-level
# @AUTHOR       A Malinowski <malinows@math.uni-goettingen.de>
# @DATE		29.08.2011
# @GLOBAL-END*******************************************************************
ZF_DEFAULT_STRING <- "RFdefault"
ZF_MODEL_PREFIX <- "RM"

isPosDef <- function(type) {
  if (is.character(type)) type <- pmatch(type, TYPENAMES, duplicates.ok=TRUE)-1
  ##  .C("isPosDef", as.integer(type))$type
  type==TcfType | type == PosDefType | type==UndefinedType
}
isVariogram <- function(type) { 
  if (is.character(type)) type <- pmatch(type, TYPENAMES, duplicates.ok=TRUE)-1
  ##  .C("isNefDef", as.integer(type))$type
  isPosDef(type) | type == VariogramType
}




###############################################################################
##                              OTHER DEFINITIONS                            ##
###############################################################################

LSQMETHODS <- c("self", "plain", "sqrt.nr", "sd.inv", "internal") 
MLMETHODS <- c("ml") # "reml", "rml1"),

par.storage <- ".RandomFields.par"
par.storage.env <- .GlobalEnv



 PL_IMPORTANT 	<- as.integer(1)
 PL_SUBIMPORTANT 	<- as.integer(2)
 PL_DETAILSUSER <- as.integer(3)## currently unused
 PL_RECURSIVE 	<- as.integer(4)
 PL_STRUCTURE 	<- as.integer(5)
 PL_ERRORS 	<- as.integer(6)

 PL_FCTN_DETAILS 	<- as.integer(7)
 PL_FCTN_SUBDETAILS 	<- as.integer(8)

 PL_COV_STRUCTURE 	<- as.integer(7)
 PL_DIRECT_SEQU 	<- as.integer(8)
 PL_DETAILS 	<- as.integer(9)
 PL_SUBDETAILS 	<- as.integer(10)
