################################################################################
adsimp <- function( ND, VRTS, NF, F, MXFS, EA, ER, KEY, partitionInfo=FALSE ){
#
# Translation of Alan Genz's matlab file adsimp.m to R.
# See below for his description of the program, input and output values,
#   contact info and copyright information.
# Translated by John Nolan, jpnolan@american.edu, 18 Aug. 2014
#
# The R version has one change:
#    Added new argument partitionInfo; if FALSE, it is ignored. If TRUE, the 
#    function returns the following information about the subdivision of the 
#    original simplices.
#       VRTS an (ND x (ND+1) x SBS) array, final partition of the simplices:
#           VRTS[,,i] contains the vertices (as columns) of subsimplex i
#       VLS an (NF x SBS) matrix of estimated integrals: VLS[,i] is
#           the estimated value of the integral over subsimplex i
#       AES an (NF x SBS) matrix of estimated abs. error: AES[,i]
#           is the estimated abs. error on subsimplex i
#       VOL a vector of length SBS: VOL[i] is the volume of
#           subsimplex i
#
# original matlab function call:
# [VL, AE, NV, FL] = adsimp( ND, VRTS, NF, F, MXFS, EA, ER, KEY )
#
#***KEYWORDS automatic multidimensional integrator,
#            n-dimensional simplex, general purpose, global adaptive
#***PURPOSE  To calculate an approximation to a vector of integrals
#              over a collection of simplices
#
#               I = I (F ,F ,...,F   ) DS
#                    S  1  2      NF
#
#            where S is a collection ND-dimensional simplices,
#            and  F = F (X ,X ,...,X  ), K = 1 : 2, ..., NF,
#                  K   K  1  2      ND
#            and try to satisfy for each component I(K) of I
#              abs( I(K) - VL(K) ) < max( EA, ER*abs(I(K)) )
#
#  EXAMPLE: n = 4; v = [ zeros(n,1) eye(n) ];
#     [r e] = adsimp(n,v,2,@(x)[x(1)^3;x(3)^4],1000,1e-5); disp([r'; e'])
#

#***DESCRIPTION Computation of integrals over simplical regions.
#            adsimp is a driver for the integration routine SMPSAD,
#            which repeatedly subdivides the region of integration and
#            estimates the integrals and the errors over the subregions
#            with greatest estimated errors until the error request
#            is met or MXFS function evaluations have been used.
#
#   ON ENTRY
#
#     ND     Integer, number of variables. 1 < ND
#     NF     Integer, number of components of the integral.
#     MXFS   Integer, maximum number of F values.
#            RULCLS is number F values for each subregion (see KEY),
#            MXFS must be >= SBS*RULCLS.
#     F      Function  for computing components of the integrand.
#            Input parameter: X, an array that defines the evaluation point.
#            If NF > 1 component, F must output an NFx1 vector.
#     EA Real requested absolute accuracy.
#     ER Real requested relative accuracy.
#     KEY    Integer, key to selected local integration rule.
#            KEY = 0 gives the user a (default) degree 7 integration rule.
#            KEY = 1 gives the user a degree 3 integration rule.
#            KEY = 2 gives the user a degree 5 integration rule.
#            KEY = 3 gives the user a degree 7 integration rule.
#            KEY = 4 gives the user a degree 9 integration rule.
#     VRTS Real array of simplex vertices for SBS simplices;
#           the coordinates of vertex J for simplex K must be VRTS(:,J,K).
#
#   ON EXIT
#
#     VL  Real array of dimension NF of integral approximations.
#     AE  Real array of dimension NF, of absolute accuracy estimates.
#     NV Integer, number of FUNSUB calls used by adsimp.
#     FL Integer.
#            FL = 0 for normal exit, when AE(K) <=  EA or
#              AE(K) <=  abs(VL(K))*ER with MXFS or less
#              function evaluations for all values of K, 1 <= K <= NF.
#            FL = 1 if MXFS was too small for adsimp to obtain
#              the required accuracy. In this case adsimp returns
#              values VL with estimated absolute accuracies AE.
#            FL = 2 if KEY < 0 or KEY > 4,
#            FL = 3 if ND < 2,
#            FL = 4 if NF < 1,
#            FL = 5 if EA < 0 and ER < 0,
#
#***AUTHOR
#
#            Alan Genz
#            Department of Mathematics
#            Washington State University
#            Pullman, WA 99164-3113, USA
#            Email: alangenz@wsu.edu
#
#
# Copyright (C) 2013, Alan Genz,  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided the following conditions are met:
#   1. Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#   2. Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in
#      the documentation and/or other materials provided with the
#      distribution.
#   3. The contributor name(s) may not be used to endorse or promote
#      products derived from this software without specific prior
#      written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#***LAST MODIFICATION 2012-06
#***end PROLOGUE adsimp

#
#   Compute RCLS, and check the input parameters.
# 
if (KEY == 0) {KEY = 3}
if (length(dim(VRTS)) != 3) stop("VRTS must be a 3 dimensional array of simplices")
SBS <- dim(VRTS)[3] # number of original simplices

### [RCLS, FL] = SMPCHC( ND, NF, MXFS, EA, ER, SBS, KEY );
b <- SMPCHC( ND, NF, MXFS, EA, ER, SBS, KEY );

if ( b$FL != 0 ) { # error in input variables
  return( list( VL = rep(NA, NF), AE = rep(Inf, NF), NV = 0, FL=b$FL, VLS=NA, AES=NA, VOL=NA  ) )
}

return( SMPSAD( ND, NF, F, MXFS, EA, ER, KEY, b$RULCLS, SBS, VRTS, partitionInfo ) ) }
#
#***end adsimp
###############################################################################
SMPCHC <- function( ND, NF, MXFS, EA, ER, SBS, KEY ){
# function [RULCLS, FL] = SMPCHC( ND, NF, MXFS, EA, ER, SBS, KEY );   
#
#***BEGIN PROLOGUE SMPCHC
#***AUTHOR
#
#            Alan Genz
#            Department of Mathematics
#            Washington State University
#            Pullman, WA 99164-3113, USA
#
#***LAST MODIFICATION 2012-06
#***PURPOSE  SMPCHC checks validity of input parameters for adsimp.
#***DESCRIPTION
#            SMPCHC computes maxSUB, RULCLS and FL as functions of
#             input parameters for adsimp, and checks the validity of
#             input parameters for adsimp.
#
#   ON ENTRY
#
#     ND   Integer, number of variables,  ND > 1.
#     NF Integer, number of components of the integral.
#     MXFS Integer, maximum number of new F calls.
#     EA Real, requested absolute accuracy.
#     ER Real, requested relative accuracy.
#     SBS Integer, initial number of simplices.
#     KEY    Integer, key to selected local integration rule.
#            KEY = 0 gives the user a (default)degree 7 integration rule.
#            KEY = 1 gives the user a degree 3 integration rule.
#            KEY = 2 gives the user a degree 5 integration rule.
#            KEY = 3 gives the user a degree 7 integration rule.
#            KEY = 4 gives the user a degree 9 integration rule.
#
#   ON EXIT
#
#     RULCLS Integer, number of function values for each subregion.
#            If
#             KEY = 0, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
#             KEY = 1, RULCLS = 2*ND+3;
#             KEY = 2, RULCLS = (ND+3)*(ND+2)/2 + 2*(ND+1);
#             KEY = 3, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
#             KEY = 4, RULCLS = (ND+5)*(ND+4)*(ND+3)*(ND+2)/24
#                               + 5*(ND+2)*(ND+1)/2 .
#     FL Integer.
#            FL = 0 for normal exit,
#            FL = 1 if MXFS < SBS*RULCLS,
#            FL = 2 if KEY < 0 or KEY > 4,
#            FL = 3 if ND < 2,
#            FL = 4 if NF < 1,
#            FL = 5 if EA < 0 and ER < 0.

#
#     Check valid KEY.
#
      FL <- 0;
      if ( (KEY < 0) | (KEY > 4) ) { FL <- 2L }
#
#     Check valid ND.
#
      if ( ND < 2 ) { FL = 3L }
#
#     Check positive NF.
#
      if ( NF < 1 ) { FL = 4L }
#
#     Check valid accuracy requests.
#
      if ( (EA < 0) & (ER < 0) ) { FL = 5 }
#
#     Compute RULCLS as a function of KEY and ND and check MXFS.
#
      if ( FL == 0 ){
         if (KEY == 0) { RULCLS <- (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1) }
         if (KEY == 1) { RULCLS <- 2*ND + 3 }
         if (KEY == 2) { RULCLS <- (ND+3)*(ND+2)/2 + 2*(ND+1) }
         if (KEY == 3) { RULCLS <- (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1) }
         if (KEY == 4) { RULCLS <- (ND+5)*(ND+4)*(ND+3)*(ND+2)/24 + 5*(ND+2)*(ND+1)/2 }
         if ( MXFS < SBS*RULCLS ) { FL <- 1 }
      }
return( list(FL=FL,RULCLS=as.integer(RULCLS)) ) }
#
#***end SMPCHC
###############################################################################
SMPSAD <- function( ND, NF, F, MXFS, EA, ER, KEY, RCLS, SBS, VRTS, partitionInfo ){
# function [VL, AE, NV, FL] = ...
#      SMPSAD( ND, NF, F, MXFS, EA, ER, KEY, RCLS, SBS, VRTS )
#
#***BEGIN PROLOGUE SMPSAD
#***KEYWORDS automatic multidimensional integrator,
#            n-dimensional simplex,
#            general purpose, global adaptive
#***AUTHOR
#
#            Alan Genz
#            Department of Mathematics
#            Washington State University
#            Pullman, WA 99164-3113, USA
#
#***LAST MODIFICATION 2012-06
#***PURPOSE  The routine calculates an approximation to a given
#            vector of definite integrals, I, over a hyper-rectangular
#            region hopefully satisfying for each component of I the
#            following claim for accuracy:
#            abs( I(K) - VL(K) ) <= max( EA, ER*abs(I(K) ) )
#***DESCRIPTION Computation of integrals over hyper-rectangular regions.
#            SMPSAD repeatedly subdivides the regions of integration
#            and estimates the integrals and the errors over the
#            subregions with  greatest estimated errors until the error
#            request is met or maxSUB subregions are stored. The regions
#            are divided into three or four equally sized parts along
#            the direction(s) with greatest absolute fourth difference.
#
#   ON ENTRY
#
#     ND     Integer, number of variables, ND > 1.
#     NF     Integer, number of components of the integral.
#     F      Function  for computing components of the integrand.
#            Input parameter: X, an array that defines the evaluation point.
#            If NF > 1 component, F must output an NFx1 vector.
#     MXFS Integer.
#            The computations proceed until further subdivision would
#            require more than MXFS F values.
#     EA Real, requested absolute accuracy.
#     ER Real, requested relative accuracy.
#     KEY    Integer, key to selected local integration rule.
#            KEY = 1 gives the user a degree 3 integration rule.
#            KEY = 2 gives the user a degree 5 integration rule.
#            KEY = 3 gives the user a degree 7 integration rule.
#            KEY = 4 gives the user a degree 9 integration rule.
#     RCLS   Integer, number of FUNSUB calls needed for each subregion.
#     SBS    Integer, number of subregions.
#     VRTS Real array of size (ND,ND+1,SBS).
#            Simplex vertices for each subregion; for subregion K vertex
#            J must have components VERTEX(I,J,K), I = 1 : 2, ..., ND.
#   ON RETURN
#
#     VL  Real array of dimension NF.
#            Approximations to all components of the integral.
#     AE  Real array of dimension NF.
#            Estimates of absolute accuracies.
#     NV Integer, number of new FUNSUB calls used by SMPSAD.
#     FL Integer.
#            FL = 0 for normal exit, when AE(K) <=  EA or
#              AE(K) <=  abs(VL(K))*ER, 1 <= K <= NF.
#            FL = 1 if MXFS was too small for SMPSAD to obtain the
#            required accuracy. In this case SMPSAD returns values
#            of VL with estimated absolute accuracies AE.
#
#***end PROLOGUE SMPSAD
#
#   Initialize for rule parameters.
#
NV <- 0; DFCOST <- 1 + 2*ND*( ND + 1 );
#
#     Initialize NV, and VL and AE arrays.
#
VL <- rep(0.0,NF); AE <- VL;
#
#        Compute weights, generators, PTS.
#
a <-  SMPRMS( ND, KEY ); FCT <- factorial(ND);
# explictly declare VLS, AES, VOL
VLS <- matrix(0.0,nrow=NF,ncol=SBS) # VLS[i,j] = estimated integral of F[i] on simplex VRTS[,,j]
AES <- matrix(0.0,nrow=NF,ncol=SBS) # AES[i,j] = estimated abs. err. in integral of F[i] on simplex VRTS[,,j]
VOL <- rep(0.0,SBS)
for (K in 1 : SBS) {
  VOL[K] <- abs( det( VRTS[ ,1:ND,K]-VRTS[ ,ND+1,K]*rep(1.0,ND) ) )/FCT;
  #
  #     Apply basic rule over each simplex.
  #
  b <- SMPRUL( ND, VRTS[ , ,K], VOL[K], NF, F, a$G, a$W, a$PTS );
  AES[ ,K] <- b$RGNERR
  VLS[ ,K] <- b$BASVAL
  VL <- VL + VLS[ ,K]; AE <- AE + AES[ ,K]; NV <- NV + RCLS;
}
FL <- max( AE > max( c(EA,ER*abs(VL))) );
#
#     End initialisation.
#

while ( (FL > 0) & (NV + DFCOST + 4*RCLS <= MXFS) ){
  #
  #     Begin loop while error is too large, and NV is not too large.
  #
  #     Adjust VL and AE.
  #
  # select simplex with biggest abs. error
  # [GM, ID] = max(AES); if NF > 1, [GM, ID] = max(GM); end
  # VL = VL - VLS(:,ID); AE = AE - AES(:,ID);
  ID <- which.max( apply(AES, 2, max) )
  VL <- VL - VLS[,ID]; AE <- AE - AES[,ID];
  #
  #     Determine NEWSBS new subregions.
  #
  #[ VRTS, NEW ] <-  SMPDFS( ND, NF, F, ID, SBS, VRTS ); VI <- VOL(ID)/NEW;\
  b <- SMPDFS( ND, NF, F, ID, SBS, VRTS );   NEW <- b$NEW; VRTS <- b$VRTS

  VI <- VOL[ID]/NEW
  #
  #     Apply basic rule, and add new contributions to VL and AE.
  #
  # explicitly extend VOL, VLS, AES
  VOL <- c( VOL, rep(0.0,NEW-1) )
  VLS <- cbind( VLS, matrix(0.0,nrow=NF, ncol=NEW-1) )
  AES <- cbind( AES, matrix(0.0,nrow=NF, ncol=NEW-1))
  for (K in c( ID, (SBS+1):(SBS+NEW-1) )) {  
    VOL[K] <- VI;
    # [VLS(:,K), AES(:,K)] <- SMPRUL( ND, VRTS(:,:,K), VI, NF, F, G, W, PTS );
    d <- SMPRUL( ND, VRTS[ , ,K], VI, NF, F, a$G, a$W, a$PTS ); 
    VLS[ ,K] <- d$BASVAL; AES[ ,K] <- d$RGNERR
    VL <- VL + VLS[ ,K]; 
    AE <- AE + AES[ ,K]; 
    NV <- NV + RCLS;
  }
  NV <- NV + DFCOST; SBS <- SBS + NEW - 1;
  #
  #     Check for error termination.
  #
  FL <- max( AE > max( c(EA,ER*abs(VL))) );
}
#
#     Compute more accurate values of VL and AE.
#
if (SBS > 1) {
   VL <- rowSums(VLS); AE <- rowSums(AES)
}
if (partitionInfo) { result <- list(VL=VL,AE=AE,NV=NV,FL=FL,VRTS=VRTS,VLS=VLS,AES=AES,VOL=VOL ) }
else { result <- list(VL=VL,AE=AE,NV=NV,FL=FL ) }

return( result ) }
#***end SMPSAD
##############################################################################
SMPDFS <- function( ND, NF, F, TOP, SBS, VRTS ) {
#function  [VRTS, NEW] =  SMPDFS( ND, NF, F, TOP, SBS, VRTS )
#
#***BEGIN PROLOGUE SMPDFS
#***PURPOSE  To compute new subregions
#***AUTHOR
#
#            Alan Genz
#            Department of Mathematics
#            Washington State University
#            Pullman, WA 99164-3113, USA
#
#***LAST MODIFICATION 2012-06
#***DESCRIPTION SMPDFS computes fourth differences along each edge
#            direction. It uses these differences to determine a
#            subdivision of the orginal subregion into either three or
#            four new subregions.
#
#   ON ENTRY
#
#   ND   Integer, number of variables.
#   NF   Integer, number of components for the vector integrand.
#   F    Function for computing the components of the integrand at a point X.
#   TOP    Integer, location in VRTS array for original subregion.
#   SBS Integer, number of subregions in VRTS BEFORE subdivision.
#   VRTS Real array of size (ND,ND+1,*), vertices of orginal
#          subregion must be in VRTS(1:ND,1:ND+1,TOP).
#
#   ON EXIT
#
#   NEW Integer, number of new subregions (3 or 4).
#   VRTS Real array of size (ND,1+ND,*).
#          The vertices of the of new subegions will be at locations
#          TOP, SBS+1 : ..., SBS+NEW-1.
#        ***** NOTE: this function increases the size of array VRTS *****
#
CUTTF <- 2; CUTTB <- 8;
#
#       Compute the differences.
#
IS <- 1; JS <- 2; DFMX <- 0; EMX <- 0; V <- VRTS[ , ,TOP];
CN <- rowMeans(V); FC <- F(CN); DFMD <- sum( abs(FC) );
#  matlab program relied on implicit definition of FRTHDF, R requires explict def.
FRTHDF <- matrix(0.0,nrow=ND,ncol=ND+1)   
for (I in 1 : ND) {
  for (J in (I+1) : (ND+1) ) {
    H <- 2*( V[ ,I]-V[ ,J] )/( 5*(ND+1) );
    EWD <- sum(abs(H)); 
    if (EWD >= EMX) { IE <- I; JE <- J; EMX <- EWD }
    DFR <- sum( abs( F(CN-2*H)+F(CN+2*H) + 6*FC - 4*(F(CN-H)+F(CN+H)) ) );
    if ((DFMD + DFR/8) == DFMD) { DFR <- 0 }
    DFR <- DFR*EWD;
    if (DFR >= DFMX) { IT <- IS; JT <- JS; DFNX <- DFMX;
      IS <- I; JS <- J; DFMX <- DFR;
    } else { if (DFR >= DFNX) { IT <- I; JT <- J; DFNX <- DFR; }
    }
    FRTHDF[I,J] <- DFR;  
  }
}
#
#     Determine whether to compute three or four new subregions.
#
if (DFNX > DFMX/CUTTF) { NEW <- 4; }
else { 
  NEW <- 3; 
  if (DFMX == 0) { IS <- IE; JS <- JE;
  } else { DFSMX <- 0;
    for (L in 1 : (ND+1) ) {
      if (( L != IS) & (L != JS) ) { IT <- min( c(L,IS,JS) ); JT <- max( c(L, IS, JS) );
	  LT <- IS + JS + L - IT - JT; DFR <-  FRTHDF[IT,LT] + FRTHDF[LT,JT];
	  if (DFR >= DFSMX) { DFSMX <- DFR; LS <- L; }
    }
    }
    DIFIL <- FRTHDF[ min(IS,LS), max(IS,LS) ];
    DIFLJ <- FRTHDF[ min(JS,LS), max(JS,LS) ];
    DFNX <- DIFIL + DIFLJ - min( DIFIL,DIFLJ );
    if ( (DFMX/CUTTB < DFNX) & (DIFIL > DIFLJ) ) { IT <- IS; IS <- JS; JS <- IT; }
  }
}
#
#     Copy vertices and volume for TOP to new subregions
#
# original matlab code was:   (L in (SBS + 1):(SBS + NEW - 1)) { VRTS[ , ,L] <- V; }
# this INCREASES THE SIZE OF VRTS in matlab, but does not work in R.  
# Replace w/ the following 4 lines:
tmp <- VRTS
VRTS <- array( 0.0, dim=c(ND, ND+1, SBS+NEW-1 ) )
VRTS[,,1:SBS] <- tmp
for (L in (SBS+1):(SBS+NEW-1)) { VRTS[,,L] <- V }

VTI <- V[ ,IS]; VTJ <- V[ ,JS];
if ( NEW == 4 ) { #     Compute four new subregions.
  VRTS[ ,JS,TOP]  <- ( VTI + VTJ )/2; VRTS[ ,IS,SBS+1] <- VTI;
  VRTS[ ,JS,SBS+1] <- ( VTI + VTJ )/2;
  VRTS[ ,IS,SBS+2] <- ( VTI + VTJ )/2; VRTS[ ,JS,SBS+2] <- VTJ;
  VRTS[ ,IS,SBS+3] <- ( VTI + VTJ )/2; VRTS[ ,JS,SBS+3] <- VTJ;
  VTI <- VRTS[ ,IT,TOP]; VTJ <- VRTS[ ,JT,TOP];
  VRTS[ ,JT,TOP]   <- ( VTI + VTJ )/2;
  VRTS[ ,IT,SBS+1] <- ( VTI + VTJ )/2; VRTS[ ,JT,SBS+1] <- VTJ;
  VTI <- VRTS[ ,IT,SBS+2]; VTJ <- VRTS[ ,JT,SBS+2];
  VRTS[ ,JT,SBS+2] <- ( VTI + VTJ )/2;
  VRTS[ ,IT,SBS+3] <- ( VTI + VTJ )/2; VRTS[ ,JT,SBS+3] <- VTJ;
} else { #     Compute three new subregions.
  VRTS[ ,JS,TOP]   <- ( 2*VTI + VTJ )/3; VRTS[ ,IS,SBS+1] <- ( 2*VTI + VTJ )/3;
  if ( DFMX/CUTTF < DFNX ) { VRTS[ ,JS,SBS+1] <- VTJ;
    VRTS[ ,IS,SBS+2] <- ( 2*VTI + VTJ )/3; VRTS[ ,JS,SBS+2] <- VTJ;
    VTJ <- VRTS[ ,JS,SBS+1]; VTL <- VRTS[ ,LS,SBS+1];
    VRTS[ ,LS,SBS+1] <- ( VTJ + VTL )/2;
    VRTS[ ,JS,SBS+2] <- ( VTJ + VTL )/2; VRTS[ ,LS,SBS+2] <- VTL;
  } else { VRTS[ ,JS,SBS+1] <- ( VTI + 2*VTJ )/3;
    VRTS[ ,IS,SBS+2] <- ( VTI + 2*VTJ )/3; VRTS[ ,JS,SBS+2] <- VTJ;
  }
}
return( list( VRTS=VRTS, NEW=NEW ) ) }
#
#***end SMPDFS
##############################################################################
SMPRUL <- function( ND, VRTS, VOL, NF, F, G, W, PTS ) {
#function [ BASVAL RGNERR ] = SMPRUL( ND, VRTS, VOL, NF, F, G, W, PTS )
#
#
#***BEGIN PROLOGUE SMPRUL
#***KEYWORDS basic numerical integration rule
#***PURPOSE  To compute basic integration rule values.
#***AUTHOR
#
#            Alan Genz
#            Department of Mathematics
#            Washington State University
#            Pullman, WA 99164-3113, USA
#            AlanGenz@wsu.edu
#
#***LAST MODIFICATION 2012-06
#***DESCRIPTION SMPRUL computes basic integration rule values for a
#            vector of integrands over a hyper-rectangular region.
#            These are estimates for the integrals. SMPRUL also computes
#            estimates for the errors.
#
#   ON ENTRY
#
#     ND    Integer, number of variables.
#     VRTS  Real array of size (ND,ND+1) of simplex vertices;
#           vertex J must have components VRTS(I,J), I = 1 : 2, ..., ND.
#     NF    Integer, number of components for the vector integrand.
#     F     Function for computing components of the integrand.
#     VOL   Volume of simplex with VRTS vertices.
#     G     matrix of generators (from SMPRMS)
#     W     matrix of weights (from SMPRMS)
#     PTS   integer vector, PTS[i]= number of F values needed for generator i
#
#   ON EXIT
#
#     BASVAL Real array of length NF, values for the basic rule for
#            each component of the integrand.
#     RGNERR Real array of length NF, error estimates for BASVAL.
#
#***end PROLOGUE SMPRUL
#

RTMN <- 1e-1;  SMALL <- 1e-12 ; ERRCOF <- 8;
#
#     Compute the rule values.
#
#[WTS, RLS] <- size(W); RULE <- zeros(NF,RLS);
WTS <- nrow(W); RLS <- ncol(W); RULE <- matrix(0.0,nrow=NF,ncol=RLS)
for (K in 1 : WTS) {
  if (PTS[K] > 0) { RULE <- RULE + VOL*SMPSMS( ND,VRTS,NF,F, G[ ,K] ) %*% W[K, ]; }
}

#
#     Scale integral values and compute the error estimates.
#
BASVAL <- rep(0.0,NF); RGNERR <- rep(0.0,NF);
for ( I in 1 : NF) { 
  BASVAL[I] <- RULE[I,1]; NMBS <-abs( BASVAL[I] ); RT <-RTMN;
#  for (K in RLS : -2 : 3) { 
  for (K in rev(seq(3,RLS,2))) {   
    NMRL <- max( abs( RULE[I,K] ), abs( RULE[I,K-1] ) );
    if ( (NMRL > SMALL*NMBS) & (K < RLS ) )  { RT <-max( c(NMRL/NMCP, RT) ) }
    RGNERR[I] <- max( c(NMRL, RGNERR[I]) ); NMCP <-NMRL;
  } 
  if ( (RT < 1) & (RLS > 3) ) { RGNERR[I] <- RT*NMCP }
  RGNERR[I] <- max( c(ERRCOF*RGNERR[I], SMALL*NMBS) );
}
#
return( list(BASVAL=BASVAL, RGNERR=RGNERR) ) }
#***end SMPRUL
#
###############################################################################
# function SYMSMS = SMPSMS( N, VERTEX, NF, F, G )
SMPSMS <- function( N, VERTEX, NF, F, G ){
#
#***BEGIN PROLOGUE SMPSMS
#***KEYWORDS fully symmetric sum
#***PURPOSE  To compute fully symmetric basic rule sums
#***AUTHOR
#
#        Alan Genz
#        Department of Mathematics
#        Washington State University
#        Pullman, WA 99164-3113, USA
#
#***LAST MODIFICATION 2012-06
#***DESCRIPTION SMPSMS computes a fully symmetric sum for a vector
#            of integrand values over a simplex. The sum is taken over
#            all permutations of the generators for the sum.
#
#   ON ENTRY
#
#   N       Integer, number of variables.
#   VERTEX  Real array of dimension (1:N,1:N+1)
#           The vertices of the simplex, one vertex per column.
#   NF      Integer, number of components for the vector integrand.
#   F       Function for computing components of the integrand at X.
#   G       Real Array of dimension (1:N+1,1).
#           The generators for the fully symmetric sum.
#
#   ON RETURN
#
#   SYMSMS  Real array of length NF, the values for the fully symmetric
#            sums for each component of the integrand.
#
#***ROUTINES CALLED: Integrand
#
#***end PROLOGUE SMPSMS
#
SYMSMS <- rep(0.0,NF); G <- sort(G,decreasing=TRUE); pr <- 1;
#
#     Compute integrand value for permutations of G
#
while (pr) { 
  SYMSMS <- SYMSMS + F(VERTEX %*% G)
  pr <- 0
#
#     Find next distinct permuation of G and loop back for value.
#     Permutations are generated in reverse lexicographic order.
#
  for (I in 2 : (N+1) ) { 
    GI <- G[I];
    if (G[I-1] > GI) { 
      IX <- I - 1;
      for (L in 1 : (IX/2) ) { 
        GL <- G[L]; if (GL <= GI) { IX <- IX - 1; }
        G[L] <- G[I-L]; G[I-L] <- GL; if (G[L] > GI) { LX <- L; }
      }
      if (G[IX] <= GI) { IX <- LX; }
      G[I] <- G[IX]; G[IX] <- GI; 
      pr <- 1; break
    }
  }
}
return( matrix(SYMSMS,nrow=NF,ncol=1) ) }
#
#***end SMPSMS
###############################################################################
# function [G, W, PTS] =  SMPRMS( N, KEY )
SMPRMS <- function( N, KEY ){
#
#***BEGIN PROLOGUE SMPRMS
#***KEYWORDS basic integration rule, degree 2*KEY+1
#***PURPOSE  To initialize a degree 2*KEY+1 basic rule and null rules.
#***AUTHOR
#
#            Alan Genz
#            Department of Mathematics
#            Washington State University
#            Pullman, WA 99164-3113, USA
#            AlanGenz@wsu.edu
#
#***LAST MODIFICATION 2012-06
#***DESCRIPTION  SMPRMS initializes a degree 2*KEY+1 rule, and
#                and max(2*KEY,2) lower degree null rules.
#
#   ON ENTRY
#
#   N    Integer, number of variables.
#   KEY    Integer, < 5 and >= 0, rule parameter.
#          If KEY > 0 a degree 2*KEY+1 rule is initialized.
#
#   ON Exit
#
#   W      Real array of weights for the basic and null rules.
#          W(1,1),...,W(WTS,1) are weights for the basic rule.
#          W(I,1),...,W(WTS,I) for I > 1 are null rule weights.
#   G      Real array of fully symmetric sum generators for the rules.
#          G(1,J), ..., G(N+1,J) are the generators for the
#          points associated with the Jth weights.
#   PTS    Integer array numbers of integrand
#          values needed for generator J.
#
#***REFERENCES
#
#  Axel Grundmann and H. M. Moller
#  "Invariant Integration Formulas for the n-Simplex by Combinatorial
#    Methods", SIAM J Numer. Anal. 15(1978), 282--290,
# and
#  A. H. Stroud
#  "A Fifth Degree Integration Formula for the n-Simplex
#  SIAM J Numer. Anal. 6(1969), 90--98,
# and
#  I. P. Mysovskikh
#  "On a cubature formula for the simplex"
#  Vopros. Vycisl. i Prikl. Mat., Tashkent 51(1978), 74--90.
#
#
#***ROUTINES CALLED NONE
#***end PROLOGUE SMPRMS

# check input values
N <- as.integer(N); KEY <- as.integer(KEY)
stopifnot( N > 1, KEY >= 0, KEY < 5 ) 
#
#
#     Initialize RLS and GMS.
#
switch(KEY, 
  { RLS  <- 3; GMS <-  2; WTS <-  3; }, # case 1
  { RLS  <- 5; GMS <-  4; WTS <-  6; }, # case 2
  { RLS  <- 7; GMS <-  7; WTS <- 11; }, # case 3
  { RLS  <- 7; GMS <- 12; WTS <- 21; if (N == 2) { GMS <- 11; WTS <- 20 } } ) # case 4

W <- matrix( 0.0, nrow=WTS, ncol=RLS); PTS <- rep( 0.0, WTS); G <- matrix( 0.0, nrow=N+1, ncol=WTS);
#
#     Compute generator, PTS and weight values for all rules.
#
NP <- N + 1; N2 <- NP*( N + 2 ); N4 <- N2*( N + 3 )*( N + 4 );
N6 <- N4*( N + 5 )*( N + 6 ); N8 <- N6*( N + 7 )*( N + 8 );
G[,1] <- 1/NP; PTS[1] <- 1;
R1 <- ( N + 4 - sqrt(15) )/( N*N + 8*N + 1 ); S1 <- 1 - N*R1; L1 <- S1 - R1;
G[1,GMS+1] <- S1; G[2:NP,GMS+1] <- R1; PTS[GMS+1] <- N + 1; IW <- RLS;
if ( KEY < 4 ){
  #
  #        Compute weights for special degree 1 rule.
  #
  W[1,IW] <- 1; IW <- IW - 1; W[GMS+1,IW] <- 1/NP; IW <- IW - 1;
}
#
#     Compute weights, generators and PTS for degree 3 rule.
#
G[1,2] <- 3/( N + 3 ); G[2:NP,2] <- 1/( N + 3 ); PTS[2] <- NP;
W[2,IW] <- ( N + 3 )^3/( 4*N2*( N + 3 ) );
if ( KEY > 1 ) { 
  IW <- IW - 1;
  #
  #        Compute weights, generators and PTS for degree 3 and 5 rules.
  #
  if ( N == 2 ) {
    #
    #           Special degree 3 rule.
    #
    L2 <- .62054648267200632589046034361711; L1 <- -sqrt( 1/2 - L2^2 );
    R1 <- ( 1 - L1 )/3; S1 <- 1 - 2*R1;
    G[1,GMS+1] <- S1; G[2:NP,GMS+1] <- R1; PTS[GMS+1] <- 3; W[GMS+1,IW] <- 1/6;
    R2 <- ( 1 - L2 )/3; S2 <- 1 - 2*R2;
    G[1,GMS+2] <- S2; G[2:NP,GMS+2] <- R2; PTS[GMS+2] <- 3; W[GMS+2,IW] <- 1/6;
  } else {
    #
    #           Degree 3 rule using Stroud points.
    #
    R2 <- ( N + 4 + sqrt(15) )/( N*N + 8*N + 1 );
    S2 <- 1 - N*R2; L2 <- S2 - R2;
    G[1,GMS+2] <- S2; G[2:NP,GMS+2] <- R2; PTS[GMS+2] <- NP;
    W[GMS+2,IW] <- ( 2/(N+3) - L1 )/( N2*(L2-L1)*L2^2 );
    W[GMS+1,IW] <- ( 2/(N+3) - L2 )/( N2*(L1-L2)*L1^2 );
  }
  IW <- IW - 1;
  #
  #        Grundmann-Moller degree 5 rule.
  #
  G[1,3] <- 5/( N + 5 ); G[2:NP,3] <- 1/( N + 5 ); PTS[3] <- NP;
  G[1,4] <- 3/( N + 5 ); G[2,4] <- 3/( N + 5 ); G[3:NP,4] <- 1/( N + 5 );
  PTS[4] <- NP*N/2; W[2,IW] <- -( N + 3 )^5/( 16*N4 );
  W[3:4,IW] <-  ( N + 5 )^5/( 16*N4*( N + 5 ) );
}

if ( KEY > 2 ) { 
  IW <- IW - 1;
  #
  #        Compute weights, generators and PTS for degree 5 and 7 rules.
  #
  #
  #        Stroud degree 5 rule.
  #
  U1 <- ( N + 7 + 2*sqrt(15) )/( N*N + 14*N - 11 );
  V1 <- ( 1 - ( N - 1 )*U1 )/2; D1 <- V1 - U1;
  G[1,GMS+3] <- V1; G[2,GMS+3] <- V1; G[3:NP,GMS+3] <- U1;
  PTS[GMS+3] <- ( ( N + 1 )*N )/2;
  U2 <- ( N + 7 - 2*sqrt(15) )/( N*N + 14*N - 11 );
  V2 <- ( 1 - ( N - 1 )*U2 )/2; D2 <- V2 - U2;
  G[1,GMS+4] <- V2; G[2,GMS+4] <- V2; G[3:NP,GMS+4] <- U2;
  PTS[GMS+4] <- ( ( N + 1 )*N )/2;
  if ( N == 2 ) {
    W[GMS+3,IW] <- ( 155 - sqrt(15) )/1200;
    W[GMS+4,IW] <- ( 155 + sqrt(15) )/1200;
    W[1,IW] <- 1 - 3*( W[GMS+3,IW] + W[GMS+4,IW] ) ;
  } else if ( N == 3 ) {
    W[GMS+1,IW] <- ( 2665 + 14*sqrt(15) )/37800;
    W[GMS+2,IW] <- ( 2665 - 14*sqrt(15) )/37800;
    W[GMS+3,IW] <- 2*15/567; PTS[GMS+4] <- 0;
  } else {
    W[GMS+1,IW] <- ( 2*(27-N)/(N+5)-L2*(13-N) )/( L1^4*(L1-L2)*N4 );
    W[GMS+2,IW] <- ( 2*(27-N)/(N+5)-L1*(13-N) )/( L2^4*(L2-L1)*N4 );
    W[GMS+3,IW]=( 2/( N + 5 ) - D2 )/( N4*( D1 - D2 )*D1^4 );
    W[GMS+4,IW]=( 2/( N + 5 ) - D1 )/( N4*( D2 - D1 )*D2^4 );
  } 
  IW <- IW - 1;
  #
  #        Grundmann-Moller degree 7 rule.
  #
  G[1,5] <- 7/( N + 7 ); G[2:NP,5] <- 1/( N + 7 ); PTS[5] <- NP;
  G[1,6] <- 5/( N + 7 ); G[2,6] <- 3/( N + 7 );
  G[3:NP,6] <- 1/( N + 7 ); PTS[6] <- NP*N;
  G[1:3,7] <- 3/( N + 7 ); if(NP>3) {G[4:NP,7] <- 1/( N + 7 ) }; PTS[7] <- NP*N*( N - 1 )/6;
  W[2  ,IW] <-  ( N + 3 )^7/( 2*64*N4*( N + 5 ) );
  W[3:4,IW] <- -( N + 5 )^7/(   64*N6 );
  W[5:7,IW] <-  ( N + 7 )^7/(   64*N6*( N + 7 ) );
}

if ( KEY == 4 ) { 
  IW <- IW - 1;
  #
  #        Compute weights, generators and PTS for degree 7, 9 rules.
  #
  #        Mysovskikh degree 7 rule.
  #  
  SG <- 1/( 23328*N6 );
  U5 <- -6^3*SG*( 52212 - N*( 6353 + N*( 1934 - N*27 ) ) );
  U6 <-  6^4*SG*(  7884 - N*( 1541 - N*9 ) );
  U7 <- -6^5*SG*(  8292 - N*( 1139 - N*3 ) )/( N + 7 );
  P0 <- -144*( 142528 + N*( 23073 - N*115 ) );
  P1 <- -12*( 6690556 + N*( 2641189 + N*( 245378 - N*1495 ) ) );
  P2 <- -16*(6503401 + N*(4020794+N*(787281+N*(47323-N*385))));
  P3 <- -( 6386660 + N*(4411997+N*(951821+N*(61659-N*665))) )*( N + 7 );
  A <- P2/( 3*P3 ); P <- A*( P1/P2 - A ); Q <- A*( 2*A*A - P1/P3 ) + P0/P3;
  R <- sqrt(-P^3); TH <- acos(-Q/(2*R))/3; R <- 2*R^(1/3); TP <- 2*pi/3;
  A1 <- -A + R*cos(TH); A2 <- -A + R*cos(TH+2*TP); A3 <- -A + R*cos(TH+TP);
  G[1,GMS+5] <- ( 1 - N*A1 )/NP; G[2:NP,GMS+5] <- ( 1 + A1 )/NP; PTS[GMS+5] <- NP;
  G[1,GMS+6] <- ( 1 - N*A2 )/NP; G[2:NP,GMS+6] <- ( 1 + A2 )/NP; PTS[GMS+6] <- NP;
  G[1,GMS+7] <- ( 1 - N*A3 )/NP; G[2:NP,GMS+7] <- ( 1 + A3 )/NP; PTS[GMS+7] <- NP;
  W[GMS+5,IW] <- ( U7-(A2+A3)*U6+A2*A3*U5 )/( A1^2-(A2+A3)*A1+A2*A3 )/A1^5;
  W[GMS+6,IW] <- ( U7-(A1+A3)*U6+A1*A3*U5 )/( A2^2-(A1+A3)*A2+A1*A3 )/A2^5;
  W[GMS+7,IW] <- ( U7-(A2+A1)*U6+A2*A1*U5 )/( A3^2-(A2+A1)*A3+A2*A1 )/A3^5;
  G[1,GMS+8] <- 4/( N + 7 ); G[2,GMS+8] <- 4/( N + 7 );
  G[3:NP,GMS+8] <- 1/( N + 7 ); PTS[GMS+8] <- NP*N/2;
  W[GMS+8,IW] <- 10*(N+7)^6/( 729*N6 );
  G[1,GMS+9] <- 11/( N + 7 )/2; G[2,GMS+9] <-  5/( N + 7 )/2;
  G[3:NP,GMS+9] <- 1/( N + 7 ); PTS[GMS+9] <- NP*N;
  W[GMS+9,IW] <- 64*(N+7)^6/( 6561*N6 );
  W[4,IW] <- W[4,IW+1]; W[7,IW] <- W[7,IW+1]; IW <- IW - 1;
  #
  #        Grundmann-Moller degree 9 rule.
  #
  G[1,8] <- 9/( N + 9 ); G[2:NP, 8] <- 1/( N + 9 ); PTS[8] <- NP;
  G[1,9] <- 7/( N + 9 ); G[2,9] <- 3/( N + 9 );
  G[3:NP, 9] <- 1/( N + 9 ); PTS[9] <- NP*N;
  G[1:2,10] <- 5/( N + 9 );
  G[3:NP,10] <- 1/( N + 9 ); PTS[10] <- NP*N/2;
  G[1,11] <- 5/( N + 9 ); G[2:3,11] <- 3/( N + 9 );
  if(NP>3) {G[4:NP,11] <- 1/( N + 9 ) }; PTS[11] <- NP*N*( N - 1 )/2;
  W[2   ,IW] <- -( N + 3 )^9/( 6*256*N6 );
  W[3:4 ,IW] <-  ( N + 5 )^9/( 2*256*N6*( N + 7 ) );
  W[5:7 ,IW] <- -( N + 7 )^9/(   256*N8 );
  W[8:11,IW] <-  ( N + 9 )^9/(   256*N8*( N + 9 ) );
  if ( N > 2 ) {
    G[1:4,12] <- 3/( N + 9 ); if (NP>4) { G[5:NP,12] <- 1/( N + 9 ) };
    PTS[12] <- NP*N*( N - 1 )*( N - 2 )/24; W[12,IW] <- W[8,IW];
  }
}
#
#     Compute unnormalized weights.
#
W[1, ] <- 1 - PTS[2:WTS] %*% W[2:WTS, ]; NB <- sum(PTS * (W[ ,1]^2));
W[ ,2:RLS] <- W[ ,2:RLS] - W[ ,1,drop=FALSE] %*% matrix(1.0,nrow=1,ncol=RLS-1); 
#
#        Orthogonalize and normalize null rules.
#
W[ ,2] = W[ ,2]*sqrt( NB/sum(PTS*W[ ,2]^2) )
for (K in 3 : RLS ) {
  W[ ,K] = W[ ,K] - W[ ,2:(K-1),drop=FALSE] %*% t(W[ ,2:(K-1),drop=FALSE]) %*% ( PTS*W[ ,K] )/NB
  W[ ,K] = W[ ,K]*sqrt( NB/sum(PTS*W[ ,K]^2) )
}

return( list( G=G, W=W, PTS=PTS ) ) }
#
#***End Smprms
#

