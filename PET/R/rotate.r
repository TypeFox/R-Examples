########################################################################
#
# Copyright Weierstrass Institute for Applied Analysis and 
#           Stochastics (WIAS) & Humboldt Universitaet zu Berlin, 
#           Institut fuer Mathematik, Germany 2006
# *********************************************************
#
# Name:          rotate.r
#                ---------------
# Author:        Joern Schulz
# Stand:         15.08.2006
#
#########################################################################

rotate <- function(A, grad=90){

# ========================================================================
#
# Description:
# 
# With this funtion is it possible to rotate a vector, matrix or an array.
# A      Has to be a vector, matrix or an array.
# grad   Can set to 0, 90, 180, 270, -90, -180 or -270.
#
# Examples:
#
# x <- c(1:5)
# rotate(x,270)
#
# A <- matrix(1:12,ncol=4,byrow=TRUE)
# rotate(A)
#
# A <- array(1:20, c(2,5,2))
# rotate(A, 180)
# ========================================================================

  gradLogical = TRUE
  if (grad != 0 && grad != 90 && grad != 180 && grad != 270 && 
      grad != -90 && grad != -180 && grad != -270 ){
    
      cat("WARNING: grad=", grad, " is not supported. The surrender value is 'NULL'.\n", sep="")
      gradLogical=FALSE
      rotA <- NULL
  }

  if (grad<0)
      grad <- 360+grad

  dimA  <- dim(A)
  ldimA <- length(dimA)

  if ( (!(is.null(dimA)) && ldimA != 2 && ldimA != 3) || is.null(A) ){
      cat("WARNING: 'A' has to be a vector, a matrix or an array.\n") 
      cat("The surrender value is 'NULL'.\n")
      gradLogical=FALSE
      rotA <- NULL
  }
    
 #########################################################################
 # Case 1: A is a vector
  if(gradLogical){
  if (is.null(dimA)){
    dimA  <- length(A)
      if(dimA>1){
        if(grad==90){
          cat("NOTE: The type of the surrender data is changed from vector to matrix. \n")
          rotA=matrix(A, nrow=dimA, ncol=1)
          rotA=matrix(rotA[dimA:1,],nrow=dimA)
        } else if (grad==180) {
          rotA=A[dimA:1]
        } else if(grad==270) {
          cat("NOTE: The type of the surrender data is changed from vector to matrix. \n")
          rotA=matrix(A, nrow=dimA, 1)
        } else
          rotA=A
      } else if(dimA==1) rotA=A
      else stop("Unsupported data typ.")

  } else if (ldimA==2){
 #########################################################################
 # Case 2: A is a matrix
    nrowA <- dimA[1]
    ncolA <- dimA[2]
    if(grad==90){
        rotA <- aperm(A,c(2,1))
        rotA[1:ncolA,] <- rotA[ncolA:1,]
    } else if (grad==180) {
        if(nrowA==1 || ncolA==1)
            rotA <- matrix(A[,ncolA:1],nrow=nrowA ,ncol=ncolA)
        else
            rotA <- A[,ncolA:1]
        rotA[1:nrowA,] <- rotA[nrowA:1,]
    } else if(grad==270) {
        rotA <- aperm(A,c(2,1))
        rotA[,1:nrowA] <- rotA[,nrowA:1]
    } else
        rotA <- A
    
  } else if (ldimA==3){
 #########################################################################
 # Case 3: A is an 3-dimensional array
    nrowA <- dimA[1]
    ncolA <- dimA[2]
    permVec <- c(2,1, 3:length(dimA))
    if(grad==90){
        rotA <- aperm(A,permVec)
        rotA <- rotA[ncolA:1,,]
    } else if (grad==180) {
        rotA <- A[,ncolA:1,]
        rotA <- rotA[nrowA:1,,]
    } else if(grad==270) {
        rotA <- aperm(A,permVec)
        rotA <- rotA[,nrowA:1,]
    } else
        rotA=A
  }
  }

  return(rotA)
}