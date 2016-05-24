# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

which.max.matrix <- function(z) {
    if (!is.matrix(z)) {
        stop("Not a matrix")
    }
    m <- nrow(z)
    n <- ncol(z)
    # take care of NAs
    ind <- which.max(z)
    iy <- trunc((ind - 1)/m) + 1
    ix <- ind - (iy - 1) * m
    return(cbind(ix, iy))
}


which.max.image <- function(obj) {
    ind.z <- which.max.matrix(obj$z)
    return(list(x = obj$x[ind.z[, 1]], y = obj$y[ind.z[, 2]], 
        z = obj$z[ind.z], ind = ind.z))
}

LKArrayShift<-  function (A, shift, periodic=FALSE) 
{
    L <- dim(A)
    B<- array( NA, L)
    N<- prod(L)
    if (any(abs(shift) > L) ) {
        stop("shift exceeds array dimensions")
    }
    dimension<- length(L)
    indexListSource<- indexListTarget<- NULL
     for( k in 1:dimension) {
     	  indexListSource<- c( indexListSource, list( c(1:L[k]) ) ) 
        if( periodic){
          tempIndex <-  (0:(L[k]-1) + shift[k]) %% L[k] +1
        }
        else{
     	  tempIndex <-  1:L[k] + shift[k]
        } 
# set indices beyond range (i.e. boundaries of lattice to NAs)
          tempIndex[  (tempIndex<1) | (tempIndex > L[k]) ] <- NA     	   
 	      indexListTarget <- c( indexListTarget, list( tempIndex ) )
     	  }  
# indices are organized as matrices with each column representing a dimension     	  
    indexSource <- as.matrix( expand.grid( indexListSource) )
    indexTarget <- as.matrix( expand.grid( indexListTarget) )
# identify NAs  not sure if this code is too cute 
    inRange <-  rowSums( is.na(indexTarget) ) == 0
    B[ indexTarget[inRange,]] <- A[ indexSource[inRange,]] 
    return( B)
}

# A<- array( 1:(2*4*3), c( 3,4,2)); shift<- c( 1,0,0); LKArrayShift( A, shift=c( 0,1,0))

LKrig.rowshift <- function(A, shift.row, shift.col) {
    mx <- dim(A)[1]
    my <- dim(A)[2]
    if (abs(shift.row) > mx) {
        stop("shift exceeds matrix dimension")
    }
    if (shift.row < 0) {
        i.target <- 1:(mx + shift.row)
        i.source <- (-shift.row + 1):mx
    }
    else {
        i.target <- (shift.row + 1):mx
        i.source <- (1:(mx - shift.row))
    }
    
    B <- matrix(NA, mx, my)
    B[i.target, ] <- A[i.source, ]
    return(B)
}

LKrig.rowshift.periodic <- function(A, shift.row) {
    mx <- dim(A)[1]
    if (abs(shift.row) > mx) {
        stop("shift exceeds matrix dimension")
    }
    if (shift.row < 0) {
        i.source <- c((-shift.row + 1):mx, 1:(-shift.row))
    }
    else {
        i.source <- c(mx - (shift.row:1) + 1, 1:(mx - shift.row))
    }
    return(A[i.source, ])
}

LKrig.shift.matrix <- function(A, shift.row = 0, shift.col = 0, 
    periodic = c(FALSE, FALSE)) {
    if (shift.row != 0) {
        if (!periodic[1]) {
            A <- LKrig.rowshift(A, shift.row = shift.row)
        }
        else {
            A <- LKrig.rowshift.periodic(A, shift.row = shift.row)
        }
    }
    if (shift.col != 0) {
        if (!periodic[2]) {
            A <- t(LKrig.rowshift(t(A), shift.row = shift.col))
        }
        else {
            A <- t(LKrig.rowshift.periodic(t(A), shift.row = shift.col))
        }
    }
    return(A)
}

repMatrix<- function(A, times=1, each=1, byrow=TRUE){
   A<- as.matrix(A)
   if( !byrow){ A<- t(A) }
   nSize<- nrow(A)*each*times
   nColumn<- ncol(A)
   B<- matrix(NA, nSize, ncol=nColumn)
   for( k in 1:nColumn){
      B[,k] <- rep( A[,k], times=times, each=each)
    }
  if( !byrow){ B<- t(B)}
   B
 }

expandMatrix0<- function( A, B){
   A<- as.matrix(A)
   B<- as.matrix(B)
   m1<- nrow( A)
   m2<- nrow( B)
   cbind(repMatrix( A,times=m2, each=1), repMatrix( B, times=1, each=m1))
 }

expandMatrix<- function( ...){
  matrices<- list(...)
  N<- length( matrices)
  tempM<- matrices[[1]]
  for( k in 2:N){
    tempM<- expandMatrix0(tempM, matrices[[k]])
  } 
  tempM
}
    
 expandMList <- function( Mlist, byrow=TRUE){
   N<- length( Mlist)
   
   for( k in 2:N){
     mtimes<- nrow( as.matrix(Mlist[[k]]))
     each<- nrow( as.matrix( Mlist[[1]]))
     for( l in 1:(k-1)){
      Mlist[[l]] <- repMatrix(Mlist[[l]], byrow=byrow,times=mtimes)
    }
     Mlist[[k]]<- repMatrix( Mlist[[k]], byrow=byrow, each =each)
   }
     Mlist
   }
      
grid2Index<- function( I, grid){
	I<- rbind( I)
	gridP<- cumprod( c( 1, c(grid)))
	nDim<- ncol( I)
	J<- rep( 0, nrow( I) )
	for( k in 1:nDim){
		J<- J + (I[,k]-1)*gridP[k] 
		}
	return( J+1)
}
   
