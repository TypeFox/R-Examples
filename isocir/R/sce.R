"sce"<-function(arg1, arg2, meanrl=1){
# DATE WRITTEN: 4 Mar 2010          LAST REVISED:  05 Jan 2012
# AUTHOR: Sandra Barragan based on the SAS routines written by Miguel A. Fernandez.
# DESCRIPTION: It computes the Sum of Circular Errors of two vectors
# REFERENCE: Mardia, K. and Jupp, P. (2000). Directional Statistics.
# SEE ALSO: CIRE.

if(!is.vector(arg1)){stop("The first argument must be a vector")}
if(all(c(!is.matrix(arg2),!is.vector(arg2)))){stop("The second argument must be a vector or a matrix")}

proveNA1 <- ifelse(complete.cases(arg1), 0, 1)
if(is.vector(arg2)){
 if(sum(proveNA1) >= 1){
  arg1 <- arg1[proveNA1 == 0]
  arg2 <- arg2[proveNA1 == 0]
  }
 proveNA2 <- rep(0,length(arg2))
 proveNA2[is.na(arg2)] <- 1 
 if(sum(proveNA2) >= 1){
  arg2 <- arg2[proveNA2 == 0]
  arg1 <- arg1[proveNA2 == 0]
  } 
 point1 <- arg1
 point2 <- arg2
 if(length(meanrl) > 1){mrls <- meanrl}
 if(length(meanrl) == 1){mrls <- rep(1,length(point1))}
 } # end if arg2 vector

if(is.matrix(arg2)){
 point2<-suppressWarnings(apply(arg2,1,mean.circular))%%(2*pi)
 if(sum(proveNA1) >= 1){
  point1 <- arg1[proveNA1 == 0]
  point2 <- point2[proveNA1 == 0]
  }
 mrls <- mrl(arg2)
 } # end if arg2 matrix

aux <- mrls*(1 - cos(point1 - point2))
SCE <- sum(aux)
return(SCE)
} # end of function