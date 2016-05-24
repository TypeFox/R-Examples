"D2Amat" <-
function(J, filter.number = 10., family = "DaubLeAsymm", OPLENGTH = 10000.,
switch = "direction", verbose = FALSE)
{
#
# Program to compute the inner product matrix of 
#  2-D discrete autocorrelation wavelets
#
####################################################
##
## This function constructs the (-3J)*(-3J) matrix A which is 
## used to correct the raw periodogram obtained by squaring the
## non-decimated wavelet decomposition of an image
##
## It is important to note that there are two alteratives available
## here. Firstly we can generate A "by direction" - as detailed in 
## Eckley, Nason and Treloar (2001). Secondly we can generate A 
##"by level" - details of which will be explained in the help file
## for this function.
####################################################
####################################################
if(verbose == TRUE) cat("Computing ipndacw\n")
#
# 
# Firstly some checks:
#
if(verbose) now <- proc.time()[1.:2.]
if(J >= 0.)
   stop("J must be negative integer")
if(J - round(J) != 0.)
   stop("J must be an integer")
#
#
if(switch != "direction" && switch != "level") #
   stop("Sorry, but switch can only take the value direction or level! Try again.")
if(switch == "direction")
   Aorig <- A2name(J = J, filter.number = filter.number, family = family, switch = switch)
#
# See if matrix already exists
#
if(switch == "level") 
   Aorig <- A2name(J = J, filter.number = filter.number, family = family, switch = switch)
#
#
#If it doesn't, then we calculate the matrix as follows.
#
if(exists(Aorig)) {
   cat("Returning precomputed version \n")
   if(verbose) {
      speed <- proc.time()[1.:2.] - now
      cat("Took ", sum(speed), " seconds \n")
   }
   return(get(Aorig, envir=DWEnv))
}
P <- D2ACW(J, filter.number = filter.number, family = family, OPLENGTH
 = OPLENGTH, switch = switch)
J <-  - J
#The no of rows of any element in P = 
#the no. of it's columns.
#
A <- #
matrix(0., nrow = (3. * J), ncol = (3. * J))
for(i in 1.:(3. * J)) {
   for(j in 1.:(3. * J)) {
      tmp1 <- P[[i]]
      tmp2 <- P[[j]]
      nri <- nrow(tmp1)
      nrj <- nrow(tmp2)
      nrmin <- min(nri, nrj)
      nri2 <- (nri - nrmin)/2.
      nrj2 <- (nrj - nrmin)/2.
      tmp3 <- tmp1[(1. + nri2):(nri - nri2), (1. + nri2):(nri - nri2)]
      tmp4 <- tmp2[(1. + nrj2):(nrj - nrj2), (1. + nrj2):(nrj - nrj2)]
      A[i, j] <- sum(tmp3 * tmp4)
   }
}
nm <- as.character(0.:( - (3. * J) + 1.))
dimnames(A) <- list(nm, nm)
assign(Aorig, A, envir=DWEnv)
A
}

