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
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2

LKrigOLD <- function(x, y = NULL, weights = rep(1, nrow(x)), 
    Z = NULL,
    fixedFunction="LKrigDefaultFixedFunction", fixedFunctionArgs=NULL, m=2,
    LKinfo = NULL,
    basisInfo=NULL, latticeInfo=NULL,
    LKgeometry="LKRectangle",              
    iseed = 123, NtrA = 20,
    use.cholesky = NULL, return.cholesky = TRUE, wPHI=NULL, return.wPHI=TRUE,
    nlevel, a.wght, alpha=NA, nu = NULL,  overlap = 2.5, V=NA,                 
    lambda = NA, sigma = NA, rho = NA, rho.object = NULL,
    normalize = TRUE, edge = FALSE, RadialBasisFunction = "WendlandFunction", 
    distance.type = "Euclidean",             
# args particular to a geometry
#e.g. NC= 10 for the LKretangle case
   ...,
#for debugging                  
                  verbose = FALSE) {
 # make sure locations are a matrix and get the rows
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- nrow(x)
    if (any(duplicated(cat.matrix(x)))) 
        stop("locations are not unique see help(LKrig) ")
    # make sure covariate is a matrix
    if (!is.null(Z)) {
        Z <- as.matrix(Z)
    }
 # check for missing values
    if (!is.null(y)) {
        if (any(is.na(y))) 
            stop("Missing values in y not allowed ")
    }
 # hard wire m argument if the default fixed function
 # is used (a low order polynomial of degree m-1 ).
    
    if( !is.null(fixedFunction)){
        if(fixedFunction == "LKrigDefaultFixedFunction" ){
            fixedFunctionArgs<- c( fixedFunctionArgs, list(m=m))
        }
    }    
 # if LKinfo is missing create it from passed arguments   
    if (is.null(LKinfo)) {
        LKinfo <- LKrigSetup(x,
                             basisInfo=basisInfo,latticeInfo=latticeInfo,
                             LKgeometry=LKgeometry,
            nlevel = nlevel, lambda = lambda, 
            sigma = sigma, rho = rho, alpha = alpha, nu = nu, 
            a.wght = a.wght, overlap = overlap, normalize = normalize, 
            RadialBasisFunction = RadialBasisFunction, 
            V = V, distance.type = distance.type, rho.object = rho.object, setupArgs= list(...))
    }
    else{
       if( !is.na(lambda)){
        LKinfo$lambda <- lambda
       }
 # Possibly set lambda from the "noise to signal" variances.
       if (!is.na(rho) & !is.na(sigma)) {
        LKinfo$lambda<- sigma^2/rho
       }
    }
#   At this point the LKinfo object should have the right lambda value. 
    lambda = LKinfo$lambda
    if( is.na(lambda)){
      stop("lambda must be specified")
    }
#    
    if( verbose) {
         print( LKinfo)
    }
 # Begin computations ....
 # weighted observation vector
    wy <- sqrt(weights) * y
 # Spatial drift matrix -- default is assumed to be linear in coordinates. (m=2)
 # and includes possible covariate(s) -- the Z matrix.
    if( !is.null(fixedFunction)){
       wT.matrix <- sqrt(weights) *
       do.call(fixedFunction, c(
                                  list(x=x, Z=Z,
                                  distance.type = LKinfo$distance.type),
                                  fixedFunctionArgs))
       nt <- ncol(wT.matrix)
       nZ <- ifelse(is.null(Z), 0, ncol(Z))
       ind.drift <- c(rep(TRUE, (nt - nZ)), rep(FALSE, nZ))
       if( verbose){
         cat("dim wT", dim( wT.matrix), fill=TRUE)
       }
     }
   else{
       nt<- 0
       ind.drift<- NULL
       nZ<- 0
       wT.matrix<- NULL
   } 
 # Matrix of sum( N1*N2) basis function (columns) evaluated at the N locations (rows)
 # and multiplied by square root of diagonal weight matrix
 # this can be a large matrix if not encoded in sparse format.
    if( is.null(wPHI) ){
       if( verbose){
         cat("Finding new wPHI", fill=TRUE)
       }
       wPHI <- diag.spam(sqrt(weights)) %*% LKrig.basis(x, LKinfo, verbose = verbose)
        if( verbose){ cat("Computed new wPHI", fill=TRUE)
                    print( dim( wPHI) )}
    }
    else{
       if( verbose){ cat("reuse wPHI", fill=TRUE)}}
 #   square root of precision matrix of the lattice process
 #   solve(t(H)%*%H) is proportional to the covariance matrix of the Markov Random Field
    Q <- LKrig.precision(LKinfo)
 # M is the regularized regression matrix that is the key to the entire algorithm:
    M <- t(wPHI) %*% wPHI + lambda * (Q)
    nzero <- length(M@entries)
    if (verbose) {
        cat("Number of nonzero elements:", nzero, fill = TRUE)
    }
 # find Cholesky square root of M
 #  This is where the heavy lifting happens!  M is in sparse format so
 #   by the overloading is a sparse cholesky decomposition.
 #  if this function has been coded efficiently this step should dominate
 #  all other computations.
 #  If  a previous sparse cholesky decoposition is passed then the
 #  pattern of sparseness is used for the decoposition.
 #  This can speed the computation as the symbolic decomposition part of the
 #  sparse Cholesky is a nontrivial step. The condition is that
 #  the current 'M' matrix  has the same sparse pattern as that
 #  which resulted in the factorization  cholesky as 'use.cholesky'
# Timing wrapper - uncomment #--#
#--# cholTime<- system.time(
    if (is.null(use.cholesky)) {
        if( verbose){ cat("new Cholesky", fill=TRUE)}
        Mc <- chol(M)       
    }
    else {
    # check that reuse object is bad
#       if( length(M@colindices) != length(use.cholesky@colindices)){
#           stop('use.cholesky@coloindices not the same length as current model')}
#       if( any( M@colindices != use.cholesky@colindices ) ){
#           stop('use.cholesky not the same sparse pattern as current model')}
       if( verbose){
         cat("reuse symbolic decomposition", fill=TRUE)
       }
       Mc <- update.spam.chol.NgPeyton( use.cholesky, M)        
    }
#--#  )
#--#  cat( "Cholesky decomp time", cholTime,  fill=TRUE)    
# partially fill object list with some components
    object <- list(x = x, y = y, weights = weights, Z = Z, nZ = nZ, 
        ind.drift = ind.drift, LKinfo = LKinfo, lambda = lambda, sigma = sigma, rho = rho)
# use Mc to find coefficients of estimate
    out1 <- LKrig.coef(Mc, wPHI, wT.matrix, wy, lambda, weights)
       if( verbose){
         cat(" d.coef", out1$d.coef, fill=TRUE)
       }
# add in components from coefficient estimates
    object <- c(object, out1)
# compute predicted values    
    fitted.values<-  (wPHI%*%out1$c.coef)/sqrt(weights)    
    if(!is.null(fixedFunction)){
        fitted.values.fixed<- (wT.matrix%*%out1$d.coef)/sqrt(weights)
        fitted.values <-  fitted.values.fixed + fitted.values
    }     
# For reference: fitted.values <- predict.LKrig(object, x, Znew = object$Z)
    residuals <- y - fitted.values
    out2 <- LKrig.lnPlike(Mc, Q, y, lambda, residuals, weights, 
                          sigma, rho)
    object <- c(object, out2)
# estimate trace by Monte Carlo if NtrA greater than zero
    if (NtrA > 0) {
        out3 <- LKrig.traceA(Mc, wPHI, wT.matrix, lambda, weights, 
                             NtrA, iseed)
     # find GCV using this trace estimate
        out3$GCV = (sum(weights * (residuals)^2)/n)/(1 - out3$trA.est/n)^2
    }
    else {
        out3 <- list(trA.est = NA, trA.SE = NA, GCV = NA)
    }
    object <- c(object, out3)  
 # the output object
 # note the ifelse switch whether to return the big cholesky decomposition
 # and/or the weighted basis matrix wPHI
    object <- c(object, list(fitted.values = fitted.values, residuals = residuals, 
        m = LKinfo$m, lambda.fixed = lambda, nonzero.entries = nzero, 
        spatialdriftorder = 2, nt = nt, eff.df = out3$trA.est,
        fixedFunction= fixedFunction, fixedFunctionArgs= fixedFunctionArgs,
        call = match.call()))
    if (return.cholesky) {
        object$Mc <- Mc
    }
    if (return.wPHI) {
        object$wPHI <- wPHI
    }
 # set the class and return.
    class(object) <- "LKrig"
    return(object)
}

