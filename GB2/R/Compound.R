# Gamma factors, which multiplied my the GB2 density give the component densities
fg.cgb2 <- function(x, shape1, scale, shape2, shape3, pl0, decomp="r"){
      # pl0 is a vector of probabilities which sums to 1
      # for ex. pl0 = c(1/3,1/3,1/3)
	  if (decomp=="r") 
         {sh <- shape3
          a0 <- shape1}
      if (decomp=="l")
         {sh <- shape2
          a0 <- -shape1}
      t <- (x/scale)^a0+1
      ncomp <- length(pl0)
      fac <- matrix(rep(0, times=length(x)*ncomp), ncol=ncomp)
      u <- qgamma(cumsum(pl0),sh)
      u1 <-c(0,u[-ncomp])     # dim = ncomp
      pq <- shape2+shape3
      fac[,1] <-pgamma(t*u[1],pq)/pgamma(u[1],sh)
      for(i in 2:ncomp){
      fac[,i] <-(pgamma(t*u[i],pq) - pgamma(t*u1[i],pq))/(pgamma(u[i],sh) - pgamma(u1[i],sh))
      }
	  dimnames(fac)[[2]] <- paste("fac",1:dim(fac)[2],sep="")
	  return(fac) 
}

dl.cgb2 <- function(x, shape1, scale, shape2, shape3, pl0, decomp="r"){
      # pl0 is a vector of probabilities which sums to 1
      # for ex. pl0 = c(1/3,1/3,1/3)
          L <- length(pl0)
	  fac <- fg.cgb2(x, shape1, scale, shape2, shape3, pl0, decomp)
	  dcl <- fac                                             
	  for (i in 1:L){
	  dcl[,i] <- dgb2(x,shape1,scale,shape2,shape3)*fac[,i]  
	  }
	  colnames(dcl) <- paste("comp",1:L,sep="")
	  return(dcl) #returns a matrix with the l-th component density in the l-th column
}

# Distribution function
library(cubature)

pl.cgb2 <- function(y, shape1, scale, shape2, shape3, pl0, decomp="r", tol=1e-05){
      ncomp <- length(pl0)
      Fl <- matrix(rep(1, times=length(y)*ncomp), ncol=ncomp)
      dimnames(Fl)[[2]] <- paste("comp",1:dim(Fl)[2],sep="")
      v <- (y<Inf)
      if (length(y[v]) !=0){
       for (i in 1:ncomp){
        integrand <- function(x) dl.cgb2(x,shape1,scale,shape2,shape3,pl0,decomp)[i]
        fun <- function(x) adaptIntegrate(integrand,0,x,tol=tol)$integral
        Fl[v,i] <- apply(as.matrix(y[v]),1,fun)
       }
      }
     vv <- (y==Inf)
     Fl[vv,] <- 1
     return(Fl)
}

# Density
dcgb2 <- function(x, shape1, scale, shape2, shape3, pl0, pl, decomp="r"){
# x: value at which the density is evaluated, could be a vector or a scalar
# shape1,scale,shape2,shape3: GB2 parameters
# vector pl0 (sums to 1) of length L, the initial proportions
# vector pl gives the fitted proportions of the L components
# if pl=pl0 the function returns the GB2 density
        ncomp <- length(pl0)
        v <- (x<Inf)        
        fl <- matrix(rep(0, times=length(x)*ncomp), ncol=ncomp)
	    fl[v,] <- dl.cgb2(x[v],shape1,scale,shape2,shape3,pl0,decomp)
	    vv <- (x==Inf)          #new 7.06.12
        fl[vv,] <- 0           #new 7.06.12
        compd <- fl%*%pl
        dimnames(compd)[[2]] <- list("compound")
        return(compd)
}


# Cumulative distribution function
pcgb2 <- function(y, shape1, scale, shape2, shape3, pl0, pl, decomp="r"){
n <- length(y)
compp <- matrix(rep(1,n),ncol=1)
v <- (y<Inf)
if (length(y[v])!=0){
integrand <- function(x) dcgb2(x,shape1,scale,shape2,shape3,pl0,pl,decomp)
fun <- function(x) adaptIntegrate(integrand,0,x,tol=1e-08)$integral
compp[v,1] <- apply(as.matrix(y[v]),1,fun)
}
vv <- (y==Inf)
compp[vv,] <- 1
return(compp)
}

prcgb2 <- function(y1, y2, shape1, scale, shape2, shape3, pl0, pl, decomp="r", tol=1e-08, debug=FALSE){
# Given  2 arguments y1 and y2, returns the probability that the GB2 r.v. Y is P(min(y1,y2) < Y < max(y1,y2))

        integrand <- function(x) dcgb2(x,shape1,scale,shape2,shape3,pl0,pl,decomp)
        M1 <- max(y1,y2)
        m1 <- min(y1,y2)
        if (m1 == M1) return(0)
        if (M1 < Inf)  {if (!debug) {
                                return(adaptIntegrate(integrand,m1,M1,tol=tol)$integral)
                                }
                        else {
                                return(adaptIntegrate(integrand,m1,M1,tol=tol))
                                }
                        }
        if (M1 == Inf) {if (!debug) {
                                return(1-adaptIntegrate(integrand,0,m1,tol=tol)$integral)
                                }
                        else {
                                return(1-adaptIntegrate(integrand,0,m1,tol=tol))
                                }
                         }
}






