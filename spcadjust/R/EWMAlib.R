################################################################
## Computing properties of EWMA charts with known parameters ##
################################################################

# Set up the transient state transition matrix with centerline at 0.

getR <- function(c,gridpoints,pobs,lambda){
  gridpointsh=floor(gridpoints/2)
  p <- c(-(gridpointsh:1)*2*c/gridpoints,0,(1:gridpointsh)*2*c/gridpoints) 
  ptarget <- c(-c,p+c/gridpoints)
  sapply(p,function(x) {res <- pobs((ptarget-(1-lambda)*x)/lambda); diff(res)})
}

ARL_EWMA_Markovapprox <- function(c,gridpoints=100,pobs,lambda){
  gridpoints <- 2*floor(gridpoints/2)+1 #Odd number
  gridpointsh <- floor(gridpoints/2)
  R <- getR(c,gridpoints,pobs,lambda)
  tryCatch(rep(1,gridpoints)%*%solve(diag(rep(1,gridpoints))-R,c(rep(0,gridpointsh),1,rep(0,gridpointsh))),error=function(e) Inf)
}

hitprob_EWMA_Markovapprox <- function(pobs,c,n,gridpoints=100,lambda){
  gridpoints <- 2*floor(gridpoints/2)+1 #Odd number
  gridpointsh <- floor(gridpoints/2)
  R <- getR(c,gridpoints,pobs,lambda)
  R <- matrix.power(R,n)
  1-rep(1,gridpoints)%*%R%*%c(rep(0,gridpointsh),1,rep(0,gridpointsh))
}




