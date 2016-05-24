## ## commented out for time (and dependency reasons)
require("DoseFinding")
if(!(require("quadprog") & require("Rsolnp")))
  stop("need packages quadprog and Rsolnp to run these tests")

## ## calculation of optimal contrast by enumerating all active sets
## allActiveSets <- function(S, mu, mult){
##   k <- length(mu)
##   CC <- cbind(-1, diag(k - 1))
##   SPa <- CC %*% S %*% t(CC)
##   muPa <- as.numeric(CC %*% mu)
##   ## generate all possible active sets
##   mat <- matrix(nrow = 2^(k-1), ncol = (k-1))
##   for(i in 1:(k-1))
##     mat[,i] <- rep(rep(c(FALSE,TRUE), each=2^(i-1)), 2^((k-1)-i))
##   val <- numeric(2^(k-1))
##   feasible <- logical(2^(k-1))
##   cont <- matrix(nrow = 2^(k-1), ncol = (k-1))
##   for(i in 1:(2^(k-1))){
##     nonzero <- mat[i,]
##     if(sum(nonzero) > 0){
##       cont[i,!nonzero] <- 0
##       cont[i,nonzero] <- solve(SPa[nonzero, nonzero]) %*% muPa[nonzero]
##       feasible[i] <- all(mult*cont[i,] >= 0)
##       contrast <- c(-sum(cont[i,]), cont[i,])
##       val[i] <- as.numeric(t(contrast)%*%mu/sqrt(t(contrast)%*%S%*%contrast))
##     }
##   }
##   if(!any(feasible))
##     return(rep(NA, k))
##   mm <- max(val[which(feasible)])
##   c(-sum(cont[val == mm,]), cont[val == mm,])  
## }


## ## helper functions
## getStand <- function(x)
##   x/sqrt(sum(x^2))
## getNCP <- function(cont, mu, S)
##   as.numeric(t(cont)%*%mu/sqrt(t(cont)%*%S%*%cont))

## set.seed(1)
## ncp1 <- ncp2 <- ncp3 <- ncp4 <- ncp5 <- numeric(1000)
## for(i in 1:1000){
##   ## simulate mean and covariance matrix
##   kk <- round(runif(1, 4, 10))
##   A <- matrix(runif(kk^2,-1,1),kk,kk)
##   S <- crossprod(A)+diag(kk)
##   mult <- sign(rnorm(1))
##   mu <- mult*sort(rnorm(kk, 1:kk, 1))
  
##   ## unconstrained solution
##   ones <- rep(1,kk)
##   unConst <- solve(S)%*%(mu - t(mu)%*%solve(S)%*%ones/(t(ones)%*%solve(S)%*%ones))
##   cont1 <- getStand(unConst)
  
##   ## function from DoseFinding package
##   cont2 <- DoseFinding:::constOptC(mu, solve(S), placAdj=FALSE,
##                                    ifelse(mult == 1, "increasing", "decreasing"))
  
##   ## alternative solution using quadratic programming
##   D <- S
##   d <- rep(0,kk)
##   tA <- rbind(rep(1, kk), 
##               mu, 
##               mult*diag(kk)*c(-1,rep(1,kk-1)))
##   A <- t(tA)
##   bvec <- c(0,1,rep(0,kk))
##   rr <- solve.QP(D, d, A, bvec, meq=2)
##   cont3 <- getStand(rr$solution)
  
##   ## using solnp
##   LB <- rep(0, kk-1)
##   UB <- rep(20, kk-1)
##   strt <- rep(1, kk-1)
##   mgetNCP <- function(x, ...){
##     cont <- c(-sum(x), x)
##     -getNCP(cont, ...)
##   }
##   res <- solnp(strt, mgetNCP, mu=mu, S=S,
##                LB=LB, UB=UB,
##                control = list(trace = 0))
##   out <- c(-sum(res$pars), res$pars)
##   cont4 <- getStand(out)

##   ## using
##   cont5 <- allActiveSets(S=S, mu=mu, mult=mult)
    
##   ## compare optimized non-centrality parameters
##   ncp1[i] <- getNCP(cont1, mu, S)
##   ncp2[i] <- getNCP(cont2, mu, S)
##   ncp3[i] <- getNCP(cont3, mu, S)
##   ncp4[i] <- getNCP(cont4, mu, S)
##   ncp5[i] <- getNCP(cont5, mu, S)  
## }
## sapply(list(ncp1, ncp2, ncp3, ncp4, ncp5), quantile)

## ## tests whether constant shapes (possible with linInt) are handled correctly
## data(biom)
## ## define shapes for which to calculate optimal contrasts
## modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
##                 linInt = rbind(c(0, 0, 0, 1), c(0, 1, 1, 1)),
##                 doses = c(0, 0.05, 0.2, 0.6, 1), placEff = 1)

## optContr(modlist, w=1, doses=c(0.05), placAdj=TRUE, type = "u")
## optContr(modlist, w=1, doses=c(0.05), placAdj=TRUE, type = "c")
## optContr(modlist, w=1, doses=c(0.05,0.5), placAdj=TRUE, type = "u")
## optContr(modlist, w=1, doses=c(0.05,0.5), placAdj=TRUE, type = "c")
## optContr(modlist, w=1, doses=c(0,0.05), placAdj=FALSE, type = "u")
## optContr(modlist, w=1, doses=c(0,0.05), placAdj=FALSE, type = "c")
## optContr(modlist, w=1, doses=c(0,0.05,0.5), placAdj=FALSE, type = "u")
## optContr(modlist, w=1, doses=c(0,0.05,0.5), placAdj=FALSE, type = "c")


## modlist2 <- Mods(linInt = rbind(c(0, 1, 1, 1), c(0,0,0,1)),
##                  doses = c(0, 0.05, 0.2, 0.6, 1), placEff = 1)

## ## all of these should throw an error
## optContr(modlist2, w=1, doses=c(0.05), placAdj=TRUE, type = "u")
## optContr(modlist2, w=1, doses=c(0.05), placAdj=TRUE, type = "c")
## optContr(modlist2, w=1, doses=c(0,0.05), placAdj=FALSE, type = "u")
## optContr(modlist2, w=1, doses=c(0,0.05), placAdj=FALSE, type = "c")
## ## these should work
## optContr(modlist2, w=1, doses=c(0.05,0.5), placAdj=TRUE, type = "u")
## optContr(modlist2, w=1, doses=c(0.05,0.5), placAdj=TRUE, type = "c")
## optContr(modlist2, w=1, doses=c(0,0.05,0.5), placAdj=FALSE, type = "u")
## optContr(modlist2, w=1, doses=c(0,0.05,0.5), placAdj=FALSE, type = "c")
