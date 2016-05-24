
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:
#  require(Rdonlp2)
#  require(Rsolnp)
#  require(Rnlminb2)
# FUNCTION:  
#  .powellTestNLP
#  .wright4TestNLP
#  .boxTestNLP 
#  .wright9TestNLP
#  .alkylationTestNLP
#  .entropyTestNLP
#  .rosensuzukiTestNLP
#  .sharperatioTestNLP
#  .rachevratioTestNLP
#  .boxmarkowitzTestNLP 
################################################################################


.powellTestNLP <-
  function()
  {
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    start <- c(-2, 2, 2, -1, -1)
    
    fun <- function(x) { exp(x[1]*x[2]*x[3]*x[4]*x[5]) }
    
    eqFun <- list(
      function(x) x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5],
      function(x) x[2]*x[3]-5*x[4]*x[5],
      function(x) x[1]*x[1]*x[1]+x[2]*x[2]*x[2])
    eqFun.bound <- c(10, 0, -1)
    
    ans.donlp = donlp2NLP(start, fun, 
                          eqFun = eqFun, eqFun.bound = eqFun.bound) 
    ans.solnp = solnpNLP(start, fun, 
                         eqFun = eqFun, eqFun.bound = eqFun.bound)
    ans.nlminb = nlminb2NLP(start, fun, 
                            eqFun = eqFun, eqFun.bound = eqFun.bound)
    
    result.par = round(rbind(
      ans.donlp$solution, ans.solnp$solution, ans.nlminb$solution), 3)
    result.fun = c(
      fun(ans.donlp$solution), fun(ans.solnp$solution), fun(ans.nlminb$solution))
    cbind(result.par, result.fun)   
  }


# ------------------------------------------------------------------------------


.wright4TestNLP <-
  function()
  {
    # A function implemented by Diethelm Wuertz
    
    # nlminb2 does not converge
    
    # FUNCTION:
    
    start <- c(1, 1, 1, 1, 1)
    
    fun <- function(x){
      (x[1]-1)^2+(x[1]-x[2])^2+(x[2]-x[3])^3+(x[3]-x[4])^4+(x[4]-x[5])^4}
    
    eqFun = list(
      function(x) x[1]+x[2]*x[2]+x[3]*x[3]*x[3],
      function(x) x[2]-x[3]*x[3]+x[4],
      function(x) x[1]*x[5] ) 
    eqFun.bound = c(2+3*sqrt(2), -2+2*sqrt(2), 2)
    
    ans.donlp = donlp2NLP(start, fun, 
                          eqFun = eqFun, eqFun.bound = eqFun.bound)
    ans.solnp = solnpNLP(start, fun, 
                         eqFun = eqFun, eqFun.bound = eqFun.bound) 
    ans.nlminb = nlminb2NLP(start, fun, 
                            eqFun = eqFun, eqFun.bound = eqFun.bound)
    
    result.par = round(rbind(
      ans.donlp$solution, ans.solnp$solution, ans.nlminb$solution), 3)
    result.fun = c(
      fun(ans.donlp$solution), fun(ans.solnp$solution), fun(ans.nlminb$solution))
    cbind(result.par, result.fun)   
  }


# ------------------------------------------------------------------------------


.boxTestNLP <-
  function()
  {
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    start <- c(1.1, 1.1, 9)
    
    fun <- function(x) { -x[1]*x[2]*x[3] }
    par.lower = rep(1, 3)
    par.upper = rep(10, 3)
    
    eqFun <- list(
      function(x) 4*x[1]*x[2]+2*x[2]*x[3]+2*x[3]*x[1] )
    eqFun.bound = 100
    
    ans.donlp = donlp2NLP(start, fun, 
                          par.lower = par.lower, par.upper = par.upper, 
                          eqFun = eqFun, eqFun.bound = eqFun.bound)  
    ans.solnp = solnpNLP(start, fun, 
                         par.lower = par.lower, par.upper = par.upper, 
                         eqFun = eqFun, eqFun.bound = eqFun.bound)        
    ans.nlminb = nlminb2NLP(start, fun, 
                            par.lower = par.lower, par.upper = par.upper, 
                            eqFun = eqFun, eqFun.bound = eqFun.bound)  
    
    result.par = round(rbind(
      ans.donlp$solution, ans.solnp$solution, ans.nlminb$solution), 3)
    result.fun = c(
      fun(ans.donlp$solution), fun(ans.solnp$solution), fun(ans.nlminb$solution))
    cbind(result.par, result.fun) 
  }


# ------------------------------------------------------------------------------ 


.wright9TestNLP <-
  function()
  {
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    start <- c(1, 1, 1, 1, 1)
    
    fun <- function(x){
      10*x[1]*x[4]-6*x[3]*x[2]*x[2]+x[2]*(x[1]*x[1]*x[1])+
        9*sin(x[5]-x[3])+x[5]^4*x[4]*x[4]*x[2]*x[2]*x[2] }
    
    ineqFun <- list(
      function(x) x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5],
      function(x) x[1]*x[1]*x[3]-x[4]*x[5],
      function(x) x[2]*x[2]*x[4]+10*x[1]*x[5]) 
    ineqFun.lower = c(-100,  -2,   5)
    ineqFun.upper = c(  20, 100, 100)
    
    ans.donlp = donlp2NLP(start, fun, 
                          ineqFun = ineqFun, 
                          ineqFun.lower = ineqFun.lower, ineqFun.upper = ineqFun.upper)  
    ans.solnp = solnpNLP(start, fun, 
                         ineqFun = ineqFun, 
                         ineqFun.lower = ineqFun.lower, ineqFun.upper = ineqFun.upper)   
    ans.nlminb = nlminb2NLP(start, fun, 
                            ineqFun = ineqFun, 
                            ineqFun.lower = ineqFun.lower, ineqFun.upper = ineqFun.upper)
    
    result.par = round(rbind(
      ans.donlp$solution, ans.solnp$solution, ans.nlminb$solution), 3)
    result.fun = c(
      fun(ans.donlp$solution), fun(ans.solnp$solution), fun(ans.nlminb$solution))
    cbind(result.par, result.fun)    
  }


# ------------------------------------------------------------------------------ 


.alkylationTestNLP <-
  function()
  {
    # A function implemented by Diethelm Wuertz
    
    # solnp does not converge
    
    # FUNCTION:
    
    start <- c(17.45, 12, 110, 30.5, 19.74, 89.2, 92.8, 8, 3.6, 145.2)
    
    fun <- function(x) { -0.63*x[4]*x[7]+50.4*x[1]+3.5*x[2]+x[3]+33.6*x[5] }
    par.lower <- c( 0,  0,   0, 10,  0, 85, 10, 3, 1, 145)
    par.upper <- c(20, 16, 120, 50, 20, 93, 95,12, 4, 162)
    
    eqFun <- list(
      function(x) 98*x[3]-0.1*x[4]*x[6]*x[9]-x[3]*x[6],
      function(x) 1000*x[2]+100*x[5]-100*x[1]*x[8],
      function(x) 122*x[4]-100*x[1]-100*x[5])
    eqFun.bound = c(0, 0, 0)
    
    ineqFun <- list(
      function(x) (1.12*x[1]+0.13167*x[1]*x[8]-0.00667*x[1]*x[8]*x[8])/x[4],
      function(x) (1.098*x[8]-0.038*x[8]*x[8]+0.325*x[6]+57.25)/x[7],
      function(x) (-0.222*x[10]+35.82)/x[9],
      function(x) (3*x[7]-133)/x[10])
    ineqFun.lower = c(  0.99,   0.99,  0.9,   0.99)
    ineqFun.upper = c(100/99, 100/99, 10/9, 100/99)
    
    ans.donlp = donlp2NLP(start, fun, 
                          par.lower = par.lower, par.upper = par.upper,
                          eqFun = eqFun, eqFun.bound = eqFun.bound,
                          ineqFun = ineqFun, 
                          ineqFun.lower = ineqFun.lower, ineqFun.upper = ineqFun.upper)    
    ans.solnp = solnpNLP(start, fun, 
                         par.lower = par.lower, par.upper = par.upper,
                         eqFun = eqFun, eqFun.bound = eqFun.bound,
                         ineqFun = ineqFun, 
                         ineqFun.lower = ineqFun.lower, ineqFun.upper = ineqFun.upper)
    ans.nlminb = nlminb2NLP(start, fun, 
                            par.lower = par.lower, par.upper = par.upper,
                            eqFun = eqFun, eqFun.bound = eqFun.bound,
                            ineqFun = ineqFun, 
                            ineqFun.lower = ineqFun.lower, ineqFun.upper = ineqFun.upper)
    
    result.par = round(rbind(
      ans.donlp$solution, ans.solnp$solution, ans.nlminb$solution), 3)
    result.fun = c(
      fun(ans.donlp$solution), fun(ans.solnp$solution), fun(ans.nlminb$solution))
    cbind(result.par, result.fun)     
  }


# ------------------------------------------------------------------------------ 


.entropyTestNLP <-
  function()
  {
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    set.seed(1953)
    start <- runif(10, 0, 1)
    
    fun <- function(x) {
      m = length(x)
      f = -sum(log(x))
      vnorm = sum((x-1)^2)^(1/2) 
      f - log(vnorm + 0.1) 
    }
    par.lower <- rep(0, 10)
    
    eqFun <- list(
      function(x) sum(x) )
    eqFun.bound <- 10
    
    ans.donlp = donlp2NLP(start, fun, 
                          par.lower = par.lower,  
                          eqFun = eqFun, eqFun.bound = eqFun.bound)     
    ans.solnp = solnpNLP(start, fun, 
                         par.lower = par.lower,  
                         eqFun = eqFun, eqFun.bound = eqFun.bound)    
    ans.nlminb = nlminb2NLP(start, fun, 
                            par.lower = par.lower,  
                            eqFun = eqFun, eqFun.bound = eqFun.bound)   
    
    result.par = round(rbind(
      ans.donlp$solution, ans.solnp$solution, ans.nlminb$solution), 3)
    result.fun = c(
      fun(ans.donlp$solution), fun(ans.solnp$solution), fun(ans.nlminb$solution))
    cbind(result.par, result.fun)        
  }


# ------------------------------------------------------------------------------


.rosensuzukiTestNLP <-
  function()
  {
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    start <- c(1, 1, 1, 1)
    
    fun <- function(x) 
      x[1]*x[1]+x[2]*x[2]+2*x[3]*x[3]+x[4]*x[4]-5*x[1]-5*x[2]-21*x[3]+7*x[4]
    
    ineqFun <- list(
      function(x) 8-x[1]*x[1]-x[2]*x[2]-x[3]*x[3]-x[4]*x[4]-x[1]+x[2]-x[3]+x[4],
      function(x) 10-x[1]*x[1]-2*x[2]*x[2]-x[3]*x[3]-2*x[4]*x[4]+x[1]+x[4],
      function(x) 5-2*x[1]*x[1]-x[2]*x[2]-x[3]*x[3]-2*x[1]+x[2]+x[4] )
    ineqFun.lower <- rep(   0, 3)
    ineqFun.upper <- rep(10, 3)
    
    ans.donlp = donlp2NLP(start, fun,  
                          ineqFun = ineqFun, 
                          ineqFun.lower = ineqFun.lower, ineqFun.upper = ineqFun.upper)     
    ans.solnp = solnpNLP(start, fun,  
                         ineqFun = ineqFun, 
                         ineqFun.lower = ineqFun.lower, ineqFun.upper = ineqFun.upper)     
    ans.nlminb = nlminb2NLP(start, fun,  
                            ineqFun = ineqFun, 
                            ineqFun.lower = ineqFun.lower, ineqFun.upper = ineqFun.upper)     
    
    result.par = round(rbind(
      ans.donlp$solution, ans.solnp$solution, ans.nlminb$solution), 3)
    result.fun = c(
      fun(ans.donlp$solution), fun(ans.solnp$solution), fun(ans.nlminb$solution))
    cbind(result.par, result.fun)   
  }


#-------------------------------------------------------------------------------


.sharperatioTestNLP <-
  function()
  { 
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    dataSet <- data("LPP2005REC", package="timeSeries", envir=environment())
    LPP2005REC <- get(dataSet, envir=environment())
    ret <- as.matrix(LPP2005REC[,2:7])
    
    Mean <- colMeans(ret)
    Cov <- cov(ret)
    
    start <- rep(1/6, times = 6)
    
    fun <- function(x) {
      return = (Mean %*% x)[[1]]
      risk = sqrt ( (t(x) %*% Cov %*% x)[[1]] )
      -return/risk }
    par.lower <- rep(0, 6)
    par.upper <- rep(1, 6)  
    
    eqFun <- list(
      function(x) sum(x) )
    eqFun.bound = 1
    
    ans.donlp <- donlp2NLP(start, fun,  
                           par.lower = par.lower, par.upper = par.upper,
                           eqFun = eqFun, eqFun.bound = eqFun.bound ) 
    ans.solnp <- solnpNLP(start, fun,  
                          par.lower = par.lower, par.upper = par.upper,
                          eqFun = eqFun, eqFun.bound = eqFun.bound)      
    ans.nlminb <- Rnlminb2::nlminb2NLP(start, fun,  
                                       par.lower = par.lower, par.upper = par.upper,
                                       eqFun = eqFun, eqFun.bound = eqFun.bound)  
    
    result.par = round(rbind(
      ans.donlp$solution, 
      ans.solnp$solution, 
      ans.nlminb$solution), 3)
    result.fun = c(
      fun(ans.donlp$solution), 
      fun(ans.solnp$solution), 
      fun(ans.nlminb$solution))
    cbind(result.par, result.fun)     
    
    result.par <- 100*round(rbind(
      ans.donlp$solution, ans.solnp$solution, ans.nlminb$solution), 4)
    result.fun <- c(
      fun(ans.donlp$solution), fun(ans.solnp$solution), fun(ans.nlminb$solution))
    cbind(result.par, result.fun)  
  }


# ------------------------------------------------------------------------------


.rachevratioTestNLP <-
  function()
  {
    # A function implemented by Diethelm Wuertz
    
    # Note:
    #    sonlp2 and nlminb do not converge!
    
    # FUNCTION:
    
    dataSet <- data("LPP2005REC", package="timeSeries", envir=environment())
    LPP2005REC <- get(dataSet, envir=environment())
    ret = as.matrix(LPP2005REC[, 2:7])
    
    Mean = colMeans(ret)
    Cov = cov(ret)
    
    set.seed(1953)
    r = runif(6)
    start <- r/sum(r)
    start <- rep(1/6, times = 6)
    
    .VaR <- function(x) { 
      quantile(x, probs = 0.05, type = 1) }
    .CVaR <- function(x) {   
      VaR = .VaR(x)
      VaR - 0.5 * mean(((VaR-ret) + abs(VaR-ret))) / 0.05 }     
    fun <- function(x) {
      port = as.vector(ret %*% x)
      ans = (-.CVaR(-port) / .CVaR(port))[[1]] 
      ans}
    par.lower = rep(0, 6)
    par.upper = rep(1, 6)  
    
    eqFun <- list(
      function(x) sum(x) )
    eqFun.bound = 1
    
    ans.donlp = donlp2NLP(start, fun,  
                          par.lower = par.lower, par.upper = par.upper,
                          eqFun = eqFun, eqFun.bound = eqFun.bound) 
    ans.solnp = solnpNLP(start, fun,  
                         par.lower = par.lower, par.upper = par.upper,
                         eqFun = eqFun, eqFun.bound = eqFun.bound)   
    ans.nlminb = nlminb2NLP(start, fun,  
                            par.lower = par.lower, par.upper = par.upper,
                            eqFun = eqFun, eqFun.bound = eqFun.bound)    
    
    result.par = 100*round(rbind(
      ans.donlp$solution, ans.solnp$solution, ans.nlminb$solution), 4)
    result.fun = c(
      fun(ans.donlp$solution), fun(ans.solnp$solution), fun(ans.nlminb$solution))
    cbind(result.par, result.fun)  
  }


#-------------------------------------------------------------------------------


.boxmarkowitzTestNLP <-
  function()
  {
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    dataSet <- data("LPP2005REC", package="timeSeries", envir=environment())
    LPP2005REC <- get(dataSet, envir=environment())
    ret <- 100 * as.matrix(LPP2005REC[, 2:7])
    
    Mean <- colMeans(ret)
    Cov <- cov(ret)
    targetReturn <- mean(Mean)
    
    start <- rep(1/6, times = 6)
    
    objective <- function(x) {
      risk <- (t(x) %*% Cov %*% x)[[1]]
      risk }
    fun <- objective
    
    par.lower <- rep(0, 6)
    par.upper <- rep(1, 6)  
    
    eqFun <- list(
      function(x) sum(x), 
      function(x) (Mean %*% x)[[1]] )
    eqFun.bound = c(1, targetReturn)
    
    ans.donlp = donlp2NLP(start, objective, 
                          par.lower = par.lower, par.upper = par.upper,
                          eqFun = eqFun, eqFun.bound = eqFun.bound)
    ans.solnp = solnpNLP(start, objective, 
                         par.lower = par.lower, par.upper = par.upper,
                         eqFun = eqFun, eqFun.bound = eqFun.bound)
    ans.nlminb = nlminb2NLP(start, objective, 
                            par.lower = par.lower, par.upper = par.upper,
                            eqFun = eqFun, eqFun.bound = eqFun.bound)
    
    result.par <- 100*round(rbind(
      ans.donlp$solution, ans.solnp$solution, ans.nlminb$solution), 4)
    result.fun <- c(
      fun(ans.donlp$solution), fun(ans.solnp$solution), fun(ans.nlminb$solution))
    cbind(result.par, result.fun)  
  }


#-------------------------------------------------------------------------------


.groupmarkowitzTestNLP <-
  function()
  {
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    dataSet <- data("LPP2005REC", package="timeSeries", envir=environment())
    LPP2005REC <- get(dataSet, envir=environment())
    ret <- 100 * as.matrix(LPP2005REC[, 2:7])
    
    Mean <- colMeans(ret)
    Cov <- cov(ret)
    targetReturn <- mean(Mean)
    
    start <- rep(1/6, times = 6)
    start[1]+start[4]
    start[2]+start[5]+start[6]
    
    # Must be feasible for nlminb !!!
    start2 <- c(0, 0, 0.3, 0.3, 0, 0.4)
    start2 <- start/sum(start2)
    sum(start2)
    start2[1]+start2[4]
    start2[2]+start2[5]+start2[6]
    
    fun <- function(x) {
      risk = (t(x) %*% Cov %*% x)[[1]]
      risk }
    par.lower <- rep(0, 6)
    par.upper <- rep(1, 6)  
    
    eqFun <- list(
      function(x) sum(x), 
      function(x) (Mean %*% x)[[1]] )
    eqFun.bound = c(1, targetReturn)
    
    ineqFun <- list(
      function(x) x[1]+x[4],
      function(x) x[2]+x[5]+x[6])
    ineqFun.lower <- c( 0.3, 0.0)
    ineqFun.upper <- c( 1.0, 0.6)
    
    ans.donlp = donlp2NLP(start, fun, 
                          par.lower = par.lower, par.upper = par.upper,
                          eqFun = eqFun, eqFun.bound = eqFun.bound,
                          ineqFun = ineqFun, 
                          ineqFun.lower = ineqFun.lower, ineqFun.upper = ineqFun.upper)
    ans.solnp = solnpNLP(start, fun, 
                         par.lower = par.lower, par.upper = par.upper,
                         eqFun = eqFun, eqFun.bound = eqFun.bound,
                         ineqFun = ineqFun, 
                         ineqFun.lower = ineqFun.lower, ineqFun.upper = ineqFun.upper)
    ans.nlminb = nlminb2NLP(start, fun, 
                            par.lower = par.lower, par.upper = par.upper,
                            eqFun = eqFun, eqFun.bound = eqFun.bound,
                            ineqFun = ineqFun, 
                            ineqFun.lower = ineqFun.lower, ineqFun.upper = ineqFun.upper)
    
    result.par = 100*round(rbind(
      ans.donlp$solution, ans.solnp$solution, ans.nlminb$solution), 4)
    result.fun = c(
      fun(ans.donlp$solution), fun(ans.solnp$solution), fun(ans.nlminb$solution))
    cbind(result.par, result.fun) 
  }


################################################################################

