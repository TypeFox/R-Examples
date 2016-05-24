## file csampling/inst/csamplingdemo.R, v 1.0.1 2013/05/14
##
##  Copyright (C) 2000-2013 Alessandra R. Brazzale 
##
##  This file is part of the "csampling" package for R.  This program 
##  is free software; you can redistribute it and/or modify it under 
##  the terms of the GNU General Public License as published by the 
##  Free Software Foundation; either version 2 of the License, or (at 
##  your option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
##  MA 02111-1307 USA or look up the web page 
##  http://www.gnu.org/copyleft/gpl.html.
##
##  Please send any comments, suggestions or errors found to:
##  Alessandra R. Brazzale, Department of Statistics, University of
##  Padova, Via C. Battisti 241/243, 35121 Padova (PD), Italy.
##  Email: alessandra.brazzale@unipd.it
##  Web: http://www.stat.unipd.it/~brazzale


## ===================================================================
##
##       CONDITIONAL SIMULATION for the log-Weibull Distribution     
##
## ===================================================================


##  This file contains the R code for running a conditional simulation
##  study similar to the one presented in Brazzale (2000, Section
##  7.3.2).  The reference model is a linear regression model with
##  log-Weibull distributed errors, six covariates and ten
##  observations; see Example 3 of DiCiccio, Field and Fraser (1990).
##  The functions needed for the conditional simulation are provided
##  by the "csampling" package included in the "hoa" package bundle.
##  The "marg" package (inference for linear nonnormal models) must
##  also be loaded.  As Metropolis-Hastings sampling is used, not all
##  the computations could entirely be automated.  Some tasks, such as
##  the choice of the candidate generation density, graphical display
##  or diagnostics, have to be tackled by the user.  The following
##  demo is intended to trace down how it might be done.

	
## References:
## ----------
##
## Brazzale, A. R. (2000) Practical Small-Sample Parametric Inference.
##   Ph.D. Thesis N. 2230, Department of Mathematics, Swiss Federal 
##   Institute of Technology Lausanne. 
##
## DiCiccio, T. J., Field, C. A. and Fraser, D. A. S. (1990) 
##   Approximations of marginal tail probabilities and inference for 
##   scalar parameters.  Biometrika, 77, 77--95.


## --> Preamble to the conditional simulation.


library(csampling) ##  load "csampling" package
library(lattice)  ##  load "lattice" package
trellis.device()  ##  or suitable graphics device


## Simulation parameters
## ---------------------
##
X <- matrix(c(0.686,0.640,0.908,0.886,0.508,0.255,0.197,0.056,0.646,0.317,
              0.566,0.632,0.130,0.480,0.669,0.930,0.869,0.204,0.961,0.321,
              0  ,  0  ,  0  ,  0  ,  0  ,  1  ,  1  ,  1  ,  1  ,  1  , 
              0  ,  0  ,  0  ,  0  ,  0  ,0.255,0.197,0.056,0.646,0.317, 
              0  ,  0  ,  0  ,  0  ,  0  ,0.930,0.869,0.204,0.961,0.321, 
              1  ,  1  ,  1  ,  1  ,  1  ,  1  ,  1  ,  1  ,  1  ,  1  ),
              ncol = 6, byrow = FALSE,
              dimnames = list(c(), c("x1","x2","x3","x4","x5","x6")))
beta.true <- c( b1 = 20, b2 = -20, b3 = 10, b4 = -2, b5 = 2, b6 = 1 )
sigma.true <- 1


## Histogram of the marginal distribution of each observation
## ----------------------------------------------------------
##
sample.lW <- c()
for(i in 1:10)
  sample.lW <- cbind( sample.lW,
                      log(rweibull(50000, shape = sigma.true, 
                                   scale = exp(X[i,]%*%beta.true))) )
##
par( mfrow=c(2,5), pty="s" )
apply( sample.lW, 2, function(x) { hist(x, xlab="", ylab="", 
                                        prob=TRUE, nclass=30, las=1) ;
                                   lines(density(x)) } )


## A number of possible sample configurations
## ------------------------------------------
##
attach(as.data.frame(X))
par( mfrow=c(4,4), pty="s" )          
for(i in 1:16)
{
  y <- sample.lW[i*50,]
  rsmObject <- rsm( y ~ x1+x2+x3+x4+x5+x6-1, family = logWeibull, 
                    control = glm.control(maxit=200) )
  plot(rep(0,10), rsmObject$resid, xlab="", ylab="", xlim=c(-1,1), 
       cex=1, pch=16, las=1)
}                                        
detach() ; rm(rsmObject, y)


## Ancillary, dependent variable, MLEs & their unscaled asymptotic 
##  covariance matrix used in the simulation study
## ---------------------------------------------------------------
##
y.obs <- sample.lW[sample(16*50,1),]  ## random choice
##
attach(as.data.frame(X))
rsmObject <- rsm( y.obs ~ x1+x2+x3+x4+x5+x6-1, family = logWeibull, 
                  control = glm.control(maxit=50) )
beta.MLE <- rsmObject$coef
sigma.MLE <- rsmObject$disp
##
anc <- rsmObject$resid
##
var.MLE <- summary(rsmObject, corr=FALSE)$cov
corr.MLE <- summary(rsmObject, corr=TRUE)$corr
detach() 


## Laplace approximation to the univariate marginal densities of the 
##  MLEs (regression coefficients + scale parameter)
## -----------------------------------------------------------------
##
attach(as.data.frame(X))
X.data <- make.sample.data(rsmObject)
X.data$coef <- beta.true    ## replace estimated values by true values
X.data$disp <- sigma.true
detach()
##
uni.dens <- list()
uni.dens$b1 <- Laplace(c, X.data, seq(  5, 35, length=30), 1)$dens 
uni.dens$b2 <- Laplace(c, X.data, seq(-35, -5, length=30), 2)$dens 
uni.dens$b3 <- Laplace(c, X.data, seq(-10, 30, length=30), 3)$dens 
uni.dens$b4 <- Laplace(c, X.data, seq(-30, 20, length=30), 4)$dens 
uni.dens$b5 <- Laplace(c, X.data, seq(-15, 15, length=30), 5)$dens 
uni.dens$b6 <- Laplace(c, X.data, seq(-25, 20, length=30), 6)$dens    
##
uni.dens$ls  <- Laplace(s, X.data, seq(0.13, 1.8, length=30))$dens
##
par( mfrow=c(3,3), pty="s" ) 
sapply( uni.dens[1:6], 
        function(x)        
          plot(x, type="l", xlab=paste("regression coefficient"), 
                            ylab="marginal density", cex=0.7, las=1) )
plot( uni.dens$ls, type="l", xlab="log-scale parameter", 
                   ylab="marginal density", cex=0.7, las=1 )


## WARNING:
## -------
##  1) Beware that sometimes NAs are generated particularly in the 
##     tails of the distribution.  These are omitted when calculating 
##     the spline interpolation returned by the "Laplace" function.  
##  2) The spline interpolation can be "unstable" especially in the 
##     tails.


## Laplace's approximation to the bivariate marginal densities of the
##  MLEs (regression coefficients + scale parameter)
## ------------------------------------------------------------------
##
bi.dens <- list()
bi.dens$b12 <- Laplace(cc, X.data, seq(  5,35,length=30), 1, 
                                    seq(-35,-5,length=30), 2)
bi.dens$b13 <- Laplace(cc, X.data, seq(  5,35,length=30), 1, 
                                    seq(-10,30,length=30), 3)
bi.dens$b14 <- Laplace(cc, X.data, seq(  5,35,length=30), 1, 
                                    seq(-30,20,length=30), 4)
bi.dens$b15 <- Laplace(cc, X.data, seq(  5,35,length=30), 1, 
                                    seq(-15,15,length=30), 5)
bi.dens$b16 <- Laplace(cc, X.data, seq(  5,35,length=30), 1, 
                                    seq(-25,20,length=30), 6)
bi.dens$b23 <- Laplace(cc, X.data, seq(-35,-5,length=30), 2, 
                                    seq(-10,30,length=30), 3)
bi.dens$b24 <- Laplace(cc, X.data, seq(-35,-5,length=30), 2, 
                                    seq(-30,20,length=30), 4)
bi.dens$b25 <- Laplace(cc, X.data, seq(-35,-5,length=30), 2, 
                                    seq(-15,15,length=30), 5)
bi.dens$b26 <- Laplace(cc, X.data, seq(-35,-5,length=30), 2, 
                                    seq(-25,20,length=30), 6)
bi.dens$b34 <- Laplace(cc, X.data, seq(-10,30,length=30), 3, 
                                    seq(-30,20,length=30), 4)
bi.dens$b35 <- Laplace(cc, X.data, seq(-10,30,length=30), 3, 
                                    seq(-15,15,length=30), 5)
bi.dens$b36 <- Laplace(cc, X.data, seq(-10,30,length=30), 3, 
                                    seq(-25,20,length=30), 6)
bi.dens$b45 <- Laplace(cc, X.data, seq(-30,20,length=30), 4, 
                                    seq(-15,15,length=30), 5)
bi.dens$b46 <- Laplace(cc, X.data, seq(-30,20,length=30), 4, 
                                    seq(-25,20,length=30), 6)
bi.dens$b56 <- Laplace(cc, X.data, seq(-15,15,length=30), 5, 
                                    seq(-25,20,length=30), 6)
##
bi.dens$b1ls <- Laplace(cs, X.data, seq(  5,35,length=30), 1, 
                                     seq(0.1, 1.6,0.3))
bi.dens$b2ls <- Laplace(cs, X.data, seq(-35,-5,length=30), 2, 
                                     seq(0.1, 1.6,0.3))
bi.dens$b3ls <- Laplace(cs, X.data, seq(-10,30,length=30), 3, 
                                     seq(0.1, 1.6,0.3))
bi.dens$b4ls <- Laplace(cs, X.data, seq(-30,20,length=30), 4, 
                                     seq(0.1, 1.6,0.3))
bi.dens$b5ls <- Laplace(cs, X.data, seq(-15,15,length=30), 5, 
                                     seq(0.1, 1.6,0.3))
bi.dens$b6ls <- Laplace(cs, X.data, seq(-25,20,length=30), 6, 
                                     seq(0.1, 1.6,0.3))
##
par( mfrow=c(4,6), pty="s" ) 
sapply( bi.dens, function(x) plot(x, nlevels=10, las=1) )


## WARNING: 
## -------
##  Beware that sometimes NAs or infinite values are generated 
##  particularly in the tails of the distribution.  These are set 
##  to 0 by the "Laplace" function. 


## ==================================================================
## NOTE on Laplace's approximation to the uni- and bivariate marginal 
##              densities of the MLEs
## ==================================================================
##
##  The plots of the uni- and bivariate marginal densities of the MLEs
##  give insight into how the 7-dimensional conditional distribution
##  from which we want to sample from is structured.  They help in
##  choosing an appropriate candidate generation density.
## 
##  In our example, a good choice is to sample from a multivariate t 
##  distribution with few degrees of freedom, say 5.  This 
##  distribution should be re-centered and re-scaled such as to 
##  maximize the acceptance rate.  Let's take the true value of the 
##  regression coefficients as corresponding multivariate location 
##  parameter.  For the scale parameter we have to re-center the 
##  candidate generation density at 
##
loc.sim <- c(beta.true, -0.7)
##
##  This distribution is further re-scaled by a multiple of the 
##  asymptotic covariance matrix of the MLEs.  This way, we hope to 
##  be able to capture as much as possible the dependence structure 
##  of the conditional distribution.
##
var.sim <- 1.8*round(var.MLE, dig=4) 
var.sim[7,] <- 1.5*var.sim[7,]  
var.sim[,7] <- 1.5*var.sim[,7]  
##
##  The scaling factors have been found empirically.  These choices 
##  should ensure an acceptance rate of about 25%.


par( mfrow=c(3,3) )
#
plot(uni.dens$b1, type="l")
lines(uni.dens$b1$x, 
      1/sqrt(var.sim[1,1])*dt((uni.dens$b1$x-loc.sim[1])/
        sqrt(var.sim[1,1]), df=3), lwd=2)
#
plot(uni.dens$b2, type="l")
lines(uni.dens$b2$x, 
      1/sqrt(var.sim[2,2])*dt((uni.dens$b2$x-loc.sim[2])/
        sqrt(var.sim[2,2]), df=3), lwd=2)
#
plot(uni.dens$b3, type="l")
lines(uni.dens$b3$x,
      1/sqrt(var.sim[3,3])*dt((uni.dens$b3$x-loc.sim[3])/
        sqrt(var.sim[3,3]), df=3), lwd=2)
#
plot(uni.dens$b4, type="l")
lines(uni.dens$b4$x, 
      1/sqrt(var.sim[4,4])*dt((uni.dens$b4$x-loc.sim[4])/
        sqrt(var.sim[4,4]), df=3), lwd=2)
#
plot(uni.dens$b5, type="l")
lines(uni.dens$b5$x, 
      1/sqrt(var.sim[5,5])*dt((uni.dens$b5$x-loc.sim[5])/
        sqrt(var.sim[5,5]), df=3), lwd=2)
#
plot(uni.dens$b6, type="l")
lines(uni.dens$b6$x, 
      1/sqrt(var.sim[6,6])*dt((uni.dens$b6$x-loc.sim[6])/
        sqrt(var.sim[6,6]), df=3), lwd=2)
#
plot(uni.dens$ls, type="l")
lines(uni.dens$ls$x, 
      1/sqrt(var.sim[7,7])*dt((uni.dens$ls$x-loc.sim[7])/
        sqrt(var.sim[7,7]), df=3), lwd=2)


par( mfrow=c(4,6), pty="s" ) 
#
plot(bi.dens$b12)
points(t(rmt(100, df=5, mm=loc.sim[c(1,2)], 
           cov=var.sim[c(1,2),c(1,2)])), pch="+", cex=0.6)
plot(bi.dens$b13)
points(t(rmt(100, df=5, mm=loc.sim[c(1,3)], 
           cov=var.sim[c(1,3),c(1,3)])), pch="+", cex=0.6)
plot(bi.dens$b14)
points(t(rmt(100, df=5, mm=loc.sim[c(1,4)], 
           cov=var.sim[c(1,4),c(1,4)])), pch="+", cex=0.6)
plot(bi.dens$b15)
points(t(rmt(100, df=5, mm=loc.sim[c(1,5)], 
           cov=var.sim[c(1,5),c(1,5)])), pch="+", cex=0.6)
plot(bi.dens$b16)
points(t(rmt(100, df=5, mm=loc.sim[c(1,6)], 
           cov=var.sim[c(1,6),c(1,6)])), pch="+", cex=0.6)
plot(bi.dens$b23)
points(t(rmt(100, df=5, mm=loc.sim[c(2,3)], 
           cov=var.sim[c(2,3),c(2,3)])), pch="+", cex=0.6)
plot(bi.dens$b24)
points(t(rmt(100, df=5, mm=loc.sim[c(2,4)], 
           cov=var.sim[c(2,4),c(2,4)])), pch="+", cex=0.6)
plot(bi.dens$b25)
points(t(rmt(100, df=5, mm=loc.sim[c(2,5)], 
           cov=var.sim[c(2,5),c(2,5)])), pch="+", cex=0.6)
plot(bi.dens$b26)
points(t(rmt(100, df=5, mm=loc.sim[c(2,6)], 
           cov=var.sim[c(2,6),c(2,6)])), pch="+", cex=0.6)
plot(bi.dens$b34)
points(t(rmt(100, df=5, mm=loc.sim[c(3,4)], 
           cov=var.sim[c(3,4),c(3,4)])), pch="+", cex=0.6)
plot(bi.dens$b35)
points(t(rmt(100, df=5, mm=loc.sim[c(3,5)], 
           cov=var.sim[c(3,5),c(3,5)])), pch="+", cex=0.6)
plot(bi.dens$b36)
points(t(rmt(100, df=5, mm=loc.sim[c(3,6)], 
           cov=var.sim[c(3,6),c(3,6)])), pch="+", cex=0.6)
plot(bi.dens$b45)
points(t(rmt(100, df=5, mm=loc.sim[c(4,5)], 
           cov=var.sim[c(4,5),c(4,5)])), pch="+", cex=0.6)
plot(bi.dens$b46)
points(t(rmt(100, df=5, mm=loc.sim[c(4,6)], 
           cov=var.sim[c(4,6),c(4,6)])), pch="+", cex=0.6)
plot(bi.dens$b56)
points(t(rmt(100, df=5, mm=loc.sim[c(5,6)], 
           cov=var.sim[c(5,6),c(5,6)])), pch="+", cex=0.6)
plot(bi.dens$b1ls)
points(t(rmt(100, df=5, mm=loc.sim[c(1,7)], 
           cov=var.sim[c(1,7),c(1,7)])), pch="+", cex=0.6)
plot(bi.dens$b2ls)
points(t(rmt(100, df=5, mm=loc.sim[c(2,7)], 
           cov=var.sim[c(2,7),c(2,7)])), pch="+", cex=0.6)
plot(bi.dens$b31ls)
points(t(rmt(100, df=5, mm=loc.sim[c(3,7)], 
           cov=var.sim[c(3,7),c(3,7)])), pch="+", cex=0.6)
plot(bi.dens$b4ls)
points(t(rmt(100, df=5, mm=loc.sim[c(4,7)], 
           cov=var.sim[c(4,7),c(4,7)])), pch="+", cex=0.6)
plot(bi.dens$b5ls)
points(t(rmt(100, df=5, mm=loc.sim[c(5,7)], 
           cov=var.sim[c(5,7),c(5,7)])), pch="+", cex=0.6)
plot(bi.dens$b6ls)
points(t(rmt(100, df=5, mm=loc.sim[c(6,7)], 
           cov=var.sim[c(6,7),c(6,7)])), pch="+", cex=0.6)


## --> Let's start with the conditional simulation.


## ==============================================================
## Sampling from the conditional distribution of the MLEs for the 
##   given value of the ancillary
## ==============================================================
##
##  The workhorse for doing this is the "rsm.sample" function
##  contained in the "csampling" package of the "hoa" bundle.  First 
##  of all we need to define a function which generates the candidate
##  values.  
##
xmpl.gen <- function(R = R, data = data, mm, cov, df = 5, ...)
{
  div <- R%/%100        
  mod <- R%%100         
  nc <- length(mm)+1    
  ret <- matrix(0, nrow=R, ncol=nc)
  if( div != 0 )
    for(i in 1:div)
    { 
      seq <- ((i-1)*100+1):(i*100)
      ret[seq,-nc] <- t(rmt(n=100, df=df, mm=mm, cov=cov))
      ret[seq,nc] <- apply(ret[seq,-nc], 1, "dmt", df=df, 
                                                   mm=mm, cov=cov)
      print(i*100)
    }
  if( mod != 0 )
  {      
    seq <- (div*100+1):R 
    ret[seq,-nc] <- t(rmt(n=mod, df=df, mm=mm, cov=cov))
    ret[seq,nc] <- apply( matrix(ret[seq,-nc], nrow=mod), 1, "dmt", 
                          df=df, mm=mm, cov=cov)               
  }    
  ret[,nc] <- exp(-ret[,nc-1])*ret[,nc]
  ret[,nc-1] <- exp(ret[,nc-1])
  ret
}

cond.sample <- rsm.sample(X.data, R = 100000, ran.gen = xmpl.gen, 
                          mm = loc.sim, cov = var.sim)


## Observations corresponding to the sampled MLEs
## ----------------------------------------------
attach(cond.sample)
y.sim <- sim[,1:6] %*% t(X) + as.vector(sim[,7]) %*% t(as.vector(anc))
detach()


## Acceptance rate
## ---------------
attach(cond.sample)
a.rate <- 1-sum(sim[-1,1]==sim[-dim(sim)[1],1])/dim(sim)[1]
detach()


## Graphical inspection of the generated Markov chain and diagnostics
## ==================================================================
##
## Chain
## -----
par( mfrow=c(3,1), pty="m", ask=TRUE )                    
apply( cond.sample$sim, 2, function(x) plot(x, type="l", las=1) )
par( ask=FALSE )


## Simulated marginal density 
##   --> to compare with Laplace approximation
## -------------------------------------------
par( mfrow=c(3,3), pty="s" )  
apply(cond.sample$sim, 2, function(x) {
                            hist(x, prob=TRUE, nclass=20, xlab="", 
                                 las=1) ; 
                            lines(density(x)) } )
hist(log(cond.sample$sim[,7]), prob=TRUE, nclass=20, xlab="", las=1)
lines(density(log(cond.sample$sim[,7])))
      

## Bivariate scatter plots --> to compare with Laplace approximation
## ----------------------- --> to compare with `corr.MLE'
sim <- cond.sample$sim[1:1000*100,]
sim[,7] <- log(sim[,7])  ## only a subset of values used
pairs(as.data.frame(sim), 
      labels=c("b1","b2","b3","b4","b5","b6","log-scale"),
      las=1, pch="+", cex=0.5)
rm(sim) 


## Autocorrelation
## ---------------
par(mfrow=c(2,1))
acf(cond.sample$sim[,1], lag=10)
acf(cond.sample$sim[,1], lag=10, type="partial")
##
acf(cond.sample$sim[,1], lag=10, plot=FALSE)
##
## --> Looks like an AR(1) process with autocorrelation ~ 0.75


## Conditional distribution of dependent variable 
##  (compared with marginal one)
## ----------------------------------------------
par( mfrow=c(4,5), pty="s" )      
apply( y.sim[,1:5], 2, function(x) {
                         hist(x, prob=TRUE, nclass=20, xlab="", 
                              ylab="", las=1)
                         lines(density(x)) })
apply( sample.lW[,1:5], 2, function(x) {
                             hist(x, prob=TRUE, nclass=20, xlab="", 
                                  ylab="", las=1)
                             lines(density(x)) })
apply( y.sim[,6:10], 2, function(x) {
                         hist(x, prob=TRUE, nclass=20, xlab="", 
                              ylab="", las=1)
                         lines(density(x)) })
apply( sample.lW[,6:10], 2, function(x) {
                             hist(x, prob=TRUE, nclass=20, xlab="",  
                                  ylab="", las=1)
                             lines(density(x)) })


## Diagnostics using CODA
## ----------------------
## 
##  CODA
##  (http://cran.r-project.org/src/contrib/Descriptions/coda.html)
##  stands for Convergence Diagnostic and Output Analysis software
##  (http://www-fis.iarc.fr/coda/) and is a set of R functions which
##  serves as an output processor for the BUGS software
##  (http://www.mrc-bsu.cam.ac.uk/bugs/welcome.shtml).  It may also 
##  be used in conjuction with Markov chain Monte Carlo output from 
##  a user's own programs, providing the output is formatted
##  appropriately.  See the CODA manual for details.
## 
##  A CODA session can be menu-driven.  All we need are two files 
##  `b1.out' and `b1.ind', the first one containing the generated 
##  chain, the second one a summary of the data.  (For the exact 
##  format, see the CODA documentation.)


library(coda) ##  load "coda" package

attach(cond.sample)
write.table(sim[,1], "b1.out", row.names = TRUE, col.names = FALSE,
            quote = FALSE)
detach()
##
string <- "b1.sim       1    100000"  ## 100000 is the length of the
	                              ## Metropolis-Hastings chain
write.table(string, "b1.ind", row.names = FALSE, col.names = FALSE,
            quote = FALSE)
rm(string)
##
##  Start the session with 
##
codamenu()
##
##  specifying as input file `b1'.  Note that as the MLEs of all
##  parameters have been generated simultaneously, you only need to
##  run the diagnostics on one chain. 
## 
##  The CODA software admits also multiple chains.


## Diagnostics using BOA
## ---------------------
## 
##  BOA (http://cran.r-project.org/src/contrib/Descriptions/boa.html)
##  stands for Bayesian Output Analysis Program
##  (http://www.public-health.uiowa.edu/boa/) and is an R program for
##  carrying out convergence diagnostics and statistical and graphical
##  analysis of Monte Carlo sampling output.  It can be used as an
##  output processor for the BUGS software
##  (http://www.mrc-bsu.cam.ac.uk/bugs/welcome.shtml) or for any other
##  program which produces sampling output.

##  A BOA session can be menu-driven.  All we need is a file `b1.txt'
##  which contains the generated chain.  (For the exact format, see 
##  the BOA documentation.)


library(boa) ##  load "boa" package

attach(cond.sample)
write.table(cbind(1:nrow(sim), sim[,1]), "b1.txt", quote = FALSE, 
            row.names = FALSE, col.names = c("iter", "b1"))
detach()
##
##  Start the session with 
##
boa.menu()
##
##  specifying as input file `b1'.  Note that as the MLEs of all
##  parameters have been generated simultaneously, you only need to
##  run the diagnostics on one chain. 
## 
##  The BOA software admits also multiple chains.


## ===============================================================
## NOTE on the simulation from the conditional distribution of the 
##  MLEs
## ===============================================================
##
##  The sample of MLEs from their conditional distribution represents
##  the basis for further computations.  In fact, once these values
##  are available, several other quantities of interest, such as the
##  distribution of test statistics, the coverage level of confidence
##  intervals and many more, may be derived from them, and their
##  properties from the ergodicity of the chain.  The following
##  portions of code show some examples.


## Distribution of pivots Q1 and Q2
## ================================
##
attach(cond.sample)
q1.sim <- ((sim[,1:6] - matrix(rep(beta.true, dim(sim)[1]), ncol=6, 
                               byrow=TRUE)) /
          matrix(sim[,7], ncol=6, nrow=dim(sim)[1]))[-(1:500),]
q2.sim <- (sim[,7]/sigma.true)[-(1:500)]       
                                               ## burn-in = 500
detach()
##
par( mfrow=c(3,3), pty="s" )
apply( q1.sim, 2, function(x) hist(x, prob=TRUE, xlab="", nclass=20, 
                                   las=1) )
hist(q2.sim, prob=TRUE, xlab="", las=1)
hist(log(q2.sim), prob=TRUE, xlab="", las=1)


## Distribution of R and R*
## ========================

## 1) parameter of interest = scale parameter
## ------------------------------------------
##
attach(cond.sample) ; attach(X.data) ; attach(as.data.frame(X))
##
R <- dim(sim)[1] - 500          ## burn-in : 500
##
MLE.s <- matrix(ncol=6, nrow=R)          ## constrained MLEs
r.sim <- r.star.sim <- u.sim <- vector(mode="numeric", length=R)
##
n <- length(anc) ; p <- dim(X)[2] ; sigma0 <- X.data$disp  
g0 <- family$g0 ; g1 <- family$g1 ; g2 <- family$g2  
df <- family$df ; k <- family$k
##
E.a <- diag(as.vector(g2(anc, df=df, k=k))) 
X.a <- cbind(X, anc)  
X.a <- t(X.a) %*% E.a %*% X.a
##
y <- y.sim[501,]          ## burn-in : 500
rsmObject  <- rsm(y ~ x1 + x2 + x3 + x4 + x5 + x6 - 1, 
                  family = logWeibull, disp = sigma0, 
                  control = glm.control(maxit=50))   
## offset : sigma0
##
MLE.s[1,] <- beta0 <- coef(rsmObject)  
beta1 <- sim[501,1:6] ; sigma1 <- sim[501,7]
Anc <- rsmObject$resid
E.A <- diag(as.vector(g2(Anc, df=df, k=k)))
H.A <- t(X) %*% E.A %*% X / (sigma0^2)
h.1 <- sum(g1(Anc, df=df, k=k)*anc)/sigma0 - n/sigma1
H.X.a <- X.a / (sigma1^2)
H.X.a[p+1,p+1] <- H.X.a[p+1,p+1] + n / (sigma1^2)
u <- h.1 * sqrt(abs(det(H.A))) /
           sqrt(abs(det(H.X.a)))
u.sim[1] <- u
r.sim[1] <- sqrt(2*( n*log(sigma0/sigma1) - 
                     sum(g0(anc, df=family$df, k=family$k)) +
                     sum(g0(Anc, df=family$df, k=family$k)) ))*       
                 sign(sigma1-sigma0) 
r.star.sim[1] <- r.sim[1] + log(abs(u/r.sim[1]))/r.sim[1]
##
for( i in 502:dim(sim)[1] )          
{
  idx <- i-500
  if( sim[idx,1] == sim[idx-1,1] ) 
  {
    r.sim[idx] <- r.sim[idx-1]
    r.star.sim[idx] <- r.star.sim[idx-1]
    MLE.s[idx,] <- MLE.s[idx-1,]
  }     
  else
  {      
    y <- y.sim[i,]
    rsmObject <- rsm(y ~ x1 + x2 + x3 + x4 + x5 + x6 - 1, 
                     family = logWeibull, disp = sigma0, 
                     control = glm.control(maxit=50))
    MLE.s[idx,] <- beta0 <- coef(rsmObject)  
    beta1 <- sim[i,1:6]  
    sigma1 <- sim[i,7]
    Anc <- rsmObject$resid
    E.A <- diag(as.vector(g2(Anc, df=df, k=k)))
    H.A <- t(X) %*% E.A %*% X / (sigma0^2)
    h.1 <- sum(g1(Anc, df=df, k=k)*anc)/sigma0 - n/sigma1
    H.X.a <- X.a / (sigma1^2)
    H.X.a[p+1,p+1] <- H.X.a[p+1,p+1] + n / (sigma1^2)
    u <- h.1 * sqrt(abs(det(H.A))) /
               sqrt(abs(det(H.X.a)))
    u.sim[idx] <- u
    r.sim[idx] <- sqrt(2*(n*log(sigma0/sigma1) - 
                       sum(g0(anc, df=family$df, k=family$k)) +
                       sum(g0(Anc, df=family$df, k=family$k)))) *
                  sign(sigma1-sigma0)  
    r.star.sim[idx] <- r.sim[idx] + 
                       log(abs(u/r.sim[idx]))/r.sim[idx]
  } 
  if(i%%100 == 0) print(i-500)
}
##
detach() ; detach() ; detach()
rm(R, rsmObject, g0, g1, g2, df, k, y, n, p, sigma0, sigma1, 
   beta0, beta1, E.a, X.a, Anc, E.A, H.A, h.1, H.X.a, u)
##
r.s.sim <- r.sim
r.star.s.sim <- r.star.sim            
u.s.sim <- u.sim
##
rm(r.sim, r.star.sim, u.sim)

## WARNING: 
## -------
##  Beware that convergence problems with "rsm" may occur, though this
##  happens quite rarely.  This is typically the case when the
##  simulated MLEs lie in the queues of the distribution, in which
##  case the observation can be omitted.

## NOTE: 
## ----
##  The statistics r and r* could have been calculated with the
##  "rsm.marg" routine of the "marg" package.  This routine, however,
##  gives the entire profile of the statistics, whereas we only need
##  them at the true parameter value.  In order to save processing
##  time we preferred not to use it and to calculate instead the
##  correction term u separately.


## 2) parameter of interest = b1
## -----------------------------
##
attach(cond.sample) ; attach(X.data) ; attach(as.data.frame(X))
##
R <- dim(sim)[1] - 500          ## burn-in : 500
##
wh <- 1           ## index of regression coefficient
MLE.b1 <- matrix(ncol=6, nrow=R)          ## constrained MLEs
r.sim <- r.star.sim <- u.sim <- vector(mode="numeric", length=R)
##
n <- length(anc) ; p <- dim(X)[2] ; b0 <- X.data$coef[1]  
g0 <- family$g0 ; g1 <- family$g1 ; g2 <- family$g2  
df <- family$df ; k <- family$k
##
E.a <- diag(as.vector(g2(anc, df=df, k=k))) 
Xa <- cbind(X, anc)  
X.a <- t(Xa) %*% E.a %*% Xa
##
y <- y.sim[501,]
rsmObject  <- rsm(y ~ offset(b0*x1) + x2 + x3 + x4 + x5 + x6 - 1, 
                  family = logWeibull, control = glm.control(maxit=50))
## offset : b0*x1
##
MLE.b1[1,] <- c(coef(rsmObject), rsmObject$disp)  
beta0 <- c(b0, coef(rsmObject))  
beta1 <- sim[501,1:6]  
b1 <- beta1[1]
sigma0 <- rsmObject$disp  
sigma1 <- sim[501,7]
Anc <- rsmObject$resid
E.A <- diag(as.vector(g2(Anc, df=df, k=k)))
H.X.a <- X.a / (sigma1^2)
H.X.a[p+1,p+1] <- H.X.a[p+1,p+1] + n / (sigma1^2)
X.A <- cbind(X[,-wh], Anc)
H.A <- t(X.A) %*% E.A %*% X.A / (sigma0^2)
H.A[p,p] <- H.A[p,p] + n/(sigma0^2)
q0.A <- g1(Anc, df=df, k=k)
h.1 <- matrix(data=0, ncol=p+1, nrow=p+1)
h.1[,-1] <- t(cbind(Xa[,-wh],Xa[,wh])) %*% E.A %*% X.A / 
              (sigma0^2)
h.1[p,p+1] <- h.1[p,p+1] + sum(anc*q0.A)/(sigma0^2)
h.1[p+1,p+1] <- h.1[p+1,p+1] + sum(X[,wh]*q0.A)/(sigma0^2)
h.1[p,1] <- 1/sigma0 * sum(anc*q0.A) - n/sigma1
h.1[p+1,1] <- 1/sigma0 * sum(X[,wh]*q0.A)
u <- abs(det(h.1))/sqrt(abs(det(H.X.a)))/
        sqrt(abs(det(H.A)))
u.sim[1] <- u
r.sim[1] <- sqrt(2*( n*log(sigma0/sigma1) -
                     sum(g0(anc, df=family$df, k=family$k)) +
                     sum(g0(Anc, df=family$df, k=family$k)) )) *
                 sign(b1-b0) 
r.star.sim[1] <- r.sim[1] + log(abs(u/r.sim[1]))/r.sim[1]
##
for( i in 502:dim(sim)[1] )          
{
  idx <- i-500
  if( sim[idx,1] == sim[idx-1,1] ) 
  {
    r.sim[idx] <- r.sim[idx-1]
    r.star.sim[idx] <- r.star.sim[idx-1]
    MLE.b1[idx,] <- MLE.b1[idx-1,]
  }  
  else
  {      
    y <- y.sim[i,]
    rsmObject  <- rsm(y ~ offset(b0*x1) + x2 + x3 + x4 + x5 + x6 - 1, 
                      family = logWeibull, 
                      control=glm.control(maxit=50))
    MLE.b1[idx,] <- c(coef(rsmObject), rsmObject$disp)  
    beta0 <- c(b0, coef(rsmObject))  
    beta1 <- sim[i,1:6]  
    b1 <- beta1[1]
    sigma0 <- rsmObject$disp  
    sigma1 <- sim[i,7]
    Anc <- ( X %*% (beta1-beta0) + sigma1*anc ) / sigma0
    E.A <- diag(as.vector(g2(Anc, df=df, k=k)))
    H.X.a <- X.a / (sigma1^2)
    H.X.a[p+1,p+1] <- H.X.a[p+1,p+1] + n / (sigma1^2)
    X.A <- cbind(X[,-wh], Anc)
    H.A <- t(X.A) %*% E.A %*% X.A / (sigma0^2)
    H.A[p,p] <- H.A[p,p] + n/(sigma0^2)
    q0.A <- g1(Anc, df=df, k=k)
    h.1 <- matrix(data=0, ncol=p+1, nrow=p+1)
    h.1[,-1] <- t(cbind(Xa[,-wh],Xa[,wh])) %*% E.A %*% X.A /
                  (sigma0^2)
    h.1[p,p+1] <- h.1[p,p+1] + sum(anc*q0.A)/(sigma0^2)
    h.1[p+1,p+1] <- h.1[p+1,p+1] + sum(X[,wh]*q0.A)/(sigma0^2)
    h.1[p,1] <- 1/sigma0 * sum(anc*q0.A) - n/sigma1
    h.1[p+1,1] <- 1/sigma0 * sum(X[,wh]*q0.A)
    u <- abs(det(h.1)) / 
             sqrt(abs(det(H.X.a)))/
            sqrt(abs(det(H.A)))
    u.sim[idx] <- u
    r.sim[idx] <- sqrt(2*( n*log(sigma0/sigma1) -
                         sum(g0(anc, df=family$df, k=family$k)) +
                         sum(g0(Anc, df=family$df, k=family$k)) )) *
                  sign(b1-b0) 
    r.star.sim[idx] <- r.sim[idx] + 
                         log(abs(u/r.sim[idx]))/r.sim[idx]
  } 
  if(i%%100 == 0) print(i-500)
}
##
detach() ; detach() ; detach()
rm(R, rsmObject, g0, g1, g2, df, k, y, n, p, sigma0, sigma1, beta0,   
   beta1, E.a, X.a, Xa, Anc, E.A, H.A, h.1, H.X.a, u, wh, b0, b1, 
   q0.A)
##
r.b1.sim <- r.sim
r.star.b1.sim <- r.star.sim            
u.b1.sim <- u.sim
##
rm(r.sim, r.star.sim, u.sim)

## WARNING: 
## -------
##  Beware that convergence problems with "rsm" may occur, though this
##  happens quite rarely.  This is typically the case when the
##  simulated MLEs lie in the queues of the distribution, in which
##  case the observation can be omitted.


## 3) and so on with the other regression coefficients ...
## -------------------------------------------------------
##
## ......


## Histogram, density estimate and normal QQ-plot of R and R*
## ----------------------------------------------------------
##
## 1) parameter of interest: scale parameter
##
par( mfrow=c(2,2), pty="s" )
##
hist(r.s.sim, prob=TRUE, xlab="", las=1)  
lines(density(r.s.sim))
curve(dnorm, min(r.s.sim), max(r.s.sim), add=TRUE, lwd=2)
##      
hist(r.star.s.sim, prob=TRUE, xlab="", las=1)  
lines(density(r.star.s.sim))
curve(dnorm, min(r.star.s.sim), max(r.star.s.sim), add=TRUE, lwd=2)
##
qqnorm(sort(r.s.sim), xlim=c(-4,4), ylim=c(-4,4), las=1, type="s",
       lwd=2) 
abline(0,1)
##
qqnorm(sort(r.star.s.sim), xlim=c(-4,4), ylim=c(-4,4), las=1, 
       type="s", lwd=2)  
abline(0,1)
##
##
## 2) parameter of interest: b1
##
par( mfrow=c(2,2), pty="s" )
##
hist(r.b1.sim, prob=TRUE, xlab="", las=1)  
lines(density(r.b1.sim))
curve(dnorm, min(r.b1.sim), max(r.b1.sim), add=TRUE, lwd=2)
##      
hist(r.star.b1.sim, prob=TRUE, xlab="", las=1)  
lines(density(r.star.b1.sim))
curve(dnorm, min(r.star.b1.sim), max(r.star.b1.sim), add=TRUE, lwd=2) 
##
qqnorm(sort(r.b1.sim), xlim=c(-4,4), ylim=c(-4,4), las=1, type="s",
       lwd=2) 
abline(0,1)
##
qqnorm(sort(r.star.b1.sim), xlim=c(-4,4), ylim=c(-4,4), las=1, 
       type="s", lwd=2)  
abline(0,1)
##
##
## 3) and so on ...


## Goodness-of-fit to Normal distribution: Shapiro-Wilk
## -----------------------------------------------------------
##  The Shapiro-Wilk test statistic requires and independent sample. 
##  In order to achieve this, we take only every 20th observation.
##
## 1) parameter of interest: scale parameter
##
wh <- length(r.s.sim) %/% 20
##
shapiro.test(r.s.sim[(0:wh)*20+1])
shapiro.test(r.star.s.sim[(0:wh)*20+1])
##
rm(wh)
##
##
## 2) parameter of interest: b1
##
wh <- length(r.b1.sim) %/% 20
##
shapiro.test(r.b1.sim[(0:wh)*20+1])
shapiro.test(r.star.b1.sim[(0:wh)*20+1])
##
rm(wh)
##
##
## 3) and so on ...


## Conditional coverage of 95% confidence intervals
## ------------------------------------------------
##
## 1) parameter of interest: scale parameter
##
ll <- length(r.s.sim)
cov.r.s <- list( low = cumsum(ifelse( r.s.sim < (-1.96), 1, 0 ))/1:ll,
                 up = cumsum(ifelse( r.s.sim < 1.96, 1, 0 ))/1:ll)
cov.rs.s <- list( low = 
                 cumsum(ifelse( r.star.s.sim < (-1.96), 1, 0 ))/1:ll,
                  up = 
                 cumsum(ifelse( r.star.s.sim < 1.96, 1, 0 ))/1:ll)
##
par( mfrow=c(2,1) , pty="m" )
plot(0, 0, ylim=c(0,1), xlim=c(0,ll), cex=1, ylab="", bty="o", las=1,
     type="n", xlab="")
polygon( c(0,ll,ll,0), c(0.05,0.05,0.95,0.95), col="yellow")
lines(1:ll, cov.r.s$low, lwd=2)
lines(1:ll, cov.r.s$up, lwd=2)
##
plot(0, 0, ylim=c(0,1), xlim=c(0,ll), cex=1, ylab="", bty="o", las=1,
     type="n", xlab="")
polygon( c(0,ll,ll,0), c(0.05,0.05,0.95,0.95), col="yellow")
lines(1:ll, cov.rs.s$low, lwd=2)
lines(1:ll, cov.rs.s$up, lwd=2)
##
##
## 2) parameter of interest: b1
##
ll <- length(r.b1.sim)
cov.r.b1 <- list( low = 
                   cumsum(ifelse( r.b1.sim < (-1.96), 1, 0 ))/1:ll,
                  up = 
                   cumsum(ifelse( r.b1.sim < 1.96, 1, 0 ))/1:ll)
cov.rs.b1 <- list( low = 
                 cumsum(ifelse( r.star.b1.sim < (-1.96), 1, 0 ))/1:ll,
                  up = 
                 cumsum(ifelse( r.star.b1.sim < 1.96, 1, 0 ))/1:ll)
##
par( mfrow=c(2,1) , pty="m" )
plot(0, 0, ylim=c(0,1), xlim=c(0,ll), cex=1, ylab="", bty="o", las=1,
     type="n", xlab="")
polygon( c(0,ll,ll,0), c(0.05,0.05,0.95,0.95), col="yellow")
lines(1:ll, cov.r.b1$low, lwd=2)
lines(1:ll, cov.r.b1$up, lwd=2)
##
plot(0, 0, ylim=c(0,1), xlim=c(0,ll), cex=1, ylab="", bty="o", las=1,
     type="n", xlab="")
polygon( c(0,ll,ll,0), c(0.05,0.05,0.95,0.95), col="yellow")
lines(1:ll, cov.rs.b1$low, lwd=2)
lines(1:ll, cov.rs.b1$up, lwd=2)
##
##
## 3) and so on ...





## --> Let's re-do it but this time by sampling from a contaminated 
##      distribution.





## ================================================================
## Sampling from a contaminated distribution for the MLEs given the
##  value of the ancillary
## ================================================================
##
##  We will consider two contaminating distributions: 
##  log-Weibull(0,10) and log-Exp(4).

eps <- 0.1		## contamination rate
ss.cont <- 10		## scale parameter for log-Weibull(0,10)
par.true <- c(beta.true, sigma.true)


## Comparison between `true' and contaminating distribution
## --------------------------------------------------------
##
## case 1: 0.9*log-Weibull(0,1) + 0.1*log-Weibull(0,10)
## ------	
par(mfrow=c(1,1))
plot(seq(-6,5,0.1),
     dweibull(exp(seq(-6,5,0.1)), shape=1, scale=1)*
       exp(seq(-6,5,0.1)), 
     type="l", xlab="", ylab="", las=1)
lines(seq(-6,5,0.1), 
      dweibull(exp(seq(-6,5,0.1)), shape=1/10, scale=1)*
        exp(seq(-6,5,0.1)))
lines(seq(-6,5,0.1), 
      0.9*dweibull(exp(seq(-6,5,0.1)), shape=1, scale=1)*
        exp(seq(-6,5,0.1))+
      0.1*dweibull(exp(seq(-6,5,0.1)), shape=1/10, scale=1)*
        exp(seq(-6,5,0.1)), 
      lwd=2)
##
##	
par(mfrow=c(2,5), pty="s")
for(i in 1:10)
{
  x.seq <- seq(X[i,]%*%par.true[c(1:6)]-10*par.true[7],
 	       X[i,]%*%par.true[c(1:6)]+3*par.true[7], length=50)
  plot(x.seq, 
       dweibull(exp(x.seq), shape=1/par.true[7], 
	        scale=exp(X[i,]%*%par.true[c(1:6)]))*exp(x.seq), 
       type="l", xlab="", ylab="", main=paste("y",i), las=1)
  lines(x.seq, 
        dweibull(exp(x.seq), shape=1/ss.cont, 
	scale=exp(X[i,]%*%par.true[c(1:6)]))*exp(x.seq))	
  lines(x.seq, 
       (1-eps)*dweibull(exp(x.seq), shape=1/par.true[7], 
                        scale=exp(X[i,]%*%par.true[c(1:6)]))*
               exp(x.seq)+
       eps*dweibull(exp(x.seq), shape=1/ss.cont, 
                    scale=exp(X[i,]%*%par.true[c(1:6)]))*
               exp(x.seq), lwd=2)
}


## case 2: 0.9*log-Weibull(0,1) + 0.1*log-Gamma(nu=1,lambda=4) 
## ------    ==> 0.9*logExp(1) + 0.1*logExp(4)
##			
par(mfrow=c(1,1))
plot(seq(-6,5,0.1), 
     dweibull(exp(seq(-6,5,0.1)), shape=1, scale=1)*
       exp(seq(-6,5,0.1)), 
     type="l", xlab="",ylab="", las=1)
lines(seq(-6,5,0.1), 1/4*dgamma(exp(seq(-6,5,0.1))/4, shape=1)*
       exp(seq(-6,5,0.1)))
lines(seq(-6,5,0.1),
	0.9*dweibull(exp(seq(-6,5,0.1)), shape=1, scale=1)*
            exp(seq(-6,5,0.1))+
	0.1*1/4*dgamma(exp(seq(-6,5,0.1))/4, shape=1)*
                exp(seq(-6,5,0.1)), lwd=2)
##
##	
par(mfrow=c(2,5), pty="s")
for(i in 1:10)
{
  x.seq <- seq(X[i,]%*%par.true[c(1:6)]-10*par.true[7],
 	       X[i,]%*%par.true[c(1:6)]+3*par.true[7], leng=50)
  plot(x.seq, 
       dweibull(exp(x.seq), shape=1/par.true[7], 
	        scale=exp(X[i,]%*%par.true[c(1:6)]))*exp(x.seq), 
       type="l", xlab="", ylab="", main=paste("y",i), las=1)
  lines(x.seq, 1/(4*exp(X[i,]%*%par.true[c(1:6)]))*
	         dgamma(exp(x.seq)/(4*exp(X[i,]%*%par.true[c(1:6)])), 
	                shape=1)*exp(x.seq))
  lines(x.seq, 
       (1-eps)*dweibull(exp(x.seq), shape=1/par.true[7], 
	                scale=exp(X[i,]%*%par.true[c(1:6)]))*
               exp(x.seq)+
       eps*1/(4*exp(X[i,]%*%par.true[c(1:6)]))*
	       dgamma(exp(x.seq)/(4*exp(X[i,]%*%par.true[c(1:6)])), 
	              shape=1)*exp(x.seq))
}


## N O T E
## -------	
##  The densities we want to simulate from do not belong to the model
##  classes considered in the "marg" package of the "hoa" bundle.  We
##  consequently have to define a new "rsm" family according to what
##  done in the demonstration file `margdemo.R' that accompanies the
##  package.


## User-defined model classes  
## --------------------------
##  Note that the expression "contam.lW" is used for both cases 
##  considered above.

	  	
## case 1: 0.9*log-Weibull(0,1) + 0.1*log-Weibull(0,10)
## ------	

contam.lW <- function()
{
  make.family.rsm("contam.lW")
}

contam.lW.distributions <- structure(
                      .Data = list(
	g0 = function(y,...)
             { 
               dens <- 0.9*exp(-exp(y)+y) + 
                         0.1*exp(-exp(y/10)+y/10)/10 
               -log(dens)
             },
        g1 = function(y,...) 
             { 
               dens <- 0.9*exp(-exp(y)+y) + 
                       0.1*exp(-exp(y/10)+y/10)/10 
               dens.1 <- 0.9*exp(-exp(y)+y)*(1-exp(y)) + 
                         0.1*exp(-exp(y/10)+y/10)*
                               (1/10-exp(y/10)/10)/10
               -dens.1/dens
             }, 
        g2 = function(y,...)  
             {
               dens <- 0.9*exp(-exp(y)+y) + 
                       0.1*exp(-exp(y/10)+y/10)/10 
               dens.1 <- 0.9*exp(-exp(y)+y)*(1-exp(y)) + 
                         0.1*exp(-exp(y/10)+y/10)*
                               (1/10-exp(y/10)/10)/10
               dens.2 <- 0.9*(exp(-exp(y)+y)*
                                (1-exp(y))^2-exp(-exp(y)+2*y)) +
                         0.1*(exp(-exp(y/10)+y/10)*
                                (1/10-exp(y/10)/10)^2 -
                                 exp(-exp(y/10)+2*y/10)/10^2)/10
               -(dens.2*dens-(dens.1)^2)/dens^2
             } ),          
                      .Dim = c(3,1),
                      .Dimnames = list(c("g0","g1","g2"), 
                                       c("contam.lW")))

dcontam.lW <- function(x)
{
     0.9*dweibull(exp(x), shape=1, scale=1)*exp(x) +
	0.1*dweibull(exp(x), shape=1/10, scale=1)*exp(x)
}

rcontam.lW <- function(n)                          
{
    if( is.na(n) )
        return(NA)
    val1 <- rweibull(n, shape=1, scale=1)
    val2 <- rweibull(n, shape=1/10, scale=1)
    alpha <- runif(n)
    log(ifelse(alpha < 0.1, val2, val1))
}

pcontam.lW <- function(q)
{
     0.9*pweibull(exp(q), shape=1, scale=1) +
	0.1*pweibull(exp(q), shape=1/10, scale=1)
}


## case 2: 0.9*log-Weibull(0,1) + 0.1*log-Gamma(nu=1,lambda=4) 
## ------    ==> 0.9*logExp(1) + 0.1*logExp(4)
			
contam.lW <- function()
{
  make.family.rsm("contam.lW")
}

contam.lW.distributions <- structure(
                      .Data = list(
	g0 = function(y,...)
             { 
               dens <- 0.9*exp(-exp(y)+y) + 0.1*exp(-exp(y)*4+y)*4 
               -log(dens)
             },
        g1 = function(y,...) 
             { 
               dens <- 0.9*exp(-exp(y)+y) + 0.1*exp(-exp(y)*4+y)*4 
               dens.1 <- 0.9*exp(-exp(y)+y)*(1-exp(y)) + 
                         0.1*exp(-exp(y)*4+y)*(1-exp(y)*4)*4
               -dens.1/dens
             }, 
        g2 = function(y,...)  
             {
               dens <- 0.9*exp(-exp(y)+y) + 0.1*exp(-exp(y)*4+y)*4 
               dens.1 <- 0.9*exp(-exp(y)+y)*(1-exp(y)) + 
                         0.1*exp(-exp(y)*4+y)*(1-exp(y)*4)*4
               dens.2 <- 0.9*(exp(-exp(y)+y)*(1-exp(y))^2-
                              exp(-exp(y)+2*y)) +
                         0.1*(exp(-exp(y)*4+y)*(1-exp(y)*4)^2 -
                             exp(-exp(y)*4+2*y)*4)*4
               ret <- -(dens.2*dens-(dens.1)^2)/dens^2
	       ret[is.na(ret)] <- 0
	       ret
             } ),          
                      .Dim = c(3,1),
                      .Dimnames = list(c("g0","g1","g2"), 
                                       c("contam.lW")))

dcontam.lW <- function(x)
{
     0.9*dweibull(exp(x), shape=1, scale=1)*exp(x) +
	0.1*dexp(exp(x), rate=4)*exp(x)
}

rcontam.lW <- function(n)                          
{
    if( is.na(n) )
        return(NA)
    val1 <- rweibull(n, shape=1, scale=1)
    val2 <- rexp(n, rate=4)
    alpha <- runif(n)
    log(ifelse(alpha < 0.1, val2, val1))
}

pcontam.lW <- function(q)
{
     0.9*pweibull(exp(q), shape=1, scale=1) +
	0.1*dexp(exp(q), rate=4)
}


##  CHOOSE THE ERROR DISTRIBUTION YOU WANT TO WORK WITH !!!


##  We mantain the same simulation parameters apart from the error 
##  distribution.


## Unscaled covariance matrix  
## --------------------------
y.obs <- rcontam.lW(10)*sigma.true + X%*%beta.true 
##
attach(as.data.frame(X))
rsmObject <- rsm(y.obs~x1+x2+x3+x4+x5+x6-1, family=contam.lW,
                 control=glm.control(maxit=200))
beta.MLE <- rsmObject$coef
sigma.MLE <- rsmObject$disp
##
anc <- rsmObject$resid
##
var.MLE <- summary(rsmObject, corr=FALSE)$cov
corr.MLE <- summary(rsmObject)$corr
detach() 


## Laplace's approximation to the univariate marginal densities of 
##  the MLEs (regression coefficients + scale parameter)
## ---------------------------------------------------------------
##
attach(as.data.frame(X))
X.data <- make.sample.data(rsmObject)
X.data$coef <- beta.true    ## replace estimated values by true values
X.data$disp <- sigma.true
detach()
##
uni.dens <- list()
##
seq.tmp <- beta.MLE[1] + qnorm(c(0.025, 0.975))*sqrt(var.MLE[1,1])
seq.tmp <- seq(seq.tmp[1], seq.tmp[2], length=30)
uni.dens$b1 <- Laplace(c, X.data, seq.tmp, 1)$dens  
##
## and so on as above ...
##
seq.tmp <- sigma.MLE[1] + qnorm(c(0.025, 0.975))*sqrt(var.MLE[7,7])
seq.tmp <- seq(seq.tmp[1], seq.tmp[2], length=30)
uni.dens$ls  <- Laplace(s, X.data, seq(0.1,3,0.3))$dens
##
rm(seq.tmp)
##
par( mfrow=c(3,3), pty="s" ) 
sapply( uni.dens[1:6], function(x)        
          plot(x, type="l", xlab=paste("regression coefficient"), 
                            ylab="marginal density", cex=0.7, las=1) )
plot( uni.dens$ls, type="l", xlab="log-scale parameter", 
      ylab="marginal density", cex=0.7, las=1 )


## WARNING:
## -------
##  1) Beware that sometimes NAs are generated particularly in the 
##     tails of the distribution.  These are omitted when calculating 
##     the spline interpolation returned by the "Laplace" function.  
##  2) The spline interpolation can be "unstable" especially in the 
##     tails.


## Laplace's approximation to the bivariate marginal densities of the 
##   MLEs (regression coefficients + scale parameter)
## ------------------------------------------------------------------
##
bi.dens <- list()
##
seq1.tmp <- beta.MLE[1] + qnorm(c(0.025, 0.975))*sqrt(var.MLE[1,1])
seq1.tmp <- seq(seq1.tmp[1], seq1.tmp[2], length=15)
seq2.tmp <- beta.MLE[2] + qnorm(c(0.025, 0.975))*sqrt(var.MLE[2,2])
seq2.tmp <- seq(seq2.tmp[1], seq2.tmp[2], length=15)
bi.dens$b12 <- Laplace(cc, X.data, seq1.tmp, 1, seq2.tmp, 2)  
##
seq1.tmp <- beta.MLE[1] + qnorm(c(0.025, 0.975))*sqrt(var.MLE[1,1])
seq1.tmp <- seq(seq1.tmp[1], seq1.tmp[2], length=15)
seq2.tmp <- beta.MLE[3] + qnorm(c(0.025, 0.975))*sqrt(var.MLE[3,3])
seq2.tmp <- seq(seq2.tmp[1], seq2.tmp[2], length=15)
bi.dens$b13 <- Laplace(cc, X.data, seq1.tmp, 1, seq2.tmp, 3)
##
## and so on as above ...
##
seq1.tmp <- beta.MLE[1] + qnorm(c(0.025, 0.975))*sqrt(var.MLE[1,1])
seq1.tmp <- seq(seq1.tmp[1], seq1.tmp[2], length=15)
seq2.tmp <- beta.MLE[6] + qnorm(c(0.025, 0.975))*sqrt(var.MLE[6,6])
seq2.tmp <- seq(seq2.tmp[1], seq2.tmp[2], length=15)
bi.dens$b56 <- Laplace(cc, X.data, seq1.tmp, 5, seq2.tmp, 6)
##
seq1.tmp <- beta.MLE[1] + qnorm(c(0.025, 0.975))*sqrt(var.MLE[1,1])
seq1.tmp <- seq(seq1.tmp[1], seq1.tmp[2], length=15)
seq2.tmp <- sigma.MLE + qnorm(c(0.025, 0.975))*sqrt(var.MLE[7,7])
seq2.tmp <- seq(seq2.tmp[1], seq2.tmp[2], length=15)
bi.dens$b1ls <- Laplace(cs, X.data, seq1.tmp, 1, seq2.tmp)
##
## and so on as above ...
##
seq1.tmp <- beta.MLE[6] + qnorm(c(0.025, 0.975))*sqrt(var.MLE[6,6])
seq1.tmp <- seq(seq1.tmp[1], seq1.tmp[2], length=15)
seq2.tmp <- sigma.MLE + qnorm(c(0.025, 0.975))*sqrt(var.MLE[7,7])
seq2.tmp <- seq(seq2.tmp[1], seq2.tmp[2], length=15)
bi.dens$b6ls <- Laplace(cs, X.data, seq1.tmp, 6, seq2.tmp)
##
rm(seq1.tmp, seq2.tmp)
##
par( mfrow=c(4,6), pty="s" ) 
sapply( bi.dens, function(x) plot(x, nlevels=10, las=1) )


## WARNING: 
## -------
##  Beware that sometimes NAs or infinite values are generated 
##  particularly in the tails of the distribution.  These are set 
##  to 0 by the "Laplace" function. 


## Conditional sampling
## --------------------
loc.sim <- c(beta.true, -0.7)
var.sim <- 1.8*round(var.MLE, dig=4)   

cond.sample <- rsm.sample(X.data, R = 100000, ran.gen = xmpl.gen, 
                          mm = loc.sim, cov = var.sim)


## Observations corresponding to the sampled MLEs
## ----------------------------------------------
attach(cond.sample)
y.sim <- sim[,1:6] %*% t(X) + as.vector(sim[,7]) %*% t(as.vector(anc))
detach()


## Acceptance rate
## ---------------
attach(cond.sample)
a.rate <- 1-sum(sim[-1,1]==sim[-dim(sim)[1],1])/dim(sim)[1]
detach()


## Graphical inspection of the generated Markov chain and diagnostics
## ==================================================================
##
## Chain
## -----
par( mfrow=c(3,1), pty="m", ask=TRUE )                    
apply( cond.sample$sim, 2, function(x) plot(x, type="l", las=1) )
par( ask=FALSE )


## Simulated marginal density 
##   --> to compare with Laplace approximation
## -------------------------------------------
par( mfrow=c(3,3), pty="s" )  
apply(cond.sample$sim, 2, function(x) {
                      hist(x, prob=TRUE, nclass=20, xlab="", las=1) ; 
                      lines(density(x)) } )
hist(log(cond.sample$sim[,7]), prob=TRUE, nclass=20, xlab="",las=1)
lines(density(log(cond.sample$sim[,7])))
      

## and so on as above ...


## Distribution of pivots Q1 and Q2
## ================================
##
## Use the same code than above ...


## Distribution of R and R*
## ========================
##
## Use the same code than above ...





## ========================== THE END ================================
