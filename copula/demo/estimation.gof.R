## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


require(copula)
options(warn=1)

### Estimation and goodness-of-fit #############################################

source(system.file("Rsource", "estim-gof-fn.R", package="copula"))
## ../inst/Rsource/estim-gof-fn.R   --> estimation.gof() etc

### setup

## Use all available estimation and GoF methods:
(estMeth <- eval(formals(enacopula)$method))
(gofTraf <- eval(formals(gofPB)$trafo.method)[-1])
(gofMeth <- eval(formals(gofCopula)$method))

set.seed(1) # set seed

n <- 128 # sample size
d <- 5 # dimension
tau <- 0.25 # Kendall's tau

### apply all procedures (to data from AMH) ####################################

simFamily <- "AMH"
cop <- getAcop(simFamily)
theta <- cop@iTau(tau) # true parameter

## start the loop
cat("\n### data from ",simFamily," (n = ",n,", d = ",d,", theta = ",
    format(theta),", tau = ", format(tau),") ###\n\n",sep="")

## As "smle" is the slowest by far, we need to leave it away here:
estM.1 <- estMeth[estMeth != "smle"]
## Hmm, but actually, we currently only recommend to use GOF for the MLE,
## and that saves CPU time, too :
RR <-
      sapply(gofTraf, simplify="array", function(gt) {
	sapply(gofMeth, simplify="array", function(gm)
	       estimation.gof(n, d=d, simFamily=simFamily, tau=tau,
			      n.bootstrap= 16, # <-- as some methods are time consuming,
			      ## please choose a larger number here, e.g., 1000,
			      ## for particular methods.
			      include.K = TRUE, estim.method.enacopula="mle",
                              estim.method.fitCopula="ml",
			      gof.trafo = gt, gof.method = gm))
    })

str(RR)
dimnames(RR)

## Now print RR
options(digits=5)

## *Not* the times here...
##RR[,c("theta_hat","tau_hat","P_value","< 0.05"),,,]
  RR[,c("theta_hat","tau_hat","P_value","< 0.05"),,]

## ... but here
## apply(RR[,c("timeEstim","timeGoF"),,,], c(3,1,2,4), mean)
   apply(RR[,c("timeEstim","timeGoF"),,],  c(3,1,2), mean)


### MLE estimation by hand (for debugging purposes) ############################

## generate the data
simFamily <- getAcop("AMH") # choose your desired family
cop <- onacopulaL(simFamily, list(theta <- simFamily@iTau(tau),1:d))
u <- rnacopula(n,cop)

## estimate the copula
copFamily <- "Joe" # family to be estimated
cop.hat <- onacopulaL(copFamily,list(NA,1:d))
mLogL <- function(theta,cop.hat,u){
    cop.hat@copula@theta <- theta
    -sum(dCopula(u, cop.hat, log=TRUE))
}
(est <- optimize(mLogL, interval=initOpt(copFamily), cop.hat=cop.hat,
                 u=u))

## evaluate the density at u for a specified parameter
theta <- 14
cop.hat@copula@theta <- theta
(log.dens <- dCopula(u, cop.hat, log=TRUE))
-sum(log.dens)

### Plots ######################################################################

if(!dev.interactive()) # e.g. when run as part of R CMD check
    pdf("demo_est-gof.pdf")
if(!exists("doPlot")) doPlot <- TRUE

### setup for plots

u <- (0:256)/256 # evaluation points

cols <- c("black","orange3","red3","darkgreen","blue") # not very light ones
labs <- c("AMH","Clayton","Frank","Gumbel","Joe")

### plots of the densities of the diagonals

d <- 5
th <- c(0.7135001, 0.5, 1.860884, 1.25, 1.25)
dDmat <- cbind(dDiag.A=dDiag(u,cop=onacopulaL("AMH", list(th[1], 1:d))),
               dDiag.C=dDiag(u,cop=onacopulaL("Clayton", list(th[2], 1:d))),
               dDiag.F=dDiag(u,cop=onacopulaL("Frank", list(th[3], 1:d))),
               dDiag.G=dDiag(u,cop=onacopulaL("Gumbel", list(th[4], 1:d))),
               dDiag.J=dDiag(u,cop=onacopulaL("Joe", list(th[5], 1:d))))

if(doPlot) {
    matplot(u, dDmat, type="l", col=cols, xlab="t", ylab="dDiag(t)")
    legend("bottomright", legend=labs, lty=1:5, col=cols, bty="n")
    ## and in log-log scale:
    matplot(u, dDmat, type="l", col=cols, xlab="t",
            log="xy", main="dDiag(t) [log-log scale]")
    legend("bottomright", legend=labs, lty=1:5, col=cols, bty="n")
}

### plots of the Kendall distribution functions

d <- 10
Kmat <- cbind(K.A=K(u,setTheta(copAMH, th[1]),d),
              K.C=K(u,setTheta(copClayton, th[2]),d),
              K.F=K(u,setTheta(copFrank, th[3]),d),
              K.G=K(u,setTheta(copGumbel, th[4]),d),
              K.J=K(u,setTheta(copJoe, th[5]),d))
head(mm <- cbind(t=u, Kmat))
tail(mm)

dK <- apply(Kmat, 2, diff)
summary(dK)
## NOTE:  AMH and Clayton have some (very slightly) negative values
## ----    <==>  K() is not increasing there (near 1)
## MM: this is  "unavoidable" because of the numerics behind ...

if(doPlot) {
    matplot(u, Kmat, type="l", col=cols, xlab="t", ylab="K(u)")
    legend("bottomright", legend=labs, lty=1:5, col=cols, bty="n")
    ## and in log-log scale:
    matplot(u, Kmat, type="l", col=cols, xlab="t",
            log="xy", main="K(u) [log-log scale]")
    legend("bottomright", legend=labs, lty=1:5, col=cols, bty="n")
}
