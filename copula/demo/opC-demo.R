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


### Demo of the two-parameter outer power Clayton copula #######################

### setup ######################################################################

library(copula)
library(bbmle)
library(lattice)
library(grid)

## specify parameters
n <- 100 # sample size
d <- 10 # dimension
thetabase <- 1 # fix thetabase
tau <- 0.5 # => psi(t) = psi_thetabase(t^(1/theta)) with (thetabase,theta) = (1,4/3) (see below)

## adjustment for initial value
h <- c(0.4,0) # h_-, h_+


### functions ##################################################################

##' Initial interval for opC
##'
##' @title Initial interval for opC
##' @param U (n x d)-matrix of simulated data
##' @param h non-negative auxiliary parameter for computing initial intervals
##' @param method "etau" via sample version of Kendall's tau (may be slow)
##'               "dmle.G" via DMLE of Gumbel (may be inaccurate)
##' @return (2 x 2)-matrix containing the initial interval [1st row: lower,
##'         2nd row: upper; 2 parameters => 2 cols]
##' @author Marius Hofert
ii.opC <- function(U, h, method=c("etau","dmle.G")){
    stopifnot(h >= 0, length(h) >= 2)
    I <- matrix(, nrow=2, ncol=2, dimnames=list(c("lower", "upper"), c("thetabase", "theta")))
    ## estimate Kendall's tau
    method <- match.arg(method)
    tau.hat <- switch(method,
                      "etau" = { # uses sample version of tau, more accurate but slower
                          tau.hat.mat <- cor(U, method="kendall")
                          mean(tau.hat.mat[upper.tri(tau.hat.mat)])
                      },
                      "dmle.G" = { # uses DMLE for Gumbel to estimate tau
                          Z <- apply(U, 1, max)
                          theta.hat.G <- log(ncol(U))/(log(length(Z))-log(sum(-log(Z))))
                          copGumbel@tau(theta.hat.G)
                      },
                      stop("wrong method:", method))
    ## compute largest value of theta (for lower right endpoint of the inital interval)
    stopifnot(tau.hat > 0)
    tau.hat.hp <- min(tau.hat+h[2], 0.995)
    th.max <- 2*tau.hat.hp/(1-tau.hat.hp)
    I[2,1] <- th.max # largest value for theta (= thetabase)
    I[1,2] <- 1 # smallest value for beta (= theta)
    ## compute smallest theta (for lower left endpoint of the inital interval)
    tau.hat.hm <- max(tau.hat-h[1], 0.005) # tau=0.005 <=> theta=0.01
    th.min <- 2*tau.hat.hm/(1-tau.hat.hm)
    I[1,1] <- th.min
    ## compute largest beta (for upper left endpoint of the inital interval)
    b.max <- 2/((2+th.min)*(1-tau.hat.hp))
    I[2,2] <- b.max
    ## result
    I
}

##' -log-likelihood
##'
##' @title -log-likelihood
##' @param thetabase parameter of the base (=Clayton) generator
##' @param theta outer power parameter
##' @param u (n x d)-matrix of simulated data
##' @return -sum(log(density))
##' @author Marius Hofert
nlogl.opC <- function(thetabase, theta, u){
    if(!is.matrix(u)) u <- rbind(u)
    if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
    cop <- opower(copClayton, thetabase)
    -sum(cop@dacopula(u, theta=theta, log=TRUE))
}

## vectorized version
nlogl.opC. <- function(theta, u) nlogl.opC(theta[1], theta=theta[2], u=u)


### estimation #################################################################

## determine theta such that tau is matched (for given thetabase)
opC <- opower(copClayton, thetabase) # outer power Clayton copula
theta <- opC@iTau(tau) # choose theta such that Kendall's tau equals tau

## define the outer power Clayton copula to be sampled and estimated
cop <- onacopulaL(opC, list(theta, 1:d))

## sample
set.seed(1000)
U <- rnacopula(n, cop)

## plot
splom2(U, cex=0.4, pscales=0, main=paste("Sample of size",n,
                              "from an outer power Clayton copula"))

## initial interval and value
I <- ii.opC(U, h)
start <- colMeans(I)

## without profiling: optim with method="L-BFGS-B"
system.time(optim(par=start, method="L-BFGS-B",
                  fn=function(x) nlogl.opC(x[1], theta=x[2], u=U),
                  lower=c(I[1,1], I[1,2]), upper=c(I[2,1], I[2,2])))

## with profiling: via mle (uses optimizer="optim" with method="L-BFGS-B")
nLL <- function(thetabase, theta) nlogl.opC(thetabase, theta, u=U)
system.time(ml <- mle(nLL, method="L-BFGS-B",
                      start=list(thetabase=mean(I[,1]), theta=mean(I[,2])),
                      lower=c(thetabase=I[1,1], theta=I[1,2]),
                      upper=c(thetabase=I[2,1], theta=I[2,2])))
summary(ml)
str(ml@details)

## with profiling: via mle2 (uses optimizer="optim" with method="L-BFGS-B")
system.time(ml2 <- mle2(nlogl.opC, data=list(u=U), method="L-BFGS-B",
                        start=list(thetabase=mean(I[,1]), theta=mean(I[,2])),
                        lower=c(thetabase=I[1,1], theta=I[1,2]),
                        upper=c(thetabase=I[2,1], theta=I[2,2])))
summary(ml2)
str(ml2@details)


### plots ######################################################################

### profile likelihood plots ###################################################

prof <- profile(ml)
if(FALSE) { ## FIXME (?)
    ## maybe this helps: https://stat.ethz.ch/pipermail/r-help/2005-July/076003.html
    ci <- confint(prof)
    ci
    plot(prof)
}

prof2 <- profile(ml2)
(ci <- confint(prof2))
plot(prof2)


### -log-likelihood plots ######################################################

## for the plots, we use the standard mathematical notation (theta, beta)
## instead of (thetabase, theta)

## build grid
m <- 20 # number of grid points = number of intervals + 1
th <- seq(I[1,1], I[2,1], length.out=m) # grid points for thetabase
beta <- seq(I[1,2], I[2,2], length.out=m) # grid points for theta
grid <- expand.grid(theta=th, beta=beta) # grid
val.grid <- apply(grid, 1, nlogl.opC., u=U) # value of the -log-likelihood on the grid


### wireframe

## plot settings
true.theta <- c(thetabase, theta)
true.val <- c(true.theta, nlogl.opC.(true.theta, u=U)) # theoretical optimum
opt <- ml@coef # optimizer-optimum
opt.val <- c(opt, nlogl.opC.(opt, u=U)) # optimizer-optimum and its value
pts <- rbind(true.val, opt.val) # points to add to wireframe plot
title <- "-log-likelihood of an outer power Clayton copula" # title
sub <- substitute(italic(n) == N ~~~  italic(d)== D ~~~
                  tau == TAU ~~~ "#{eval}:" ~ NIT,
                  list(N=n, D=d, TAU= tau, NIT= ml@details$counts[[1]]))
## lattice bug:
sub <- as.expression(sub)

## wireframe
wireframe(val.grid~grid[,1]*grid[,2], screen=list(z=70, x=-55), zoom=0.95,
          xlab=expression(italic(theta)), ylab=expression(italic(beta)),
          zlab = list(as.expression(-log~L * group("(",list(theta,beta),")")), rot=90),
          main=title, sub=sub, pts=pts, scales=list(col=1, arrows=FALSE),
          par.settings=list(axis.line=list(col="transparent"),
          clip=list(panel="off")), zlim=c(min(val.grid, pts[,3]),
                                   max(val.grid, pts[,3])), aspect=1,
          panel.3d.wireframe = function(x,y,z,xlim,ylim,zlim,xlim.scaled,
          ylim.scaled,zlim.scaled,pts,...){
              panel.3dwire(x=x, y=y, z=z, xlim=xlim, ylim=ylim, zlim=zlim,
                           xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
                           zlim.scaled=zlim.scaled, alpha.regions=0.8, ...)
              panel.3dscatter(x=pts[,1], y=pts[,2], z=pts[,3],
                              xlim=xlim, ylim=ylim, zlim=zlim,
                              xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
                              zlim.scaled=zlim.scaled, type="p", col=c("red","blue"),
                              pch=c(3,4), lex=2, cex=1.4, .scale=TRUE, ...)
          }, key=list(x=0.64, y=1.01, points=list(pch=c(3,4), col=c("red","blue"),
                                      lwd=2, cex=1.4),
             text=list(c("True value", "Optimum of optimizer")), padding.text=3,
             cex=1, align=TRUE, transparent=TRUE))


### levelplot

## plot settings
xlim. <- c(min(grid[,1]),max(grid[,1]))
ylim. <- c(min(grid[,2]),max(grid[,2]))
xeps <- (xlim.[2] - xlim.[1]) * 0.04
yeps <- (ylim.[2] - ylim.[1]) * 0.04
cols <- adjustcolor(colorRampPalette(c("darkgreen", "green", "orange", "yellow"),
                                     space="Lab")(100), 0.8)

## levelplot
levelplot(val.grid~grid[,1]*grid[,2],
          par.settings=list(layout.heights=list(main=3, sub=2),
          regions=list(col=cols)),
          xlim=c(xlim.[1]-xeps, xlim.[2]+xeps),
          ylim=c(ylim.[1]-yeps, ylim.[2]+yeps),
          xlab=expression(italic(theta)), ylab=expression(italic(beta)),
          main=title, sub=sub, pts=pts, aspect=1,
          scales=list(alternating=c(1,1), tck=c(1,0)), contour=TRUE,
          panel=function(x, y, z, pts, ...){
              panel.levelplot(x=x, y=y, z=z, ...)
              grid.points(x=pts[1,1], y=pts[1,2], pch=3,
                          gp=gpar(lwd=2, col="red")) # + true value
              grid.points(x=pts[2,1], y=pts[2,2], pch=4,
                          gp=gpar(lwd=2, col="blue")) # x optimum
          }, key=list(x=0.18, y=1.08, points=list(pch=c(3,4), col=c("red","blue"),
                                      lwd=2, cex=1.4),
             columns=2, text=list(c("True value", "Optimum of optimizer")),
             align=TRUE, transparent=TRUE))
