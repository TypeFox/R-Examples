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
set.seed(1)
n <- 2000
N <- 100 # number of variates used for the Kolmogorov-Smirnov test

## maximal deviation for deciding if sample versions of Kendall's tau are
## "close enough" to population versions; NB: depends on 'n'
eps.tau <- 0.06

doPlots <- (Sys.getenv("USER") == "maechler")


### 3d check functions #########################################################

##' correlation check function and run time measuring
##'
##' @title Check correlation matrix and measure run times of 3d fully nested Archimedean copulas
##' @param n number of variates to be drawn
##' @param th0 theta0
##' @param th1 theta1
##' @param cop acopula
##' @return a list containing run times for V0 and V01 and Kendall's taus
##' @author Marius Hofert, Martin Maechler
corCheck <- function(n, th0,th1, cop) {
    mat <- matrix(0,nrow = n,ncol = 3)
    V0time <- system.time(V0 <- cop@V0(n,th0))
    V01time <- system.time(V01 <- cop@V01(V0,th0,th1))
    mat <- cbind(runif(n),
                 exp(-V0*cop@iPsi(cop@psi(rexp(n)/V01,th1),th0)),
                 exp(-V0*cop@iPsi(cop@psi(rexp(n)/V01,th1),th0)))
    mat[,] <- cop@psi(-log(mat[,])/V0,th0)
    list(V0time,V01time, name = cop@name, cor = cor(mat, method = "kendall"))
}

##' create output
prt.tau.diff <- function(c1, c2){
    stopifnot(is.matrix(c1))
    delta.c <- 1000 * abs(c1 - c2)
    cat(sprintf("Max & Mean distance * 1000 to true pairwise Kendall's taus: %7.1f %7.1f\n",
                max(delta.c), mean(delta.c)))
    invisible()
}

##' create output
corCheckout <- function(x, trCorr, famName = x$name) {
    cat(sprintf("Time [ms] V0  for '%s': %7.1f\n", famName, 1000*x[[1]][1]))
    cat(sprintf("Time [ms] V01 for '%s': %7.1f\n", famName, 1000*x[[2]][1]))
    prt.tau.diff(x[["cor"]], trCorr) ; cat("\n")
}

##' Function implementing the chi^2 test
##'
##' @title The Chi-square test
##' @param n [integer] sample size
##' @param N [integer] number of replications
##' @param cop outer_nacopula to generate from
##' @param nInt positive integer: the number of intervals used for each grid
##' dimension
##' @return an "chiSqChk_cop" object; just a list(...) and a print method
##' @author Marius Hofert, Martin Maechler
chiSq_check_cop <- function(n,N,cop,nInt, verbose = interactive()){
    copName <- deparse(substitute(cop)) # copula name
    d <- dim(cop) # copula dimension
    stopifnot(is.numeric(d), d >= 1, is.numeric(nInt), nInt >= 1)
    pts <- (1:nInt)/nInt # (upper) division points [lower = upper - h
                                        # = {0, ..., 1-h}; h=1/nInt]
    mygrid <- do.call("expand.grid", rep.int(list(pts), d)) # build grid
    m <- nInt^d # == grid length == nrow(mygrid)
    v.cube <- nInt^(0:(d-1))

    ## build a function that returns the number of the cube in which each row
    ## of U falls
    cube <- function(U, pts) {
        di <- dim(U)
        intervals <- array(cut(U, breaks = pts, include.lowest = TRUE,
                               labels = FALSE), # find "interval number" for
					# each component of U; these numbers
					# are in {NA,1,2,3,...,nInt-1}
                           dim = di)
        intervals[is.na(intervals)] <- 0 # NAs correspond to smallest interval
        as.vector(intervals %*% v.cube) + 1
    }

    ## determine the expected number of observations in each cube
    prob_up <- function(u) {
        ## probability mass in cube with upper corner 'u'
        mesh <- 1/nInt
        l <- u - mesh
        prob(cop, l, u)
    }
    masscube <- apply(mygrid, 1, prob_up)
    E_nobs <- n * masscube # expected number of observations in each cube

    ## now simulate data, count observations in each cube, and compute test
    ## statistic
    k <- 0
    CPU <- system.time({
        T <- replicate(N, {
            if(verbose) cat(sprintf("%2d%1s",{k <<- k+1}, if(k %% 20)""
            else "\n"))
            U <- rnacopula(n,cop) # generate data
            cubenumbers <- cube(U,pts) # for each row vector of U, find the
                                        # number of the cube in which the
					# vector falls
            nobs <- tabulate(cubenumbers, nbins = m) # number of observations
                                        # in each cube
            sum((nobs - E_nobs)^2 / E_nobs) # chi^2 test statistic
        }); if(verbose) cat("done\n")
    })[1]

    structure(class = "chiSqChk_cop",
              list( ## compute the result of the Kolmogorov-Smirnov test based
                   ## on the N realizations of the chi^2 test statistics:
                   ks = ks.test(T, "pchisq", df = m-1),
                   T = T, CPU = CPU,
                   n=n, N=N, copName = copName, m = m,
                   ## percentage of cubes that fulfill the rule of thumb:
                   percentrot = (sum(E_nobs >= 5)/m)*100
                   ))
}

##' a print method for this class:
print.chiSqChk_cop <- function(x, ...) {
    stopifnot(is.list(x), all(c("ks","T","CPU","n","N") %in% names(x)),
	      is.numeric(pv <- x$ks[[2]]))
    cat(sprintf("%s (3d)NAcopula (n=%d):\n  %s (N=%d): %s\n  ",
		x$copName, x$n,
                "P-value of the chi-square test", x$N, format.pval(pv)),
        sprintf("Percentage fulfilling chi^2 rule of thumb: %4.1f\n",
                x$percentrot),
        sprintf("Time (user) needed = c(N,n; cop) = %8.1f [ms]\n",
		1000 * x$CPU), sep="")
    if(pv < 0.05) {
	if(pv < 0.01)
	    cat("\n*************** P-value < 0.01 <<<<<<<<<<<<<<<<<<<<<<<<\n",
		"\n*************** ============== <<<<<<<<<<<<<<<<<<<<<<<<\n\n")
	else cat("\n*** > > > P-value < 0.05 <<<<<<<<<<<<<<<<<<<<<<<<\n\n")
	stopifnot(pv > 0.001)
    }
    invisible(x)
}

##' compute the probability to fall in a cube with
##' lower point l and upper point u for d=3
probin3dcube <- function(cop,l,u) {
    pCopula(u, cop)+
        - pCopula(c(l[1],u[2],u[3]), cop)+
            - pCopula(c(u[1],l[2],u[3]), cop)+
                - pCopula(c(u[1],u[2],l[3]), cop)+
                    + pCopula(c(l[1],l[2],u[3]), cop)+
                        + pCopula(c(l[1],u[2],l[3]), cop)+
                            + pCopula(c(u[1],l[2],l[3]), cop)+
                                - pCopula(l, cop)
}


### 3d examples ################################################################

### AMH ########################################################################

theta0 <- 0.7135 # tau_{12}=tau_{13}=0.2, tau_{23}=0.3
theta1 <- 0.9430

## check 1
corCheckAMH <- corCheck(n,theta0,theta1,copAMH)
trCorr <- rbind(c(1,0.2,0.2),
                c(0.2,1,0.3),
                c(0.2,0.3,1))
corCheckout(corCheckAMH,trCorr)
stopifnot(max(abs(corCheckAMH[["cor"]]-trCorr)) < eps.tau)

## check 2
AMH3d <-
    new("outer_nacopula", copula = setTheta(copAMH, theta0),
        comp = as.integer( 1 ),
        childCops = list(new("nacopula",
        copula = setTheta(copAMH, theta1),
        comp = as.integer(c(2,3)))) # no childCops
        )

## constructor forms of the above:
rr <- onacopula("A",   C(0.7135, 1, list(C(0.943, 2:3, NULL))))
r0 <- onacopula("A",   C(0.7135, 1,      C(0.943, 2:3, NULL)))
r1 <- onacopula("A",   C(0.7135, 1,      C(0.943, 2:3, )))
r2 <- onacopula("AMH", C(0.7135, 1,      C(0.943, 2:3  )))
stopifnot(identical(AMH3d, rr), identical(rr, r0),
          identical(r0, r1), identical(r1, r2))

## check
(chkAMH <- chiSq_check_cop(n,N,AMH3d,5))

## check probability
l <- c(.1, .05, .3)
u <- c(.4, .7,  .6)
stopifnot(all.equal(print(  prob(AMH3d,l,u)),
		    probin3dcube(AMH3d,l,u), tolerance=1e-14))

### Clayton ####################################################################

theta0 <- 0.5 # tau_{12}=tau_{13}=0.2, tau_{23}=0.5
theta1 <- 2

## check 1
corCheckClayton <- corCheck(n,theta0,theta1,copClayton)
trCorr <- rbind(c(1,0.2,0.2),
                c(0.2,1,0.5),
                c(0.2,0.5,1))
corCheckout(corCheckClayton,trCorr)
stopifnot(max(abs(corCheckClayton[["cor"]]-trCorr)) < eps.tau)

## check 2
Clayton3d <- onacopula("Clayton", C(theta0, 1, C(theta1, 2:3)))
(chkClayton <- chiSq_check_cop(512,100,Clayton3d,5))

## check probability
stopifnot(all.equal(print(  prob(Clayton3d,l,u)),
		    probin3dcube(Clayton3d,l,u), tolerance=1e-14))

### Frank ######################################################################

theta0 <- 1.8609 # tau_{12}=tau_{13}=0.2, tau_{23}=0.5
theta1 <- 5.7363

## check 1
corCheckFrank <- corCheck(n,theta0,theta1,copFrank)
corCheckout(corCheckFrank,trCorr)
stopifnot(max(abs(corCheckFrank[["cor"]]-trCorr)) < eps.tau)

## check 2
Frank3d <- onacopula("F", C(theta0, 1, C(theta1, 2:3)))
(chkFrank <- chiSq_check_cop(n,N,Frank3d,5))

## check probability
stopifnot(all.equal(print(  prob(Frank3d,l,u)),
		    probin3dcube(Frank3d,l,u), tolerance=1e-14))

### Gumbel #####################################################################

theta0 <- 1.25
theta1 <- 2    #--> tau_{12}=tau_{13}=0.2, tau_{23}=0.5
trCorr <- rbind(c(1,0.2,0.2),
                c(0.2,1,0.5),
                c(0.2,0.5,1))
## check 1
corCheckGumbel <- corCheck(n,theta0,theta1,copGumbel)
corCheckout(corCheckGumbel,trCorr)
stopifnot(max(abs(corCheckGumbel[["cor"]]-trCorr)) < eps.tau)

## check 2
Gumbel3d <- onacopula("Gumbel", C(theta0, 1, C(theta1, 2:3)))
(chkGumbel <- chiSq_check_cop(n,N,Gumbel3d,5))

## check probability
stopifnot(all.equal(print(  prob(Gumbel3d,l,u)),
		    probin3dcube(Gumbel3d,l,u), tolerance=1e-14))

### Joe ########################################################################

theta0 <- 1.4438#tau_{12}=tau_{13}=0.2, tau_{23}=0.5
theta1 <- 2.8562

## check 1
corCheckJoe <- corCheck(n,theta0,theta1,copJoe)
corCheckout(corCheckJoe,trCorr)
stopifnot(max(abs(corCheckJoe[["cor"]]-trCorr)) < eps.tau)

## check 2
Joe3d <- onacopula("J", C(theta0, 1, C(theta1, 2:3)))
(chkJoe <- chiSq_check_cop(n,N,Joe3d,5))

## check probability
stopifnot(all.equal(print(  prob(Joe3d,l,u)),
		    probin3dcube(Joe3d,l,u), tolerance=1e-14))

### Examples that check pnacopula() and rnacopula() ############################

## generate output for the examples
prt.stats <- function(c1,c2, rt) {
    cat("Time [ms] for generating", n,
        "vectors of variates:  ", round(1000*rt[1],1), "\n")
    prt.tau.diff(c1, c2) ; cat("\n")
}

### 3d Ali-Mikhail-Haq copula example ##########################################

c3 <- onacopula("A", C(0.7135, 1, list(C(0.943, 2:3))))

## basic check
d <- dim(c3)
stopifnot(d == 3,
	  allComp(c3) == 1:3,
	  allComp(c3@childCops[[1]]) == 2:3)

## test pCopula(., <nacopula>)  {was pnacopula()}
u <- c(.3, .4, .5)
## with function:
v <- pCopula(u, c3)
## by hand
psi <- function(t,theta) { (1-theta)/(exp(t)-theta) }
iPsi <- function(t,theta) { log((1-theta*(1-t))/t) }
th0 <- 0.7135
th1 <- 0.9430
level1 <- psi(iPsi(u[2],th1) + iPsi(u[3],th1), th1)
level0 <- psi(iPsi(u[1],th0) + iPsi(level1, th0), th0)
stopifnot(all.equal(v, level0, tolerance = 1e-14))

## test rnacopula()
rt <- system.time(rC3 <- rnacopula(n,c3))
C3 <- cor(rC3,method = "kendall")
trCorr <- rbind(c(1,0.2,0.2),
                c(0.2,1,0.3),
                c(0.2,0.3,1)) # tau_{12}=tau_{13}=0.2, tau_{23}=0.3
stopifnot(is.numeric(rC3), is.matrix(rC3),
	  dim(rC3) == c(n, 3),max(abs(C3-trCorr)) < eps.tau)
prt.stats(C3,trCorr,rt)
if(doPlots) {
    stopifnot(require("KernSmooth"))## for smoothScatter()
    pairs(rC3, panel = function(...) { par(new = TRUE); smoothScatter(...) })
}

### 2d Clayton copula example ##################################################

c2 <- onacopula("Clayton", C(0.5, c(1,2))) # no childCops
## or simply c2 <- onacopula("Clayton", C(0.5, 1:2))

## basic check
d <- dim(c2)
stopifnot(d ==  2,
	  allComp(c2) == 1:2)

## test pCopula()
v <- pCopula(c(.3, .4), c2)
stopifnot(all.equal(v,
                    local( { u1 <- .3; u2 <- .4
                             (u1^(-1/2)+u2^(-1/2)-1)^(-2) }),
                    tolerance = 1e-14))

## test rnacopula()
racopula <- copula:::racopula
set.seed(17) ; rt  <- system.time(rC2 <- rnacopula(n,c2))
set.seed(17) ; rt. <- system.time(rc2 <- racopula (n, c2@copula, d=2))
stopifnot(identical(rC2, rc2))

C2 <- cor(rC2,method = "kendall")
trCorr <- rbind(c(1,0.2),
                c(0.2,1)) # tau_{12}=0.2
stopifnot(is.numeric(rC2), is.matrix(rC2),
	  dim(rC2) == c(n, 2), max(abs(C2-trCorr)) < eps.tau)
prt.stats(C2,trCorr,rt)
if(doPlots)
    smoothScatter(rC2)

### 3d Clayton copula example ##################################################

c3 <- onacopula("C", C(0.5, 1, C(2., c(2,3))))

## basic check
d <- dim(c3)
stopifnot(d == 3,
	  allComp(c3) == 1:3,
	  allComp(c3@childCops[[1]]) == 2:3)

## test pCopula()
v <- pCopula(c(.3, .4, .5), c3)
stopifnot(all.equal(v,
                    local( { u1 <- .3; u2 <- .4; u3 <- .5
                             1/((1/u2^2 +1/u3^2 -1)^(1/4) -1 +1/sqrt(u1))^2 }),
                    tolerance = 1e-14))

## test rnacopula()
rt <- system.time(rC3 <- rnacopula(n,c3))
C3 <- cor(rC3,method = "kendall")
trCorr <- matrix(c(1,0.2,0.2,0.2,1,0.5,0.2,0.5,1),nrow = 3,byrow = TRUE)
                                        # tau_{12}=tau_{13}=0.2, tau_{23}=0.5
stopifnot(is.numeric(rC3), is.matrix(rC3),
	  dim(rC3) == c(n, 3),max(abs(C3-trCorr)) < eps.tau)
prt.stats(C3,trCorr,rt)

if(doPlots)
    pairs(rC3, panel = function(...) { par(new = TRUE); smoothScatter(...) })

### 9d Clayton copula example ##################################################

c9 <- onacopula("Clayton", C(0.5, c(3,6,1),
			     C(2., c(9,2,7,5),
			       C(3., c(8,4)))))
c9Lis <- list(0.5, c(3,6,1),
              list(list(2., c(9,2,7,5),
                        list(list(3., c(8,4))))))
## consistency onacopula()  <->  onacopulaL() :
stopifnot(identical(c9, onacopulaL("Clayton", c9Lis)))


## basic check
d <- dim(c9)
stopifnot(d == 9,
          allComp(c9) == c(3,6,1,9,2,7,5,8,4),
          allComp(c9@childCops[[1]]) == c(9,2,7,5,8,4),
          allComp(c9@childCops[[1]]@childCops[[1]]) == c(8,4))

## test pCopula()
u <- seq(0.1,0.9,by = 0.1)
v <- pCopula(u, c9)
## by hand
psi <- function(t,theta) { (1+t)^(-1/theta) }
iPsi <- function(t,theta) { t^(-theta) - 1 }
th0 <- 0.5
th1 <- 2
th2 <- 3
level2 <- psi(iPsi(u[8],th2) + iPsi(u[4],th2), th2)
level1 <- psi(iPsi(u[9],th1)+
              iPsi(u[2],th1)+
              iPsi(u[7],th1)+
              iPsi(u[5],th1) +
              iPsi(level2, th1), th1)
level0 <- psi(iPsi(u[3],th0)+
              iPsi(u[6],th0)+
              iPsi(u[1],th0)+
              iPsi(level1, th0), th0)
stopifnot(all.equal(v, level0, tolerance = 1e-14))

## test rnacopula()
rt <- system.time(rC9 <- rnacopula(n,c9))
C9 <- cor(rC9,method = "kendall")

## Theoretical values:
## (11,12,13,14,15,16,17,18,19)=(1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2)
##    (22,23,24,25,26,27,28,29)=(1,0.2,0.5,0.5,0.2,0.5,0.5,0.5)
##       (33,34,35,36,37,38,39)=(1,0.2,0.2,0.2,0.2,0.2,0.2)
##          (44,45,46,47,48,49)=(1,0.5,0.2,0.5,0.6,0.5)
##             (55,56,57,58,59)=(1,0.2,0.5,0.5,0.5)
##                (66,67,68,69)=(1,0.2,0.2,0.2)
##                   (77,78,79)=(1,0.5,0.5)
##                      (88,89)=(1,0.5)

C9.true <- rbind(c(1. ,rep(0.2,8)),
                 c(0.2,1. ,0.2,0.5,0.5,0.2, rep(0.5,3)),
                 c(0.2,0.2,1. , rep(0.2,6)),
                 c(0.2,0.5,0.2,1. ,0.5,0.2,0.5,0.6,0.5),
                 c(0.2,0.5,0.2,0.5,1. ,0.2, rep(0.5,3)),
                 c(rep(0.2,5),         1. , rep(0.2,3)),
                 c(0.2,0.5,0.2,0.5,0.5,0.2,1. ,0.5,0.5),
                 c(0.2,0.5,0.2,0.6,0.5,0.2,0.5,1. ,0.5),
                 c(0.2,0.5,0.2,0.5,0.5,0.2,0.5,0.5,1. ))
stopifnot(dim(rC9) == c(n, 9),
          max(abs(C9-C9.true)) < eps.tau)
prt.stats(C9,C9.true,rt)
if(doPlots && dev.interactive(orNone=TRUE)) # "large"
    pairs(rC9, gap = .1, pch = 20, cex = 0.2, col = rgb(.2,.1,.7, alpha = .5),
          main = paste(n," vectors of a ",d,
          "-dimensional nested Clayton copula",sep = ""))

### 25d Clayton ==============================================
c25 <- onacopula("Clayton", C(0.5, 17,
                               list(C(2,   20:18),
                                    C(2.5, c(25,23, 8:12)),
                                    C(2.25,c(24,21), C(4, 3:5)),
                                    C(1.5, c(22,15:16), C(1.7, 1:2)),
                                    C(3,   c(6:7, 14:13)))))
stopifnot(identical(sort(allComp(c25)), 1:25))
c25

### 125d Clayton copula example ################################################

c125 <- onacopula("Clayton", C(0.5, , # no direct components
                               list(C(2,  1:10),
                                    C(3, 11:40),
                                    C(2, 41:60),
                                    C(2, 61:85),
                                    C(3, 86:105),
                                    C(2,106:125))))
c125Lis <- list(0.5, integer(0), # <- could use NULL
                list(list(2,  1:10),
                     list(3, 11:40),
                     list(2, 41:60),
                     list(2, 61:85),
                     list(3, 86:105),
                     list(2,106:125)))
## consistency onacopula()  <->  onacopulaL() :
stopifnot(identical(c125, onacopulaL("C", c125Lis)))



## basic check
d <- dim(c125)
stopifnot(d == 125,
	  allComp(c125) == 1:125,
	  allComp(c125@childCops[[1]]) == 1:10,
	  allComp(c125@childCops[[2]]) == 11:40,
	  allComp(c125@childCops[[3]]) == 41:60,
	  allComp(c125@childCops[[4]]) == 61:85,
	  allComp(c125@childCops[[5]]) == 86:105,
	  allComp(c125@childCops[[6]]) == 106:125
	  )

## test rnacopula()
rt <- system.time(rC125 <- rnacopula(n,c125))
stopifnot(is.numeric(rC125), is.matrix(rC125), dim(rC125) == c(n, 125))
cat("Time elapsed for generating ",n," vectors of variates:\n",sep = "")
rt

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
