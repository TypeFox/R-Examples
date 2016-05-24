HW.original <- function (data, cod, Nsim=500) {
 A <- as.double(0.0)
 B <- as.double(0.0)
 # FOR UNBIASED ESTIMATES (OF THE LAMBDA'S) SET A=B=ZERO. OTHERWISE,
 # PLOTTING-POSITION ESTIMATORS ARE USED, BASED ON THE PLOTTING POSITION
 # (J+A)/(N+B) FOR THE J'TH SMALLEST OF N OBSERVATIONS. FOR EXAMPLE,
 # A=-0.35D0 AND B=0.0D0 YIELDS THE ESTIMATORS RECOMMENDED BY
 # HOSKING ET AL. (1985, TECHNOMETRICS) FOR THE GEV DISTRIBUTION.

 reglmom <- REGLMR(data, cod)$RMOM
 parkappa <- PEL(reglmom, distr="KAP")
 out <- REGTST(data, cod, A=A, B=B, Nsim)
 out$RPARA <- parkappa 
 
 class(out) <- c("HWorig")
 return(out)
}


#.First.lib <- function(libname, pkgname) {
# library.dynam("nsRFA", pkgname, libname)
#}
.onLoad <- function(libname, pkgname){
 # do whatever needs to be done when the package is loaded
 # some people use it to bombard users with 
 # messages using 
 #packageStartupMessage( "my package is so cool" )
 #packageStartupMessage( "so I will print these lines each time you load it")
 library.dynam("nsRFA", pkgname, libname)
}


# ------------------------------------------------------------------- #

print.HWorig <- function(x, ...) {
 # x = object of classHWorig

 cat("\n\nRESULTS OBTAINED WITH THE FORTRAN ROUTINE OF HOSKING\n\n")

 names <- x$NAMES
 len <- x$LEN
 xmom <- x$XMOM
 d <- x$D
 rlm <- x$RMOM[2:4]
 names(rlm) <- c("L-CV", "L-SKEW", "L-KURT")
 if (length(d) > 1) {
  tabella01 <- data.frame(cbind(len,t(xmom)[,2:4],d), row.names=names)
 }
 else {
  tabella01 <- c(len,t(xmom)[,2:4],d)
 }
 names(tabella01) <- c("n", "L-CV", "L-SKEW", "L-KURT", "D(I)")
 print(tabella01)

 cat("\n\nREGIONAL L-MOMENTS:\n")
 print(rlm)

 cat("\n\nPARAMETERS OF REGIONAL KAPPA DISTRIBUTION:\n")
 print(x$RPARA)


 cat("\n\n***** HETEROGENEITY MEASURES *****\n")
 cat("\nNUMBER OF SIMULATIONS = ", x$NSIM, "\n", sep="")

 cat("\nOBSERVED S.D. OF GROUP L-CV     = ", x$VOBS[1], "\n", sep="")
 cat("SIM. MEAN OF S.D. OF GROUP L-CV = ", x$VBAR[1], "\n", sep="")
 cat("SIM. S.D. OF S.D. OF GROUP L-CV = ", x$VSD[1], "\n", sep="")
 cat("STANDARDIZED TEST VALUE H(1)    = ", x$H[1], "\n", sep="")

 cat("\nOBSERVED AVE. OF L-CV / L-SKEW DISTANCE  = ", x$VOBS[2], "\n", sep="")
 cat("SIM. MEAN OF AVE. L-CV / L-SKEW DISTANCE = ", x$VBAR[2], "\n", sep="")
 cat("SIM. S.D. OF AVE. L-CV / L-SKEW DISTANCE = ", x$VSD[2], "\n", sep="")
 cat("STANDARDIZED TEST VALUE H(2)             = ", x$H[2], "\n", sep="")

 cat("\nOBSERVED AVE. OF L-SKEW/L-KURT DISTANCE  = ", x$VOBS[3], "\n", sep="")
 cat("SIM. MEAN OF AVE. L-SKEW/L-KURT DISTANCE = ", x$VBAR[3], "\n", sep="")
 cat("SIM. S.D. OF AVE. L-SKEW/L-KURT DISTANCE = ", x$VSD[3], "\n", sep="")
 cat("STANDARDIZED TEST VALUE H(3)             = ", x$H[3], "\n", sep="")


 cat("\n\n***** GOODNESS-OF-FIT MEASURES *****\n")
 cat("\nNUMBER OF SIMULATIONS = ", x$NSIM, "\n", sep="")

 cat("\nGEN. LOGISTIC        Z VALUE = ", x$Z[1], "\n", sep="")
 cat("GEN. EXTREME VALUE   Z VALUE = ", x$Z[2], "\n", sep="")
 cat("GEN. NORMAL          Z VALUE = ", x$Z[3], "\n", sep="")
 cat("PEARSON TYPE III     Z VALUE = ", x$Z[4], "\n", sep="")
 cat("GEN. PARETO          Z VALUE = ", x$Z[5], "\n", sep="")

 cat("\nPARAMETER ESTIMATES FOR DISTRIBUTIONS\n")
 cat("GEN. LOGISTIC        =", x$PARA[1:3,1], "\n", sep=" ")
 cat("GEN. EXTREME VALUE   =", x$PARA[1:3,2], "\n", sep=" ")
 cat("GEN. NORMAL          =", x$PARA[1:3,3], "\n", sep=" ")
 cat("PEARSON TYPE III     =", x$PARA[1:3,4], "\n", sep=" ")
 cat("GEN. PARETO          =", x$PARA[1:3,5], "\n", sep=" ")
 cat("WAKEBY               =", x$PARA[1:5,6], "\n", sep=" ")
}


plot.HWorig <- function(x, interactive=TRUE, ...) {
 lcv <- x$XMOM[2,]
 lca <- x$XMOM[3,]
 lkur <- x$XMOM[4,]
 n <- x$LEN
 crit <- criticalD()
 Dcrit <- approx(crit[,1], crit[,2], n, method="linear")$y
 Dcrit[n > 15] <- 3
 Dcrit[n < 5] <- 1.333
 outlrs <- ifelse(x$D > Dcrit, 2, 1)
 rlm <- x$RMOM[3:4]
 plot(lca, lcv, col=outlrs, pch=outlrs, ...)
 legend("topleft","D > crit D",pch=2,col=2)
 if(interactive) {par(ask = interactive())}
 plot(lca, lkur, col=outlrs, pch=outlrs, ...)
 Lmoment.ratio.diagram(); points(rlm[1],rlm[2], col = 2, pch=19)
}







# ------------------------------------------------------------------------------ #

LMR <- function (PARA, distr="EXP") {
 if ((distr == "GAM") || (distr == "PE3")) {NMOM <- as.integer(4)}
 else {NMOM <- as.integer(6)}
 XMOM <- as.double(rep(1,NMOM))
 PARA <- as.double(PARA)
 
 if(distr == "EXP") {out <- .Fortran("LMREXP", PARA=PARA, XMOM=XMOM, NMOM=NMOM, PACKAGE="nsRFA")}
 else if (distr == "GAM") {out <- .Fortran("LMRGAM", PARA=PARA, XMOM=XMOM, NMOM=NMOM, PACKAGE="nsRFA")}
 else if (distr == "GEV") {out <- .Fortran("LMRGEV", PARA=PARA, XMOM=XMOM, NMOM=NMOM, PACKAGE="nsRFA")}
 else if (distr == "GLO") {out <- .Fortran("LMRGLO", PARA=PARA, XMOM=XMOM, NMOM=NMOM, PACKAGE="nsRFA")}
 else if (distr == "GNO") {out <- .Fortran("LMRGNO", PARA=PARA, XMOM=XMOM, NMOM=NMOM, PACKAGE="nsRFA")}
 else if (distr == "GPA") {out <- .Fortran("LMRGPA", PARA=PARA, XMOM=XMOM, NMOM=NMOM, PACKAGE="nsRFA")}
 else if (distr == "GUM") {out <- .Fortran("LMRGUM", PARA=PARA, XMOM=XMOM, NMOM=NMOM, PACKAGE="nsRFA")}
 else if (distr == "KAP") {out <- .Fortran("LMRKAP", PARA=PARA, XMOM=XMOM, NMOM=NMOM, PACKAGE="nsRFA")}
 else if (distr == "NOR") {out <- .Fortran("LMRNOR", PARA=PARA, XMOM=XMOM, NMOM=NMOM, PACKAGE="nsRFA")}
 else if (distr == "PE3") {out <- .Fortran("LMRPE3", PARA=PARA, XMOM=XMOM, NMOM=NMOM, PACKAGE="nsRFA")}
 else if (distr == "WAK") {out <- .Fortran("LMRWAK", PARA=PARA, XMOM=XMOM, NMOM=NMOM, PACKAGE="nsRFA")}

 output <- out$XMOM
 if ((distr == "GAM") || (distr == "PE3")) {names(output) <- c("lambda1","lambda2","tau3","tau4")}
 else {names(output) <- c("lambda1","lambda2","tau3","tau4","tau5","tau6")}
 return(output)
}


PEL <- function (XMOM, distr="EXP") {
 XMOM <- as.double(XMOM)
 PARA <- as.double(rep(0,5))
 IFAIL <- as.integer(0)

 if(distr == "EXP") {out <- .Fortran("PELEXP", XMOM=XMOM, PARA=PARA, PACKAGE="nsRFA")}
 else if (distr == "GAM") {out <- .Fortran("PELGAM", XMOM=XMOM, PARA=PARA, PACKAGE="nsRFA")}
 else if (distr == "GEV") {out <- .Fortran("PELGEV", XMOM=XMOM, PARA=PARA, PACKAGE="nsRFA")}
 else if (distr == "GLO") {out <- .Fortran("PELGLO", XMOM=XMOM, PARA=PARA, PACKAGE="nsRFA")}
 else if (distr == "GNO") {out <- .Fortran("PELGNO", XMOM=XMOM, PARA=PARA, PACKAGE="nsRFA")}
 else if (distr == "GPA") {out <- .Fortran("PELGPA", XMOM=XMOM, PARA=PARA, PACKAGE="nsRFA")}
 else if (distr == "GUM") {out <- .Fortran("PELGUM", XMOM=XMOM, PARA=PARA, PACKAGE="nsRFA")}
 else if (distr == "KAP") {out <- .Fortran("PELKAP", XMOM=XMOM, PARA=PARA, IFAIL=IFAIL, PACKAGE="nsRFA")}
 else if (distr == "NOR") {out <- .Fortran("PELNOR", XMOM=XMOM, PARA=PARA, PACKAGE="nsRFA")}
 else if (distr == "PE3") {out <- .Fortran("PELPE3", XMOM=XMOM, PARA=PARA, PACKAGE="nsRFA")}
 else if (distr == "WAK") {out <- .Fortran("PELWAK", XMOM=XMOM, PARA=PARA, IFAIL=IFAIL, PACKAGE="nsRFA")}

 output <- out$PARA
 if(distr == "EXP") {output <- output[c(1:2)]; names(output) <- c("xi","alfa")}
 else if (distr == "GAM") {output <- output[c(1:2)]; names(output) <- c("alfa","beta")}
 else if (distr == "GEV") {output <- output[c(1:3)]; names(output) <- c("xi","alfa","k")}
 else if (distr == "GLO") {output <- output[c(1:3)]; names(output) <- c("xi","alfa","k")}
 else if (distr == "GNO") {output <- output[c(1:3)]; names(output) <- c("xi","alfa","k")}
 else if (distr == "GPA") {output <- output[c(1:3)]; names(output) <- c("xi","alfa","k")}
 else if (distr == "GUM") {output <- output[c(1:2)]; names(output) <- c("xi","alfa")}
 else if (distr == "KAP") {output <- output[c(1:4)]; names(output) <- c("xi","alfa","k","h")}
 else if (distr == "NOR") {output <- output[c(1:2)]; names(output) <- c("mu","sigma")}
 else if (distr == "PE3") {output <- output[c(1:3)]; names(output) <- c("mu","sigma","gamm")}
 else if (distr == "WAK") {names(output) <- c("xi","alfa","beta","gamm","delta")}
 return(output)
}



# --------------------------------------------------------------------------- #

SAMLMR <- function (X, A=0, B=0) {
 X <- as.double(sort(X))
 N <- as.integer(length(X))
 NMOM <- as.integer(6)
 XMOM <- as.double(rep(1,NMOM))
 A <- as.double(A)
 B <- as.double(B)

 out <- .Fortran("SAMLMR", X=X, N=N, XMOM=XMOM, NMOM=NMOM, A=A, B=B, PACKAGE="nsRFA")

 output <- out$XMOM
 names(output) <- c("l1","l2","t3","t4","t5","t6")
 return(output)
}


SAMLMU <- function (X) {
 X <- as.double(sort(X))
 N <- as.integer(length(X))
 NMOM <- as.integer(6)
 XMOM <- as.double(rep(1,NMOM))

 out <- .Fortran("SAMLMU", X=X, N=N, XMOM=XMOM, NMOM=NMOM, PACKAGE="nsRFA")

 output <- out$XMOM
 names(output) <- c("l1","l2","t3","t4","t5","t6")
 return(output)
}


SAMPWM <- function (X, A=0, B=0) {
 X <- as.double(sort(X))
 N <- as.integer(length(X))
 NMOM <- as.integer(6)
 XMOM <- as.double(rep(1,NMOM))
 A <- as.double(A)
 B <- as.double(B)
 KIND <- 2

 out <- .Fortran("SAMPWM", X=X, N=N, XMOM=XMOM, NMOM=NMOM, A=A, B=B, KIND=KIND, PACKAGE="nsRFA")

 output <- out$XMOM
 names(output) <- c("b0","b1","b2","b3","b4","b5")
 return(output)
}


# -------------------------------------------------------------------------------------- #

REGLMR <- function (data, cod) {
 NSITE <- as.integer(length(unique(cod)))
 ni <- tapply(data,cod, length)
 NMOM <- as.integer(6)
 NXMOM <- NMOM
 XMOM <- array(sapply(split(data,cod),SAMLMR), dim=c(NMOM,NSITE))
 WEIGHT <- as.double(ni)
 RMOM <- as.double(rep(0,NMOM))

 out <- .Fortran("REGLMR", NSITE=NSITE, NMOM=NMOM, NXMOM=NXMOM, XMOM=XMOM, WEIGHT=WEIGHT, RMOM=RMOM, PACKAGE="nsRFA")

 output <- out
 return(output)
}


# ------------------------------------------------------------------------------------ #

REGTST <- function(data, cod, A=0, B=0, Nsim=500) {
 # x = data
 # cod = codes
 #
 # PARAMETERS OF ROUTINE:
 # NSITES * INPUT* NUMBER OF SITES IN REGION
 # NAMES * INPUT* CHARACTER*12 ARRAY OF LENGTH NSITES. SITE NAMES.
 # LEN * INPUT* ARRAY OF LENGTH NSITES. RECORD LENGTHS AT EACH SITE.
 # XMOM * INPUT* ARRAY OF DIMENSION (5,NSITES). ARRAY CONTAINING
 #   THE FIRST 5 SAMPLE L-MOMENTS FOR EACH SITE, IN THE
 #   ORDER MEAN, L-CV, L-SKEWNESS, L-KURTOSIS, T-5, I.E
 #   XMOM(I,J) CONTAINS THE I'TH L-MOMENT FOR SITE J.
 #   N.B. XMOM(2,.) CONTAINS L-CV, NOT THE USUAL L-2!
 # A * INPUT* ) PARAMETERS OF
 # B * INPUT* ) PLOTTING POSITION.
 #   NOTE: A AND B SHOULD BE THE SAME AS THE VALUES USED
 #   TO CALCULATE THE MOMENTS IN THE XMOM ARRAY.
 # SEED * INPUT* SEED FOR RANDOM NUMBER GENERATOR. SHOULD BE A WHOLE
 #   NUMBER IN THE RANGE 2D0 TO 2147483647D0.
 # NSIM * INPUT* NUMBER OF SIMULATED WORLDS FOR HETEROGENEITY AND
 #   GOODNESS-OF-FIT TESTS.
 #   NOTE: NSIM=0 WILL FORCE RETURN AT COMPLETION OF
 #   OUTLIER TEST. NSIM=1 WILL SUPPRESS CALCULATION OF
 #   H AND Z STATISTICS, BUT PARAMETER AND QUANTILE
 #   ESTIMATES WILL BE FOUND.
 # NPROB * INPUT* NUMBER OF QUANTILES TO BE CALCULATED
 # PROB * INPUT* ARRAY OF LENGTH NPROB. PROBABILITIES FOR WHICH
 #   QUANTILES ARE TO BE CALCULATED.
 # KPRINT * INPUT* OUTPUT FLAG. SHOULD BE SET TO
 #   0 TO SUPPRESS OUTPUT
 #   1 TO PRINT OUTPUT
 # KOUT * INPUT* CHANNEL TO WHICH OUTPUT IS DIRECTED
 # RMOM *OUTPUT* ARRAY OF LENGTH 5. ON EXIT, CONTAINS THE REGIONAL
 #   WEIGHTED AVERAGE L-MOMENT RATIOS.
 # D *OUTPUT* ARRAY OF LENGTH NSITES. ON EXIT, CONTAINS THE
 #   DISCORDANCY MEASURE (D STATISTIC) FOR EACH SITE.
 # VOBS *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE REGIONAL
 #   OBSERVED VALUES OF 3 HETEROGENEITY STATISTICS:
 #   (1) WEIGHTED S.D. OF L-CVS;
 #   (2) AVERAGE OF L-CV/L-SKEW DISTANCES;
 #   (3) AVERAGE OF L-SKEW/L-KURTOSIS DISTANCES.
 # VBAR *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE MEAN OF THE
 #   SIMULATED VALUES OF THE 3 HETEROGENEITY STATISTICS.
 # VSD *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE S.D. OF THE
 #   SIMULATED VALUES OF THE 3 HETEROGENEITY STATISTICS.
 # H *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS HETEROGENEITY
 #   MEASURES (H STATISTICS), I.E. H=(VOBS-VBAR)/VSD.
 # Z *OUTPUT* ARRAY OF LENGTH 5. ON EXIT, CONTAINS GOODNESS-OF-FIT
 #   MEASURES (Z STATISTICS) FOR 5 DISTRIBUTIONS:
 #   (1) GEN. LOGISTIC, (2) GEN. EXTREME VALUE,
 #   (3) GEN. NORMAL, (4) PEARSON TYPE III,
 #   (5) GEN. PARETO.
 # PARA *OUTPUT* ARRAY OF DIMENSION (5,6). ON EXIT, IF NSIM.GE.1,
 #   CONTAINS PARAMETERS OF GROWTH CURVES FITTED BY THE
 #   ABOVE 5 DISTRIBUTIONS, PLUS WAKEBY.
 # RPARA *OUTPUT* PARAMETERS OF REGIONAL KAPPA DISTRIBUTION

 x <- data

 # Input
 NSITES <- as.integer(length(unique(cod)))
 NAMES <- unique(cod)
 LEN <- as.integer(tapply(x, cod, length))
 NMOM <- as.integer(5)
 XMOM <- array(sapply(split(x,cod),SAMLMR, A=A, B=B)[c(1:NMOM),], dim=c(NMOM,NSITES))
 XMOM[2,] <- XMOM[2,]/XMOM[1,]
 A <- as.double(A)  # -0.35
 B <- as.double(B)
 # FOR UNBIASED ESTIMATES (OF THE LAMBDA'S) SET A=B=ZERO. OTHERWISE,
 # PLOTTING-POSITION ESTIMATORS ARE USED, BASED ON THE PLOTTING POSITION
 # (J+A)/(N+B) FOR THE J'TH SMALLEST OF N OBSERVATIONS. FOR EXAMPLE,
 # A=-0.35D0 AND B=0.0D0 YIELDS THE ESTIMATORS RECOMMENDED BY
 # HOSKING ET AL. (1985, TECHNOMETRICS) FOR THE GEV DISTRIBUTION.
 #SEED <- as.integer(ifelse(exists(".Random.seed"), .Random.seed[1], 56389))
 SEED <- 619145091
 NSIM <- as.integer(Nsim)
 PROB <- as.double(c(.01, .02, .05, .1, .2, .5, .9, .95, .99, .999))
 NPROB <- as.integer(length(PROB))
 KPRINT <- as.integer(0)
 KOUT <- as.integer(6)   # ???

 # Output initialization
 RMOM <- double(length=NMOM)   # as.numeric(rep(1, 5))
 D <- double(length=NSITES)   # rep(1, NSITES)
 VOBS <- double(length=3)   # rep(1, 3)
 VBAR <- double(length=3)   # rep(1, 3)
 VSD <- double(length=3)   # rep(1, 3)
 H <- double(length=3)   # rep(1, 3)
 Z <- double(length=5)   # rep(1, 5)
 PARA <- matrix(1, nrow=5, ncol=6)

 out <- .Fortran("REGTST", NSITES=NSITES, NAMES=NAMES, LEN=LEN, XMOM=XMOM, A=A, B=B,
                           SEED=SEED, NSIM=NSIM, NPROB=NPROB, PROB=PROB, KPRINT=KPRINT,
                           KOUT=KOUT, RMOM=RMOM, D=D, VOBS=VOBS, VBAR=VBAR, VSD=VSD, H=H, Z=Z,
                           PARA=PARA, PACKAGE="nsRFA")
 return(out)
}
