#' Evalute ACME Reduction Factor R*
#'
#' Calculates R*, the reduction factor (inverse of inflation factor for mortality
#' estimates), based on parameter estimates and maximum number of previous
#' search intervals to consider.
#'
#'@param Iij Search interval length (days)
#'@param arabt 5-element vector: alpha and rho parameters for the Weibull
#' persistence distribution, a and b parameters for the Exponential search
#' proficiency distribution, and bt as the bleed-through parameter
#'@param kmax number of intervals to use in calculation - includes current
#' interval and number of look-back intervals. Minimum is 1 (only current
#' interval).
#'@param v logical. Verbose flag - see Value for what is reported             
#' 
#'@return If verbose, returns upper bound on truncation error (total and as a
#' fraction of R*), R* calculated from only current interval, and expected
#'  fraction of "old" carcasses discovered. If not verbose, returns only R*.
#' 
#
Rst <- function(Iij=7, arabt=c(alp=0.4695,rho=0.0809,
                a=1.0322,b=0.0706,bt=0.9573), kmax=5, v=FALSE) {
  #Compute R* from Wolpert (2015) Eqn (7a,b)
  
    alp <- arabt[1]; a <- arabt[3];
    rho <- arabt[2]; b <- arabt[4]; bt <- arabt[5];

    pars <- c(a=a, bI=(bI<-b*Iij), alp=alp, rI=(rI<-rho*Iij), bt=bt);
    kmn  <- recur(kmax);
    Q    <- numeric(nQ<-dim(kmn)[1]);
    for(i in 1:nQ) {
      kk <- kmn[i,"k"]; mm <- kmn[i,"m"]; nn <- kmn[i,"n"];
      Q[i] <- - bt^kk * (-1)^mm * simp(Fs,, kmn[i,], pars);
    } # Truncation error bound:
    bnd <- (bt^kmax) * Fs(0, c(k=kmax,m=1,n=kmax), pars) *
      min((1-exp(-a))^kmax, exp(-a-kmax*bI)/max(bI,(1-bt*exp(-bI))));
    R   <- sum(Q);
    if(v) { return(c(Rstar=R, err=bnd, frc=bnd/R, T0=Q[1], old=1-Q[1]/R)); }
    else  { return(list(Rstar=R, T0 = Q[1])); }
}
