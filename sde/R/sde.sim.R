sde.sim <- function (t0 = 0, T = 1, X0 = 1, N = 100, delta, drift, sigma, 
    drift.x, sigma.x, drift.xx, sigma.xx, drift.t, method = c("euler", 
        "milstein", "KPS", "milstein2", "cdist","ozaki","shoji","EA"), 
		alpha = 0.5, eta = 0.5, pred.corr = T, rcdist = NULL, theta = NULL,
	model = c("CIR", "VAS", "OU", "BS"),
	k1, k2, phi, max.psi = 1000, rh, A, M=1) 
{
	method <- match.arg(method)
	if(!missing(model)){
	 model <- match.arg(model)
	 method <- "model"
    }

    x0 <- rep(X0,M)[1:M]

	if (missing(drift)){
     if (method == "cdist" || !missing(model))
      drift <- expression(NULL)
     else  
        stop("please specify al least the drift coefficient of the SDE")
    }
	
    if (missing(sigma)) 
        sigma <- expression(1)
    if (!is.expression(drift) || !is.expression(sigma)) 
        stop("coefficients must be expressions in `t' and `x'")
    if (pred.corr == F) {
        alpha <- 0
        eta <- 0
        sigma.x <- NULL
    }
    needs.sx <- FALSE
    needs.dx <- FALSE
    needs.sxx <- FALSE
    needs.dxx <- FALSE
    needs.dt <- FALSE    

	if (method == "cdist" && is.null(rcdist)) 
        stop("please provide a random number generator `rcdist'")
    if (method == "milstein") 
        needs.sx = TRUE
    if ((method == "euler" && pred.corr == T)) 
        needs.sx = TRUE
    if (method == "KPS" || method == "milstein2") {
        needs.sx <- TRUE
        needs.dx <- TRUE
        needs.sxx <- TRUE
        needs.dxx <- TRUE
    }
	if(method == "ozaki" || method == "shoji" || method == "EA")
	    needs.dx <- TRUE
	if(method == "shoji"){
	    needs.dxx <- TRUE
		needs.dt <- TRUE
	}
    if (needs.sx && missing(sigma.x)) {
        message("sigma.x not provided, attempting symbolic derivation.\n")
        sigma.x <- D(sigma, "x")
    }
    if (needs.dx && missing(drift.x)) {
        message("drift.x not provided, attempting symbolic derivation.\n")
        drift.x <- D(drift, "x")
    }
    if (needs.dxx && missing(drift.xx)) {
        message("drift.xx not provided, attempting symbolic derivation.\n")
        drift.xx <- D(D(drift, "x"), "x")
    }
    if (needs.sxx && missing(sigma.xx)) {
        message("sigma.xx not provided, attempting symbolic derivation.\n")
        sigma.xx <- D(D(sigma, "x"), "x")
    }
    if (needs.dt && missing(drift.t)) {
        message("drift.t not provided, attempting symbolic derivation.\n")
        drift.t <- D(drift, "t")
    }
    d1 <- function(t, x) eval(drift)
    d1.x <- function(t, x) eval(drift.x)
    d1.xx <- function(t, x) eval(drift.x)
    d1.t <- function(t, x) eval(drift.t)
    s1 <- function(t, x) eval(sigma)
    s1.x <- function(t, x) eval(sigma.x)
    s1.xx <- function(t, x) eval(sigma.xx)

    if (t0 < 0 || T < 0) 
        stop("please use positive times!")
    if (missing(delta)) {
        t <- seq(t0, T, length = N + 1)
    }
    else {
        t <- c(t0, t0 + cumsum(rep(delta, N)))
        T <- t[N + 1]
        message(sprintf("\nT set to = %f\n", T))
    }
    Dt <- (T - t0)/N


	if(method == "model"){
      if(is.null(theta)) 
        stop("please provide a vector of parameters for the model")
     if(model == "CIR")
 	    X <- sde.sim.cdist(x0, t0, Dt, N, M, rcCIR, theta)
	 if(model == "OU")
	    X <- sde.sim.cdist(x0, t0, Dt, N, M, rcOU, theta)
	 if(model == "BS")
	    X <- sde.sim.cdist(x0, t0, Dt, N, M, rcBS, theta)
	  }
	 
     if (method == "EA")	 
      X <- sde.sim.ea(X0, t0, Dt, N, d1, d1.x, k1, k2, phi, max.psi, rh, A)

     if (method == "cdist"){
      if(is.null(theta)) 
        stop("please provide a vector of parameters for `rcdist'")
      else
        X <- sde.sim.cdist(x0, t0, Dt, N, M, rcdist, theta)
     }
	
     if (method == "ozaki"){ 
	    vd <- all.vars(drift)
	    vs <- all.vars(sigma)		
		if((length(vd)!=1) || (length(vs)>0))
		 stop("drift must depend on `x' and volatility must be constant")
  		if((length(vd) == 1) && (vd != "x"))
		 stop("drift must depend on `x'")
		X <- sde.sim.ozaki(x0, t0, Dt, N, M, d1, d1.x, s1)
	 }
     if (method == "shoji"){ 
	    vd <- all.vars(drift)
	    vs <- all.vars(sigma)		
		if(length(vd)>2 || length(vd)<1  || length(vs)>0)
		 stop("drift must depend on `x' and/or `t' and volatility must be constant")
  		if((length(vd) == 1) && (vd != "x"))
		 stop("drift must depend at least on `x'")
		X <- sde.sim.shoji(x0, t0, Dt, N, M, d1, d1.x, d1.xx, d1.t, s1)
	 }
     if (method == "euler") 
        X <- sde.sim.euler(x0, t0, Dt, N, M, d1, s1, s1.x, alpha, 
            eta, pred.corr)
     if (method == "milstein") 
        X <- sde.sim.milstein(x0, t0, Dt, N, M, d1, s1, s1.x)
     if (method == "milstein2") 
        X <- sde.sim.milstein2(x0, t0, Dt, N, M, d1, d1.x, d1.xx, 
            s1, s1.x, s1.xx)
     if (method == "KPS") {
        Sigma <- matrix(c(Dt, 0.5 * Dt^2, 0.5 * Dt^2, 1/3 * Dt^3), 
            2, 2)
        tmp <- mvrnorm(N*M, c(0, 0), Sigma)
        Z <- tmp[, 1]
        U <- tmp[, 2]
        X <- sde.sim.KPS(x0, t0, Dt, N, M, d1, d1.x, d1.xx, s1, 
            s1.x, s1.xx, Z, U)
     }

    nm <- "X"
    nm <- if(M>1)
	 paste("X",1:M,sep="")
    X <- ts(X, start = t0, deltat = Dt, names=nm)
    invisible(X)
}
