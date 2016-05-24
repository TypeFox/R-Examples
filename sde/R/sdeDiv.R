sdeDiv <- 
function (X, theta1, theta0, phi= expression( -log(x) ), C.phi, K.phi, b, s, b.x, s.x, s.xx, B, B.x, H, S, guess, 
    ...) 
{
    if(missing(theta0))
     stop("please specify null hypotesis `theta0'")

    if(!missing(theta1) & !missing(theta0)){
      if(length(theta1) != length(theta0))
       stop("initial values and null hypothesis of different dimensions")
    }
    e1 <- new.env()
	
    if(missing(C.phi)){
	 d1 <- D(phi, "x")
     e1$x <- 1
	 C.phi <- eval(d1,e1)
	}

    if(missing(K.phi)){
	 d2 <- D(D(phi, "x"),"x")
     e1$x <- 1
	 K.phi <- eval(d2,e1)
	}

    n <- length(X)
    DELTA <- deltat(X)
    if (missing(theta1) || is.null(theta1)) {
        if (missing(guess)) 
            stop("cannot estimate the model. Specify initial guess values for theta1")
 
	    Euler.lik <- function(theta)
         -sum( dcEuler( x=X[2:n], t=DELTA, x0=X[1:(n-1)], t0=0, theta=theta, 
	               d=function(t,x,theta) b(x,theta), s=function(t,x,theta) s(x,theta),
                   log=TRUE) ) 

        cat("estimating the model...\n")
        est <- optim(guess, Euler.lik, ...)
        theta1 <- est$par		
        cat("\nestimated parameters\n")
		cat(theta1)
		cat("\n")

    }

    if (missing(B)) 
        B <- function(x, theta) b(x, theta)/s(x, theta) - 0.5 * 
            s.x(x, theta)
    if (missing(B.x) && !(missing(b.x) || missing(s.x) || missing(s.xx))) {
        B.x <- function(x, theta) (b.x(x, theta)/s(x, theta) - 
            b(x, theta) * s.x(x, theta)/(s(x, theta)^2) - 0.5 * 
            s.xx(x, theta))
    }
    else {
        stop("error")
    }
    C1 <- function(x,theta) (B(x, theta)^2)/3 + 0.5 * B.x(x, theta) * 
        s(x, theta)
    g <- function(x, y, theta) -0.5 * (C1(x, theta) + C1(y, theta) + B(x, theta) * 
        B(y, theta)/3)
    h <- function(x,theta) B(x, theta)/s(x, theta)
    if (missing(H)) 
        H <- function(x, y, theta) 
 		  sapply(1:length(x), function(i) integrate(h, x[i], y[i], theta)$value) 
		 
    s1 <- function(x, theta) 1/s(x, theta)
    if (missing(S)) 
        S <- function(x, y, theta) 
		  sapply(1:length(x), function(i) integrate(s1, x[i], y[i], theta)$value)
		 
    u <- function(x, y, theta) (-0.5 * log(2 * pi * DELTA) - 
        log(s(y, theta)) - S(x, y, theta)^2/(2 * DELTA) + H(x, y, theta) + 
        DELTA * g(x, y, theta))

    DC.lik <- function(theta)
      sum(u(X[1:(n - 1)], X[2:n], theta))

	u1 <- DC.lik(theta1)

	u0 <- DC.lik(theta0)
 	
	val <- u1 - u0 
    switched <- FALSE 
	if(val>0) {
	 val <- -val
     switched <- TRUE
    }
	U <- exp(val)
    	
	LRT <- -2 * val

    DIV <- NA

    n.par <- length(theta0)
	e1$x <- U
	DIV <- eval(phi, e1)
    p.DIV <- NA
    p.LRT <- pchisq(LRT, df=n.par, lower.tail = FALSE)
    est.pval <- FALSE
    chi <- NULL
	coef <- 1
	maxsim <- 500000
	if(switched)
	 coef <- -1
    if(K.phi != 0){
	 chi <- rchisq(maxsim, df=n.par)
     if(C.phi==0){
		T <- 0.5*K.phi*chi^2
	 } else {
	  T <- coef*0.5*(C.phi*chi) + 0.5*(C.phi+K.phi)*chi^2
	 } 
	 p.DIV <- sum(T>DIV)/maxsim
	 est.pval <- TRUE
	}
    
	if(C.phi!=0 & K.phi==0)
	  p.DIV <- pchisq( coef*DIV/C.phi, df=n.par, lower.tail = FALSE)

    cat("\nTesting H0 against H1\n\nH0: ")
	cat(theta0)
	cat("\nH1: ")
	cat(theta1)


	cat(sprintf("\n\nDivergence statistic: %.5g (p-value=%.5g)\nLikelihood ratio test statistic: %.5g (p-value=%.5g)\n\n",DIV, p.DIV, LRT, p.LRT))
    cat("\n(p-value of phi-divergence approximated by simulation)\n")
	
	invisible(
	list(theta1=theta1, theta0=theta0, DIV=DIV, p.DIV=p.DIV, C.phi=C.phi, K.phi=K.phi, 
	 LRT=LRT, p.LRT=p.LRT, est.pval = est.pval))
}


