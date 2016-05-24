balloon <-
function(stim, totaltime, acc, par=list(),verbose=TRUE){

	#require(deSolve, quietly=TRUE)
	
	if(missing(par)){
	  if(verbose==TRUE){
	    warning("Default parameter values are used. See ?balloon.fnc for more details.")
	  }
	  par <- list()
	  par$kappa <- 2
	  par$tau1 <- 3
	  par$tauf <- 4
	  par$taum <- 4
	  par$f1 <- 1.5
	  par$deltat <- 1
	  par$n <- 3
	  par$E0 <- 0.4
	  par$V0 <- 0.03
	  par$a1 <- 3.4
	  par$a2 <- 1.0
	  par$tauMTT <- 3
	  par$tau <- 20
	  par$alpha <- 0.4
	}

	par$deltatf <- 0
	par$deltatm <- par$deltat - par$deltatf
	par$m1 <- (par$f1-1)/par$n +1

	t <- seq(acc, totaltime, acc)
	it <- c(1:(totaltime/acc))

	par.I <- c(par$tau1, par$kappa)
	inhib <- function(t,y,p){
		yd1 <- (1/p[1])*(p[2]*stim[t] - (p[2]+1)*y[1])
		list(c(yd1))
	}
	I <- ode(c(0), it, inhib, par.I, method=rkMethod("ode45"))[,2]
	N <- stim - I
	F <- 1 + convolve((par$f1-1)*gammaHRF(t-par$deltatf,FWHM=par$tauf,verbose=FALSE),rev(N), type="o")
	M <- 1 + convolve((par$m1-1)*gammaHRF(t-par$deltatm,FWHM=par$taum,verbose=FALSE),rev(N), type="o")
	E <- par$E0*M/F
	par.balloon <- c(par$E0, par$tauMTT, par$tau, par$alpha)
	de.balloon <- function(t,y,p){
		dv <- (F[t]-y[1]^(1/p[4]))/(p[2]+p[3])
		dq <- 1/p[2]*(F[t]*E[t]/p[1] - y[2]/y[1]*(y[1]^(1/p[4]) + p[3]/(p[2] + p[3])*(F[t] - y[1]^(1/p[4]))))
		list(c(dv,dq))
	}
	balloon <- ode(c(1,1),it,de.balloon,par.balloon,method=rkMethod("ode45"))
	V <- balloon[,2]
	Q <- balloon[,3]
	bold <- par$V0*(par$a1*(1-Q)-par$a2*(1-V))

	return(bold)
}

