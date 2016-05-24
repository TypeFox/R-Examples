####  Program (package) with the t-walk implementation in R
####  see http://www.cimat.mx/~jac/twalk/
####  Author J. Andres Christen




###############################################################
#### Some auxiliary functions and constants:

IntProd <- function(x) { sum(x*x)  } ## square of the norm.
DotProd <- function( x, y) { sum( x*y )  }  ##dot product




########## h1 function for the "traverse" kernel 
Simh1 <- function( dim, pphi, x, xp, beta)
{	
	phi <- (runif(dim) < pphi) 

	rt <- NULL 
	for (i in 1:dim)
		if (phi[i])
			rt <- append( rt, xp[i] + beta*(xp[i] - x[i]))
		else
			rt <- append( rt, x[i]) 
			
	list( rt=rt, nphi=sum(phi)) 
}

### Simulation of the beta parameter for kernel h1
Simfbeta <- function(at)
{
	
	if (runif(1) < (at-1)/(2*at))
		exp(1/(at + 1)*log(runif(1)))
	else
		exp(1/(1 - at)*log(runif(1))) 
}




########## h function for the "walk" kernel 
Simh2 <- function( dim, pphi, aw, x, xp)
{

	 u <- runif(dim) 
	 phi <- (runif(dim) < pphi) 

	 z <- (aw/(1+aw))*(aw*u^2 + 2*u -1) 
	 z <- z*phi 
	 list( rt=x + (x - xp)*z, nphi=sum(phi)) 
	
}



########## h function for the "blow" kernel 
Simh3 <- function( dim, pphi, x, xp)
{
	phi <- (runif(dim) < pphi) 
	
	sigma <- max(phi*abs(xp - x)) 

	x + sigma*rnorm(dim)*phi 

	list( rt=xp*phi + sigma*rnorm(dim)*phi + x*(1-phi), nphi=sum(phi), phi=phi) 

}

## -log of g3, the density of h_b.
G3U <- function( nphi, phi, h, x, xp)
{	sigma <- max(phi*abs(xp - x)) ##different sigma for the proposal and the reverse, but same phi
	if (nphi > 0)
		(nphi/2)*log(2*pi) + nphi*log(sigma) + 0.5*IntProd(h - xp)/(sigma^2)
	else
		0 
}



########## h function for the "hop" kernel 
Simh4 <- function( dim, pphi, x, xp)
{
	phi <- (runif(dim) < pphi) 
	
	sigma <- max(phi*abs(xp - x))/3 
	
	rt <- NULL 
	for (i in 1:dim)
		if (phi[i])
			rt <- append( rt, x[i] + sigma*rnorm(1))
		else
			rt <- append( rt, x[i]) 
			
	list( rt=rt, nphi=sum(phi), phi=phi) 

}
log2pi <- log(2*pi); log3 <- log(3)
## -log of g4, the density of h_h.
G4U <- function( nphi, phi, h, x, xp)
{	sigma <- max(phi*abs(xp - x))/3 ##different sigma for the proposal and the reverse, but same phi
	if (nphi > 0)
		(nphi/2)*log2pi - nphi*log3 + nphi*log(sigma) + 0.5*9*IntProd((h - x))/(sigma^2)
	else
		0 
}





############ This is the twalk implementation.  It requires the "energy" of the
############ objective function ie. U = -log f and its support.  Also the initial
############ values x0 and x0p.  Tr is the number of iterations required.
############ The dimension is dim, defaults to the global variable n
############ in the calling NAMESPACE for backward compatibility. 
############ The rest of the parameters are for ploting dynamically the trajectories
############ of the twalk for objectives of dim=2.  Otherwise set
############ PlotObj=FALSE.  See examples.R for details.
Runtwalk <- function( Tr, dim = length(x0), Obj, Supp, x0, xp0, PlotObj=FALSE, PlotLogPost=TRUE, dynty="b", pathcol="grey", add=FALSE,
	at=6, aw=1.5, pphi=min( dim, 4)/dim, F1=0.4918, F2=F1+0.4918, F3=F2+0.0082, ...) {
	
	
	## Initial values
	x <- x0 
	xp <- xp0 
	
	if (Supp(x) && Supp(xp))
	{ 
		cat("This is the twalk for R.\n") 
		cat("Evaluating the objective density at initial values.")  
		flush(stdout())   
		U <- Obj(x, ...) 
		Up <- Obj(xp, ...) 

		acc <- 0 
		cat("\nOpening", 2+2*dim, "X", Tr+1, "matrix to save output.") 
		flush(stdout())   
		rec <- matrix( 0, ncol=2+2*dim, nrow=Tr+1)   ##to save U, Up, x and xp
		rec[1,] <- append( U, append(Up, append(x, xp))) 

		recacc <- matrix( 0, ncol=2, nrow=Tr) 
		
		flush(stdout())   
		if (is.function(PlotObj))
		{
			cat("\nInitializing graphics.") 	
			PlotObj(add=add) 
			Plpoints <- TRUE 
			PlotLogPost <- FALSE
		}
		else
			Plpoints <- FALSE
	}
	else
	{
		cat(paste("Initial values out of support,\n  x=", x,"\n xp=", xp)) 
		Tr <- 0 
	}
	
	if (any(abs(x0 -xp0) <= 0))
	{
		cat(paste("Not all entries of initial values different,\n  x=", x,"\n xp=", xp)) 
		
		Tr <- 0 
	}
	

	if (Plpoints || PlotLogPost)
		cat("  Sampling, check graphics divice.\n")
	else
		cat("  Sampling (no graphics mode).\n")

	
	acc <- 0
			
	for (it in 1:Tr)
	{

		if (((it %% 1000) == 0)  && PlotLogPost) {
	 		plot( max(1,it-2000):it, -rec[max(1,it-2000):it,1], type="l",
	 			xlab="Iteration", ylab="LogPost", main=paste("dim=", dim) 
) 
		}
			
#########


		move <- OneMove( dim=dim, Obj=Obj, Supp=Supp, x, U, xp, Up,
			at=at, aw=aw, pphi=pphi, F1=F1, F2=F2, F3=F3, ...)


########

		#cat( "move:", move$y, move$yp, move$x, move$xp)
		
		if (runif(1) < move$A)
		{ ## accepted
			if (Plpoints)
			{
				if (dynty == "p")
					points( move$y[1], move$y[2], pch=".", col="blue")
				else {
					lines( c( x[1], move$y[1]), c( x[2], move$y[2]), col=pathcol) 
				}
			}
			
			recacc[it,] <- append( move$funh, move$nphi/dim)  
			acc <-  acc + move$nphi/dim 
			
			x <- move$y 
			U <- move$propU 
			
			xp <- move$yp 
			Up <- move$propUp 

		}
		else
			recacc[it,] <- append( move$funh, 0) 
			
		
		rec[it+1,] <- c( U, Up, x, xp) 
		
	}
	
		
	if (is.function(PlotObj))
		PlotObj(add=TRUE) 
	
	list( dim=dim, Tr=Tr, acc=acc, Us=rec[,1], Ups=rec[,2],
		output=rec[,2+(1:dim)], outputp=rec[,2+dim+(1:dim)], recacc=recacc) 
}	


OneMove <- function( dim, Obj, Supp, x, U, xp, Up,
	at=6, aw=1.5, pphi=min( dim, 4)/dim, F1=0.4918, F2=0.9836, F3=0.9918, ...) {

		ker <- runif(1)  ##choose a kernel
		
		
		if (ker < F1) ## the t-walk, kerlnel h1: traverse
		{	
			dir <- runif(1) 
			funh <- 1 

			if ((0 <= dir) && (dir < 0.5))
			{	
				beta <- Simfbeta(at) 

				tmp <- Simh1( dim, pphi, xp, x, beta) 
				yp <- tmp$rt
				nphi <- tmp$nphi
				
				y  <- x 
				propU <- U 
				
				if (Supp(yp))
				{
					propUp <- Obj(yp, ...) 
					## The propolsal is simetric

					if (nphi == 0)
						A <- 1 ###Nothing moved
					else
						A <- exp((U - propU) + (Up - propUp) +  (nphi-2)*log(beta)) 
				}
				else {
					propUp <- NULL
					A <- 0  ##out of support, not accepted
				}
			}

			if ((0.5 <= dir) && (dir < 1.0))
			{
				beta <- Simfbeta(at) 
				
				tmp <- Simh1( dim, pphi, x, xp, beta) 
				y <- tmp$rt
				nphi <- tmp$nphi

				yp  <- xp 
				propUp <- Up 

				if (Supp(y))
				{
					propU <- Obj(y, ...) 
					## The propolsal is simetric

					if (nphi == 0)
						A <- 1 ###Nothing moved
					else
						A <- exp((U - propU) + (Up - propUp) +  (nphi-2)*log(beta)) 
				}
				else {
					propU <- NULL
					A <- 0  ##out of support, not accepted
				}
			}
		}


		
		if ((F1 <= ker) && (ker < F2)) ## the t-walk, kerlnel h2: walk
		{	
			dir <- runif(1) 
			funh <- 2 

			if ((0 <= dir) && (dir < 0.5))  ## x as pivot
			{
				tmp <- Simh2( dim, pphi, aw, xp, x) 
				yp <- tmp$rt
				nphi <- tmp$nphi

				y  <- x 
				propU <- U 

				if ((Supp(yp)) && (all(abs(yp - y) > 0)))
				{
					propUp <- Obj(yp, ...) 
					A <- exp((U - propU) + (Up - propUp))  
				}
				else{
					propUp <- NULL
					A <- 0  ##out of support, not accepted
				}
			}

			if ((0.5 <= dir) && (dir < 1.0))  ## xp as pivot
			{
				tmp <- Simh2( dim, pphi, aw, x, xp) 
				y <- tmp$rt
				nphi <- tmp$nphi

				yp  <- xp 
				propUp <- Up 

				if ((Supp(y)) && (all(abs(yp - y) > 0)))
				{
					propU <- Obj(y, ...) 
					A <- exp((U - propU) + (Up - propUp))  
				}
				else{
					propU <- NULL
					A <- 0  ##out of support, not accepted
				}
			}
		}


		if ((F2 <= ker) && (ker < F3)) ## the t-walk, kernel h3: blow
		{	
			dir <- runif(1) 
			funh <- 3 

			if ((0 <= dir) && (dir < 0.5))  ## x as pivot
			{
				tmp <- Simh3( dim, pphi, xp, x) 
				yp <- tmp$rt
				nphi <- tmp$nphi
				phi <- tmp$phi
				
				y  <- x 
				propU <- U 
				if ((Supp(yp)) && all(yp != x))
				{
					propUp <- Obj(yp, ...) 
					W1 <- G3U( nphi, phi, yp, xp,  x) 
					W2 <- G3U( nphi, phi, xp, yp,  x)  
					A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))
				}
				else{
					propUp <- NULL
					A <- 0  ##out of support, not accepted
				}
			}

			if ((0.5 <= dir) && (dir < 1.0))  ## xp as pivot
			{
				tmp <- Simh3( dim, pphi, x, xp) 
				y <- tmp$rt
				nphi <- tmp$nphi
				phi <- tmp$phi

				yp  <- xp 
				propUp <- Up 

				if ((Supp(y)) && all(y != xp))
				{
					propU <- Obj(y, ...) 
					W1 <- G3U( nphi, phi, y,  x, xp) 
					W2 <- G3U( nphi, phi, x,  y, xp) 
					A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))
				}
				else{
					propU <- NULL
					A <- 0  ##out of support, not accepted
				}
			}

		}
		

		if (F3 <= ker) ## the t-walk, kernel h4: hop
		{	
			dir <- runif(1) 
			funh <- 4 

			if ((0 <= dir) && (dir < 0.5))  ## x as pivot
			{
				tmp <- Simh4( dim, pphi, xp, x) 
				yp <- tmp$rt
				nphi <- tmp$nphi
				phi <- tmp$phi

				y  <- x 
				propU <- U 

				if ((Supp(yp)) && all(yp != x))
				{
					propUp <- Obj(yp, ...) 
					W1 <- G4U( nphi, phi, yp, xp,  x) 
					W2 <- G4U( nphi, phi, xp, yp,  x) 
					A <- exp((U - propU) + (Up - propUp) +  (W1 - W2)) 
				}
				else{
					propUp <- NULL
					A <- 0  ##out of support, not accepted
				}
			}

			if ((0.5 <= dir) && (dir < 1.0))  ## xp as pivot
			{
				tmp <- Simh4( dim, pphi, x, xp) 
				y <- tmp$rt
				nphi <- tmp$nphi
				phi <- tmp$phi

				yp  <- xp 
				propUp <- Up 

				if ((Supp(y)) && all(y != xp))
				{
					propU <- Obj(y, ...) 
					W1 <- G4U( nphi, phi, y,  x, xp) 
					W2 <- G4U( nphi, phi, x,  y, xp) 
					A <- exp((U - propU) + (Up - propUp) +  (W1 - W2)) 
				}
				else{
					propU <- NULL
					A <- 0  ##out of support, not accepted
				}
			}
		
		}
				

		if (is.nan(A))  #### debugging line
		{
			cat("Rtwalk: ERROR, in evaluating the objective.  Value returned by objective function:", propU)
		}
		
		list( y=y, propU=propU, yp=yp, propUp=propUp, A=A, funh=funh, nphi=nphi)		
}




######### Simple basic functions to analyse the output of Runtwalk, commonly saved in info.

######### To plot a time series of the log of the posterios
PlotLogObj <- function(info, from=0, to=length(info$Us))
{
	plot( from:(to-1), -info$Us[(from+1):to], type="l",
		ylab="Log of Objective", xlab="Iteration", main="") 
}


######### To plot a histogram of any parameter
PlotHist <- function( info, par=1, from=0, xlab=paste("Parameter", par), main="", ...)
{
	if (info$dim == 1)
		hist( info$output[from:(info$Tr)], xlab=xlab, main=main, ...)
	else
		hist( info$output[from:(info$Tr), par], xlab=xlab, main=main, ...)
}

SaveOutput <- function( info, file, pars=1:(info$dim), from=1, to=info$Tr,
	row.names=FALSE, col.names=paste("X", pars), ...) {

	if (info$dim == 1)
		write.table( info$output[from:(info$Tr)], file=file,
			row.names=row.names, col.names=col.names, ...)
	else
		write.table( info$output[from:(info$Tr), pars], file=file,
			row.names=row.names, col.names=col.names, ...)
}


############ To plot an outline of the output:
####  Calculate IAT and acceptance ratios
Ana <- function(info, from=1, to=info$Tr, par=0, file="")
{
	sel <- from:to 
		
	accrt <- rep(0, 4) 
	for (h in 1:4)
	{
		selh <- which(info$recacc[sel,1] == h) 
		accrt[h] <- sum(info$recacc[selh,2])/length(selh) 
	}

	
	#### No plots

	Tint <- IAT( info, par=par)  ### defined below
	
	itmap = which(-info$Us == max(-info$Us))[1] 
	
	
	cat( "Ratio of moved coodinates per it=\n",
		 accrt[1], accrt[2], accrt[3], accrt[4],
		 "\ndim=", info$dim, "AcceptanceRatio=", info$acc/info$Tr,
		"MAPlogPost=", -info$Us[itmap], "IAT=", Tint, "IAT/dim=", Tint/info$dim,"\n\n") 

	if (file != "")
	 cat(file=file, info$dim, 
		 accrt[1], accrt[2], accrt[3], accrt[4], info$acc/info$Tr,
		-info$Us[itmap], Tint/info$dim,"\n") 
	
}
	

### Plot time series of parameters

TS <- function(info, pars=1:(info$dim), from=1, to=info$Tr, prime=FALSE)
{	
	sel <- from:to 
	
	if (length(pars) <= 10)
	{
		if (info$dim == 1)
			if (!prime)
				plot(as.ts(as.matrix(info$output[sel])),  main="x")
			else
				plot(as.ts(as.matrix(info$outputp[sel])), main="xp")
		else
			if (!prime)
				plot(as.ts(as.matrix(info$output[sel, pars])),  main="x")
			else
				plot(as.ts(as.matrix(info$outputp[sel, pars])), main="xp") 
	}
	else
		cat("Cannot print time series for more than 10 parameters, select a subset with arg. pars\n\n") 
	
}





###### These functions are for calculating and plotting autcorrelations and
###### Integrated Autocorrelation Times



GetAutoCorr <- function( info, par=0, from=1, to=info$Tr, lag=30*info$dim) {
	
	if (par>0)
		cor( info$output[from:(to-lag), par], info$output[(from+lag):to, par])
	else
		cor( info$Us[from:(to-lag)], info$Us[(from+lag):to]) 
}
	
	

GetAutoCov <- function( dt, lags)
{
	n <- length(dt) 
	
	aut <- rep( 0.0, length(lags)) 
	mu <- mean(dt) 
	
	for (i in 1:length(lags))
	{
		lg <- lags[i]
		aut[i] <- sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n 
		#cat( i, lg, aut[i], "\n")
	}
	
	aut 
}


##### A much better, although slower, way to calculate the
##### Integrated Autocorrelation Time (not plots though)
IAT <- function( info, par=0, from=1, to=info$Tr) {

	##-lag/log(GetAutoCorr( info, lag, par=par, from=from, to=to)) 
	
	
	#### we get the desired time series, the parameter and the from - to selection
	if (par>0) {
		if (info$dim > 1)
			dt <-  info$output[from:to, par]
		else
			dt <-  info$output[from:to]
	}
	else
		dt <-  info$Us[from:to] 

	n <- to-from 	
	mu <- mean(dt)  ### with its mean and variance
	s2 <- var(dt) 
	
	### The maximum lag is half the sample size
	maxlag <- max( 3, floor(n/2)) 
	
	#### The gammas are sums of two consecutive autocovariances
	Ga <- rep(0,2)  ## two consecutive gammas
	
	lg <- 0 
	Ga[1] <- s2  #sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n 
	lg <- 1 
	Ga[1] <- Ga[1] + sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n 
	
	m <- 1 
	lg <- 2*m 
	Ga[2] <- sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/(n-lg) 
	lg <- 2*m+1 
	Ga[2] <- Ga[2] + sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/(n-lg) 

	IAT <- Ga[1]/s2  ### Add the autocorrelations
	
	
	### RULE: while Gamma stays positive and decreasing
	while  ((Ga[2] > 0.0) & (Ga[2] < Ga[1])) {
		m <- m+1 
		if (2*m+1 > maxlag) {
			cat("Not enough data, maxlag=", maxlag, "\n") 
			break 
		}
		Ga[1] <- Ga[2] 
		
		lg <- 2*m 
		Ga[2] <- sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n 
		lg <- 2*m+1 
		Ga[2] <- Ga[2] + sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n 
		
		IAT <- IAT + Ga[1]/s2 
	}
	 
	IAT <- -1 + 2*IAT   ##Calculates the IAT from the gammas
	
	#cat("IAT: IAT=", IAT, ", last lag=", 2*m+1, ", last Gamma=", Ga[1], "\n") 

	IAT 
	
}
	




 














