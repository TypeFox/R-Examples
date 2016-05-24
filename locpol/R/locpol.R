#	kernels.R
#	- useful kernels...
#	NOTES:
#	- list should increase.. 
#	- 'gaussK' is very slowly to do computations, better a compact supported kernel...
#	ERRORS:
	
#
#	Kernels (ok)
#

# Gaussian kernel
gaussK <- function(x) dnorm(x,0,1)
## OK
attr(gaussK,"RK") <- 1/( 2 * sqrt(pi) )
attr(gaussK,"RdK") <- 1/( 4 * sqrt(pi) )
attr(gaussK,"mu0K") <- 1
attr(gaussK,"mu2K") <- 1
attr(gaussK,"K4") <- 0.1994711
attr(gaussK,"dom") <- c(-Inf,Inf)

# Gaussian kernel(locfit version)
gaussKlf <- function(x) exp(-(2.25 * x)^2);	## va ahí el signo '-'
## OK
attr(gaussKlf,"RK") <- 0.557029	
attr(gaussKlf,"RdK") <- 2.81996
attr(gaussKlf,"mu0K") <- 0.787757
attr(gaussKlf,"mu2K") <- 0.0778032
attr(gaussKlf,"K4") <- 0.2444259
attr(gaussKlf,"dom") <- c(-Inf,Inf)

# Squared kernel
SqK <- function(x) ifelse( abs(x) <= 1 , 1 , 0 );
## OK
attr(SqK,"RK") <- 2.0	
attr(SqK,"RdK") <- 0.0
attr(SqK,"mu0K") <-	2.0
attr(SqK,"mu2K") <- 0.666667
attr(SqK,"K4") <- NA

# Triangle kernel
TrianK <- function(x) ifelse( abs(x) <= 1 , (1 - abs(x)) , 0 );
## OK
attr(TrianK,"RK") <- 0.666667	
attr(TrianK,"RdK") <- 2
attr(TrianK,"mu0K") <-	1.
attr(TrianK,"mu2K") <- 0.166667
attr(TrianK,"K4") <- 0.4793729
attr(TrianK,"dom") <- c(-1,1)
 	
# Epanechnikov kernel
EpaK <- function(x) ifelse( abs(x) <= 1 , 3/4*(1 - x^2) , 0 );
## OK
attr(EpaK,"RK") <- 0.6
attr(EpaK,"RdK") <- 1.5		
attr(EpaK,"mu0K") <- 1.
attr(EpaK,"mu2K") <- 0.2
attr(EpaK,"K4") <- 0.4337657
attr(EpaK,"dom") <- c(-1,1)

# Epanechnikov kernel
Epa2K <- function(x) ifelse( abs(x) <= 1 , 3/(4*sqrt(5))*(1 - x^2/5) , 0 );
## OK
attr(Epa2K,"RK") <- 0.1968
attr(Epa2K,"RdK") <- 0.012
attr(Epa2K,"mu0K") <- 0.626099
attr(Epa2K,"mu2K") <- 0.196774
attr(Epa2K,"K4") <- NA
attr(Epa2K,"dom") <- c(-1,1)

# Quartic kernel
QuartK <- function(x) ifelse( abs(x) <= 1 , 15/16*(1 - x^2)^2 , 0 );
## OK
attr(QuartK,"RK") <- 0.714286
attr(QuartK,"RdK") <- 2.14286
attr(QuartK,"mu0K") <- 1.
attr(QuartK,"mu2K") <- 0.142857
attr(QuartK,"K4") <- 0.5164128
attr(QuartK,"dom") <- c(-1,1)

biweigK <- function(x) ifelse( abs(x) <= 1 , (1 - x^2)^2 , 0 );
## OK
attr(biweigK,"RK") <- 0.812698;
attr(biweigK,"RdK") <- 2.4381;
attr(biweigK,"mu0K") <- 1.06667;
attr(biweigK,"mu2K") <- 0.152381;
attr(biweigK,"K4") <- 0.6685162
attr(biweigK,"dom") <- c(-1,1)

# Triweigth kernel
TriweigK <- function(x) ifelse( abs(x) <= 1 , 35/32*(1 - x^2)^3 , 0 );
## OK
attr(TriweigK,"RK") <- 0.815851
attr(TriweigK,"RdK") <- 3.18182
attr(TriweigK,"mu0K") <- 1.
attr(TriweigK,"mu2K") <- 0.111111
attr(TriweigK,"K4") <- 0.5879012
attr(TriweigK,"dom") <- c(-1,1)

tricubK <- function(x) ifelse( abs(x) <= 1 , 70/81*(1 - abs(x)^3)^3 , 0 );
## OK
attr(tricubK,"RK") <- 0.708502;
attr(tricubK,"RdK") <- 2.24599;
attr(tricubK,"mu0K") <- 1;
attr(tricubK,"mu2K") <- 0.144033;
attr(tricubK,"K4") <- 0.5879012
attr(tricubK,"dom") <- c(-1,1)

tricubKlf <- function(x) ifelse( abs(x) <= 1 , (1 - abs(x)^3)^3 , 0 );
## OK, (locfit version)
attr(tricubKlf,"RK") <- 0.94867;
attr(tricubKlf,"RdK") <- 3.00733
attr(tricubKlf,"mu0K") <- 1.15714;
attr(tricubKlf,"mu2K") <- 0.166667;
attr(tricubKlf,"K4") <- NA
attr(tricubKlf,"dom") <- c(-1,1)

# cosin kernel
CosK <- function(x) ifelse( abs(x) <= 1 , pi/4*cos(pi*x/2) , 0 );
## OK
attr(CosK,"RK") <- 0.6168
attr(CosK,"RdK") <- 1.52202		
attr(CosK,"mu0K") <- 1.	
attr(CosK,"mu2K") <- 0.189431
attr(CosK,"K4") <- 0.4464387
attr(CosK,"dom") <- c(-1,1)


.kernelList <- c( "gaussK", "EpaK", "Epa2K", "TrianK", "QuartK", "biweigK", 
                  "TriweigK", "tricubK","CosK",  "SqK" )


#	R(K) wrapper
RK <- function(K) return( attr(K,"RK") )

#	R(K') wrapper
RdK <- function(K) return( attr(K,"RdK") )

#	mu2K wrapper
mu2K <- function(K) return( attr(K,"mu2K") )

#	mu2K wrapper
mu0K <- function(K) return( attr(K,"mu0K") ) 

#	K4 wrapper
K4 <- function(K) return( attr(K,"K4") )

#	dom wrapper
dom <- function(K) return( attr(K,"dom") ) 



##
##	Construcción del kernel evquivalente
##
equivKernel <- function(kernel,nu,deg,
			lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],subdivisions=25)
##OK
#	p = local polinomial order.
#	nu = which of the equivalent kernels do you want.
#	kernel = kernel.
#	NOTES:
#		Al menos para p=1,1 funciona
{
	mu <- array( dim=c(deg+1,deg+1) )
	for(i in 0:deg)
		for(j in 0:deg)
			mu[i+1,j+1] <- computeMu(i+j,kernel,lower,upper,subdivisions)
	invMu <- solve(mu)
	res <- function(t)
	{
		return( (invMu %*% outer(0:deg,t,function(x,y) y^x))[nu+1,]*kernel(t) )
		##return( (invMu %*% outer(nu+1,t,function(x,y) y^x))*kernel(t) )
	}
	return( res )
}

##
##	Construcción C_{\nu p}(K)
##
cteNuK <- function(nu,p,kernel,
			lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],subdivisions=25)
##
#	nu = ;
#	p = ;
#	kernel = ;
#	NOTES:
#		Al menos para EpaK, y gaussK y TriweigK con p=0,1,... la cosa funciona...
#		Es decir  funciona...
{
	
	if( p==0 || p==1 ){
		a <- RK(kernel)
		b <- mu2K(kernel)
		return( ((p+1)^2*(2*nu+1)*a/(2*(p+1-nu)*b^2))^(1/(2*p+3)) )
	}else{
		f <- 1
		for(i in 2:(p+1)) 
			f <- f * i
		eqKern <- equivKernel(kernel,nu,p,lower,upper,subdivisions)
		a <- computeRK(eqKern,lower,upper,subdivisions)
		b <- computeMu(p+1,eqKern,lower,upper,subdivisions)
		return( (f^2*(2*nu+1)*a/(2*(p+1-nu)*b^2))^(1/(2*p+3)) )
	}
}

##
##	Construcción adj_{\nu p}(K), para ajustar el bandwidth a la estimación de 
##	derivadas.
##
adjNuK <- function(nu,p,kernel,
			lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],subdivisions=25)
##
#	nu = ;
#	p = ;
#	kernel = ;
#	RETURNS
#		Just 'cteNuK(nu,p,kernel)/cteNuK(0,p,kernel)', useful for derivative 
#	bandwidth computations.
{
	a <- cteNuK(nu,p,kernel,lower,upper,subdivisions)
	b <- cteNuK(0,p,kernel,lower,upper,subdivisions)
	return( a/b )
}


##
##	Some kernel miscelanea functions
##
Kconvol <- function(kernel,
			lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],subdivisions=25)
##
#	RETURNS
#		given Kernel Convolution function. Computed by Num. Int.
#	kernel = Kernel use in computation
#	order =	currently onder two is the only one computed...
#	NOTES:
#		Use 'integrate' from R, seems to be quick and accurate!!!
{
	ff <- function(x) integrate(
			function(u) return( kernel(u) * kernel(x-u) ),
			lower,upper,subdivisions=subdivisions)$value
	res <- function(x) sapply(x,ff)
	return( res )
}


computeRK <- function(kernel,lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],subdivisions=25)
##
#
#	kernel = Kernel use in computation
#	lower = lower limit of interval.
#	upper = upper limit of interval.
#	subdivisions = number of subdivisions used to compute integral.
#		Computes RK for the given kernel.
{
	f <- function(u) kernel(u)^2
	return( integrate(f,lower,upper,subdivisions=subdivisions)$value )
}

computeK4 <- function(kernel,
			lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],subdivisions=25)
##
#
#	kernel = Kernel use in computation
#	lower = lower limit of interval.
#	upper = upper limit of interval.
#	subdivisions = number of subdivisions used to compute integral.
#		Computes K4 for the given kernel, the autoconvolution fo the
#	autoconvolution for the given kernel .
{
	f <- Kconvol(kernel,lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],subdivisions=25)
	ff <- function(x) f(x)^2
	return( integrate(ff,2*lower,2*upper,subdivisions=subdivisions)$value )
}



# computeRdK <- function(kernel,lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],subdivisions=25)
# {
# 	dd <- 
# 	f <- function(u) kernel(u)^2
# 	return( integrate(f,lower,upper,subdivisions)$value )
# }

computeMu0 <- function(kernel,
				lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],subdivisions=25)
##
#
#	kernel = Kernel use in computation
#	lower = lower limit of interval.
#	upper = upper limit of interval.
#	subdivisions = number of subdivisions used to compute integral.
#		Computes mu_0 for the given kernel.
{
	return( integrate(kernel,lower,upper,subdivisions=subdivisions)$value )
}

computeMu <- function(i,kernel,
				lower=dom(kernel)[[1]],upper=dom(kernel)[[2]],subdivisions=25)
##
#
#	kernel = Kernel use in computation
#	lower = lower limit of interval.
#	upper = upper limit of interval.
#	subdivisions = number of subdivisions used to compute integral.
#		Computes mu_i for the given kernel.
{
	f <- function(u) u^i*kernel(u)
	return( integrate(f,lower,upper,subdivisions=subdivisions)$value )
}

# locpolClass.R
#  local polynomial class
# NOTES:
#  -- Bias estimations !! ??
#  -- Keep track of errors... So you'll be able to give always
#	a suggestion in case errors appears
#  -- Imporve error messages !!
# ERRORS:




locpol <- function(formula, data, weig=rep(1,nrow(data)),
                   bw=NULL, kernel=EpaK, deg=1, xeval=NULL,xevalLen=100)
##	
{
    ##  checking
    stopifnot(nrow(data)==length(weig))
    ## compute result
    res <- list()
    res$bw <- bw
    res$KName <- match.call()
    res$kernel <- kernel
    res$deg <- deg
    res$xeval <- xeval
    ## get info from formula
    res$mf <- model.frame(formula,data)
    datCla <- attr(attr(res$mf, "terms"),"dataClasses") 
    varNames <- names(datCla)[datCla=="numeric"] 
    stopifnot(length(varNames)==2) 
    res$Y <- varNames[1] 
    res$X <- varNames[2] 
    ##	sort x's
    xo <- order(res$mf[,res$X])
    res$mf <- res$mf[xo,]
    res$weig <- weig[xo]
    ## xeval
    if( is.null(xeval) )
        res$xeval <- seq(min(res$mf[,res$X]),max(res$mf[,res$X]),len=xevalLen)
      else 
        res$xeval <- sort(xeval)
    if( is.null(res$bw) ) 
        res$bw <- regCVBwSelC(data[,res$X],data[,res$Y],res$deg,
                    res$kernel,res$weig)
    ## regression estimation
    res$lpFit <- locPolSmootherC(res$mf[,res$X], res$mf[,res$Y], res$xeval, 
                    res$bw, res$deg, res$kernel, DET = TRUE, res$weig )
    names(res$lpFit)[] <- c(res$X,res$Y,paste(res$Y,1:deg,sep=""),"xDen")
    res$lpFit$xDen <- res$lpFit$xDen^(1/(deg+1))/(nrow(data)*res$bw)
    ## CI comp.(##Should depned on nu, to be able to)
    nu <- 0
    res$CIwidth <- computeRK(equivKernel(kernel,nu,deg),
                    lower=dom(res$kernel)[[1]], upper=dom(res$kernel)[[2]], 
                    subdivisions = 25) * factorial(nu)^2 
    res$CIwidth <- res$CIwidth / ( nrow(data)*res$bw )
    ## residuals
    res$residuals <- res$mf[,res$Y]-locLinSmootherC(res$mf[,res$X],
                res$mf[,res$Y], res$mf[,res$X],res$bw, res$kernel, 
                res$weig)$beta0 
    ## variance estimation
    res$lpFit$var <- locCteSmootherC(res$mf[,res$X], res$residuals^2, 
                res$xeval,1.2*res$bw, res$kernel, res$weig)$beta0                 
    ## setupclass
    class(res) <- "locpol"
    return(res)
}


residuals.locpol <- function(object,...)
##	
{
    return( object$residuals )
}


fitted.locpol <- function(object,deg=0,...)
##	
{
    stopifnot(object$deg>=deg)
    return( object$lpFit[,2+deg] )
}


summary.locpol <- function(object,...)
##	
{
    cat("\nKernel = \n\t")
    print( body(object$kernel) )
    cat("\n")
    print(  data.frame(n=nrow(object$mf), deg=object$deg, bw=object$bw,
            ase=mean(resid(object)^2), row.names=" ") )
}


print.locpol <- function(x,...)
##	
{
    summary.locpol(x)
}


confInterval <- function(x)
##
{
    plot(x$mf[,x$X],x$mf[,x$Y],pch="+",main="95% Conf. Int. for x-Points",
			xlab=x$X,ylab=x$Y )
    points(x$lpFit[,x$X],x$lpFit[,x$Y],type="l",col="darkgreen")
    dev <- sqrt(x$CIwidth * x$lpFit$var/x$lpFit$xDen)
    points(x$lpFit[,x$X],x$lpFit[,x$Y]+2*dev,type="l",col="green")
    points(x$lpFit[,x$X],x$lpFit[,x$Y]-2*dev,type="l",col="green")
}


plot.locpol <- function(x,...)
##	
{
    par(ask=TRUE,cex=.6)
    plot(x$mf[,x$X],x$mf[,x$Y],pch="+",main="Data and Regres.",
        xlab=x$X,ylab=x$Y )
    points(x$lpFit[,x$X],x$lpFit[,x$Y],type="l",col="blue")
    ym <- max(x$lpFit$xDen)
    plot(x$lpFit[,c(x$X,"xDen")], main="X dens.",type="l",ylim=c(0,1.25*ym), 
        xlab=x$X,ylab="den" )
    ym <- max(x$lpFit$var)
    plot(x$lpFit[,c(x$X,"var")], main="Var.",type="l",ylim=c(0,1.25*ym),
        xlab=x$X,ylab="var" )
    confInterval(x)
    par(ask=FALSE)
}

#
#	locpol.R
#		Interface para el polinomio local en C
#	NOTAS:
#	- Añadir todos los kernels...
#	- Ojito, devuelve 0 si det(X^TWX)=0!!
#	ERRORES:


##  old s3 classes stuff
# .First.lib <- function(lib, pkg) {
#   library.dynam("locpol", pkg, lib)
# }


##  s4 classes dim lib load stuff
.onLoad <- function(lib, pkg) {
  library.dynam("locpol", pkg, lib)
}


.maxEvalPts <- 5000
.lokestOptInt <- c(0.0005,1.5)

selKernel <- function(kernel)
{
	if( RK(kernel)==0.6 )	## Epanechnikov kernel(ok R)
		return(1)
	else if ( RK(kernel)==0.1968 )	## 2nd. Epanechnikov kernel(ok R)
		return(2)
	else if ( RK(kernel)==0.666667 )	## Triangle kernel(ok R)
		return(3)
	else if ( RK(kernel)==0.714286 )	## Quartic kernel(ok R)
		return(4)
	else if ( RK(kernel)==0.812698 )	## biweight kernel(ok R)
		return(5)
	else if ( RK(kernel)==0.815851 )	## Triweigth kernel(ok R)
		return(6)
	else if ( RK(kernel)==0.708502 )	## tricube kernel(ok R)
		return(7)
	else if ( RK(kernel)==0.6168 )		## cosin kernel(ok R)
		return(9)	
  else if ( RK(kernel)==2.0 )		    ## Square kernel(ok R)
		return(10)
	else 	## gaussian kernel
		return(0)	
}

locWeightsEval <- function(lpweig,y)
##OK
#
# 	lpweig = an x-data frame.
# 	y = eval. points.
{
	stopifnot(	ncol(lpweig) == length(y), nrow(lpweig)<.maxEvalPts )
    return( lpweig %*% y )
}

locWeightsEvalC <- function(lpweig,y)
##OK
#
# 	lpweig = an x-data frame.
# 	y = eval. points.
{
	stopifnot(	ncol(lpweig) == length(y), nrow(lpweig)<.maxEvalPts )
	res <- .C("locWeightsEval",
		as.double(lpweig),as.integer(nrow(lpweig)),
		as.double(y),as.integer(length(y)),
		res=double(nrow(lpweig)), PACKAGE="locpol"
		)
	#rownames(res$res) <- rownames(lpweig)
    return( res$res )
}

simpleSmootherC <- function(x,y,xeval,bw,kernel,weig = rep(1,length(y)))
##OK
#
# 	x, y = x and y vectors.
# 	xeval = eval. points.
# 	bw = bandwidth.
# 	weig = vector of weigths for observations.
# 	kernel = kernel used in estimation.
#		An interface to C function 'simpleSmoother'.
#		Computes m(x)f(x)mu_0.
{
	stopifnot(length(x)==length(y), length(y)==length(weig), bw > 0)
	Ktype <- selKernel(kernel)
	reg <- .C("simpleSmoother",
		as.double(xeval),as.integer(length(xeval)),
		as.double(x),as.double(y),as.double(weig),
		as.integer(length(y)),as.double(bw), as.integer(Ktype),
		res = double(length(xeval)), PACKAGE="locpol"
		)$res
	res <- data.frame(x = xeval, reg = reg)
    return( res )
}


PRDenEstC <- function(x,xeval,bw,kernel,weig = rep(1,length(x)))
##OK
#
# 	x = an x vector.
# 	xeval = eval. points.
# 	bw = bandwidth.
# 	weig = vector of weigths for observations.
# 	kernel = kernel used in estimation.
#		An interface to C function 'parzenRossen'.
#		Computes the density of x
{
	stopifnot(length(x)==length(weig), bw > 0)
	Ktype <- selKernel(kernel)
	den <- .C("parzenRossen",
    	as.double(xeval),as.integer(length(xeval)),
     	as.double(x),as.double(weig),as.integer(length(x)),
		  as.double(bw), as.integer(Ktype),
     	res = double(length(xeval)), PACKAGE="locpol"
		)$res
	res <- data.frame(x = xeval, den = den)
    return( res )
}


locCteSmootherC <- function(x,y,xeval,bw,kernel,weig = rep(1,length(y)))
##OK
#
# 	x, y = x and y vectors.
# 	xeval = eval. points.
# 	bw = bandwidth.
# 	weig = vector of weigths for observations.
# 	kernel = kernel used in estimation.
#		An interface to C function 'locCteSmoother'.
#		Computes m(x), Nadaraya-Watson estimator, as a by--product we 
#		got the marginal density for X.
{
	stopifnot(length(x)==length(y), length(y)==length(weig), bw > 0,
				length(xeval)<.maxEvalPts)
	Ktype <- selKernel(kernel)
	res <- .C("locCteSmoother",
		as.double(xeval),as.integer(length(xeval)),
		as.double(x),as.double(y),as.double(weig),
		as.integer(length(y)),as.double(bw), as.integer(Ktype),
		den=double(length(xeval)),res=double(length(xeval)), PACKAGE="locpol"
		)
	res <- data.frame(x = xeval, beta0 = res$res, den=res$den)
    return( res )
}

locCteWeightsC <- function(x,xeval,bw,kernel,weig = rep(1,length(x)))
##OK
#
# 	x = x vector.
# 	xeval = eval. points.
# 	bw = bandwidth.
# 	weig = vector of weigths for observations.
# 	kernel = kernel used in estimation.
#		An interface to C function 'locCteSmoother'.
#		Computes m(x), Nadaraya-Watson estimator, as a by--product we 
#		got the marginal density for X.
{
	stopifnot(length(x)==length(weig), bw > 0, length(xeval)<.maxEvalPts)
	Ktype <- selKernel(kernel)
	res <- .C("locCteWeights",
		as.double(xeval),as.integer(length(xeval)),
		as.double(x),as.double(weig),as.integer(length(x)),
		as.double(bw), as.integer(Ktype),
		den=double(length(xeval)),res=mat.or.vec(length(xeval),length(x)),
		PACKAGE="locpol"
		)
	rownames(res$res) <- xeval
	colnames(res$res) <- 1:length(x)
    return( list(den=res$den,locWeig=res$res) )
}

locLinSmootherC <- function(x,y,xeval,bw,kernel,weig = rep(1,length(y)))
##OK
#
# 	x, y = x and y vectors.
# 	xeval = eval. points.
# 	bw = bandwidth.
# 	weig = vector of weigths for observations.
# 	kernel = kernel used in estimation.
#		An interface to C function 'locLinSmoother'.
#		Computes m(x), Local linear estimator, as a by--product we 
#		got f(x)^2, the squared marginal density for X.
{
	stopifnot(length(x)==length(y), length(y)==length(weig), bw > 0,
				length(xeval)<.maxEvalPts)
	Ktype <- selKernel(kernel)
	res <- .C("locLinSmoother",
		as.double(xeval),as.integer(length(xeval)),
		as.double(x),as.double(y),as.double(weig),
		as.integer(length(y)),
		as.double(bw), as.integer(Ktype),den=double(length(xeval)),
		beta0=double(length(xeval)),beta1=double(length(xeval)), 
		PACKAGE="locpol"
		)
	res <- data.frame(x=xeval, beta0=res$beta0, beta1=res$beta1, den=res$den) 
    return( res )
}

locLinWeightsC <- function(x,xeval,bw,kernel,weig = rep(1,length(x)))
##OK
#
# 	x = x vector.
# 	xeval = eval. points.
# 	bw = bandwidth.
# 	weig = vector of weigths for observations.
# 	kernel = kernel used in estimation.
#		An interface to C function 'locCteSmoother'.
#		Computes m(x), Nadaraya-Watson estimator, as a by--product we 
#		got the marginal density for X.
{
	stopifnot(length(x)==length(weig), bw > 0, length(xeval)<.maxEvalPts)
	Ktype <- selKernel(kernel)
	res <- .C("locLinWeights",
		as.double(xeval),as.integer(length(xeval)),
		as.double(x),as.double(weig),
		as.integer(length(x)),
		as.double(bw), as.integer(Ktype),
		den=double(length(xeval)),res=mat.or.vec(length(xeval),length(x)), 
		PACKAGE="locpol"
		)
	rownames(res$res) <- xeval
	colnames(res$res) <- 1:length(x)
    return( list(den=res$den,locWeig=res$res) )
}

locCuadSmootherC <- function(x,y,xeval,bw,kernel,weig = rep(1,length(y)))
##OK
#
# 	x, y = x and y vectors.
# 	xeval = eval. points.
# 	bw = bandwidth.
# 	weig = vector of weigths for observations.
# 	kernel = kernel used in estimation.
#		An interface to C function 'locCuadSmoother'.
#		Computes m(x), Local cuadratic estimator, as a by--product we 
#		got f(x)^3, the squared marginal density for X.
{
	stopifnot(length(x)==length(y), length(y)==length(weig), bw > 0,
				length(xeval)<.maxEvalPts)
	Ktype <- selKernel(kernel)
	res <- .C("locCuadSmoother",
		as.double(xeval),as.integer(length(xeval)),
		as.double(x),as.double(y),as.double(weig),
		as.integer(length(y)),as.double(bw), as.integer(Ktype),
		den=double(length(xeval)), beta0=double(length(xeval)),
		beta1=double(length(xeval)), beta2=double(length(xeval)), 
		PACKAGE="locpol"
		)
	res <- data.frame(x = xeval, beta0=res$beta0, beta1=res$beta1, 
						beta2=res$beta2, den=res$den) 
    return( res )
}

locPolSmootherC <- function(x,y,xeval,bw,deg,kernel,DET=FALSE,
												weig=rep(1,length(y)))
##OK
#
# 	x, y = x and y vectors.
# 	xeval = eval. points.
# 	bw = bandwidth.
# 	weig = vector of weigths for observations.
# 	kernel = kernel used in estimation.
#	DET = flag to compute the determinant of t(X)WX(X local desing matrix) 
#		An interface to C function 'locPolSmootherC'.
#		Computes m(x), Local pol. estimator, as a by--product we 
#		got f(x)^(p+1), the squared marginal density for X.
{
	stopifnot(length(x)==length(y), length(y)==length(weig), bw > 0,
				length(xeval)< .maxEvalPts,0<=deg && deg<10)
	Ktype <- selKernel(kernel)
	res <- .C("locPolSmoother",
		as.double(xeval), as.integer(length(xeval)),
		as.double(x), as.double(y), as.double(weig),
		as.integer(length(y)), as.double(bw), as.integer(deg), 
		as.integer(Ktype), as.integer(DET),
		den=double(length(xeval)), 
		beta=mat.or.vec(length(xeval),(deg+1)), PACKAGE="locpol"
		)
	res <- data.frame(x=xeval, beta=matrix(res$beta,ncol=deg+1),res$den) 
	colnames(res) <- c("x",paste("beta",0:deg,sep=""),"den")
	if( !DET ) res$den <- NULL
    return( res )
}

simpleSqSmootherC <- function(x,y,xeval,bw,kernel)
##OK
#
# 	x, y = x and y vectors.
# 	xeval = eval. points.
# 	bw = bandwidth.
# 	kernel = kernel used in estimation.
#		An interface to C function 'simpleSmoother'. It returns
#	E[y^2|x]f(x)R(K), can be useful to estimate variance...
{
	stopifnot(length(x)==length(y), bw > 0)
	Ktype <- selKernel(kernel)
	reg <- .C("simpleSqSmoother",
		as.double(xeval),as.integer(length(xeval)),
		as.double(x),as.double(y),
		as.integer(length(y)),
		as.double(bw), as.integer(Ktype),
     	res = double(length(xeval)), PACKAGE="locpol")$res
	res <- data.frame(x = xeval, reg = reg)
    return( res )
}

## 
## 	Leave One out computations.
## 
looLocPolSmootherC <- function(x,y,bw,deg,kernel,weig=rep(1,length(y)),DET=FALSE)
##OK
#	x, y = x and y vectors.
# 	xeval = eval. points.
# 	bw = bandwidth.
# 	weig = vector of weigths for observations.
# 	kernel = kernel used in estimation.
#	DET = flag to compute the determinant of t(X)WX(X local desing matrix) 
#		An interface to C function 'looLocPolSmootherC'. Estimation of 
#		m(x_i) i=1,...n without ith observation (x_i,y_i).
#		Computes m(x), Local pol. estimator, as a by--product we 
#		got f(x)^(p+1), the squared marginal density for X.
{
	stopifnot(length(x)==length(y), length(y)==length(weig), bw > 0,
				length(x)<.maxEvalPts,0<=deg && deg<10)
	Ktype <- selKernel(kernel)
	res <- .C("looLocPolSmoother",
		as.double(x), as.double(y), as.double(weig),
		as.integer(length(y)), as.double(bw), as.integer(deg), 
		as.integer(Ktype), as.integer(DET),
		den=double(length(y)), beta=mat.or.vec(length(y),(deg+1)), 
		PACKAGE="locpol"
		)
	res <- data.frame(x=x, beta=matrix(res$beta,ncol=deg+1),res$den) 
	colnames(res) <- c("x",paste("beta",0:deg,sep=""),"den")
    return( res )
}

##
##	Bandwidth selection procedures...
##

denCVBwSelC <- function(x,kernel=gaussK,weig=rep(1,length(x)),
							interval=.lokestOptInt)
##OK
#
# 	x = x vector
# 	xeval = eval. points.
# 	weig = vector of weigths for observations.
# 	kernel = kernel used in estimation.
#		Computes Cross Validation bandwidth selector for the 
#	density estimator...	
{
	stopifnot(length(x)==length(weig))
	Ktype <- selKernel(kernel)
	cvFunc <-	function(h) .C("denCVBwEval",as.double(h),
		as.double(x),as.double(weig),as.integer(length(x)),
		as.integer(Ktype),res = double(1), PACKAGE="locpol")$res
    return( optimise(cvFunc,interval)$minimum )
}

regCVBwSelC <- function(x,y,deg,kernel=gaussK,weig=rep(1,length(y)),
							interval=.lokestOptInt)
##OK
#
# 	dat = an x-data frame.
# 	xeval = eval. points.
# 	weig = vector of weigths for observations.
# 	kernel = kernel used in estimation.
#		Computes Cross Validation bandwidth selector for the 
#	regression function estimator...
{
	stopifnot(length(x)==length(y), length(y)==length(weig), 0<=deg && deg<10 )
	Ktype <- selKernel(kernel)
	cvFunc <-	function(h) .C("regCVBwEvalB",as.double(h),
     	as.double(x), as.double(y), as.double(weig),
		as.integer(length(y)), as.integer(deg),
		as.integer(Ktype),res = double(1), PACKAGE="locpol")$res
    return( optimise(cvFunc,interval)$minimum )
}



# locPolWeights.R
#
# NOTAS:
# ERRORES:


locPolWeights <- function (x, xeval, deg, bw, kernel, 
                           weig=rep(1,length(x))) 
# Suggestion from munevvere@hotmail.com 
{
  stopifnot( length(x)==length(weig), bw>0, deg>=0 )
  den <- array( dim=length(xeval),dimnames=list(x=1:length(xeval)) )
  res <- array( dim=c(length(xeval),deg+1,length(x)),
                dimnames=list(x=1:length(xeval),deg=0:deg,xData=1:length(x)))
  for(i in 1:length(xeval)) 
  {
    xx <- xeval[i]
    dMat <- (x-xx)/bw
    xx <- outer(dMat, 0:deg, function(a, b) a^b)
    w <- diag(kernel(dMat)*weig)
    aux <- t(xx) %*% w
    sMat <- aux %*% xx
    den[i] <- det(sMat)
    res[i,,] <- solve(sMat,aux)/bw^(0:deg)
  }
  invisible( list(den=den*bw^deg,locWeig=res[,1,],allWeig=res) ) 
}





# mvNoparEst.R
#   Multivariate density and regression.
# NOTAS:
#   . Check normalization constant !!
# ERRORS:


##
## some bivariate kernels, 
##i
epaK2d <- function(x) (2/pi) * (1-x[[1]]^2 -x[[2]]^2)*( (x[[1]]^2+x[[2]]^2)< 1 )
gauK2d <- function(x) (1/(2*pi)) * exp(-0.5*(x[[1]]^2 + x[[2]]^2) )   ##check this, not sure !!!


##
##  ... a bandwidth selector ??
##
##myDR <- function(x) diff(range(x))
mayBeBwSel <- function(X,prop=.45)
##
# X = data(data.frame or matrix), only first two cols are used.
{ 
  ## uses 'prop' of the range of points.
  return( diag(prop*c(diff(range(X[,1])), diff(range(X[,2]))) ))
}


##
## density estimation
##
bivDens <- function(X,weig=rep(1,nrow(X)),K=epaK2d,H=mayBeBwSel(X) )
##
# X = data(data.frame or matrix), only first two cols are used.
# weig = weigths for each point.
# K = bivariate kernel.
# H = bandwidht matrix, by default
{
  stopifnot( length(weig)==nrow(X) )
  weig <- weig/sum(weig)
  ## compute kernel dens. est. for each data point. 
  denEst <- function(xa)
  {
    we <- apply(X[,1:2],1,function(x)K( (x-xa) %*% solve(H)))
    return( sum(we*weig) )
  }
  ## apply it to each data point
  res <- list( X=X[,1:2], H=H, estFun=denEst )
  class(res) <- "bivNpEst"
  return(res)
}


##
## Regression
##
bivReg <- function(X,Y,weig=rep(1,nrow(X)),K=epaK2d,H=mayBeBwSel(X) )
##
# X = X data(data.frame or matrix), only first two cols are used.
# Y = Y data,a vector. 
# weig = weigths for each point.
# K = bivariate kernel.
# H = bandwidht matrix.
{
  stopifnot( length(weig)==nrow(X), length(Y)==nrow(X))
  weig <- weig/sum(weig)
  ## compute NW est. for each data point. 
  regEst <- function(xa)
  {
    we <- apply(X[,1:2],1,function(x)K( (x-xa) %*% solve(H)))
    return( sum(we*weig*Y)/sum(we*weig) )
  }
  ## apply it to each data point
  res <- list( X=X[,1:2], Y=Y, H=H, estFun=regEst )
  class(res) <- "bivNpEst"
  return(res)
}


predict.bivNpEst <- function(object,newdata=NULL,...)
##
# object = "bivNpEst" object
# newdata = new data where the density should be predicted 
{
  if( is.null(newdata) ) 
    newdata <- object$X
  else{ 
   stopifnot( is.data.frame(newdata),all(names(object$X) %in% names(newdata)) ) 
  }
  return( apply(newdata[,names(object$X)],1,function(xa) object$estFun(xa)) )
}


plotBivNpEstOpts <- list(  
  pathLen=10, phi=30,  theta=15, r=sqrt(3), 
  ticktype="detailed", nticks=3, main="NP est." 
  )


plot.bivNpEst <- function(x,...)
##
# x = a 'bivNpEst' object
# ... = 'persp' options plus 'plRng' and 'pathLen':
#         plRng = matrix with the range for x and y in each row. 
#         pathLen = number of points in the grid 
{
  ## default options
  mc <- match.call()
  parIdx <- na.omit( match(names(mc)[-1],names(plotBivNpEstOpts)) )
  if( length(parIdx)>0 ) 
    plotBivNpEstOpts <- plotBivNpEstOpts[-parIdx] 
  for( i in names(plotBivNpEstOpts) )
      mc[[i]] <- plotBivNpEstOpts[[i]]
  if( is.null(mc$plRng) ) 
      plRng <- rbind(range(x$X[,1]),range(x$X[,2]))
    else{
      plRng <- eval(mc$plRng)
      mc$plRng <- NULL
  }
  pathLen <- eval(mc$pathLen)
  mc$pathLen <- NULL
  if( is.null(mc$xlab) ) 
      mc$xlab <- names(x$X)[1]
  if( is.null(mc$ylab) ) 
      mc$ylab <- names(x$X)[2]
  if( is.null(mc$zlab) ) 
      mc$zlab <- ifelse( is.null(x$Y),"den","reg")
  ##plot 
  mc[[1]] <- persp
  mc[["x"]] <- seq(from=plRng[[1,1]],to=plRng[[1,2]],len=pathLen)
  mc[["y"]]<- seq(from=plRng[[2,1]],to=plRng[[2,2]],len=pathLen)
  dd <- expand.grid(mc[["x"]],mc[["y"]])
  names(dd) <- names(x$X)
  mc[["z"]]<- predict(x, dd)
  dd$den <- mc[["z"]]
  mc[["z"]] <- matrix(mc[["z"]],nrow=pathLen)
  eval(mc)
  return(dd)
}









#
#	thumbBw.R
#		Rule of thumb bandwith selector for loc. pol.
#	NOTES:
#		Pesos no implementados...
#	ERRORS:

##source('kernels.R')

compDerEst <- function(x,y,p,weig=rep(1,length(y)))
##OK
#	x = x data.
#	y = y data.
#	p = loc. pol. degree.
#	weig = weigths used in the estimation
{
	##	compute global pol. regression
	xnam <- paste("x^", 2:(p+3), sep="")
	xnam <- paste("I(",xnam, ")")
	fmla <- as.formula(paste("y ~ 1 + x + ", paste(xnam, collapse= "+")))
	lmFit <- lm(fmla,data=data.frame(x,y),weights=weig)
	##	compute (p+1)-der 
	cp <- cumprod(1:(p+3))
	coef <- coefficients(lmFit)
	der <- coef[[p+2]]*cp[[p+1]]+coef[[p+3]]*cp[[p+2]]*x+
				coef[[p+4]]*cp[[p+3]]*x^2/2
	##der <- 
	res <- data.frame(x,y,res=residuals(lmFit)*weig,der)
	return( res )
}



thumbBw <- function(x,y,deg,kernel,weig=rep(1,length(y)))
##OK
#	x = x data.
#	y = y data.
#	deg = loc. pol. degree.
#	kernel = kenrel used in estimation.
#	weig = weigths used in the estimation
{
	k <- 3
	##	compute reg. residuals and derivatives
	rd <- compDerEst(x,y,deg,weig)
	denom <- sum( rd$der^2 )	
	numer <- mean( rd$res^2 )
	##constante kernel, eta, p
	cte <- cteNuK(0,deg,kernel,lower=dom(kernel)[[1]],
					upper=dom(kernel)[[2]],subdivisions=100)
	##valor de h
	res <- cte * (numer/denom)^(1/(2*deg+k))
  return( res )
}


pluginBw <- function(x,y,deg,kernel,weig=rep(1,length(y)))
##OK
#	x = x data.
#	y = y data.
#	deg = loc. pol. degree.
#	kernel = kenrel used in estimation.
#	weig = weigths used in the estimation
#	Only valid for odd degrees,  base on Fand&Gijbels book.
{
	stopifnot(deg%%2==1) 
	##	thumb bandwithd
	thBw <- thumbBw(x, y, deg, kernel)
	## compute loo res.
	regComp <- looLocPolSmootherC(x, y, thBw, deg, kernel, weig)
	res <- (y-regComp$beta0)*weig
	numer <- mean( res^2 )
	## compute loo der.
	thBwBis <- thumbBw(x, y, deg+2, kernel)
	derBw <- adjNuK(deg+1,deg+2,kernel,lower=dom(kernel)[[1]],
					upper=dom(kernel)[[2]],subdivisions=50) * thBwBis			
	regCompBis <- looLocPolSmootherC(x, y, derBw, deg+2, kernel, weig)
	cp <- cumprod(1:(deg+2))
	der <- regCompBis[,1+deg+2]/cp[[deg+1]]
	denom <- sum( der^2 )
	##constante kernel, eta, p
	cte <- cteNuK(0,deg,kernel,subdivisions=100)
	##valor de h
	res <- cte * (numer/denom)^(1/(2*deg+3))
	return(res)
}


