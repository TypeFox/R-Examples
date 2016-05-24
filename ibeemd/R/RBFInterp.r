# Created by Mao-Gui Hu (humg@lreis.ac.cn), 2013.10.
# This is the R translation of the MATLAB impletation of "Scattered Data 
#  Interpolation and Approximation using Radial Base Functions" by Alex Chirokov in
#  "http://www.mathworks.com/matlabcentral/fileexchange/10056-scattered-data-interpolation-and-approximation-using-radial-base-functions", 2006

# fit RBF model
rbfcreate <- function(x, y, rbffun='linear', rbfconst=NULL, rbfsmooth=0)
{
	if(!is(x,"matrix"))
	{
		if(is(x,"vector")) x <- matrix(x, ncol=1)
		else x <- as.matrix(x)
	}
	if(!is(y,"matrix"))
	{
		if(is(y,"vector")) y <- matrix(y, ncol=1)
		else y <- as.matrix(y)
	}
	
	# approx. average distance between the nodes
	if(is.null(rbfconst))
	{
		rbfconst <- (prod(apply(x,2,max)-apply(x,2,min))/nrow(x))^(1/ncol(x))
	}
	
	phi <- switch(tolower(rbffun),
			linear = rbfphi_linear,
			cubic = rbfphi_cubic,
			gaussian = rbfphi_gaussian,
			multiquadric = rbfphi_multiquadric,
			thinplate = rbfphi_thinplate)

	rbfAssemble <- function()
	{
		A <- phi(as.matrix(dist(x)), rbfconst)
		diag(A) <- diag(A)-rbfsmooth;
		
		# Polynomial part
		P <- cbind(rep(1,nrow(x)), x)
		A <- rbind(cbind(A,P), cbind(t(P),matrix(0,nrow=ncol(x)+1,ncol=ncol(x)+1)))
	}

	A <- rbfAssemble()
	b <- rbind(y, matrix(0,nrow=ncol(x)+1,ncol=1))
	rbfcoeff <- solve(A,b)
	
	param <- list(x=x, y=y, phi=phi, rbfconst=rbfconst, rbfcoeff=rbfcoeff)
	return(param)
}

rbfinterp <- function(param, xx)
{
	if(!is(xx,"matrix"))
	{
		if(is(xx,"vector")) xx <- matrix(xx, ncol=1)
		else xx <- as.matrix(xx)
	}
	
	nodes <- param$x
	phi <- param$phi
	rbfconst <- param$rbfconst
	rbfcoeff <- param$rbfcoeff	
	
	m <- ncol(nodes)
	n <- nrow(nodes)
	f <- matrix(NA, nrow=nrow(xx), ncol=1)
	if(nrow(nodes) < 1000 || ncol(nodes) != 2)
	{
		for(i in 1:nrow(xx))
		{
			r <- as.matrix(dist(rbind(xx[i,], nodes)))[,1][-1]
			s <- rbfcoeff[n+1] + sum(rbfcoeff[1:n]*phi(r,rbfconst))	
			s <- s+sum(rbfcoeff[(1:m)+n+1]*xx[i,])	# linear part
			f[i,1] <- s
		}
	}
	else
	{
		for(i in 1:nrow(xx))
		{
			#r <- as.matrix(dist(rbind(xx[i,], nodes)))[,1][-1]
			r <- sp::spDistsN1(nodes, xx[i,1,drop=FALSE])
			s <- rbfcoeff[n+1] + sum(rbfcoeff[1:n]*phi(r,rbfconst))	
			s <- s+sum(rbfcoeff[(1:m)+n+1]*xx[i,])	# linear part
			f[i,1] <- s
		}
	}
	
	return(f)
}

rbfcheck <- function(param)
{
	nodes <- param$x
	y <- param$y
	s <- rbfinterp(param, nodes)
	
	cat('RBF Check\n')
	cat(sprintf('max|y - yi| = %e \n', max(abs(s-y))))
}

#**************************************************************************
# Radial Base Functions
#************************************************************************** 
rbfphi_linear <- function(r, const)
{
	u <- r
	return(u)
}

rbfphi_cubic <- function(r, const)
{
	u <- r^3
	return(u)
}

rbfphi_gaussian <- function(r, const)
{
	u <- exp(-0.5*r^2/const^2)
	return(u)
}

rbfphi_multiquadric <- function(r, const)
{
	u <- sqrt(1+r^2/const^2)
	return(u)
}

rbfphi_thinplate <- function(r, const)
{
	u <- r^2*log(r+1)
	return(u)
}

rbfphi_inversemultiquadric <- function(r, const)
{
	u <- 1/sqrt(1+r^2/const^2)
	return(u)
}

rbfphi_inversequadratic <- function(r, const)
{
	u <- 1/(1+r^2/const^2)
	return(u)
}
