#########################################
### Auxiliary functions for penalized ###
### multicategory logistic regression ###
#########################################

####################################################
# The negative log-likelihood and its derivatives: #
####################################################

Lambda.only <- function(XX,Y,Theta)
{
	L <- dim(Theta)[2]
	tmp0 <- XX%*%Theta
	tmp0 <- tmp0 - apply(tmp0,MARGIN=1,FUN=max)%*%t(rep(1,L))
	Lambda <- sum(log(rowSums(exp(tmp0)))) -
		sum(diag(1,L,L)[Y,]*tmp0)
	return(Lambda)
}

Lambda.full <- function(XX,Y,Theta)
{
	# Auxiliary objects:
	dd <- dim(Theta)[1]
	L <- dim(Theta)[2]
	n <- dim(XX)[1]

	tmp0 <- XX%*%Theta
	tmp0 <- tmp0 -
		apply(tmp0,MARGIN=1,FUN=max)%*%t(rep(1,L))
	tmp1 <- exp(tmp0)
	tmp2 <- rowSums(tmp1)
	PM <- tmp1 / tcrossprod(tmp2, rep(1,L))
	
	# Computation of functional Lambda(Theta):
	Lambda <- sum(log(tmp2)) -
		sum(diag(1,L,L)[Y,]*tmp0)

	# Computation of gradient G(Theta):
	IL <- diag(rep(1,L))
	G <- as.vector(crossprod(XX, PM - IL[Y,]))
	
	# Computation of Hessian H(Theta):
	H <- matrix(0,nrow=L*dd,ncol=L*dd)

	for (y in 1:L)
	{
		tmpi <- ((y-1)*dd + 1):(y*dd)
		for (z in y:L)
		{
			tmpj <- ((z-1)*dd + 1):(z*dd)
			if (y == z)
			{
				H[tmpi,tmpj] <- crossprod(XX,
                                     (tcrossprod(PM[,y]-PM[,y]*PM[,z],rep(1,dd))*XX))
			}
			else
			{
				tmp.mat <- crossprod(XX,
                                            (tcrossprod(PM[,y]*PM[,z], rep(1,dd))*XX))
				H[tmpi,tmpj] <- - tmp.mat
				H[tmpj,tmpi] <- - tmp.mat
			}
		}
	}
	
	return(list(Lambda=Lambda,G=G,H=H,PM=PM))
}


#####################
# Regularization 0: #
#####################

R0.only <- function(Theta)
{
	R <- sum(rowSums(Theta)^2)/2
	return(R)
}

R0.full <- function(Theta)
{
	# Auxiliary objects:
	dd <- dim(Theta)[1]
	L <- dim(Theta)[2]
	
	# Computation of the functional R(Theta):
	R <- sum(rowSums(Theta)^2)/2
	
	# Computation of its Hessian and gradient:
	HR <- kronecker(matrix(1,L,L),diag(1,dd,dd))
	GR <- HR %*% as.vector(Theta)
	
	return(list(R=R,GR=GR,HR=HR))
}


###########################
# Regularization 0 and 1: #
###########################

R1.only <- function(Theta,tau,eps=2^(-10))
{
	R <- sum(rowSums(Theta)^2)/2 +
		sum(tau *
			sqrt(eps^2 + rowSums(Theta^2)))
	return(R)
}

R1.full <- function(Theta,tau,eps=2^(-10))
{
	# Auxiliary objects:
	dd <- dim(Theta)[1]
	L <- dim(Theta)[2]
	# Vector of smoothed norms of Theta's rows:
	rho.v <- sqrt(eps^2 + rowSums(Theta^2))
	# Matrix of smoothly normalized rows of Theta:
	g.m <- Theta / tcrossprod(rho.v, rep(1,L))
	
	# Computation of the functional R(Theta):
	R <- sum(rowSums(Theta)^2)/2 +
			sum(tau*rho.v)
	
	# Computation of its gradient:
	HR <- kronecker(matrix(1,L,L),diag(1,dd,dd))
	# == contribution of regularization 0 to Hessian.
	GR <- HR %*% as.vector(Theta) +
		as.vector(tcrossprod(tau, rep(1,L)) * g.m)
	
	# Computation of its Hessian:
	for (y in 1:L)
	{
		tmpi <- ((y-1)*dd+1):(y*dd)
		for (z in y:L)
		{
			tmpj <- ((z-1)*dd+1):(z*dd)
			if (y == z)
			{
				HR[tmpi,tmpj] <- HR[tmpi,tmpj] +
					diag(tau*(1 - g.m[,y]*g.m[,z])/rho.v)
			}
			else
			{
				HR[tmpi,tmpj] <- HR[tmpi,tmpj] +
					- diag(tau*(g.m[,y]*g.m[,z])/rho.v)
				HR[tmpj,tmpi] <- HR[tmpi,tmpj]
			}
		}
	}
	
	return(list(R=R,GR=GR,HR=HR))
}


###########################
# Regularization 0 and 2: #
###########################

R2.only <- function(Theta,tau,eps=2^(-10))
{
	sigma <- 1 - (tau > 0)
	R <- sum(sigma *
			rowSums(Theta)^2)/2 +
		sum(tau *
			rowSums(sqrt(eps^2+Theta^2)))
	return(R)
}

R2.full <- function(Theta,tau,eps=2^(-10))
{
	# Auxiliary objects:
	dd <- dim(Theta)[1]
	L <- dim(Theta)[2]
	sigma <- 1 - (tau > 0)
	rho.m <- sqrt(eps^2 + Theta^2)
	
	# Computation of the functional R(Theta):
	R <- sum(sigma *
			rowSums(Theta)^2)/2 +
		sum(diag(tau)%*%rho.m)
	
	# Computation of its gradient:
	HR <- kronecker(matrix(1,L,L),diag(sigma))
	# == contribution of regularization 0 to Hessian. 
	GR <- HR %*% as.vector(Theta) +
		as.vector(diag(tau) %*% (Theta/rho.m))
	
	# Computation of its Hessian:
	HR <- HR +
		diag(as.vector(diag(tau) %*% (eps^2/rho.m^3)))
	
	return(list(R=R,GR=GR,HR=HR))
}


#####################################
# Augmented negative log-likelihood #
# for usual logistic regression:    #
#####################################

LR0.only <- function(XX,Y,Theta)
{
	LR <- Lambda.only(XX,Y,Theta) +
		R0.only(Theta)
	return(LR)
}

LR0.full <- function(XX,Y,Theta)
{
	Lf <- Lambda.full(XX,Y,Theta)
	Rf <- R0.full(Theta)
	LR <- Lf$Lambda + Rf$R
	GLR <- Lf$G + Rf$GR
	HLR <- Lf$H + Rf$HR
	return(list(LR=LR,GLR=GLR,HLR=HLR,PM=Lf$PM))
}


#######################################
# Penalized negative log-likelihoods: #
#######################################

LR1.only <- function(XX,Y,Theta,tau)
{
	LR <- Lambda.only(XX,Y,Theta) +
		R1.only(Theta,tau)
	return(LR)
}

LR1.full <- function(XX,Y,Theta,tau)
{
	Lf <- Lambda.full(XX,Y,Theta)
	Rf <- R1.full(Theta,tau)
	LR <- Lf$Lambda + Rf$R
	GLR <- Lf$G + Rf$GR
	HLR <- Lf$H + Rf$HR
	return(list(LR=LR,GLR=GLR,HLR=HLR,PM=Lf$PM))
}

LR2.only <- function(XX,Y,Theta,tau)
{
	LR <- Lambda.only(XX,Y,Theta) +
		R2.only(Theta,tau)
	return(LR)
}

LR2.full <- function(XX,Y,Theta,tau)
{
	Lf <- Lambda.full(XX,Y,Theta)
	Rf <- R2.full(Theta,tau)
	LR <- Lf$Lambda + Rf$R
	GLR <- Lf$G + Rf$GR
	HLR <- Lf$H + Rf$HR
	return(list(LR=LR,GLR=GLR,HLR=HLR,PM=Lf$PM))
}



####################
# find optimal tau #
####################

find.tau.opt <- function(Xa,Ya,theta,tau.o=10, delta=2, tau.max=80, tau.min=1,                                             pen.method=c("vectors","simple","none"), a0=a0,b0=b0)
  {
    stopifnot(delta > 1)

    S <- stab(Xa,Ya,theta,tau.o=tau.o,pen.method, a0=a0,b0=b0)
    tau.tmp <- tau.o * delta
    if(tau.tmp <= tau.max)
      {
        S.new <- stab(Xa,Ya,theta,tau.o=tau.tmp,pen.method, a0=a0,b0=b0)
      } else
    {
      S.new <- S
    }
    if(S > S.new && tau.o * delta <= tau.max)
      {
        
        while(S > S.new && tau.tmp * delta <= tau.max)
          {
            S <- S.new
            tau.tmp <- tau.tmp * delta
            S.new <- stab(Xa,Ya,theta,tau.o=tau.tmp,pen.method, a0=a0,b0=b0)
          }
        if(S <= S.new)
          {
            tau.tmp <- tau.tmp / delta
          }
        return(tau.tmp)
      } else
    {
      if(tau.o / delta >= tau.min)
        {
          tau.tmp <- tau.o / delta
          S.new <- stab(Xa,Ya,theta,tau.o=tau.tmp,pen.method, a0=a0,b0=b0)
        } else {
          return(tau.o)
        }
      while(S > S.new && tau.tmp / delta >= tau.min)
        {      
          S <- S.new
          tau.tmp <- tau.tmp / delta
          S.new <- stab(Xa,Ya,theta,tau.o=tau.tmp,pen.method, a0=a0,b0=b0)
        }
      if(S <= S.new)
        {
          tau.tmp <- tau.tmp * delta
        }
      return(tau.tmp)
    }
  }


stab <- function(Xa,Ya,theta,tau.o,pen.method=c("vectors","simple","none"), a0=a0,b0=b0)
  {
    T <- rep(0, length(Ya))
    for(i in which(Ya != theta)) {
      Ytmp <- Ya
      Ytmp[i] <- theta
      tmp <- penlogreg(Xa, Ytmp, tau.o=tau.o, a0=a0, b0=b0)
      T[i] <- tmp$PM[i, theta]
    }
    S <- sum(T)
    return(S)
  }
