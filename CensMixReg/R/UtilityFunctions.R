###########  Density of Normal and t with location - scale ###########

dNormal <- function(cc, y, mu, sigma2 = 1)
{
 densN        <- vector(mode = "numeric", length = length(y))
 densN[cc==0] <- dnorm(y[cc==0], mu[cc==0], sqrt(sigma2))
 densN[cc==1] <- 1-pnorm((y[cc==1] - mu[cc==1])/sqrt(sigma2))
 if(length(which(densN == 0)) > 0) densN[which(densN == 0)] <- .Machine$double.xmin
 return(densN)
}

dT <- function(cc, y, mu, sigma2 = 1,nu=4){
  densN        <- vector(mode = "numeric", length = length(y))
  aux          <- (y-mu)/sqrt(sigma2)
  densN[cc==0] <- dt(aux[cc==0],nu)/sqrt(sigma2)
  densN[cc==1] <- 1-pt(aux[cc==1],nu)
  if(length(which(densN == 0)) > 0) densN[which(densN == 0)] <- .Machine$double.xmin
  return(densN)
}


###########    Densities of Normal Mixtures   ##################

d.mixedN <- function(cc, x, pi1, mu, sigma2){
    # x: data vector
    # mu[,] a matrix
    # other parameters must be of the vector type c() of g dimensions (number of mixtures)
    g    <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dNormal(cc, x, mu[, j], sigma2[j])
    return(dens)
  }


d.mixedT <- function(cc, x, pi1, mu, sigma2, nu)
{
 # x: data vector
 ## mu[,] a matrix
 # other parameters must be of the vector type c() of g dimensions (numbe of mixtures)
 g    <- length(pi1)
 dens <- 0
 for(j in 1:g) dens <- dens + pi1[j]*dT(cc, x, mu[, j], sigma2[j],nu)
 return(dens)
}
##########     End  the densities Normal and t   ############
################################################################


MomN <- function(mu,sigma2,a)
{
 z   <- (a-mu)/sqrt(sigma2)
 aux <-(1-pnorm(z))
 if(length(which(aux == 0)) > 0) aux[which(aux == 0)] <- .Machine$double.xmin
 Ey  <- mu+sqrt(sigma2)*dnorm(z)/aux
 Ey2 <- mu^2+sigma2+(dnorm(z)/aux)*(mu+a)*sqrt(sigma2)
 return(list(Ey=Ey,Ey2=Ey2))
}

MomT <- function(mu,sigma2,nu,a)
{
 a1  <- (a-mu)/sqrt(sigma2)
 aux <- 1-pt(a1,nu)
 if(length(which(aux == 0)) > 0){ aux[which(aux == 0)] <- .Machine$double.xmin }
 G1  <- 0.5*(gamma((nu-1)/2)*nu^(nu/2))/(aux*gamma(nu/2)*gamma(1/2))
 Ex  <- G1*(nu+a1^2)^(-(nu-1)/2)
 Ey  <- mu+sqrt(sigma2)*Ex
 Ex2 <- nu/(nu-2)+G1*(a1*(nu+a1^2)^(-(nu-1)/2))
 Ey2 <- mu^2+sigma2*Ex2+2*mu*sqrt(sigma2)*Ex
 return(list(Ey=Ey,Ey2=Ey2))
}

cdfNI  <- function(x,mu,sigma2,nu,type)
{
  resp <- matrix(0,length(x),1)

	if(type=="Normal"){ resp<-pnorm(x,mu,sqrt(sigma2)) }

	if(type=="T")
  {
		z    <- (x-mu)/sqrt(sigma2)
		resp <- pt(z,df=nu)
	}
  return(resp)
}

MomNIT <- function(mu,sigma2,nu,a,type)
{
	if(type=="Normal")
  {
		z   <- (a-mu)/sqrt(sigma2)
    aux <- (1-pnorm(z))
    if(length(which(aux == 0)) > 0) aux[which(aux == 0)] <- .Machine$double.xmin
		Ey  <- mu+sqrt(sigma2)*dnorm(z)/aux
 		Ey2 <- mu^2+sigma2+(dnorm(z)/aux)*(mu+a)*sqrt(sigma2)
	}

	if(type=="T")
  {
	 a1   <- (a-mu)/sqrt(sigma2)
   aux1 <- (1-cdfNI(a1,0,1,nu,type))
   if(length(which(aux1 == 0)) > 0) {  aux1[which(aux1 == 0)] <- .Machine$double.xmin }
	 Aux  <- (1-cdfNI(a1,0,nu/(nu-2),nu-2,type))/aux1
	 G1   <- 0.5*(gamma((nu-1)/2)*nu^(nu/2))/(aux1*gamma(nu/2)*gamma(1/2))
	 Ex   <- G1*(nu+a1^2)^(-(nu-1)/2)
	 Ey   <- mu+sqrt(sigma2)*Ex
	 Ex2  <- Aux*(nu/(nu-2))+G1*(a1*(nu+a1^2)^(-(nu-1)/2))
	 Ey2  <- mu^2+sigma2*Ex2+2*mu*sqrt(sigma2)*Ex
	}

	return(list(Ey=Ey,Ey2=Ey2))
}


CalMom <- function(mu,sigma2,nu,a,type)
{
	n<-length(mu)

	if(type=="Normal")
  {
		Mom <- MomNIT(mu=mu,sigma2=sigma2,a=a,type=type)
		Ey0 <- rep(1,n)
		Ey1 <- Mom$Ey
		Ey2 <- Mom$Ey2
	}

	if(type=="T")
  {
	 sigma2a <- sigma2*nu/(nu+2)
   aux0    <- (1-cdfNI(a,mu,sigma2,nu,type))
   if(length(which(aux0 == 0)) > 0) aux0[which(aux0 == 0)] <- .Machine$double.xmin
	 aux1    <- (1-cdfNI(a,mu,sigma2a,nu+2,type))/aux0
	 aux2    <- gamma((nu+1)/2)*gamma((nu+2)/2)*(nu+1)/(nu*gamma(nu/2)*gamma((nu+3)/2))
	 Ey0     <- aux1*aux2*rep(1,n)
	 Mom     <- MomNIT(mu,sigma2a,nu+2,a,type)
	 Ey1     <- aux1*aux2*Mom$Ey
	 Ey2     <- aux1*aux2*Mom$Ey2
	}

	return(list(Ey0=Ey0,Ey1=Ey1,Ey2=Ey2))
}


#To obtain the observed information matrix
imm.fmcr = function(y,x1,cc,model)
{
 #if((class(model) != "T") && (class(model) != "Normal")) stop(paste("Family",class(model),"not recognized.",sep=" "))
 n <- length(y)
 p <- ncol(x1)-1
 g <- length(model$pii)

 AbetasIMM <- model$Abetas
 medjIMM   <- model$medj
 sigma2IMM <- model$sigma2
 nuIMM     <- model$nu
 piiIMM    <- model$pii

 beta0IMM  <- AbetasIMM[1]
 betasIMM  <- as.matrix(AbetasIMM[2:(p+1)])   # parameters of regression dimension "p"
 varphi    <- rep(0,g) # beta0 + mu_j
 muIMM = mu1IMM <- matrix(0,n,g)
 x         <-  as.matrix(x1[,2:(p+1)])

 for (k in 1:g)
 {
  varphi[k] <- beta0IMM+medjIMM[k]
  mu1IMM[,k]<- x%*%betasIMM
  muIMM[,k] <- mu1IMM[,k]+varphi[k]
 }

 tal        <- matrix(0,n,g)
 Sibeta     <- matrix(0,p,n)
 Sibeta0    <- matrix(0,1,n)
 Simu       <- matrix(0,2,n)
 Sisigma    <- matrix(0,2,n)

 for(j in 1:g)
 {
  d        <- ((y-muIMM[,j])^2)/sigma2IMM[j]
  u0       <- (nuIMM+1)/(nuIMM+d)
  u1       <- y*(nuIMM+1)/(nuIMM+d)
  u2       <- y^2*(nuIMM+1)/(nuIMM+d)

  d1       <- dT(cc, y, muIMM[,j], sigma2IMM[j],nuIMM)
  d2       <- d.mixedT(cc, y, piiIMM, muIMM, sigma2IMM,nuIMM)
  tal[,j]  <- d1*piiIMM[j]/d2

  Sibeta      <- Sibeta + t(x)%*%diag(tal[,j]*c(u1)/sigma2IMM[j]) - t(x)%*%diag(tal[,j]*c(u0)*(mu1IMM[,j] + varphi[j])/sigma2IMM[j])
  Sibeta0     <- Sibeta0 + (1/sigma2IMM[j])*tal[,j]*c(u1) - tal[,j]*c(u0)*(mu1IMM[,j] + varphi[j])
  Simu[j,]    <- (1/sigma2IMM[j])*(tal[,j]*c(u1) - tal[,j]*c(u0)*(mu1IMM[,j] + varphi[j]))
  Sisigma[j,] <- -0.5*(1/sigma2IMM[j]^2)*tal[,j]*(sigma2IMM[j] - c(u2) + 2*c(u1)*(mu1IMM[,j] + varphi[j]) - c(u0)*(mu1IMM[,j] + varphi[j])^2)
 }

 Sip <- matrix(0,1,n)
 for(j in 1:(g-1))
   Sip[j,] <- (1/piiIMM[j])*tal[,j]
 Sip <- Sip - 1


 Isum <- matrix(0,7 + length(AbetasIMM)-2,7 + length(AbetasIMM)-2)
 si   <- matrix(0,n,7 + length(AbetasIMM)-2)
 for(i in 1:n)
 {
  si[i,] <- c(Sibeta0[,i],Sibeta[,i],Simu[,i],Sisigma[,i],Sip[,i])
  Isum   <- Isum + cbind(si[i,])%*%si[i,]
 }

 return(list(IM=Isum,class=class(model)))
}






