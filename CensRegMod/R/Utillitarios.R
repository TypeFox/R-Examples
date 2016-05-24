
cdfNI<-function(x,mu,sigma2,nu,type)
{
  
  resp<-matrix(0,length(x),1)

	if(type=="Normal")
  {
		resp<-pnorm(x,mu,sqrt(sigma2))
  }

	if(type=="T")
  {
		z=(x-mu)/sqrt(sigma2)
		resp=pt(z,df=nu)
	}
  return(resp)
}

MomNIT<-function(mu,sigma2,nu,a,type){

	if(type=="Normal")
  {
		z<-(a-mu)/sqrt(sigma2)
		Ey<-mu+sqrt(sigma2)*dnorm(z)/(1-pnorm(z))
		Ey2<-mu^2+sigma2+(dnorm(z)/(1-pnorm(z)))*(mu+a)*sqrt(sigma2)
	}

	if(type=="T")
  {
		a1<- (a-mu)/sqrt(sigma2)
		Aux <- (1-cdfNI(a1,0,nu/(nu-2),nu-2,type))/(1-cdfNI(a1,0,1,nu,type))
		G1<- 0.5*(gamma((nu-1)/2)*nu^(nu/2))/((1-cdfNI(a1,0,1,nu,type))*gamma(nu/2)*gamma(1/2))
		Ex<- G1*(nu+a1^2)^(-(nu-1)/2)
		Ey<- mu+sqrt(sigma2)*Ex
		Ex2<- Aux*(nu/(nu-2))+G1*(a1*(nu+a1^2)^(-(nu-1)/2))
		Ey2<- mu^2+sigma2*Ex2+2*mu*sqrt(sigma2)*Ex
	}

	return(list(Ey=Ey,Ey2=Ey2))
}

CalMom<-function(mu,sigma2,nu,a,type)
{
	n<-length(mu)

	if(type=="Normal")
  {
		Mom<-MomNIT(mu=mu,sigma2=sigma2,a=a,type=type)
		Ey0<-rep(1,n)
		Ey1<-Mom$Ey
		Ey2<-Mom$Ey2
	}

	if(type=="T")
  {
		sigma2a<- sigma2*nu/(nu+2)
		aux1<-(1-cdfNI(a,mu,sigma2a,nu+2,type))/(1-cdfNI(a,mu,sigma2,nu,type))
		aux2 <- gamma((nu+1)/2)*gamma((nu+2)/2)*(nu+1)/(nu*gamma(nu/2)*gamma((nu+3)/2)  )
		Ey0<- aux1*aux2*rep(1,n)
		Mom<- MomNIT(mu,sigma2a,nu+2,a,type)
		Ey1<- aux1*aux2*Mom$Ey
		Ey2<- aux1*aux2*Mom$Ey2
	}
  
	return(list(Ey0=Ey0,Ey1=Ey1,Ey2=Ey2))
}

CalMomMIt<-function(mu,sigma2,nu,a)
{
	n<-length(mu)
	sigma2a<- sigma2*nu/(nu+4)
	kv <- gamma((nu+1)/2)*gamma((nu+4)/2)*(nu+1)*(nu+1)/(nu*nu*gamma(nu/2)*gamma((nu+5)/2)  )
	aux1<- (1-cdfNI(a,mu,sigma2a,nu+4,type="T"))/(1-cdfNI(a,mu,sigma2,nu,type="T"))
	Ey0<- kv*aux1*rep(1,n)
	Mom<- MomNIT(mu,sigma2a,nu+4,a,type="T")
	Ey1<- kv*aux1*Mom$Ey
	Ey2<- kv*aux1*Mom$Ey2
	return(list(Ey0=Ey0,Ey1=Ey1,Ey2=Ey2))
}

HessianaQ<-function(cc,x,y,beta,sigma2,u0,u1,u2)
{
 	p <- ncol(x)
	n <- nrow(x)    
  mu<-x%*%beta

  suma1<-matrix(0,p,p)
  suma2<-0
  suma3<-matrix(0,p,1)

  for(i in 1:n)
  {
    suma1<-suma1+(-1/sigma2)*(u0[i]*((x[i,])%*%t(x[i,])))
    suma2<-suma2+((0.5/(sigma2^2))-((1/(sigma2^3))*(u2[i]-2*u1[i]*x[i,]%*%beta+((x[i,]%*%beta)^2)*u0[i]))    )
    suma3<-suma3+((1/(sigma2^2))*((u0[i]*((x[i,])%*%t(x[i,])%*%beta))-u1[i]*(x[i,])))
  }
  
  derbeta<-suma1
  dersigma2<-suma2
	derbetasigma2<-suma3

	MatrizQ<-matrix(0,nrow=(p+1),ncol=(p+1))
  MatrizQ[1:p,1:p]<-derbeta
  MatrizQ[p+1,1:p]<-t(derbetasigma2)
  MatrizQ[1:p,p+1]<-derbetasigma2
  MatrizQ[p+1,p+1]<-dersigma2
    
  MatrizQbeta<-derbeta
  MatrizQsigma2<-dersigma2

  obj.out <- list(MatrizQ = MatrizQ, MatrizQbeta = MatrizQbeta, MatrizQsigma2 = MatrizQsigma2)

  return(obj.out)
}

Q.function<-function(cc,x,y,beta,sigma2,u0,u1,u2)
{
	p <- ncol(x)
	n <- nrow(x)
  mu<-x%*%beta
  suma<-0

  for(i in 1:n)
  {
    suma <- suma + ((-0.5*log(2*pi)) - (0.5*log(sigma2)) - ((0.5/sigma2)*(u2[i]-2*u1[i]*x[i,]%*%beta+((x[i,]%*%beta)^2)*u0[i])) )
  }
	Q.theta<-suma
	
	return(Q.theta)
}


QD<-function(cc,x,y,beta,sigma2,u0,u1,u2,MatrizQ)
{

 	p <- ncol(x)
	n <- nrow(x)
  mu<-x%*%beta
            
  theta<-c(beta,sigma2)

  Q.theta<-Q.function(cc,x,y,beta,sigma2,u0,u1,u2)
  Q.theta<-rep(Q.theta,n)
  Q.thetai<-rep(0,n)

  for (i in 1:n)
  {
    derthetai<-rep(0,p+1)
    thetai<-rep(0,p+1+1)

    derbeta<-(1/sigma2)*(x[i,]*u1[i]-u0[i]*((x[i,])%*%t(x[i,])%*%beta))
    dersigma2<-  -(0.5/sigma2^2) - ((0.5/sigma2)*(u2[i]-2*u1[i]*x[i,]%*%beta+((x[i,]%*%beta)^2)*u0[i]))

    derthetai<-matrix(c(t(derbeta),dersigma2),nrow=(p+1),ncol=1)
    thetai<-theta-(solve((-1*MatrizQ))%*%(derthetai))

    betai<-matrix(c(thetai[1:p]),p,1)
    sigma2i<-as.numeric(thetai[p+1])

    Q.thetai[i]<- Q.function(cc,x,y,betai,sigma2i,u0,u1,u2)
  }
  
  QD<-abs(2*(Q.theta-Q.thetai))

  return(QD)
}


CD<-function(cc,x,y,beta,sigma2,u0,u1,u2,MatrizQ,MatrizQbeta,MatrizQsigma2)
{

	p <- ncol(x)
	n <- nrow(x)
      
  mu<-x%*%beta
  theta<-c(t(beta),sigma2)

  GD<-rep(0,n)
	GDbeta<-rep(0,n);GDsigma2<-rep(0,n)

  for (i in 1:n)
  {
	  derthetai<-rep(0,p+1)
    thetai<-rep(0,p+1)

    derbeta<-(1/sigma2)*(x[i,]*u1[i]-u0[i]*((x[i,])%*%t(x[i,])%*%beta))
    dersigma2<-( -(0.5/sigma2) + ((0.5/sigma2^2)*(u2[i]-2*u1[i]*x[i,]%*%beta+((x[i,]%*%beta)^2)*u0[i])) )

    derthetai<-matrix(c(t(derbeta),dersigma2),nrow=(p+1),ncol=1)
    thetai<-theta-(solve((-1*MatrizQ))%*%(derthetai))

    thetabeta<- theta[1:p]; thetabetai<- thetai[1:p]
    thetasigma2<-theta[p+1]; thetasigma2i<-thetai[p+1]

   	GD[i]<-t(theta-thetai)%*%(-1*MatrizQ)%*%(theta-thetai)
    GDbeta[i]<-t(thetabeta-thetabetai)%*%(-1*MatrizQbeta)%*%(thetabeta-thetabetai)
    GDsigma2[i]<-t(thetasigma2-thetasigma2i)%*%(-1*MatrizQsigma2)%*%(thetasigma2-thetasigma2i)
  }

  obj.out<-list(GD = GD, GDbeta = GDbeta ,GDsigma2 = GDsigma2)

  return(obj.out)
}

If<-function(cc,x,y,beta,sigma2,u0,u1,u2,MatrizQ,Perturb)
{
  p <- ncol(x)
	n <- nrow(x)
  mu<-x%*%beta          
  theta<-c(t(beta),sigma2)
 	mdelta<-matrix(0,(p+1),n)

  for (i in 1:n )
  {
	  if(Perturb==1)
    {
  	  delta1<-(1/sigma2)*(x[i,]*u1[i]-u0[i]*((x[i,])%*%t(x[i,])%*%beta)) #derbetaw
    	delta2<-( -(0.5/sigma2) + ((0.5/sigma2^2)*(u2[i]-2*u1[i]*x[i,]%*%beta+((x[i,]%*%beta)^2)*u0[i])) ) #dersigma2w
    	mdelta[,i]<-c(delta1,delta2)
    }

    if(Perturb==2)
    {
		  delta1<-(1/sigma2)*(x[i,]*u1[i]-u0[i]*((x[i,])%*%t(x[i,])%*%beta)) #derbetaw
    	delta2<-((0.5/sigma2^2)*(u2[i]-2*u1[i]*x[i,]%*%beta+((x[i,]%*%beta)^2)*u0[i])) #dersigma2w
    	mdelta[,i]<-c(delta1,delta2)
   	}
  }

  If<-t(mdelta)%*%solve(MatrizQ)%*%mdelta
  return(If)
}

IM.CensGrad <- function(cc,x,y,BETA,SIGMA2,nu,type)
{      
  x <- as.matrix(x)
  betas <- BETA
  sigma2 <- SIGMA2
  mu <- x%*%(betas)
  n <- length(y)
  
  if(type=="Normal")
  {
    u0 <- rep(1,n)
    u1 <- y
    u2 <- y^2
  }
  if(type=="T")
  {
    d <- ((y-mu)^2)/sigma2
    u0 <- (nu+1)/(nu+d)
    u1 <- y*(nu+1)/(nu+d) 
    u2 <- y^2*(nu+1)/(nu+d)
  }
  
  if(sum(cc)>0)
  {
    aux <- CalMom(mu,sigma2,nu,y,type)
    u1[cc==1]<- aux$Ey1[cc==1]
    u2[cc==1]<- aux$Ey2[cc==1]
    u0[cc==1]<- aux$Ey0[cc==1]  
  }
  
  p <- length(betas)
  se <- matrix(0,1,p)   
  AnsBet <- c()
  AnsSig <- c()
  MIE <- matrix(0,p+1,p+1)   
  A <- matrix(0,1,p+1)   
  for(i in 1:n)
  {
    AnsBet <- (x[i,]*u1[i]-u0[i]*x[i,]*mu[i])/sigma2
    AnsSig <- 0.5*((-1/sigma2) + (1/sigma2^2)*(u2[i]-2*u1[i]*mu[i]+ u0[i]*mu[i]^2))
    A <- c(AnsBet,AnsSig)
    MIE1 <- A%*%t(A)
    ind <- lower.tri(MIE1)
    MIE1[ind] <- t(MIE1)[ind]
    MIE <- MIE1 + MIE
  }
  se <- sqrt(diag(solve(MIE)))
  return(se)
}

