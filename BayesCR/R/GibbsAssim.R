Gibbs <- function(y,x,cc,distribution,cens,M,burnin,n.thin,n.chains){

  n.iter <- M
  cad <- 0
  set1<- set <- c()
	n <- length(y)
	p <- ncol(x)
  Cont <- ceiling(M/10)
  M <- Cont*10

	if(cens=="right")
	{
		lower <- y
		upper <- rep(Inf,n)
	}
	if(cens=="left")
	{
		lower <- rep(-Inf,n)
		upper <- y
	}
  y1 <- y
  
  beta <- matrix(0,ncol=p,nrow=M*n.chains)
  delta <- matrix(0,ncol=1,nrow=M*n.chains)
  tau <- matrix(0,ncol=1,nrow=M*n.chains)
  q <- 1
  nu <- matrix(0,ncol=1,nrow=M*n.chains)
  u <- matrix(0,ncol=1,nrow=n)
  t <- matrix(0,ncol=1,nrow=n)
  gama <- matrix(0,ncol=1,nrow=M*n.chains)
  media<- matrix(0,nrow=M*n.chains,ncol=n)

	mu.beta <- rep(0,p)
	sigma.beta <- 10*diag(p)
	mu.delta <- 0
	sigma.delta <- 10
	a.tau <- 2.1
	b.tau <- 3
	if(distribution=="ST")
	{
		a.gamma <- 0.02
		b.gamma <- 0.49
	}
	if(distribution=="SSL")
	{
		a.gamma <- 0.02
		b.gamma <- 0.9
	}
	
  cat('% of iterations \n')
  
  while(cad< n.chains)
  {
    
    cad <- cad+1
    Lin <- n.iter*(cad-1)
    W <- seq(Cont+Lin,n.iter*cad,Cont)
    z <- 1
        
    reg <- lm(y1[cc==0] ~ x[cc==0,1:p] - 1)
    beta[1+Lin,]<- as.vector(coefficients(reg),mode="numeric")
    m1 <- mean(y1[cc==0])
    m2 <- mean((y1[cc==0]-m1)^2)
    m3 <- mean((y1[cc==0]-m1)^3)
    t <- rnorm(n,mean=10,sd=10)
    
    if(distribution=="SN")
    {
      u <- rep(1,n)
      k1 <- k2 <- k3 <- 1
    }
    if(distribution=="ST")
    {
      nu[1+Lin] <- 4.7 
      u <- rgamma(n,shape=15,rate=0.5) #shape=.1 rate=.01
      k1 <- sqrt(nu[1+Lin]/2)*gamma((nu[1+Lin]-1)/2)/gamma(nu[1+Lin]/2)
      k2 <- (nu[1+Lin]/2)*gamma((nu[1+Lin]-2)/2)/gamma(nu[1+Lin]/2)
      k3 <- ((nu[1+Lin]/2)^(3/2))*gamma((nu[1+Lin]-3)/2)/gamma(nu[1+Lin]/2)
      gama[1+Lin] <- (a.gamma + b.gamma)/2 #0.25
    }
    if(distribution=="SSL")
    {
      nu[1+Lin] <- 4 
      u <- runif(n)
      k1 <- 2*nu[1+Lin]/(2*nu[1+Lin]-1)
      k2 <- 2*nu[1+Lin]/(2*nu[1+Lin]-2) 
      k3 <- 2*nu[1+Lin]/(2*nu[1+Lin]-3) 
      gama[1+Lin] <- (a.gamma + b.gamma)/2
    }
    f.ini <- function(w){
      s2.ini <- m2/(k2 - 2*k1^2*w^2/pi)
      m.ini <- m1 - k1*sqrt(2/pi)*s2.ini*w
      p1 <- sqrt(s2.ini)*w*sqrt(2/pi)
      resp <- m3 - (m.ini^3 + p1*(3*m.ini^2*k1 - k3 - 3*s2.ini) + 3*m.ini*k2*s2.ini + 2*p1^3 )
      return(resp)
    }
    d.ini <- uniroot.all(f=f.ini,interval=c(-1,1),tol=.Machine$double.eps,maxiter=10000000)
    f <- f.ini(d.ini)
    d.ini <- d.ini[abs(f.ini(d.ini))==min(abs(f))]
    sigma2.ini <- m2/(k2 - 2*k1^2*d.ini^2/pi)
    delta[1+Lin] <- d.ini*sqrt(sigma2.ini)
    tau[1+Lin] <- (1-d.ini^2)*sigma2.ini
    
    y <- y1
    
    cat("\n")
    Fil <- 2 + Lin
    Col <- n.iter*cad
    
	  for(k in Fil:Col)
	  {
		  if(distribution=="SN")
		  {
			  k1 <- 1
		  }
		  if(distribution=="ST")
		  {
			  k1 <- sqrt(nu[k-1]/2)*gamma((nu[k-1]-1)/2)/gamma(nu[k-1]/2)
		  }
		  if(distribution=="SSL")
		  {
			  k1 <- nu[k-1]/(nu[k-1]-0.5)
		  }
		  b1 <- -sqrt(2/pi)*k1

		  mu <- x%*%beta[k-1,]

		  mean.t <- (delta[k-1]/(delta[k-1]^2 + tau[k-1]))*(y-mu+b1*tau[k-1]/delta[k-1])
		  var.t <- tau[k-1]/(u*(delta[k-1]^2+tau[k-1]))
		  t <- rtrunc(n,"norm",a=b1,b=Inf,mean=mean.t,sd=sqrt(var.t))

		  mean.y <- mu + delta[k-1]*t
		  var.y <- tau[k-1]/u
		  for(i in 1:n)
		  {
			  if(cc[i]==1)
			  {
				  y[i] <- rtrunc(1,"norm",a=lower[i],b=upper[i],mean=mean.y[i],sd=sqrt(var.y[i]))
			  } 
		  }   

		  x.ast <- sqrt(diag(u))%*%x
		  y.ast <- y*sqrt(u)
		  t.ast <- t*sqrt(u)	

		  var.beta <- solve((t(x.ast)%*%x.ast/tau[k-1]) + solve(sigma.beta) )
		  ind <- lower.tri(var.beta)
      var.beta[ind] <- t(var.beta)[ind] 
		  mean.beta <- var.beta%*%(solve(sigma.beta)%*%mu.beta + (t(x.ast)%*%y.ast/tau[k-1]) - (delta[k-1]*t(x.ast)%*%t.ast/tau[k-1]))			
		  beta[k,] <- rmvnorm(1,mean=mean.beta,sigma=var.beta,method="chol")

		  mu <- x%*%beta[k,]

		  var.delta <- 1/(sum(u*t*t)/tau[k-1] + 1/sigma.delta)
		  mean.delta <- var.delta*(mu.delta/sigma.delta + (sum(u*t*(y-mu)))/tau[k-1])
		  delta[k] <- rnorm(1,mean=mean.delta,sd=sqrt(var.delta))

		  shape.tau <- a.tau+0.5*n
		  rate.tau <- b.tau+0.5*sum(u*(y-mu-delta[k]*t)^2)
		  tau[k] <- 1/rgamma(1,shape=shape.tau,rate=rate.tau)
      
		  media[k,] <- mu - delta[k]*sqrt(2/pi)*k1 
    
		  if(distribution=="ST")
		  {
			  A <- (y-mu-delta[k]*t)^2/tau[k] + (t-b1)^2
      	for(i in 1:n)
			  {
				  u[i] <- rgamma(1,shape=(0.5*nu[k-1]+1), rate=((nu[k-1]+A[i])/2) )			
			  }
        
			  gama[k] <- rtrunc(1,"gamma",a=a.gamma,b=b.gamma,shape=2,scale=(1/nu[k-1]))
       	mh <- MHnuST(nu[k-1],u,gama[k])
			  nu[k] <- mh$last 
          
		  }
		  if(distribution=="SSL")
		  {	
			  A <- (y-mu-delta[k]*t)^2/tau[k] + (t-b1)^2
			  u <- rtrunc(n,"gamma",a=0.001,b=1,shape=(nu[k-1]+1),scale=1/(0.5*A))
			  gama[k] <- rtrunc(1,"gamma",a=a.gamma,b=b.gamma,shape=2,scale=(1/nu[k-1])) 
			  nu[k] <- rtrunc(1,"gamma",a=1.001,b=Inf,shape=(n+1),scale=1/(gama[k]-sum(log(u))))		
		  }
		  
		  if(k==W[z])
		  {
		    z <- z +1
		    BayCens(Fil,Col,k,cad)
		  }
	  }
	  
    cat('\r')
    initial <- burnin 
    set1 <- seq(initial+n.thin +n.iter*(cad-1) ,n.iter*cad,n.thin)
    set <- c(set,set1)  
    
  }
  
  lambda <- delta/sqrt(tau)
  sigma2 <- tau*(1+lambda^2)

	if(distribution=="SN")
	{
	  resposta <- list(beta=beta[set,],sigma2=sigma2[set],lambda=lambda[set],mu=media[set,])
	}
	else
	{
		resposta <- list(beta=beta[set,],sigma2=sigma2[set],lambda=lambda[set],nu=nu[set],mu=media[set,])
	}
	return(resposta)
}