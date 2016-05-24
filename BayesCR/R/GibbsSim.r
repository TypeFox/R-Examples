GibbsTruncSMN<-function(cc,y,x,n.iter,n.thin,burnin,type="normal", cens="1",prior="Jeffreys",hyper=hyper, n.chains=n.chains)
{
  
  n <- length(y)
  
  if(cens=="left"){cens <- "1"}
  if(cens=="right"){cens <- "2"}
  
  
  if(cens=="2")
  {
    lower <- y
    upper <- rep(Inf,n)
  }
  if(cens=="1")
  {
    lower <- rep(-Inf,n)
    upper <- y
  }
    
  cad <- 0
  x <- as.matrix(x)
  p <- ncol(x)                                            
  n <- nrow(x)                                                
  sigma2<- lambda<- Aux <- V <- rho <- matrix(0,n.iter*n.chains,1)
  media<- matrix(0,n.iter*n.chains,n)
  beta<- matrix(0,n.iter*n.chains,p)
  y1<-y
   
  mu0<-matrix(0,p,1)
  Sigma0<-1000*diag(p)
  alfa1<-0.1
  alfa2<-0.01
  a <-  r1 <- 0.02 
  b <-  s1 <- 0.5     
  cat('% of iterations \n')
  
  if(type=="T")
  {
    U  <-matrix(0,n.iter*n.chains,n)
    nu <-matrix(0,n.iter*n.chains,1)
    set1<- set <- c()
    while(cad< n.chains)
    {
      cad <- cad +1
      Cont <- ceiling(n.iter/10)
      n.iter <- Cont*10
      Lin <- n.iter*(cad-1)
      W <- seq(Cont+Lin,n.iter*cad,Cont)
      reg <- lm(y ~ x[,2:p])
      z <- 1 
      lambda[1+ Lin] <- 1
      beta[1+ Lin,]<- as.vector(coefficients(reg),mode="numeric")
      sigma2[1+ Lin] <- sum((y-x%*%(beta[1+Lin,]))^2)/(n-p)
      media[1+ Lin,]<-x%*%(beta[1+Lin,])
      U[1+Lin,]<-1/rgamma(1,1)       
      nu[1+Lin]<- 6
      cat("\n")
      Fil <- 2 + Lin
      Col <- n.iter*cad
      for (j in Fil:Col)
      {
        for(i in 1:n)
        {
          if(cc[i]==1)
          {
            y[i] <- rtrunc(1,"norm",a=lower[i],b=upper[i],mean=media[j-1,i],sd=sqrt((U[j-1,i])^(-1)*sigma2[j-1]))
          } 
        }
        
        sigma2[j] <- 1/rgamma(1, shape=alfa1+n/2, rate = (sum(U[j-1,]*(y-media[j-1,])^2)/2+alfa2))                       ### Densid de sigma2/Resto   
        xast <-sqrt(U[j-1,])*x
        yast <-sqrt(U[j-1,])*y
        SigmaA <-solve(solve(Sigma0)+t(xast)%*%(xast)/sigma2[j])
        muA <- SigmaA%*%(solve(Sigma0)%*%mu0+t(xast)%*%yast/sigma2[j])
        beta[j,]<-rmvnorm(1,muA,SigmaA)                                                                           ### Densid de beta/Resto   
        media[j,]<-x%*%(beta[j,])
        meany2<-(y-media[j,])^2/sigma2[j-1]
        par2gam<- (meany2+nu[j-1])/2
        U[j,]<-rgamma(n,nu[j-1]/2+0.5,par2gam)                                                                  ### Densid de U/Resto    
        Aux[j] <- runif(1,0,1)
        V[j] <- pgamma(a,nu[j-1]) + (pgamma(b,2,nu[j-1])- pgamma(a,2,nu[j-1]))*Aux[j]
        lambda[j] <- qgamma(V[j],2,nu[j-1])                                                                                                                
        nu[j] <- MHnu(nu[j-1],U[j,],lambda[j],prior,hyper)                                                     ### Densid de Nu/Resto    

        if(j==W[z])
        {
          z <- z +1
          BayCens(Fil,Col,j,cad)
        }
      }
      cat('\r')
      initial <- burnin 
      set1 <- seq(initial+n.thin +n.iter*(cad-1) ,n.iter*cad,n.thin)
      set <- c(set,set1)
    }
    return(list(beta=beta[set,],sigma2=sigma2[set],U=U[set,],nu=nu[set],mu=media[set,]))
  }
  
  if(type=="Normal")
  {
    set1<- set <- c()
    while(cad< n.chains)
    {
      cad <- cad +1
      Cont <- ceiling(n.iter/10)
      n.iter <- Cont*10
      Lin <- n.iter*(cad-1)
      W <- seq(Cont+Lin,n.iter*cad,Cont)
      reg <- lm(y ~ x[,2:p])
      z <- 1 
      beta[1+ Lin,]<- as.vector(coefficients(reg),mode="numeric")
      sigma2[1+ Lin] <- sum((y-x%*%(beta[1+Lin,]))^2)/(n-p)
      media[1+ Lin,]<-x%*%(beta[1+Lin,])
      cat("\n")
      Fil <- 2 + Lin
      Col <- n.iter*cad
      for (j in Fil:Col)
      {
        for(i in 1:n)
        {
          if(cc[i]==1)
          {
            y[i] <- rtrunc(1,"norm",a=lower[i],b=upper[i],mean=media[j-1,i],sd = sqrt(sigma2[j-1]))
          } 
        }
        sigma2[j]<- 1/rgamma(1, alfa1+n/2, rate = (sum((y-media[j-1,])^2)/2+alfa2))                                        ### Densid de sigma2/Resto    
        xast <- x
   	    yast <- y
        SigmaA <- solve(solve(Sigma0)+t(xast)%*%(xast)/sigma2[j])
        muA <- SigmaA%*%(solve(Sigma0)%*%mu0+t(xast)%*%yast/sigma2[j])
        beta[j,] <- rmvnorm(1,muA,SigmaA)                                                                                    ### Densidade de beta/Resto
        media[j,] <- x%*%(beta[j,])

        if(j==W[z])
        {
          z <- z +1
          BayCens(Fil,Col,j,cad)
        }
      }
      cat('\r')
      initial <- burnin 
      set1 <- seq(initial+n.thin +n.iter*(cad-1) ,n.iter*cad,n.thin)
      set <- c(set,set1)    
    }
    return(list(beta=beta[set,],sigma2=sigma2[set],mu=media[set,]))
  }
  
  if(type=="Slash")                     
  {
  U <-matrix(0,n.iter*n.chains,n)
  nu <-matrix(0,n.iter*n.chains,1)
  set1<- set <- c()
  while(cad< n.chains)
  {
    r1 <- 0.01                             
    s1 <- 1
    cad <- cad +1
    Cont <- ceiling(n.iter/10)
    n.iter <- Cont*10
    Lin <- n.iter*(cad-1)
    W <- seq(Cont+Lin,n.iter*cad,Cont)
    W 
    reg <- lm(y ~ x[,2:p])
    z <- 1 
    beta[1+ Lin,]<- as.vector(coefficients(reg),mode="numeric")
    sigma2[1+ Lin] <- sum((y-x%*%(beta[1+Lin,]))^2)/(n-p)
    media[1+ Lin,]<- x%*%(beta[1+Lin,])
    U[1+Lin,] <- rbeta(n,1,1)       
    nu[1+Lin] <- 1.5
    cat("\n")
    Fil <- 2 + Lin
    Col <- n.iter*cad
    for (j in Fil:Col)
    {
      for(i in 1:n)
      {
        if(cc[i]==1)
        {
          y[i] <- rtrunc(1,"norm",a=lower[i],b=upper[i],mean=media[j-1,i],sd = sqrt((U[j-1,i])^(-1)*sigma2[j-1]))
        } 
      }
      
      sigma2[j] <- 1/rgamma(1, alfa1+n/2, rate = (sum(U[j-1,]*(y-media[j-1,])^2)/2+alfa2))                               ### Densid de sigma2/Resto   
      xast <- sqrt(U[j-1,])*x                     
      yast <- sqrt(U[j-1,])*y                                                                    
 	    SigmaA <- solve(solve(Sigma0)+t(xast)%*%(xast)/sigma2[j])
      muA <- SigmaA%*%(solve(Sigma0)%*%mu0+t(xast)%*%yast/sigma2[j])
      beta[j,] <- rmvnorm(1,muA,SigmaA)                                                                                   ### Densid de beta/Resto   
      media[j,] <- x%*%(beta[j,])      
      meany2 <- (y-media[j,])^2/sigma2[j-1]
      par2gam <- (meany2)/2                                                                                            
  	  U[j,] <- rtrunc(n, "gamma", a =0.001 , b =1 , shape=nu[j-1]+0.5, scale=1/par2gam)                                    ### Densid de U/Resto             
      rho[j] <- rtrunc(1, "gamma", a = r1, b =s1 , shape=2, scale=1/nu[j-1])
      nu[j] <- rgamma(1,shape= n+1,rate = rho[j]-sum(log(U[j,])))                                                     ### Densidade de Nu/Resto    

      if(j==W[z])
      {
        z <- z +1
        BayCens(Fil,Col,j,cad)
      }
    }
    cat('\r')
    initial <- burnin 
    set1 <- seq(initial+n.thin +n.iter*(cad-1) ,n.iter*cad,n.thin)
    set <- c(set,set1)   
    }
    return(list(beta=beta[set,],sigma2=sigma2[set],U=U[set,],nu=nu[set],mu=media[set,]))
  }

  if(type=="NormalC")
  {
    U <-matrix(0,n.iter*n.chains,n)
    nu <-rho <- rho1 <- matrix(0,n.iter*n.chains,1)
    p1 <- p2 <- ptot <- A <- matrix(0,n.iter*n.chains,n)
    set1<- set <- c()
    
    while(cad< n.chains)
    {
      AuxCN <- matrix(0,n,1)
      v0 <- s0 <- 2 
      v1 <- s1 <- 2
      cad <- cad +1 
      Cont <- ceiling(n.iter/10)
      n.iter <- Cont*10
      Lin <- n.iter*(cad-1)
      W <- seq(Cont+Lin,n.iter*cad,Cont)
      W 
      reg <- lm(y ~ x[,2:p])
      z <- 1 
      beta[1+ Lin,]<- as.vector(coefficients(reg),mode="numeric")
      sigma2[1+ Lin] <- sum((y-x%*%(beta[1+Lin,]))^2)/(n-p)
      media[1+ Lin,]<- x%*%(beta[1+Lin,])
      U[1+Lin,] <- 1
      rho[1+Lin] <- nu[1+Lin] <- 0.4
      cat("\n")
      Fil <- 2 + Lin
      Col <- n.iter*cad
      cont  <- c()       
      for (j in Fil:Col)
      {
        for(i in 1:n)
        {
          if(cc[i]==1)
          {
            y[i] <- rtrunc(1,"norm",a=lower[i],b=upper[i],mean=media[j-1,i],sd = sqrt((U[j-1,i])^(-1)*sigma2[j-1]))
          } 
        }
        sigma2[j] <- 1/rgamma(1, alfa1+n/2, rate = (sum(U[j-1,]*(y-media[j-1,])^2)/2+alfa2))                               ### Densid de sigma2/Resto   
        xast <- sqrt(U[j-1,])*x
        yast <- sqrt(U[j-1,])*y                                                                    
 	      SigmaA <- solve(solve(Sigma0)+t(xast)%*%(xast)/sigma2[j])
        muA <- SigmaA%*%(solve(Sigma0)%*%mu0+t(xast)%*%yast/sigma2[j])
        beta[j,] <- rmvnorm(1,muA,SigmaA)                                                                                   ### Densid de beta/Resto   
        media[j,] <- x%*%(beta[j,])
        meany2 <- (y-media[j,])^2/sigma2[j-1]
        par2gam <- (meany2)/2

        p1[j,] <- nu[j-1]*sqrt(rho[j-1])*exp(-0.5*rho[j-1]*meany2)
        p2[j,] <- (1-nu[j-1])*exp(-0.5*meany2)
        ptot[j,] <- p1[j,]+p2[j,]
        A[j,] <- p1[j,]/ptot[j,]
        AuxCN <-  rbinom(n=n, size=1, prob=A[j,])
        U[j,] <- rho[j-1]*AuxCN + (1-AuxCN)                                                                                ### Densid de U/Resto          
        cont <- sum(U[j,]==rho[j-1])
        nu[j] <- rbeta(1,v0+cont,v1+n-cont)                                                                                 ### Densid de Nu/Resto

        rho1[j] <- MHrhoCN(rho[j-1],meany2,cc,sigma2[j],c(nu[j],rho[j-1]),s0,s1,cens)
        rho[j] <- rho1[j]/(1 + rho1[j])  
 
        if(j==W[z])
        {
          z <- z +1
          BayCens(Fil,Col,j,cad)
        }
      }
      cat('\r')
      initial <- burnin 
      set1 <- seq(initial+n.thin +n.iter*(cad-1) ,n.iter*cad,n.thin)
      set <- c(set,set1)   
    }
    return(list(beta=beta[set,],sigma2=sigma2[set],U=U[set,],nu=nu[set],rho=rho[set],mu=media[set,]))
  }
}


