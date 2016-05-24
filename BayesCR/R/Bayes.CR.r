Bayes.CR <- function(cc, x,y,cens="left",dist="Normal",criteria="FALSE",influence="FALSE",spacing="NULL",prior=NULL,hyper=NULL,n.thin=10,burnin=100,n.iter=2000,n.chains=2,chain="TRUE")
{

  if(dist=="SSL")
  {
    print("Calculations may take a while due to numerical integrations. You may want to disable options criteria and/or influence or specify a big number to spacing")
  }
  
  type <- dist
  if(cens==1){cens <- "left"}
  if(cens==2){cens <- "right"}
  
	namesx <- ('x1')
  if(ncol(as.matrix(x))>1){
  	for(i in 2:ncol(as.matrix(x))){namesx <- cbind(namesx, paste("x",i,"     ",sep=""))}
  }
	if(n.chains > 4) stop("The number of chains must be less than 5")
  if(ncol(as.matrix(y)) > 1) stop("Only univariate linear regression supported!")
  if(ncol(as.matrix(cc)) > 1) stop("Only univariate linear regression supported!")
	if( length(y) != nrow(as.matrix(x)) ) stop("X variable does not have the same number of lines than y")
	if( (length(x) == 0) | (length(y) == 0) ) stop("All parameters must be provided.")
	if( (dist != "T") && (dist != "Normal") && (dist != "Slash") && (dist != "NormalC")  && (dist != "SN")  && (dist != "ST")  && (dist != "SSL")) stop("Distribution family not supported. Check documentation!")
  if( (cens != "left") && (cens != "right")) stop("Censored type not supported. 'left' for left censoring and 'right' for right censoring.")
  if ((burnin >= n.iter)|| burnin < 1) 
  {
  stop("Invalid burnin number")
  }
  if (!is.numeric(n.iter) || n.iter < 1) 
  {
  stop("Invalid number of iterations")
  }
  if (n.thin >= n.iter) 
  {
  stop("Invalid number of lag (n.thin) for posterior sample")
  }
  
  if(type == "T")
  {
    if(length(prior) == 0) stop("The nu prior must be provided.")
    if( (prior != "Exp") && (prior != "Jeffreys") && (prior != "Hierar") && (prior != "Unif")) stop("Prior distribution not supported. Check documentation.")
    if(prior == "Exp")
    { 
      if(length(hyper) == 0) stop("The hyper parameter must be provided.")
      if(!is.numeric(hyper)) stop("Invalid value of hyper parameter")
      if(hyper<=0) stop("The hyper parameter must be positive.")  
    }
  }

  M <- n.iter
	p <- ncol(as.matrix(x))
  namespar <- colnames(x) 
  colx <- ncol(as.matrix(x))
  if(length(namespar)==0)namespar <- namesx[1:colx]
  
  
  if(dist == "Normal")
  {
    out <- GibbsTruncSMN(cc,y,x,n.iter,n.thin,burnin,type=type, cens=cens, prior, hyper, n.chains=n.chains)
    
    betas.f <- as.matrix(apply(out$beta,2,mean))
    sigma2.f <- mean(out$sigma2)
    param <- rbind(betas.f,sigma2.f)
    
    se.betas <- as.matrix(apply(out$beta,2,sd))
    se.sigma2 <- sd(out$sigma2)
    se <- rbind(se.betas,se.sigma2)
    
    HPD <- matrix(0,nrow=p,ncol=2)
    for (i in 1:p)
    {
      HPD[i,]<- hpd(out$beta[,i],alpha=0.05)  
    }
    HPDS2  <- hpd(out$sigma2,alpha=0.05)
    HPDTot <- rbind(HPD,HPDS2)
    
    ver<-LogVerosCens(cc,y,apply(out$mu,2,mean),mean(out$sigma2),0,1,type=dist, cens=cens)
    
    if(n.chains>1)
    {
      RBeta <- Rhat1(out$beta,n.iter,burnin,n.chains,n.thin) 
      Rsigma <- Rhat1(out$sigma2,n.iter,burnin,n.chains,n.thin) 
      RhatFin <- rbind(RBeta,Rsigma)
      
      paramT <- round(cbind(param, se, HPDTot,RhatFin),digits=5)
      dimnames(paramT) <- list(c(namespar,expression(sigma^2)),c("Mean", "Sd", " HPD(95%)","","Rhat"))
    }else{
      paramT <- round(cbind(param, se, HPDTot),digits=5)
      dimnames(paramT) <- list(c(namespar,expression(sigma^2)),c("Mean", "Sd", " HPD(95%)",""))
    }
  }
  if((dist == "T") | (dist == "Slash") ) 
  {
    out <- GibbsTruncSMN(cc,y,x,n.iter,n.thin,burnin,type=type, cens=cens, prior, hyper, n.chains=n.chains)
    
    betas.f <- as.matrix(apply(out$beta,2,mean))
    sigma2.f <- mean(out$sigma2)
    nu.f <- mean(out$nu)
    param <- rbind(betas.f,sigma2.f,nu.f)
    
    se.betas <- as.matrix(apply(out$beta,2,sd))
    se.sigma2 <- sd(out$sigma2)
    se.nu <- sd(out$nu)
    se <- rbind(se.betas,se.sigma2,se.nu)
    
    HPD <- matrix(0,nrow=p,ncol=2)
    for (i in 1:p)
    {
      HPD[i,]<- hpd(out$beta[,i],alpha=0.05)  
    }
    HPDS2  <- hpd(out$sigma2,alpha=0.05)
    HPDnu  <- hpd(out$nu,alpha=0.05) 
    HPDTot <- rbind(HPD,HPDS2,HPDnu)
    
    mu.f <- mean(out$mu)
    ver<-LogVerosCens(cc,y,apply(out$mu,2,mean),mean(out$sigma2),0,mean(out$nu),type=dist, cens=cens)
    
    
    if(n.chains>1)
    {
      RBeta <- Rhat1(out$beta,n.iter,burnin,n.chains,n.thin) 
      Rsigma <- Rhat1(out$sigma2,n.iter,burnin,n.chains,n.thin)
      Rnu <- Rhat1(out$nu,n.iter,burnin,n.chains,n.thin)
      RhatFin <- rbind(RBeta,Rsigma,Rnu)
      
      paramT <- round(cbind(param, se, HPDTot,RhatFin),digits=5)
      dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(nu)),c("Mean", "Sd", " HPD(95%)","", "Rhat"))
    }else{
      paramT <- round(cbind(param, se, HPDTot),digits=5)
      dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(nu)),c("Mean", "Sd", " HPD(95%)",""))
    }
  }
  if(dist=="NormalC")
  {
    out <- GibbsTruncSMN(cc,y,x,n.iter,n.thin,burnin,type=type, cens=cens, prior, hyper, n.chains=n.chains)
   
    betas.f <- as.matrix(apply(out$beta,2,mean))
    sigma2.f <- mean(out$sigma2)
    nu.f <- mean(out$nu)
    rho.f <- mean(out$rho)
    param <- rbind(betas.f,sigma2.f,nu.f,rho.f)
    
    se.betas <- as.matrix(apply(out$beta,2,sd))
    se.sigma2 <- sd(out$sigma2)
    se.nu <- sd(out$nu)
    se.rho <- sd(out$rho)
    se <- rbind(se.betas,se.sigma2,se.nu,se.rho)
    
    HPD <- matrix(0,nrow=p,ncol=2)
    for (i in 1:p)
    {
      HPD[i,]<- hpd(out$beta[,i],alpha=0.05)  
    }
    HPDS2  <- hpd(out$sigma2,alpha=0.05)
    HPDnu  <- hpd(out$nu,alpha=0.05)
    HPDrho <- hpd(out$rho,alpha=0.05)
    HPDTot <- rbind(HPD,HPDS2,HPDnu,HPDrho)
    
    ver<-LogVerosCens(cc,y,apply(out$mu,2,mean),mean(out$sigma2),0,c(mean(out$nu),mean(out$rho)),type=dist, cens=cens)
    
    if(n.chains>1)
    {
      RBeta <- Rhat1(out$beta,n.iter,burnin,n.chains,n.thin) 
      Rsigma <- Rhat1(out$sigma2,n.iter,burnin,n.chains,n.thin)
      Rrho <- Rhat1(out$rho,n.iter,burnin,n.chains,n.thin)
      Rnu <- Rhat1(out$nu,n.iter,burnin,n.chains,n.thin)
      RhatFin <- rbind(RBeta,Rsigma,Rnu,Rrho)
      
      paramT <- round(cbind(param, se, HPDTot,RhatFin),digits=5)
      dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(nu),expression(rho)),c("Mean", "Sd", " HPD(95%)","", "Rhat"))
    }else{
      paramT <- round(cbind(param, se, HPDTot),digits=5)
      dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(nu),expression(rho)),c("Mean", "Sd", " HPD(95%)",""))
    }
  }
  
	if(dist == "SN")
	{
	  out <- Gibbs(y,x,cc,dist,cens,n.iter,burnin,n.thin,n.chains)
    
	  betas.f <- as.matrix(apply(out$beta,2,mean))
	  sigma2.f <- mean(out$sigma2)
	  lambda.f <- mean(out$lambda)
	  param <- rbind(betas.f,sigma2.f,lambda.f)
    
	  se.betas <- as.matrix(apply(out$beta,2,sd))
	  se.sigma2 <- sd(out$sigma2)
	  se.lambda <- sd(out$lambda)
	  se <- rbind(se.betas,se.sigma2,se.lambda)
    
	  HPD <- matrix(0,nrow=p,ncol=2)
	  for (i in 1:p)
	  {
	    HPD[i,]<- hpd(out$beta[,i],alpha=0.05)  
	  }
	  HPDS2  <- hpd(out$sigma2,alpha=0.05)
	  HPDlambda  <- hpd(out$lambda,alpha=0.05)
	  HPDTot <- rbind(HPD,HPDS2,HPDlambda)
	  
	  ver<-LogVerosCens(cc,y,apply(out$mu,2,mean),mean(out$sigma2),mean(out$lambda),0,type="SN",cens=cens)
	  
	  if(n.chains>1)
	  {
	    RBeta <- Rhat1(out$beta,n.iter,burnin,n.chains,n.thin) 
	    Rsigma <- Rhat1(out$sigma2,n.iter,burnin,n.chains,n.thin)
	    Rlambda <- Rhat1(out$lambda,n.iter,burnin,n.chains,n.thin)
	    RhatFin <- rbind(RBeta,Rsigma,Rlambda) 
      
	    paramT <- round(cbind(param, se, HPDTot,RhatFin),digits=5)
	    dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(lambda)),c("Mean", "Sd", " HPD(95%)","", "Rhat"))
	    
	  }else{
	    paramT <- round(cbind(param, se, HPDTot),digits=5)
	    dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(lambda)),c("Mean", "Sd", " HPD(95%)",""))
	  }
	}
  
  if((dist == "ST") | (dist == "SSL"))
  {
    out <- Gibbs(y,x,cc,dist,cens,n.iter,burnin,n.thin,n.chains)  
    betas.f <- as.matrix(apply(out$beta,2,mean))
    sigma2.f <- mean(out$sigma2)
    lambda.f <- mean(out$lambda)
    nu.f <- mean(out$nu)
    param <- rbind(betas.f,sigma2.f,lambda.f,nu.f)
    
    se.betas <- as.matrix(apply(out$beta,2,sd))
    se.sigma2 <- sd(out$sigma2)
    se.lambda <- sd(out$lambda)
    se.nu <- sd(out$nu)
    se <- rbind(se.betas,se.sigma2,se.lambda,se.nu)
    
    HPD <- matrix(0,nrow=p,ncol=2)
    for (i in 1:p)
    {
      HPD[i,]<- hpd(out$beta[,i],alpha=0.05)  
    }
    HPDS2  <- hpd(out$sigma2,alpha=0.05)
    HPDlambda  <- hpd(out$lambda,alpha=0.05) 
    HPDnu <- hpd(out$nu,alpha=0.05)
    HPDTot <- rbind(HPD,HPDS2,HPDlambda,HPDnu)
    
    ver <- LogVerosCens(cc,y,apply(out$mu,2,mean),mean(out$sigma2),mean(out$lambda),mean(out$nu),type=dist,cens=cens)
    
    if(n.chains>1)
    {
      RBeta <- Rhat1(out$beta,n.iter,burnin,n.chains,n.thin) 
      Rsigma <- Rhat1(out$sigma2,n.iter,burnin,n.chains,n.thin)
      Rlambda <- Rhat1(out$lambda,n.iter,burnin,n.chains,n.thin)
      Rnu <- Rhat1(out$nu,n.iter,burnin,n.chains,n.thin)
      RhatFin <- rbind(RBeta,Rsigma,Rlambda,Rnu)
      
      paramT <- round(cbind(param, se, HPDTot,RhatFin),digits=5)
      dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(lambda),expression(nu)),c("Mean", "Sd", " HPD(95%)","", "Rhat"))
    }else{
      paramT <- round(cbind(param, se, HPDTot),digits=5)
      dimnames(paramT) <- list(c(namespar,expression(sigma^2),expression(lambda),expression(nu)),c("Mean", "Sd", " HPD(95%)",""))
    }
  }   
  
  logver <- sum(log(ver))
  
  cat('\n') 
  cat('-------------------------------------------------------------\n')
  cat('Posterior mean(Mean), standard deviation(Sd) and HPD interval\n')
  cat('-------------------------------------------------------------\n')
  print(paramT)
  cat('-------------------------------------------------------------\n')
  cat('Log-likelihood')
  print(logver)
  cat('-------------------------------------------------------------\n')
  cat('\r \n')
  x <- as.matrix(x)
  param <- ncol(x) 
  if(criteria=="TRUE")
  {
    if(influence=="FALSE")
    {
      crit <- criterios(cc,y,espac=spacing,cadeia=out,type=type, cens=cens, p=length(param),influence="FALSE") 
    }else{
      crit <- criterios(cc,y,espac=spacing,cadeia=out,type=type, cens=cens, p=length(param),influence="TRUE")
      KL <- crit$KL
      JDist <- crit$JDist
      LDist <- crit$LDist
      ChiDist <- crit$ChiDist
    }
    critFin <- c(crit$CPO, crit$DIC, crit$EAIC, crit$EBIC,crit$WAIC1,crit$WAIC2)
    critFin <- round(t(as.matrix(critFin)),digits=3)
    dimnames(critFin) <- list(c("Value"),c("LPML", "DIC", "EAIC","EBIC","WAIC1","WAIC2")) 
    cat('\n') 
    cat('Model selection criteria\n')
    cat('-------------------------------------------------------------\n')
    print(critFin)
    cat('-------------------------------------------------------------\n')
    cat('\r \n')
  }else{
    if(influence=="TRUE")
    {
      crit <- criterios(cc,y,espac=spacing,cadeia=out,type=type, cens=cens, p=length(param),influence="TRUE")
      KL <- crit$KL
      JDist <- crit$JDist
      LDist <- crit$LDist
      ChiDist <- crit$ChiDist
    }
  }
 
  if(influence=="TRUE")
  {
    if(chain=="TRUE")
    {
      if(dist=="Normal")
      {
        return(list(beta=out$beta[,],sigma2=out$sigma2,KL=KL,JDist=JDist,LDist=LDist,ChiDist=ChiDist))   
      }
      if((dist == "T") | (dist == "Slash"))
      {
        return(list(beta=out$beta[,],sigma2=out$sigma2,nu=out$nu,KL=KL,JDist=JDist,LDist=LDist,ChiDist=ChiDist))
      }
      if(dist=="NormalC")
      {
        return(list(beta=out$beta[,],sigma2=out$sigma2,nu=out$nu,rho=out$rho,KL=KL,JDist=JDist,LDist=LDist,ChiDist=ChiDist))   
      }
      if((dist == "ST") | (dist == "SSL"))
      {
        return(list(beta=out$beta[,],sigma2=out$sigma2,lambda=out$lambda,nu=out$nu,KL=KL,JDist=JDist,LDist=LDist,ChiDist=ChiDist))
      }
      if(dist == "SN")
      {
        return(list(beta=out$beta[,],sigma2=out$sigma2,lambda=out$lambda,KL=KL,JDist=JDist,LDist=LDist,ChiDist=ChiDist)) 
      }
    }else{
      return(list(KL=KL,JDist=JDist,LDist=LDist,ChiDist=ChiDist))  
    }  
  }else{
    if(chain=="TRUE")
    {
      if(dist=="Normal")
      {
        return(list(beta=out$beta[,],sigma2=out$sigma2))   
      }
      if((dist == "T") | (dist == "Slash"))
      {
        return(list(beta=out$beta[,],sigma2=out$sigma2,nu=out$nu))
      }
      if(dist=="NormalC")
      {
        return(list(beta=out$beta[,],sigma2=out$sigma2,nu=out$nu,rho=out$rho))   
      }
      if((dist == "ST") | (dist == "SSL"))
      {
        return(list(beta=out$beta[,],sigma2=out$sigma2,lambda=out$lambda,nu=out$nu))
      }
      if(dist == "SN")
      {
        return(list(beta=out$beta[,],sigma2=out$sigma2,lambda=out$lambda)) 
      }
    }
  }      
}




