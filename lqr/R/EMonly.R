########################################################################
#EM ALGORITHM
########################################################################

EM <- function(y,x,p=0.5,dist = "normal",nu="",gama="",precision = 0.000001,envelope=FALSE){
  
  n <- length(y)
  d = dim(x)[2]
  beta<-solve(t(x)%*%x)%*%t(x)%*%y
  sigma2<-sum((y-x%*%beta)^2)/(n-d)
  teta <- c(beta,sigma2)
  beta.old <- beta
  sigma2.old <- sigma2
  nu.old = nu
  gama.old = gama
  if(dist == ""){dist = "normal"}
  if(nu.old=="" && dist == "t"){nu=4}
  if(nu.old=="" && dist == "slash"){nu=2}
  if(nu.old=="" && dist == "cont"){nu=0.1}
  if(gama.old=="" && dist == "cont"){gama=0.1}
  E = rep(1,n)
  pvec = rep(1,n)
  criterio <- 1
  count <- 0
    
  while(criterio > precision){
    count <- count + 1
    ### E-step: calculando ui, tui###
    
    if(dist == "t"){
      alpha1 = (nu+1)/2
      alpha2 = (nu+1)/2
      beta1  = (nu + ((4*(1-p)^2)/sigma2)*(y[y-x%*%beta<=0]-x[y-x%*%beta<=0,]%*%beta)^2)/2
      beta2  = (nu + ((4*p^2)/sigma2)*((y[y-x%*%beta>0]-x[y-x%*%beta>0,]%*%beta)^2))/2
      
      E[y-x%*%beta<=0] = alpha1/beta1
      E[y-x%*%beta>0]  = alpha2/beta2
    }
    
    if(dist == "laplace"){
      Gi1    = (((y[y-x%*%beta<=0]-x[y-x%*%beta<=0,]%*%beta)/sqrt(sigma2))^2)
      Gi2    = (((y[y-x%*%beta>0]-x[y-x%*%beta>0,]%*%beta)/sqrt(sigma2))^2)
      
      #require(ghyp)
      E[y-x%*%beta<=0] = Egig(lambda = 1/2, chi = 2*Gi1*(1-p)^2, psi = 1/2, func = "1/x")
      E[y-x%*%beta>0]  = Egig(lambda = 1/2, chi = 2*Gi2*(p^2), psi = 1/2, func = "1/x")
    }
    
    if(dist == "slash"){
      Gi1    = ((y[y-x%*%beta<=0]-x[y-x%*%beta<=0,]%*%beta)/sqrt(sigma2))^2
      Gi2    = ((y[y-x%*%beta>0]-x[y-x%*%beta>0,]%*%beta)/sqrt(sigma2))^2
      aux1   = nu + (3/2)
      aux2   = 2*Gi1*(1-p)^2
      aux3   = nu + (1/2)
      aux4   = 2*Gi2*(p)^2
      
      E[y-x%*%beta<=0] = (aux3/aux2)*pgamma(q = 1,shape = aux1,rate = aux2)/pgamma(q = 1,shape = aux3,rate = aux2)
      E[y-x%*%beta>0]  = (aux3/aux4)*pgamma(q = 1,shape = aux1,rate = aux4)/pgamma(q = 1,shape = aux3,rate = aux4)
    }
    
    if(dist == "cn"){
      mean1 =  x[y-x%*%beta<=0,]%*%beta
      mean2 =  x[y-x%*%beta>0,]%*%beta
      var1  = sigma2/(4*(1-p)^2)
      var2  = sigma2/(4*(p)^2)
      
      a    = nu*p*dnorm(x = y[y-x%*%beta<=0],mean = mean1,sd = sqrt(var1/gama))
      b    = (1-nu)*p*dnorm(x = y[y-x%*%beta<=0],mean = mean1,sd = sqrt(var1))
      
      c    = nu*(1-p)*dnorm(x = y[y-x%*%beta>0],mean = mean2,sd = sqrt(var2/gama))
      d    = (1-nu)*(1-p)*dnorm(x = y[y-x%*%beta>0],mean = mean2,sd = sqrt(var2))
      
      E[y-x%*%beta<=0] = ((gama*a)+b)/(a+b)
      E[y-x%*%beta>0]  = ((gama*c)+d)/(c+d)
    }
    
    ### M-step: atualizar mu, sigma2 ###
    
    pvec[y-x%*%beta<=0] = (1-p)^2
    pvec[y-x%*%beta>0] = (p)^2
    
    den = t(x)%*%diag(c(E*pvec))%*%x
    num = t(x)%*%diag(c(E*pvec))%*%y
    beta = solve(den)%*%(num)
    sigma2 = 4*(c(c(E*pvec)%*%((y-x%*%beta)^2)))/n
    
    ### Updating nu parameters
    
    if(nu.old=="" && dist == "t"){
      optnu  = optimize(f = loglikT,interval = c(2,100),maximum = T,x=y,mu=x%*%beta,sigma = sqrt(sigma2),p=p)
      loglik = optnu$objective
      nu     = optnu$maximum
    }
    
    if(nu.old=="" && dist == "slash"){
      optnu  = optimize(f = loglikSl,interval = c(0,10),maximum = T,x=y,mu=x%*%beta,sigma = sqrt(sigma2),p=p)
      loglik = optnu$objective
      nu     = optnu$maximum
    }
    
    if(nu.old!="" && gama.old =="" && dist == "cont"){
      optgamma  = optimize(f = loglikNC,interval = c(0,1),maximum = T,x=y,mu=x%*%beta,sigma = sqrt(sigma2),nu=nu,p=p)
      loglik = optgamma$objective
      gama  = optgamma$maximum
    }
    
    if(nu.old=="" && gama.old !="" && dist == "cont"){
      optnu  = optimize(f = loglikNC,interval = c(0,1),maximum = T,x=y,mu=x%*%beta,sigma = sqrt(sigma2),gama=0.1,p=p)
      loglik = optnu$objective
      nu     = optnu$maximum
    }

    if(nu.old=="" && gama.old=="" && dist == "cont"){
      opt  = optim(par = c(nu,gama),fn = AUXloglikNC,control = list(fnscale = -1),method = "L-BFGS-B",
                   lower = c(0.001,0.001),upper = c(0.999,0.999),x=y,mu=x%*%beta,
                   sigma = sqrt(sigma2),p=p)
      loglik = opt$value
      nu     = opt$par[1]
      gama  = opt$par[2]
    }
    
    param <- teta
    teta <- c(beta,sqrt(sigma2))
    criterio <- c(sqrt((teta-param)%*%(teta-param)))
  }
  
  #### Computing loglikelihood
  
  if(dist == "normal"){
    loglik =loglikN(y,x%*%teta[1:d],sigma = teta[d+1],p=p)
  }
  
  if(nu.old!="" && dist == "t"){
    loglik =loglikT(y,x%*%teta[1:d],sigma = teta[d+1],nu=nu,p=p)
  }
  
  if(dist == "laplace"){
    loglik =loglikL(y,c(x%*%teta[1:d]),sigma = teta[d+1],p=p)
  }
  
  if(nu.old!="" && dist == "slash"){
    loglik =loglikSl(y,x%*%teta[1:d],sigma = teta[d+1],nu=nu,p=p)
  }
  
  if(nu.old!="" && dist == "cont"){
    loglik =loglikNC(y,x%*%teta[1:d],sigma = teta[d+1],nu=nu,gama = gama,p=p)
  }
  
  npar = length(c(teta))
  AIC  = -2*loglik +2*npar
  BIC  = -2*loglik +log(n)*npar
  HQ   = -2*loglik +2*log(log(n))*npar
  
  ###### Computing SE's
  
  #Summarized score (sum_i)
  #   derbeta  = (2/sigma2)*t(x)%*%diag(c(E*pvec))%*%(y - x%*%beta)
  #   dersigma = c((2/sigma2)*t((y - x%*%beta))%*%diag(c(E*pvec))%*%(y - x%*%beta)) - n
  #   grad     = rbind(derbeta,dersigma)
  
  #Individual score (by i)
  derbetai  = (2/sigma2)*diag(c(E*pvec*(y - x%*%beta)))%*%x
  dersigmai = c((2/sigma2)*diag(c(E*pvec*(y - x%*%beta)))%*%(y - x%*%beta)) - 1
  gradi     = cbind(derbetai,dersigmai)
  
  IEI = 0
  for(i in 1:n)
  {
    IEI = gradi[i,]%*%t(gradi[i,]) + IEI
  }
  
  SE         = sqrt(diag(solve(IEI)))
  table      = data.frame(beta,SE[1:d],beta/SE[1:d],2*pnorm(abs(beta/SE[1:d]),lower.tail = F))
  asteriscos = apply(X = table[4],MARGIN = 1,FUN = defast)
  table      = data.frame(round(table,5),asteriscos)
  rownames(table) = paste("beta",1:d)
  colnames(table) = c("Estimate","Std. Error","z value","Pr(>|z|)","")
  
   ########ENVELOPES: NORMAL                       
  
  if(envelope==TRUE && dist == "normal"){
    n <-length(y)
    tau=p
    muc<- (y-x%*%teta[1:d])
    Ind<- (muc<0)+0  
    d2s<- (muc/(teta[d+1]))*(tau-Ind)
    d2s=sort(d2s)
    
    xq2 <- qnorm(seq(from = 0.525,to = 0.975,length.out = n),sd = 0.5)
    
    Xsim<-matrix(0,200,n)
    for(i in 1:200){
      Xsim[i,]<-abs(rnorm(n,sd = 0.5))
    }
    
    Xsim2<-apply(Xsim,1,sort)
    d21<-matrix(0,n,1)
    d22<-matrix(0,n,1)
    for(i in 1:n){
      d21[i]  <- quantile(Xsim2[i,],0.05)
      d22[i]  <- quantile(Xsim2[i,],0.95)
    }
    
    d2med <-apply(Xsim2,1,mean)
    
    fy <- range(d2s,d21,d22)
    plot(xq2,d2s,xlab = expression(paste("Theoretical ",HN(0.5), " quantiles")), 
         ylim = fy,ylab="Sample values and envelope",pch=20,cex=0.8,type="n",
         main = paste("Normal model for p =",p,"quantile",sep=" "))
    
    polygon(c(rev(xq2), xq2), c(rev(d22), d21), col = 'grey80', border = NA) 
    
    points(xq2,d2s,pch=20,cex=0.8)
    
    lines(xq2,d21,type="l",ylim=fy,xlab="",ylab="",lwd=0.5)
    lines(xq2,d2med,type="l",ylim=fy,xlab="",ylab="")
    lines(xq2,d22,type="l",ylim=fy,xlab="",ylab="",lwd=0.5)
    
  }
  
  ######## ENVELOPES: T                            
  
  if(envelope==TRUE && dist == "t"){
    muc<- (y-x%*%teta[1:d])
    Ind<- (muc<0)+0  
    d2s<- (muc/(teta[d+1]))*(p-Ind)
    d2s=sort(d2s)
    
    xq2 <- qt(seq(from = 0.525,to = 0.975,length.out = n),df = nu)*0.5
    
    Xsim<-matrix(0,200,n)
    for(i in 1:200){
      Xsim[i,]<-abs(rt(n,df = nu))*0.5
    }
    
    Xsim2<-apply(Xsim,1,sort)
    d21<-matrix(0,n,1)
    d22<-matrix(0,n,1)
    for(i in 1:n){
      d21[i]  <- quantile(Xsim2[i,],0.05)
      d22[i]  <- quantile(Xsim2[i,],0.95)
    }
    
    d2med <-apply(Xsim2,1,mean)
    
    fy <- range(d2s,d21,d22)
    plot(xq2,d2s,xlab = expression(paste("Theoretical ",HT(0.5,nu), " quantiles")), 
         ylim = fy,ylab="Sample values and envelope",pch=20,cex=0.8,type="n",
         main = paste("Student's t model for p =",p,"quantile",sep=" "))
    
    polygon(c(rev(xq2), xq2), c(rev(d22), d21), col = 'grey80', border = NA) 
    
    points(xq2,d2s,pch=20,cex=0.8)
    
    lines(xq2,d21,type="l",ylim=fy,xlab="",ylab="",lwd=0.5)
    lines(xq2,d2med,type="l",ylim=fy,xlab="",ylab="")
    lines(xq2,d22,type="l",ylim=fy,xlab="",ylab="",lwd=0.5)
  }
  
  ###### ENVELOPES: LAPLACE
  
  if(envelope==TRUE && dist == "laplace"){
    n <-length(y)
    tau = p
    muc<- (y-x%*%teta[1:d])
    Ind<- (muc<0)+0  
    d2s<- (muc/(teta[d+1]))*(tau-Ind)
    d2s=sort(d2s)
    
    xq2 <- qexp(ppoints(n), 2)
    
    Xsim<-matrix(0,200,n)
    for(i in 1:200){
      Xsim[i,]<-rexp(n,2)
    }
    
    Xsim2<-apply(Xsim,1,sort)
    d21<-matrix(0,n,1)
    d22<-matrix(0,n,1)
    for(i in 1:n){
      d21[i]  <- quantile(Xsim2[i,],0.05)
      d22[i]  <- quantile(Xsim2[i,],0.95)
    }
    
    d2med <-apply(Xsim2,1,mean)
    
    fy <- range(d2s,d21,d22)
    plot(xq2,d2s,xlab = expression(paste("Theoretical ",exp(0.5)," quantiles")),
    ylim = fy,ylab="Sample values and envelope",pch=20,cex=0.8,type="n",
    main = paste("Laplace model for p =",p,"quantile",sep=" "))

    polygon(c(rev(xq2), xq2), c(rev(d22), d21), col = 'grey80', border = NA) 
    
    points(xq2,d2s,pch=20,cex=0.8)
    
    lines(xq2,d21,type="l",ylim=fy,xlab="",ylab="",lwd=0.5)
    lines(xq2,d2med,type="l",ylim=fy,xlab="",ylab="")
    lines(xq2,d22,type="l",ylim=fy,xlab="",ylab="",lwd=0.5)
  }
  
  ######## ENVELOPES: SLASH                           
  
  if(envelope==TRUE && dist == "slash"){
    muc<- (y-x%*%teta[1:d])
    Ind<- (muc<0)+0  
    d2s<- (muc/(teta[d+1]))*(p-Ind)
    d2s=sort(d2s)
    
    #COMPUTE QUANTILES
    #require(spatstat)
    a1<-ewcdf(gendistSLK(n=50000,dist = "slash",nu=nu))    #empricial *weighted* cdf and quantile function
    xq2 <- quantile(a1,ppoints(n))     #calls quantile.ecdf()
    
    Xsim<-matrix(0,200,n)
    for(i in 1:200){
      Xsim[i,]<-gendistSLK(n = n,dist = "slash",nu=nu)
    }
    
    Xsim2<-apply(Xsim,1,sort)
    d21<-matrix(0,n,1)
    d22<-matrix(0,n,1)
    for(i in 1:n){
      d21[i]  <- quantile(Xsim2[i,],0.05)
      d22[i]  <- quantile(Xsim2[i,],0.95)
    }
    
    d2med <-apply(Xsim2,1,mean)
    
    fy <- range(d2s,d21,d22)
    plot(xq2,d2s,xlab = expression(paste("Theoretical ",HSl(0.5,nu), " quantiles")), 
         ylim = fy,ylab="Sample values and envelope",pch=20,cex=0.8,type="n",
         main = paste("Slash model for p =",p,"quantile",sep=" "))
    
    polygon(c(rev(xq2), xq2), c(rev(d22), d21), col = 'grey80', border = NA) 
    
    points(xq2,d2s,pch=20,cex=0.8)
    
    lines(xq2,d21,type="l",ylim=fy,xlab="",ylab="",lwd=0.5)
    lines(xq2,d2med,type="l",ylim=fy,xlab="",ylab="")
    lines(xq2,d22,type="l",ylim=fy,xlab="",ylab="",lwd=0.5)
  }
  
  ######## ENVELOPES: CONTAMINATED NORMAL                           
  
  if(envelope==TRUE && dist == "cont"){
    muc<- (y-x%*%teta[1:d])
    Ind<- (muc<0)+0  
    d2s<- (muc/(teta[d+1]))*(p-Ind)
    d2s=sort(d2s)
    
    #COMPUTE QUANTILES
    #require(spatstat)
    sec = seq(from = 0,to = 6,length.out = 100000)
    dd = densdist(di = sec,dens = densNC,nu=nu,gama=gama)
    a1<-ewcdf(x = sec,weights = dd/sum(dd))
    xq2 <- quantile(a1,ppoints(n))     #calls quantile.ecdf()
    
    Xsim<-matrix(0,200,n)
    for(i in 1:200){
      Xsim[i,]<-gendistSLK(n = n,dist = "cont",nu=nu,gama=gama)
    }
    
    Xsim2<-apply(Xsim,1,sort)
    d21<-matrix(0,n,1)
    d22<-matrix(0,n,1)
    for(i in 1:n){
      d21[i]  <- quantile(Xsim2[i,],0.05)
      d22[i]  <- quantile(Xsim2[i,],0.95)
    }
    
    d2med <-apply(Xsim2,1,mean)
    
    fy <- range(d2s,d21,d22)
    plot(xq2,d2s,xlab = expression(paste("Theoretical ",HCN(0.5,nu,gamma)," quantiles")), 
         ylim = fy,ylab="Sample values and envelope",pch=20,cex=0.8,type="n",
         main = paste("Contaminated Normal model for p =",p,"quantile",sep=" "))
    
    polygon(c(rev(xq2), xq2), c(rev(d22), d21), col = 'grey80', border = NA) 
    
    points(xq2,d2s,pch=20,cex=0.8)
    
    lines(xq2,d21,type="l",ylim=fy,xlab="",ylab="",lwd=0.5)
    lines(xq2,d2med,type="l",ylim=fy,xlab="",ylab="")
    lines(xq2,d22,type="l",ylim=fy,xlab="",ylab="",lwd=0.5)
  }
  
  #fitted and residuals values
  fitted.values    = c(x%*%teta[1:d])
  residuals        = y - fitted.values
   
  # OUTPUT
  if(dist == "normal" || dist == "laplace"){
    return(list(iter = count,criterio = criterio,theta = teta,table = table,SE = SE,loglik = loglik,AIC=AIC,BIC=BIC,HQ=HQ,fitted.values = fitted.values,residuals=residuals))
  }
  
  if(dist == "t" || dist == "slash"){
    return(list(iter =count,criterio = criterio,theta = teta, nu=nu,table = table,SE = SE,loglik = loglik,AIC=AIC,BIC=BIC,HQ=HQ,fitted.values = fitted.values,residuals=residuals))
  }
  
  if(dist == "cont"){
    return(list(iter =count,criterio = criterio,theta = teta, nu=nu,gamma = gama,table = table,SE = SE,loglik = loglik,AIC=AIC,BIC=BIC,HQ=HQ,fitted.values = fitted.values,residuals=residuals))
  }
}