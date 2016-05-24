mmeta <-
  function(data, rhow, type, k, method) {
    
    nn.reml <-
      function(y1, s1, y2, s2, rhow) {
        S=cbind(s1^2, rhow*s1*s2, s2^2)
        myreml= mvmeta(cbind(y1, y2), S, method="reml")
        coef = coef(myreml)[1,]
        vcov = vcov(myreml)
        rhob <- cov2cor(myreml$Psi)[2,1]    
        myresults = list(coefficients = coef, vcov= vcov, rhob=rhob)
        return(myresults)
      }
    
    nn.cl <-
      function(y1, s1, y2, s2) {
        m=length(y1)
        estim.pseudo=rep(NA, 4)
        my.pseudo1= mvmeta(y1,s1^2,method="reml")
        my.pseudo2= mvmeta(y2,s2^2,method="reml")
        estim.pseudo[c(1:2)]=c(coef(my.pseudo1), my.pseudo1$Psi)
        estim.pseudo[c(3:4)]=c(coef(my.pseudo2), my.pseudo2$Psi)
        estim.pseudo = matrix(estim.pseudo, nrow=1)
        colnames(estim.pseudo)=c("beta1", "tau1.2", "beta2", "tau2.2")
        Score1.beta = (y1-coef(my.pseudo1))/(s1^2+my.pseudo1$Psi)
        Score2.beta = (y2-coef(my.pseudo2))/(s2^2+my.pseudo2$Psi)
        Score1.tau2 = -1/(2*(s1^2+my.pseudo1$Psi)) + (y1-coef(my.pseudo1))^2/(s1^2+my.pseudo1$Psi)^2/2
        Score2.tau2 = -1/(2*(s2^2+my.pseudo2$Psi)) + (y2-coef(my.pseudo2))^2/(s2^2+my.pseudo2$Psi)^2/2
        Score1 = rbind(Score1.beta, Score1.tau2); Score2 = rbind(Score2.beta, Score2.tau2)
        I11.hat = Score1%*%t(Score1)/m; I22.hat = Score2%*%t(Score2)/m
        I12.hat = Score1%*%t(Score2)/m
        myoff.diag = solve(I11.hat,tol=1e-100)%*%I12.hat%*%solve(I22.hat,tol=1e-100)/m
        myupper = cbind(solve(m*I11.hat,tol=1e-100), (myoff.diag))
        mylower = cbind(t(myoff.diag), solve(m*I22.hat,tol=1e-100))
        myV=rbind(myupper,mylower)
        colnames(myV)=c("beta1", "tau1.2", "beta2", "tau2.2")
        rownames(myV)=c("beta1", "tau1.2", "beta2", "tau2.2")
        
        myresults = list(coefficients = estim.pseudo, vcov= myV)
        return(myresults)
      }
    
    
    nn.mom= function(y1, s1, y2, s2){
      v1 = s1^2
      v2 = s2^2
      w1 = 1/(v1)
      w2 = 1/(v2)
      y1.weight = sum(w1*y1)/sum(w1)
      y2.weight = sum(w2*y2)/sum(w2)
      n1 = sum(1-1*(v1 > 10^4))  
      n2 = sum(1-1*(v2 > 10^4)) 
      Q1 = sum(w1*(y1-y1.weight)^2)
      Q2 = sum(w2*(y2-y2.weight)^2)  	   
      tau1.2.hat = max(0, (Q1-(n1-1))/(sum(w1)-sum(w1^2)/sum(w1)))
      tau2.2.hat = max(0, (Q2-(n2-1))/(sum(w2)-sum(w2^2)/sum(w2)))	   
      ## variance estimate:
      w1.star = 1/(v1 + tau1.2.hat)
      w2.star = 1/(v2 + tau2.2.hat)	
      beta1.hat = sum(w1.star*y1)/sum(w1.star)
      beta2.hat = sum(w2.star*y2)/sum(w2.star)	
      var.beta1.hat = 1/sum(w1.star)
      var.beta2.hat = 1/sum(w2.star) 
      mycov.beta = sum((w1.star/sum(w1.star))*(w2.star/sum(w2.star))*(y1 - beta1.hat)*(y2 - beta2.hat))
      beta.hat = cbind(beta1=beta1.hat,beta2=beta2.hat)
      myV = matrix(c(var.beta1.hat,mycov.beta,mycov.beta,var.beta2.hat),nrow = 2, byrow = T)
      colnames(myV)=c("beta1", "beta2")
      rownames(myV)=c("beta1", "beta2")
        myresults = list(coefficients = beta.hat, vcov= myV)
        return(myresults)
    }
    
    
    log.Riley.Lik.reml.func = function(mygamma, y1, s1, y2, s2){
      m = length(y1)
      beta1 = mygamma[1]; beta2= mygamma[2]
      log.psi1.2 = mygamma[3]; log.psi2.2 = mygamma[4]; omega = mygamma[5]
      psi1.2 = exp(log.psi1.2); psi2.2 = exp(log.psi2.2)
      rho = 2*plogis(omega)-1
      s1.2 = s1^2
      s2.2 = s2^2
      Sigma = matrix(0, nrow=(2*m), ncol=(2*m))
      for(i in 1:m){
        cov.y1.y2 = rho*sqrt((psi1.2+s1.2[i])*(psi2.2+s2.2[i]))
        Sigma.i = matrix(c(s1.2[i]+psi1.2, cov.y1.y2, cov.y1.y2, s2.2[i]+psi2.2), nrow=2)
        Sigma[((i-1)*2+c(1:2)), ((i-1)*2+c(1:2))] = Sigma.i
      }
      y1.temp = rep(y1, each=2); y2.temp = rep(y2, each=2)
      y1.temp2 = y1.temp%*%diag(rep(c(1,0), times=m))
      y2.temp2 = y2.temp%*%diag(rep(c(0,1), times=m))
      Y = as.vector(y1.temp2+y2.temp2)
      x.design=cbind(rep(c(1,0), times=m), rep(c(0,1),times=m))
      mybeta=x.design%*%c(beta1,beta2)
      mylik = -0.5*((m-2)*log(2*pi)-log(det(t(x.design)%*%x.design)) +
                      log(det(Sigma))+log(det(t(x.design)%*%solve(Sigma,tol=1e-100)%*%x.design)) +
                      t(Y-mybeta)%*%solve(Sigma,tol=1e-100)%*%(Y-mybeta))
      mylik = as.numeric(mylik)
      return(mylik)
    }
    
    log.Riley.Lik.reml.i.func=function(mygamma, y1.i, s1.i, y2.i, s2.i){
      result=log.Riley.Lik.reml.func(mygamma, y1.i, s1.i, y2.i, s2.i)
      return(result)
    }
    
    robustV.ftn=function(mygamma, y1, s1, y2, s2){
      m = length(y1)
      SigmaB=array(NA, c(5,5,m))
      SigmaA=array(NA, c(5,5,m))
      for (i in 1:m){
        y1.i=y1[i]
        s1.i=s1[i]
        y2.i=y2[i]
        s2.i=s2[i]
        myB=jacobian(log.Riley.Lik.reml.i.func, mygamma, y1.i=y1.i, s1.i=s1.i, y2.i=y2.i, s2.i=s2.i)
        SigmaB[,,i]=t(myB)%*%myB
        myA=hessian(log.Riley.Lik.reml.i.func, mygamma, y1.i=y1.i, s1.i=s1.i, y2.i=y2.i, s2.i=s2.i)
        SigmaA[,,i]=-myA
      }
      E.A=apply(SigmaA, c(1:2), mean) # sensitivity matrix
      E.B=apply(SigmaB, c(1:2), mean) # variability matrix
      myV=solve(E.A)%*%E.B%*%solve(E.A)
      return(myV)
    }
    
    nn.rs=function(y1, s1, y2, s2){
      mygamma.init = c(0,0,0,0,0)
      myresults = optim(mygamma.init, log.Riley.Lik.reml.func, y1=y1, s1=s1, y2=y2, s2=s2, hessian=TRUE,
                        control = list(fnscale=-1,maxit=1000))
      beta1=myresults$par[1]
      tau1.2=exp(myresults$par[3])
      beta2=myresults$par[2]
      tau2.2=exp(myresults$par[4])
      rhos=2*plogis(myresults$par[5])-1
      par=cbind(beta1, tau1.2, beta2, tau2.2, rhos)
      myV=robustV.ftn(myresults$par, y1=y1, s1=s1, y2=y2, s2=s2)/length(y1)
      colnames(myV)=c("beta1", "log.tau1.2", "beta2", "log.tau2.2", "omega")
      rownames(myV)=c("beta1", "log.tau1.2", "beta2", "log.tau2.2", "omega")
      
      return(list(coefficients=par, vcov=myV))
    }
    
    
    myLik.indep.log=function(mypar, mydat) {
      a1.temp <- mypar[1]; b1.temp <- mypar[2]
      a2.temp <- mypar[3]; b2.temp <- mypar[4]
      
      a1 <- exp(a1.temp); b1 <- exp(b1.temp)
      a2 <- exp(a2.temp); b2 <- exp(b2.temp)
      
      temp1 <- (lgamma(a1+mydat$y1) + lgamma(b1+mydat$n1-mydat$y1)
                + lgamma(a2+mydat$y2) + lgamma(b2+mydat$n2-mydat$y2)
                + lgamma(a1+b1) + lgamma(a2+b2))
      temp2 <- (lgamma(a1) + lgamma(b1) + lgamma(a2) + lgamma(b2)
                + lgamma(a1+b1+mydat$n1) + lgamma(a2+b2+mydat$n2))
      
      myLogLik <- sum(temp1 - temp2)
      return(myLogLik)
    }
    
    par.cal=function(mypar) {
      a1 <- exp(mypar[1]); b1 <- exp(mypar[2])
      a2 <- exp(mypar[3]); b2 <- exp(mypar[4])
      eta <- mypar[5]
      cc <- sqrt(a1*a2*b1*b2)/sqrt((a1+b1+1)*(a2+b2+1))
      upper.bound <- cc/max(a1*b2, a2*b1)
      lower.bound <- -cc/max(a1*a2, b1*b2)
      expit.eta= exp(eta)/(1+exp(eta))
      rho <- (upper.bound-lower.bound)*expit.eta + lower.bound
      return(c(a1,b1,a2,b2,rho,eta))
    }
    
    bb.cl=function(y1,n1,y2,n2){
      nstudy=length(y1)
      
      initial.val.gen <- function(y, n) {
        BBfit <- betabin(cbind(y, n-y)~1, ~1, data=data.frame(y=y,n=n))
        expit.BB <- exp(as.numeric(BBfit@param[1]))/(1+exp(as.numeric(BBfit@param[1])))
        
        a.ini <- expit.BB*(1/as.numeric(BBfit@param[2])-1)
        b.ini <- (1/as.numeric(BBfit@param[2])-1)*(1-expit.BB)
        return(list(a=a.ini, b=b.ini))
      }
      
      mle.CL <- function(y1=y1,n1=n1,y2=y2,n2=n2) {
        init.val <- rep(0, 5)
        fit1 <- initial.val.gen(y1, n1)
        fit2 <- initial.val.gen(y2, n2)
        init.val[1] <- log(fit1$a);
        init.val[2] <- log(fit1$b)
        init.val[3] <- log(fit2$a);
        init.val[4] <- log(fit2$b)
        
        MLE.inde.log <- optim(init.val[1:4], myLik.indep.log, method = "L-BFGS-B",
                              lower=rep(-20,4), upper=rep(20,4),
                              control = list(fnscale=-1,maxit=1000),
                              hessian = T, mydat=list(y1=y1,n1=n1,y2=y2,n2=n2))
        
        mypar<- par.cal(MLE.inde.log$par)
        rho<-0;
        eta<-NA;
        hessian.log<-MLE.inde.log$hessian
        colnames(hessian.log)<-c("loga1","logb1","loga2","logb2")
        rownames(hessian.log)<-c("loga1","logb1","loga2","logb2")
        conv=MLE.inde.log$convergence
        
        a1 <- mypar[1]; b1 <- mypar[2];
        a2 <- mypar[3]; b2 <- mypar[4];
        prior.MLE<-c(a1, b1, a2, b2, rho ,eta)
        return(list(MLE=prior.MLE,hessian.log=hessian.log,conv=conv))
      }
      
      
      myLik.indep.vector=function(mypar, mydat) {
        a1.temp <- mypar[1]; b1.temp <- mypar[2]
        a2.temp <- mypar[3]; b2.temp <- mypar[4]
        
        a1 <- exp(a1.temp); b1 <- exp(b1.temp)
        a2 <- exp(a2.temp); b2 <- exp(b2.temp)
        
        temp1 <- (lgamma(a1+mydat$y1) + lgamma(b1+mydat$n1-mydat$y1)
                  + lgamma(a2+mydat$y2) + lgamma(b2+mydat$n2-mydat$y2)
                  + lgamma(a1+b1) + lgamma(a2+b2))
        temp2 <- (lgamma(a1) + lgamma(b1) + lgamma(a2) + lgamma(b2)
                  + lgamma(a1+b1+mydat$n1) + lgamma(a2+b2+mydat$n2))
        
        myLogLik <- (temp1 - temp2)
        return(myLogLik)
      }
      
      sandwich.var.cal=function(mypar, myhessian.indep, mydat){
        npar=length(mypar)
        delta=1e-8
        myD = matrix(NA, nrow=npar, ncol=nstudy)
        for(i in 1:npar){
          par.for=par.back=mypar
          par.for[i]=par.for[i]+delta
          par.back[i]=par.back[i]-delta
          myD[i,] = (myLik.indep.vector(par.for, mydat)-myLik.indep.vector(par.back,mydat))/(2*delta)
        }
        score.squared = myD%*%t(myD)
        inv.hessian = solve(-myhessian.indep)
        results=inv.hessian%*% score.squared %*%inv.hessian
        return(results)
      }
      
      out=mle.CL(y1=y1,n1=n1,y2=y2,n2=n2)
      mle=out$MLE
      hessian.log=out$hessian.log
      a1=mle[1];b1=mle[2]
      a2=mle[3];b2=mle[4]
      mypar=log(c(a1,b1,a2,b2))
      mydat=list(y1=y1,n1=n1,y2=y2,n2=n2)
      covar.sandwich=sandwich.var.cal(mypar, hessian.log, mydat)
      myD <- matrix(c(-1, 1, 1, -1), nrow=1)
      myVar.log <- as.numeric(myD %*% covar.sandwich %*% t(myD))
      logOR <- log((a2/b2)/(a1/b1))
      OR_CI=c(exp(logOR-1.96*sqrt(myVar.log)), exp(logOR+1.96*sqrt(myVar.log)))
      #return(list(logOR=logOR, OR_CI=OR_CI, se=sqrt(myVar.log),hessian.log=hessian.log,conv=out$conv,mle=mle))
      return(list(logOR=logOR, OR_CI=OR_CI, se=sqrt(myVar.log),hessian.log=hessian.log,conv=out$conv,mle=mle))
    }
    
    
    logit = function(p) log(p/(1-p))
    expit = function(b) exp(b)/(1+exp(b))
    
    dlogitnorm = function(p, mu, tau2){
      (2*pi*tau2)^(-1/2)*exp(-(qlogis(p)-mu)^2/(2*tau2))/(p*(1-p))
    }
    
    integrand1<-function(n,y, mypar2, p.grid){
      dbinom(y, size=n, prob=p.grid)*dlogitnorm(p.grid, mypar2[1], mypar2[2])
    }
    
    integrand2<-function(n,y,mypar2, p.grid){
      dbinom(y, size=n, prob=p.grid)*dlogitnorm(p.grid, mypar2[1], mypar2[2])*2*(qlogis(p.grid)-mypar2[1])/(2*mypar2[2])
    }
    
    integrand3<-function(n,y,mypar2, p.grid){
      dbinom(y, size=n, prob=p.grid)*dlogitnorm(p.grid, mypar2[1], mypar2[2])*(-1/(2*mypar2[2])+(qlogis(p.grid)-mypar2[1])^2/(2*mypar2[2]^2))
    }
    
    mystandize = function(MM){
      diag(1/sqrt(diag(MM)))%*%MM%*%diag((1/sqrt(diag(MM))))
    }
    
    bn.cl=function(y1, n1, y2, n2){
      id=c(1:length(y1))
      mydat1 = data.frame(id,n1,y1,n2,y2)
      names(mydat1) = c("id","Nd","SeY","Nn","SpY")
      
      m = nrow(mydat1) 
      
      estim = estim.orig = mbse.orig = matrix(NA, nrow=4, ncol=1)
      mySandwich = matrix(NA,nrow=4,ncol=4)
      
      fit.SeY = glmmML(cbind(SeY, Nd-SeY)~1, data = mydat1, cluster = id, method= "ghq", n.points=8, 	control=list(epsilon=1e-08 , maxit=50000,trace=FALSE))
      fit.SpY = glmmML(cbind(SpY, Nn-SpY)~1, data = mydat1, cluster = id, method= "ghq", n.points=8, control=list(epsilon=1e-08 , maxit=50000,trace=FALSE))
      
      estim[c(1: 2), 1] = c(as.numeric(fit.SeY$coefficients), (fit.SeY$sigma)^2)
      estim[c(3: 4), 1] = c(as.numeric(fit.SpY$coefficients), (fit.SpY$sigma)^2)
      
      estim.orig[c(1: 2), 1] = c(plogis(estim[1, 1]), estim[2, 1])
      estim.orig[c(3: 4), 1] = c(plogis(estim[3, 1]), estim[4, 1])
      
      rownames(estim)=c("SeY", "SeY.S.2", "SpY", "SpY.S.2")
      rownames(estim.orig)=c("SeY", "SeY.S.2", "SpY", "SpY.S.2")
      
      I11.hat = solve(m*fit.SeY$variance)
      I22.hat = solve(m*fit.SpY$variance)
      
      int1 = int2  = int3 = rep(NA, length=m)
      est1 = estim[c(1: 2), 1] ##plug in the point estimate based on sensitivity
      for(i in 1: m){
        int1[i] = integrate(integrand1, lower=0, upper=1, n=mydat1$Nd[i], mypar2=est1, y=mydat1$SeY[i])$value
        int2[i] = integrate(integrand2, lower=0, upper=1, n=mydat1$Nd[i], mypar2=est1, y=mydat1$SeY[i])$value
        int3[i] = integrate(integrand3, lower=0, upper=1, n=mydat1$Nd[i], mypar2=est1, y=mydat1$SeY[i])$value
      }
      
      deriveSeY.mu = int2/int1
      deriveSeY.tau2 = int3/int1
      B.SeY = cbind(deriveSeY.mu, deriveSeY.tau2)
      
      int1 = int2  = int3 = rep(NA, length=m)
      est2 = estim[c(3: 4), 1] ##plug in the point estimate based on specificity
      for(i in 1: m){
        int1[i] = integrate(integrand1, lower=0, upper=1, n=mydat1$Nn[i], mypar2=est2, y=mydat1$SpY[i])$value
        int2[i] = integrate(integrand2, lower=0, upper=1, n=mydat1$Nn[i], mypar2=est2, y=mydat1$SpY[i])$value
        int3[i] = integrate(integrand3, lower=0, upper=1, n=mydat1$Nn[i], mypar2=est2, y=mydat1$SpY[i])$value
      }
      
      deriveSpY.mu = int2/int1
      deriveSpY.tau2 = int3/int1
      B.SpY = cbind(deriveSpY.mu, deriveSpY.tau2)
      
      I12.hat = t(B.SeY)%*%B.SpY/m
      
      myoff.diag = solve(I11.hat)%*%I12.hat%*%solve(I22.hat)/m
      myupper = cbind(fit.SeY$variance, myoff.diag)
      mylower = cbind(t(myoff.diag), fit.SpY$variance)
      myV = rbind(myupper, mylower)
      colnames(myV)=c("SeY", "SeY.S.2", "SpY", "SpY.S.2")
      rownames(myV)=c("SeY", "SeY.S.2", "SpY", "SpY.S.2")
      
      return(list(coefficients=estim, par.orig=estim.orig, vcov=myV))
    }
    
    
    score.Li.func = function(PiY,SeY,SpY,n,n1,n0,PiY.est,SeY.est,SpY.est){
      beta0 = PiY.est[1]; beta1 = SeY.est[1]; beta2 = SpY.est[1]
      mu0 = expit(beta0); mu1 = expit(beta1); mu2 = expit(beta2)
      phi0 = PiY.est[2]; phi1 = SeY.est[2]; phi2 = SpY.est[2]
      theta0 = phi0/(1-phi0); theta1 = phi1/(1-phi1); theta2 = phi2/(1-phi2)
      k01 = c(0:(PiY-1)); k02 = c(0:(n-PiY-1)); k03 = c(0:(n-1))
      k11 = c(0:(SeY-1)); k12 = c(0:(n1-SeY-1)); k13 = c(0:(n1-1))
      k21 = c(0:(SpY-1)); k22 = c(0:(n0-SpY-1)); k23 = c(0:(n0-1))
      # notice for outbound
      if(PiY == 0) k01 = NULL; if(PiY == n) k02 = NULL
      if(SeY == 0) k11 = NULL; if(SeY == n1) k12 = NULL
      if(SpY == 0) k21 = NULL; if(SpY == n0) k22 = NULL
      Dbeta0 = Dbeta1 = Dbeta2 = Dphi0 = Dphi1 = Dphi2 = 0
      Dphi0 = (sum(k01/(mu0+k01*theta0))+sum(k02/(1-mu0+k02*theta0))- sum(k03/(1+k03*theta0)))/((1-phi0)^2)
      Dphi1 = (sum(k11/(mu1+k11*theta1))+sum(k12/(1-mu1+k12*theta1))-sum(k13/(1+k13*theta1)))/((1-phi1)^2)
      Dphi2 = (sum(k21/(mu2+k21*theta2))+sum(k22/(1-mu2+k22*theta2))-sum(k23/(1+k23*theta2)))/((1-phi2)^2)
      Dbeta0 = sum(mu0*(1-mu0)/(mu0+k01*theta0))-sum(mu0*(1-mu0)/(1-mu0+k02*theta0))
      Dbeta1 = sum(mu1*(1-mu1)/(mu1+k11*theta1))-sum(mu1*(1-mu1)/(1-mu1+k12*theta1))
      Dbeta2 = sum(mu2*(1-mu2)/(mu2+k21*theta2))-sum(mu2*(1-mu2)/(1-mu2+k22*theta2))
      score.PiY = as.vector(c(Dbeta0,Dphi0))
      score.SeY = as.vector(c(Dbeta1,Dphi1))
      score.SpY = as.vector(c(Dbeta2,Dphi2))
      return(list(score.PiY=score.PiY, score.SeY=score.SeY, score.SpY=score.SpY))
    }
    
    
    score.Li.all.func= function(SeY,SpY,n1,n0,SeY.est,SpY.est){
      beta1 = SeY.est[1] ; beta2 = SpY.est[1]
      mu1 = expit(beta1) ; mu2 = expit(beta2)
      phi1 = SeY.est[2] ; phi2 =SpY.est[2]
      theta1 = phi1/(1-phi1) ; theta2 = phi2/(1-phi2)
      k11 = c(0:(SeY-1)); k12 = c(0:(n1-SeY-1)); k13 = c(0:(n1-1))
      k21 = c(0:(SpY-1)); k22 = c(0:(n0-SpY-1)); k23 = c(0:(n0-1))
      # notice for outbound
      if(SeY == 0) k11 = NULL ; if(SeY == n1) k12 = NULL
      if(SpY == 0) k21 = NULL ; if(SpY == n0) k22 = NULL
      Dbeta1 = Dbeta2 = Dphi1 = Dphi2 = 0
      Dphi1 = (sum(k11/(mu1+k11*theta1))+sum(k12/(1-mu1+k12*theta1))- sum(k13/(1+k13*theta1)))/((1-phi1)^2)
      Dphi2 = (sum(k21/(mu2+k21*theta2))+sum(k22/(1-mu2+k22*theta2))- sum(k23/(1+k23*theta2)))/((1-phi2)^2)
      Dbeta1 = sum(mu1*(1-mu1)/(mu1+k11*theta1))- sum(mu1*(1-mu1)/(1-mu1+k12*theta1))
      Dbeta2 = sum(mu2*(1-mu2)/(mu2+k21*theta2))-sum(mu2*(1-mu2)/(1-mu2+k22*theta2))
      score.SeY = as.vector(c(Dbeta1,Dphi1))
      score.SpY = as.vector(c(Dbeta2,Dphi2))
      return(list(score.SeY=score.SeY, score.SpY=score.SpY))
    }
    
    
    Fisher.inf = function(PiY.del,SeY.del,SpY.del,n.del,n1.del, n0.del, SeY,SpY,n1,n0,PiY.est,SeY.est, SpY.est){
      m = length(SeY) # the length of all study, including case control study & cohort study
      m2 = length(PiY.del) # the length of the cohort study only
      I00.array = I01.array = I02.array = array(NA,c(2,2,m2))
      I11.array = I12.array = I22.array = array(NA,c(2,2,m))
      for(i in 1:m2){
        temp = score.Li.func (PiY.del[i], SeY.del[i],SpY.del[i],n.del[i], n1.del[i],n0.del[i],PiY.est, SeY.est, SpY.est)
        score.PiY = temp$score.PiY
        score.SeY = temp$score.SeY
        score.SpY = temp$score.SpY
        I00.array[ , , i] = score.PiY%*%t(score.PiY)
        I01.array[ , , i] = score.PiY%*%t(score.SeY)
        I02.array[ , , i] = score.PiY%*%t(score.SpY)
      }
      for(i in 1:m){
        temp = score.Li.all.func (SeY[i],SpY[i],n1[i],n0[i],SeY.est,SpY.est)
        score.SeY = temp$score.SeY
        score.SpY = temp$score.SpY
        I11.array[ , , i] = score.SeY%*%t(score.SeY)
        I12.array[ , , i] = score.SeY%*%t(score.SpY)
        I22.array[ , , i] = score.SpY%*%t(score.SpY)
      }
      I00 = apply(I00.array, c(1,2), mean)
      I01 = apply(I01.array, c(1,2), mean)
      I02 = apply(I02.array, c(1,2), mean)
      I11 = apply(I11.array, c(1,2), mean)
      I12 = apply(I12.array, c(1,2), mean)
      I22 = apply(I22.array, c(1,2), mean)
      return(list(I00=I00, I01=I01, I02=I02, I11=I11, I12=I12, I22=I22))
    }

    tb.cl=function(n, PiY, SeY, n1, SpY, n2){
      m = length(SeY) # the length of all study, including case control study & cohort study
      m2 = sum(is.na(PiY)!=1) # the length of the cohort study only
      
      estim = estim.orig = mbse.orig = matrix(NA,nrow=6,ncol=1)
      mySandwich = matrix(NA,nrow=6,ncol=6)
      mydat1 = data.frame(n,PiY,n1,SeY,n2, SpY)
      names(mydat1) = c("n","PiY", "n1","SeY","n0","SpY")
      mydat2 = mydat1[is.na(mydat1$PiY)==FALSE, ] ## dataset used for estimating prevalence
      
      fit.PiY = betabin(cbind(PiY, n-PiY)~1, ~1, hessian=T, data=mydat2, link="logit")
      fit.SeY = betabin(cbind(SeY, n1-SeY)~1, ~1, hessian=T, data=mydat1, link="logit")
      fit.SpY = betabin(cbind(SpY, n0-SpY)~1, ~1, hessian=T, data=mydat1, link="logit")
      
      estim[c(1:2), 1] = c(as.numeric(fit.PiY@param))
      estim[c(3:4), 1] = c(as.numeric(fit.SeY@param))
      estim[c(5:6), 1] = c(as.numeric(fit.SpY@param))
      
      estim.orig[c(1: 2), 1] = c(expit(estim[1, 1]), estim[2, 1])
      estim.orig[c(3: 4), 1] = c(expit(estim[3, 1]), estim[4, 1])
      estim.orig[c(5: 6), 1] = c(expit(estim[5, 1]), estim[6, 1])
      
      rownames(estim)=c("PiY1", "PiY2", "SeY1", "SeY2", "SpY1", "SpY2")
      rownames(estim.orig)=c("PiY1", "PiY2", "SeY1", "SeY2", "SpY1", "SpY2")
      
      fit.PiY.var = matrix(c(fit.PiY@varparam[1,1], fit.PiY@varparam[1,2], fit.PiY@varparam[2,1],fit.PiY@varparam[2,2]), nrow=2,byrow=TRUE)
      fit.SeY.var = matrix(c(fit.SeY@varparam[1,1], fit.SeY@varparam[1,2], fit.SeY@varparam[2,1],fit.SeY@varparam[2,2]), nrow=2,byrow=TRUE)
      fit.SpY.var = matrix(c(fit.SpY@varparam[1,1], fit.SpY@varparam[1,2], fit.SpY@varparam[2,1],fit.SpY@varparam[2,2]), nrow=2,byrow=TRUE)
      
      I00.hat = solve(m2*fit.PiY.var)
      I11.hat = solve(m*fit.SeY.var)
      I22.hat = solve(m*fit.SpY.var)
      estim.PiY = estim[c(1:2), 1]
      estim.SeY = estim[c(3:4), 1]
      estim.SpY = estim[c(5:6), 1]
      
      temp = Fisher.inf(mydat2$PiY,mydat2$SeY,mydat2$SpY,mydat2$n,mydat2$n1,mydat2$n0, mydat1$SeY,mydat1$SpY,mydat1$n1,mydat1$n0,estim.PiY, estim.SeY, estim.SpY)
      
      I01.hat = temp$I01 ; I02.hat = temp$I02 ; I12.hat = temp$I12
      myoff.diag01 = solve(I00.hat)%*%I01.hat%*%solve(I11.hat)/m
      myoff.diag02 = solve(I00.hat)%*%I02.hat%*%solve(I22.hat)/m
      myoff.diag12 = solve(I11.hat)%*%I12.hat%*%solve(I22.hat)/m
      myupper = cbind((fit.PiY.var), sqrt(m/m2)*(myoff.diag01), sqrt(m/m2)*(myoff.diag02))
      mymiddle = cbind(sqrt(m/m2)*t(myoff.diag01), (fit.SeY.var), (myoff.diag12))
      mylower = cbind(sqrt(m/m2)*t(myoff.diag02), t(myoff.diag12), (fit.SpY.var))
      myV = rbind(myupper, mymiddle, mylower)
      colnames(myV)=c("PiY1", "PiY2", "SeY1", "SeY2", "SpY1", "SpY2")
      rownames(myV)=c("PiY1", "PiY2", "SeY1", "SeY2", "SpY1", "SpY2")
      
      mu0.mbse = (estim.orig[1, 1]*(1-estim.orig[1, 1]))*sqrt(diag(myV)[1])
      mu1.mbse = (estim.orig[3, 1]*(1-estim.orig[3, 1]))*sqrt(diag(myV)[3])
      mu2.mbse = (estim.orig[5, 1]*(1-estim.orig[5, 1]))*sqrt(diag(myV)[5])
      s1.2.mbse = sqrt(diag(myV)[2])
      s2.2.mbse = sqrt(diag(myV)[4])
      s3.2.mbse = sqrt(diag(myV)[6])
      mbse.orig = c(mu0.mbse, s1.2.mbse, mu1.mbse, s2.2.mbse, mu2.mbse, s3.2.mbse)
      return(list(coefficients=estim, estim.orig=estim.orig, vcov=myV))
    }
    
    
    
    tn.cl=function(n, PiY2, SeY1, SeY2, SpY1, SpY2, Nd, Nn){
      
      estim = estim.orig = mbse.orig = matrix(NA, nrow=6, ncol=1)
      
      m1=length(SeY1) 
      m2=length(SeY2)
      m=m1+m2
      
      PiY = c(rep(NA, length=m1),PiY2)
      SeY = c(SeY1, SeY2)
      SpY = c(SpY1, SpY2)
      
      id=c(1:m)
      mydat = data.frame(id, rep(n, length=m), PiY, Nd ,SeY,Nn, SpY)
      names(mydat) = c("id", "Nt","PiY", "Nd","SeY","Nn","SpY") ## codebook for mydat: study id, number of total subjects, number of disease ppl among total subjects, number of diease ppl, number of being postive test among disease ppl, number of disease-free ppl, number of being negative test among disease-free ppl
      
      mydat2=na.omit(mydat)
      
      fit.PiY = glmmML(cbind(PiY, Nt-PiY)~1, data = mydat2, cluster = id, method= "ghq", n.points=8, 	control=list(epsilon=1e-08 , maxit=50000,trace=FALSE))
      fit.SeY = glmmML(cbind(SeY, Nd-SeY)~1, data = mydat, cluster = id, method= "ghq", n.points=8, 	control=list(epsilon=1e-08 , maxit=50000,trace=FALSE))
      fit.SpY = glmmML(cbind(SpY, Nn-SpY)~1, data = mydat, cluster = id, method= "ghq", n.points=8, control=list(epsilon=1e-08 , maxit=50000,trace=FALSE))
      
      estim[c(1: 2),1] = c(as.numeric(fit.PiY$coefficients), (fit.PiY$sigma)^2)
      estim[c(3: 4),1] = c(as.numeric(fit.SeY$coefficients), (fit.SeY$sigma)^2)
      estim[c(5: 6),1] = c(as.numeric(fit.SpY$coefficients), (fit.SpY$sigma)^2)
      
      estim.orig[c(1: 2),1] = c(plogis(estim[1,1]), estim[2,1])
      estim.orig[c(3: 4),1] = c(plogis(estim[3,1]), estim[4,1])
      estim.orig[c(5: 6),1] = c(plogis(estim[5,1]), estim[6,1])
      rownames(estim)=c("PiY1", "PiY2", "SeY1", "SeY2", "SpY1", "SpY2")
      rownames(estim.orig)=c("PiY1", "PiY2", "SeY1", "SeY2", "SpY1", "SpY2")
      
      I00.hat = solve(m2*fit.PiY$variance)
      I11.hat = solve(m*fit.SeY$variance)
      I22.hat = solve(m*fit.SpY$variance)
      
      int1 = int2  = int3 = rep(NA, length=m2)
      est0 = estim[c(1: 2),1] ##plug in the point estimate based on disease prevalence
      for(i in 1: m2){
        int1[i] = integrate(integrand1, lower=0, upper=1, n=mydat2$Nt[i], mypar2=est0, y=mydat2$PiY[i])$value
        int2[i] = integrate(integrand2, lower=0, upper=1, n=mydat2$Nt[i], mypar2=est0, y=mydat2$PiY[i])$value
        int3[i] = integrate(integrand3, lower=0, upper=1, n=mydat2$Nt[i], mypar2=est0, y=mydat2$PiY[i])$value
      }
      ##
      derivePiY.mu = int2/int1
      derivePiY.tau2 = int3/int1
      B.PiY = cbind(derivePiY.mu, derivePiY.tau2)
      
      int1 = int2  = int3 = rep(NA, length=m)
      est1 = estim[c(3: 4),1] ##plug in the point estimate based on sensitivity
      for(i in 1: m){
        int1[i] = integrate(integrand1, lower=0, upper=1, n=mydat$Nd[i], mypar2=est1, y=mydat$SeY[i])$value
        int2[i] = integrate(integrand2, lower=0, upper=1, n=mydat$Nd[i], mypar2=est1, y=mydat$SeY[i])$value
        int3[i] = integrate(integrand3, lower=0, upper=1, n=mydat$Nd[i], mypar2=est1, y=mydat$SeY[i])$value
      }
      deriveSeY.mu = int2/int1
      deriveSeY.tau2 = int3/int1
      B.SeY = cbind(deriveSeY.mu, deriveSeY.tau2)
      B.SeY.del = B.SeY[-c(1:m1), ] ## delete the first m1 studies
      
      int1 = int2  = int3 = rep(NA, length=m)
      est2 = estim[c(5: 6),1] ##plug in the point estimate based on specificity
      for(i in 1: m){
        int1[i] = integrate(integrand1, lower=0, upper=1, n=mydat$Nn[i], mypar2=est2, y=mydat$SpY[i])$value
        int2[i] = integrate(integrand2, lower=0, upper=1, n=mydat$Nn[i], mypar2=est2, y=mydat$SpY[i])$value
        int3[i] = integrate(integrand3, lower=0, upper=1, n=mydat$Nn[i], mypar2=est2, y=mydat$SpY[i])$value
      }
      deriveSpY.mu = int2/int1
      deriveSpY.tau2 = int3/int1
      B.SpY = cbind(deriveSpY.mu, deriveSpY.tau2)
      B.SpY.del = B.SpY[-c(1:m1), ] ## delete the first m1 studies 
      
      I01.hat = t(B.PiY)%*%B.SeY.del/m2 
      I02.hat = t(B.PiY)%*%B.SpY.del/m2 
      I12.hat = t(B.SeY)%*%B.SpY/m
      
      myoff.diag01 = solve(I00.hat)%*%I01.hat%*%solve(I11.hat)/m
      myoff.diag02 = solve(I00.hat)%*%I02.hat%*%solve(I22.hat)/m
      myoff.diag12 = solve(I11.hat)%*%I12.hat%*%solve(I22.hat)/m
      myupper = cbind(fit.PiY$variance, sqrt(m/m2)*(myoff.diag01), sqrt(m/m2)*(myoff.diag02))
      mymiddle = cbind(sqrt(m/m2)*t(myoff.diag01), fit.SeY$variance, (myoff.diag12))
      mylower = cbind(sqrt(m/m2)*t(myoff.diag02), t(myoff.diag12), fit.SpY$variance)
      myV = rbind(myupper, mymiddle, mylower)
      colnames(myV)=c("PiY1", "PiY2", "SeY1", "SeY2", "SpY1", "SpY2")
      rownames(myV)=c("PiY1", "PiY2", "SeY1", "SeY2", "SpY1", "SpY2")
      
      return(list(coefficients=estim, estim.orig=estim.orig, vcov=myV))
    }
    
    if (missing(data)){data=NULL}
    if (is.null(data)){
      stop("The dataset must be specified.")
    }   
    
    if (missing(k)){k=NULL}
    if (is.null(k)){
      stop("The number of outcomes must be specified.")
    }  
    
    if (missing(type)){type=NULL}
    if (is.null(type)){
      stop('The type of outcome must be specified. Please input "continuous" or "binary". ')
    }      
    
    if (missing(method)){method=NULL}
    if (is.null(method)){
      stop("A method must be specified. ")
    }
    if (missing(rhow)){rhow=NULL}
    if ((is.null(rhow))&(method=="nn.reml")){
      stop('The within-study correlations are required for the input method "nn.reml".')
    }
    if (k>2){
      stop('The method for MMA with more than 2 outcomes are currently under development. ')
    }
    
    if (type=="continuous"){
      y1=data$y1
      s1=data$s1
      y2=data$y2
      s2=data$s2
      
      if(method=="nn.cl"){
        fit=nn.cl(y1, s1, y2, s2)}
      
      if(method=="nn.reml"){
        fit=nn.reml(y1, s1, y2, s2, rhow)}
      
      if(method=="nn.mom"){
        fit=nn.mom(y1, s1, y2, s2)}
      
      if(method=="nn.rs"){
        fit=nn.rs(y1, s1, y2, s2)}
      
    
    if ((method%in%c("nn.reml", "nn.cl", "nn.mom", "nn.rs"))!=1){
      stop('The input method is not available for continuous outcomes. Please choose from
           "nn.reml", "nn.cl", "nn.mom", or "nn.rs". ')
    }
    }
if(type=="binary"){
  if (method=="bb.cl"){
    y1=data$y1
    n1=data$n1
    y2=data$y2
    n2=data$n2
    fit=bb.cl(y1,n1,y2,n2)
  }
  if (method=="bn.cl"){
    y1=data$y1
    n1=data$n1
    y2=data$y2
    n2=data$n2
    fit=bn.cl(y1,n1,y2,n2)
  }
  if (method=="tb.cl"){
    n=data$n
    PiY=data$PiY
    SeY=data$SeY
    n1=data$n1
    SpY=data$SpY
    n2=data$n2       
    fit=tb.cl(n, PiY, SeY, n1, SpY, n2) 
  }
  if (method=="tn.cl"){
    n=data$n
    PiY2=data$PiY2
    SeY1=data$SeY1
    SeY2=data$SeY2
    SpY1=data$SpY1
    SpY2=data$SpY2
    Nd=data$Nd
    Nn=data$Nn
    fit=tn.cl(n, PiY2, SeY1, SeY2, SpY1, SpY2, Nd, Nn)
  }
  if ((method%in%c("bb.cl", "bn.cl", "tb.cl", "tn.cl"))!=1){
    stop('The input method is not available for binary outcomes. Please choose from
         "bb.cl", "bn.cl", "tb.cl", or "tn.cl". ')
  }
  }

res=list(type=type, k=k, method=method, coefficients=fit$coefficients, vcov=fit$vcov)
if(method=="bb.cl"){
res=list(type=type, k=k, method=method, logOR=fit$logOR, OR_CI=fit$OR_CI, se=fit$se,hessian.log=fit$hessian.log,conv=fit$conv,mle=fit$mle)
}

class(res) = c("mmeta")
return(res)
}
