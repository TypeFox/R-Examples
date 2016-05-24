algEMbssn<-function(ti,alpha,beta,delta,loglik=F,accuracy = 1e-8,show.envelope="FALSE",iter.max=1000)
{
 start.time <- Sys.time() #alpha=alpha0; beta=beta0; delta=delta0

 betamax<-function(beta,ti,alpha,delta,hi,hi2)
 {
  n  <- length(ti)
  at <- (1/alpha)*(sqrt(ti/beta)-sqrt(beta/ti))
  Q  <- n*log(alpha)+(n/2)*log(beta)-sum(log(ti+beta))+(n/2)*log(1-delta^2)+(1/(2*(1-delta^2)))*sum(at^2)-(delta/(1-delta^2))*sum(hi*at)+(delta^2/(2 * (1-delta^2)))*sum(hi2)
  return(Q)
 }

 n         <- length(ti)
 tetha     <- c(alpha,beta,delta)
 beta0     <- beta

 criterion <- sqrt(tetha %*% tetha)
 count     <- 0

 while(criterion>accuracy && (count <= iter.max))
 {
  count     <- count+1

  # E-step: computing hi, hi2
   at       <- (1/alpha)*(sqrt(ti/beta)-sqrt(beta/ti))
   hi       <- (delta*at)+(dnorm((delta*at)/sqrt(1-delta^2))/pnorm((delta*at)/sqrt(1-delta^2)))*sqrt(1-delta^2)
   hi2      <- (delta*at)^2+(1-delta^2)+(dnorm((delta*at)/sqrt(1-delta^2))/pnorm((delta*at)/sqrt(1-delta^2)))*sqrt(1-delta^2)*(delta*at)

  # M-step: updating alpha, beta, and delta
   psi      <- (sqrt(ti/beta)-sqrt(beta/ti))
   th       <- sum(hi2)/n
   alpha    <- sqrt((1/n)*sum(psi^2)+(1-th)*(sum(hi*psi)/(n*th))^2)
   delta    <- (1/alpha)*(sum(hi*psi)/(n*th))

  beta      <- nlminb(start=beta0,betamax,gradient=NULL,hessian=NULL,ti,alpha,delta,hi,hi2,scale=1,control=list(),lower=1e-3,upper=Inf)$par

  par       <- tetha
  tetha     <- c(alpha,beta,delta)
  criterion <- sqrt((tetha-par)%*%(tetha-par))

  lambda    <- delta/sqrt(1-delta^2)

 }

 if(loglik == T)
 {
  lk        <- sum(log(dbssn(ti,alpha,beta,lambda)))
  result    <- list(alpha=alpha, beta=beta, lambda=lambda, delta=delta,loglik=lk,convergence=criterion < accuracy,iter=count,n=length(ti))
 } else {
  result    <- list(alpha=alpha, beta=beta, lambda=lambda, delta=delta,convergence=criterion < accuracy,iter=count,n= length(ti))
 }
 EP              <- sqrt(diag(solve(Infmatrix(ti,alpha,beta,lambda))))
 parameters <- rbind(alpha,beta,lambda)
 table      <- data.frame(parameters,EP,parameters/EP,2*pnorm(abs(parameters/EP),lower.tail = F))

 p          <- length(tetha)
 loglik     <- sum(log(dbssn(ti, alpha, beta, lambda)))
 AIC        <- (-loglik/n)+(p/n)
 BIC        <- (-loglik/n)+((p/2)*(log(n)/n))
 HQC        <- (-loglik/n)+(p*log(log(n))/n)

 rownames(table) <- c("alpha","beta","lambda")
 colnames(table) <- c("Estimate","Std. Error","z value","Pr(>|z|)")

 end.time   <- Sys.time()
 time.taken <- end.time - start.time

 if(show.envelope == TRUE)
 {
  d2s       <- ((1/alpha)*(sqrt(ti/beta) - sqrt(beta/ti)))^2
  d2s       <- sort(d2s)
  xq2       <- qchisq(ppoints(n), 1)

  Xsim      <- matrix(0,100,n)
  for(i in 1:100)
  {
   Xsim[i,] <-rchisq(n, 1)
  }

  Xsim2     <- apply(Xsim,1,sort)
  d21       <- matrix(0,n,1)
  d22       <- matrix(0,n,1)

  for(i in 1:n)
  {
   d21[i]   <- quantile(Xsim2[i,],0.05)
   d22[i]   <- quantile(Xsim2[i,],0.95)
  }

  d2med     <- apply(Xsim2,1,mean)
  fy        <- range(d2s,d21,d22)
  plot(xq2,d2s,xlab = expression(paste("Theorical quantile")),
        ylab="Empirical quantile and generated envelope",pch=20,ylim=fy)
  par(new=T); plot(xq2,d21,type="l",ylim=fy,xlab="",ylab="")
  par(new=T); plot(xq2,d2med,type="l",ylim=fy,xlab="",ylab="")
  par(new=T); plot(xq2,d22,type="l",ylim=fy,xlab="",ylab="")
 }

 result     <- list(alpha=alpha,betat=beta,lambda=lambda,iter=count,criterion = criterion, n=length(t),EP=EP,table=table,loglik=loglik, AIC=AIC, BIC=BIC, HQC=HQC, time = time.taken, convergence = criterion<accuracy, iter.max=iter.max)

 obj.out    <- list(result = result)
 class(obj.out)  =  "bssn"
 return(obj.out)
}

#algEMbssn(x,alpha,beta1,delta,loglik=T,show.envelope="FALSE")












