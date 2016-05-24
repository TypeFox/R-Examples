lnorm.ml <-
function(dd){
#    ML estimates for lognormal sample with non-detects 
#   REVISE arguments to optim() as required by
#    R 2.6.0 Change cont TO control=cont on lines 43 and 48 
# USAGE: mlndln( dd )
# ARGUMENT: matrix dd with x[i] in column 1 and det[i] in col 2
#     x[i] is positive lognormal data 
#     det[i]=0 for non-detect ; 1 for detect
# NOTE:
#     y= log(x) is normal with mean mu and standard deviation sigma
#     E(X) = exp( mu + 0.5*sig2)= exp(logE) where sig2 = sigma^2
#     m is number of detects and Conver is convergence check
# VALUE: ML estimates of following in 2 by 6 matrix format:
#     mu    sigma      logE     sig2     -2Log(L)     Conver
#   se.mu  se.sigma  se.logE  se.sig2   cov(mu,sig)   m
# REFERENCE: Cohen, A.C (1991) Truncated and Censored Samples
#            Marcel Decker, New York        
#  REQUIREs:  ndln() ndln2() loglikelihood functions for optim()        
#             see R help file for details on optim() and dlnorm() 
#   
ndln <- function(p=est,xd){
# -log likelihood functions for optim() for mu and sig
 mu<-p[1]; sig<-p[2]; x<-xd[,1]
xx<-ifelse(xd[,2]==1,dlnorm(x,mu,sig,log=TRUE) , plnorm(x,mu,sig,log.p=TRUE))
 -sum(xx)
}
ndln2 <- function(p=est,xd){
# -log likelihood functions for optim() for E(x) and sig^2
 mu<-p[1] - 0.5*p[2]; sig<-sqrt(p[2]);x<-xd[,1]
xx<-ifelse(xd[,2]==1,dlnorm(x,mu,sig,log=TRUE) , plnorm(x,mu,sig,log.p=TRUE))
  -sum(xx)
}                                                        

   m <- sum(dd[,2])   #  number of non-detects
   n<- length(dd[,1])   #  samp[le size   
#  initial estimate of mu and sig (sigma)   
      yt <- ifelse(dd[,2]==0,dd[,1]/2,dd[,1] )
      est <- c( mean(log(yt)), sd(log(yt)) )
# ML estimates  mu and sig using derivative free method
est <- optim(est,ndln, method = c("Nelder-Mead"),xd=dd )$par

cont <- list(parscale=abs(est),trace=0)
opt1 <- optim(est,ndln ,NULL,  method ="L-BFGS-B",lower=c(-Inf,0.0),
          upper=c(Inf,Inf),control=cont, hessian=TRUE,xd=dd )

conv1 <- opt1$conv          # convergenc check from optim()
mle <- opt1$par             # ML estimate of mu and sig
vcm <- solve(opt1$hessian)
semle <- sqrt(diag( vcm ))  # standard Errors of mu and sig
cov <- vcm[1,2] #  covariace(mu,sig) needed for Tolerance bound
#
# ML estimate of logE(x)  and sig2 (sigmma^2)
#
 est[1] <- mle[1] + 0.5*mle[2]^2
 est[2] <- mle[2]^2
cont <- list(parscale=abs(est))
opt2 <- optim(est,ndln2 ,NULL,  method ="L-BFGS-B",lower=c(-Inf,0.0),
          upper=c(Inf,Inf),control=cont, hessian=TRUE,xd=dd )
#  next line adds ML estimate of logE  sig2 -2Log(L) and Conver
#  If Conver is not equal to 0 CHECK RESULTs--- see optim() help

mle <- c(mle,opt2$par,2*opt2$value,opt2$conv+conv1 )
se <- c(semle,sqrt(diag(solve(opt2$hessian))),vcm[1,2],m )


out<-list("mu"=mle[1],"sigma"=mle[2],"logEX"=mle[3],"Sigmasq"=mle[4],
"se.mu"=se[1],"se.sigma"=se[2],"se.logEX"=se[3],"se.Sigmasq"=se[4],
"cov.musig"=cov,"m"=m,"n"=n,"m2log(L)"=mle[5],"convergence"=mle[6] )

out
}

