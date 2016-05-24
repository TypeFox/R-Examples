##BD_calc.R: generating function for BDI model (so, "BD_calc" would be better as
##  "BDI_calc").

# Generating function: Joint birth counts (r) and ending state (s)
#
# lambda - birth rate
# mu - death rate
# nu - immigration rate (=lambda for TKF91)
# X0 - starting state

add.generator = function(r,s,t,lambda,mu,nu,X0) {
	x1 = (lambda + mu - sqrt((lambda+mu)^2 - 4*lambda*mu*r )) / (2*lambda*r)
	x2 = (lambda + mu + sqrt((lambda+mu)^2 - 4*lambda*mu*r )) / (2*lambda*r)
	exp.term = exp(-lambda*r*(x2-x1)*t)
	ratio.term = (s - x1) / (s - x2)
	H0 = (x1 - x2*ratio.term*exp.term) / (1 - ratio.term*exp.term)
	H1 = ((x2 - x1) / ( (ratio.term*exp.term - 1)*(s-x2) ))^(nu/lambda) * exp(-nu*(1-r*x1)*t)
 	(H0)^(X0) * H1  
}

#
# Laplace transform of pseudo-CDF: Joint particle-time-average (r) and ending state (s)
#
# lambda - birth rate
# mu - death rate
# nu - immigration rate (= lambda for TFK91)
# X0 - starting state


hold.generator <- function(w,s,t,lambda,mu,nu,X0){timeave.laplace(r=w,s,t,lambda,mu,nu,X0)};

### same as timeave.laplace
## hold.generator2 <-function(w,s,t,lambda,mu,nu,X0){
##   addremhold.generator(u=1,v=1,w=w,s=s,t=t,lambda=lambda,mu=mu,nu=nu,X0=X0)
## }



#
# Generating function: Joint death counts (r) and ending state (s)
#
# lambda - birth rate
# mu - death rate
# nu - immigration rate (=lambda for TKF91)
# X0 - starting state

rem.generator = function(r,s,t,lambda,mu,nu,X0) {
	x1 = (lambda + mu - sqrt((lambda+mu)^2 - 4*lambda*mu*r )) / (2*lambda)
	x2 = (lambda + mu + sqrt((lambda+mu)^2 - 4*lambda*mu*r )) / (2*lambda)
	exp.term = exp(-lambda*(x2-x1)*t)
	ratio.term = (s - x1) / (s - x2)
	H0 = (x1 - x2*ratio.term*exp.term) / (1 - ratio.term*exp.term)
	H1 = ((x2 - x1) / ( (ratio.term*exp.term - 1)*(s-x2) ))^(nu/lambda) * exp(-nu*(1-x1)*t)
 	(H0)^(X0) * H1  
}




###old version. backup.
#Generating function: product add,rem
#addrem.generator <- function( u, v, s, t, X0, lambda, mu, nu){
#  a1 <- (lambda + mu - sqrt((lambda+mu)^2 - 4*lambda*mu*u*v) ) / (2*lambda*u);
#  a2 <- (lambda + mu + sqrt((lambda+mu)^2 - 4*lambda*mu*u*v) ) / (2*lambda*u);
#  H(a1=a1,a2=a2, exparg1=a1, exparg2=a2, r=u,s,t,X0,lambda,nu);  #Like "Hplus" but modified alphai
#}



###############One generator that yields all of them.
###Note: some of the generating functions are left as they are rather than being
### recoded to come from this one (since they're the same things).
#u=N+, v=N-, w=Rt=holdtime
addremhold.generator <- function( u, v, w, s, t, X0, lambda, mu, nu){
  a1 <- (lambda + mu + w - sqrt((lambda+mu+w)^2 - 4*lambda*mu*u*v) ) / (2*lambda*u);
  a2 <- (lambda + mu + w + sqrt((lambda+mu+w)^2 - 4*lambda*mu*u*v) ) / (2*lambda*u);
  H(a1=a1,a2=a2, exparg1=a1, exparg2=a2, r=u,s,t,X0,lambda,nu);  #Like "Hplus" but modified alphai
}


remhold.generator <- function( v, w, s, t, X0, lambda, mu, nu){
  addremhold.generator( u=1, v, w, s, t, X0, lambda, mu, nu)
}

addhold.generator <- function( u, w, s, t, X0, lambda, mu, nu){
  addremhold.generator( u=u, v=1, w, s, t, X0, lambda, mu, nu)
}

addrem.generator <- function( u, v, s, t, X0, lambda, mu, nu){
  addremhold.generator(u=u,v=v,w=0,s,t,X0,lambda,mu,nu)
}

##mean number events, conditional only on starting point, not ending.
add.uncond.mean.one <-function(t,X0,lambda,mu,nu,
                               delta=0.001,r=4){
  a.gen <- function(u){addremhold.generator(u=u,v=1,w=0,s=1, t=t,X0=X0,lambda=lambda,mu=mu,nu=nu)}
  return(genD(a.gen,x=1,method.args=list(d=delta,eps=delta,r=r))$D[1])
}

rem.uncond.mean.one <-function(t,X0,lambda,mu,nu,
                               delta=0.001,r=4){
  r.gen <- function(v){addremhold.generator(u=1,v=v,w=0,s=1, t=t,X0=X0,lambda=lambda,mu=mu,nu=nu)}
  return(genD(r.gen,x=1,method.args=list(d=delta,eps=delta,r=r))$D[1])
}


hold.uncond.mean.one <-function(t,X0,lambda,mu,nu,
                               delta=0.001,r=4){
  h.gen <- function(w){addremhold.generator(u=1,v=1,w=w,s=1, t=t,X0=X0,lambda=lambda,mu=mu,nu=nu)}
  return(-genD(h.gen,x=0,method.args=list(d=delta,eps=delta,r=r))$D[1])
}









####################################################
######## All the below *.cond.mean*.one are now in BD_calc_helpers.R via the Curry function
################################################################

## addrem.cond.mean.one <- function(t,lambda,mu,nu=0,X0=1,Xt,
##                                  trans.prob=NULL, joint.mean=NULL,
##                                  delta=0.001,n=1024, r=4, prec.tol=1e-12, prec.fail.stop=TRUE){
##   if (is.null(joint.mean))
##     joint.mean <- addrem.joint.mean.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,
##                                         delta=delta,n=n, r=r)
##   if(is.null(trans.prob))
##     trans.prob <- process.prob.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,n=n)
##   prec.error.handler(joint.mean,trans.prob,prec.tol,prec.fail.stop,
##                      fnid="remhold.cond.mean.one");
## }


## addhold.cond.mean.one <- function(t,lambda,mu,nu=0,X0=1,Xt,
##                                   trans.prob=NULL,
##                                   delta=0.001,n=1024,r=4, prec.tol=1e-12, prec.fail.stop=TRUE){
##   joint.mean <- addhold.joint.mean.one(t,lambda,mu,nu=nu,X0=X0,Xt=Xt,delta=delta,
##                                        n=n, r=r);
##   if(is.null(trans.prob))
##     trans.prob = process.prob.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,n=n);
##   prec.error.handler(joint.mean,trans.prob,prec.tol,prec.fail.stop,
##                      fnid="addhold.cond.mean.one");
## }

## remhold.cond.mean.one <- function(t,lambda,mu,nu=0,X0=1,Xt,
##                                   trans.prob=NULL,
##                                   delta=0.001,n=1024,  r=4, prec.tol=1e-12, prec.fail.stop=TRUE){
##   joint.mean <- remhold.joint.mean.one(t,lambda,mu,nu=nu,X0=X0,Xt=Xt,
##                                        delta=delta,n=n,r=r);
##   if(is.null(trans.prob))
##     trans.prob = process.prob.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,n=n);
##   prec.error.handler(joint.mean,trans.prob,prec.tol,prec.fail.stop,
##                      fnid="remhold.cond.mean.one");
## ##   if (joint.mean < 10^-13) {
## ##     stop("remhold.cond.mean.one: joint mean too small.");
## ##     #return(NA)
## ##   }
## ##   else
## ##     return(joint.mean/trans.prob)
## }


## hold.cond.meanSq.one <- function(t,lambda,mu,nu=0,X0=1,Xt,
##                                  trans.prob=NULL,
##                                  n=1024,delta=1e-4,r=4,
##                                  prec.tol=1e-12, prec.fail.stop=TRUE){
##     joint.meanSq <- hold.joint.meanSq.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,
##                                       Xt=Xt,
##                                       n=n, delta=delta, r=r);
##   if(is.null(trans.prob))
##     trans.prob = process.prob.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,n=n)
##   prec.error.handler(joint.meanSq,trans.prob,prec.tol,prec.fail.stop,
##                      fnid="hold.cond.meanSq.one");
##   ##   if (joint.mean<10^-13) {
##   ##     stop("hold.cond.meanSq.one: joint mean too small.");
## ##                                         #return(NA)
## ##   }
## ##   else return(joint.mean/trans.prob);
## }

## rem.cond.meanSq.one <- function(t,lambda,mu,nu=0,X0=1,Xt,
##                                 trans.prob=NULL, joint.mean=NULL,
##                                 delta=0.001,n=1024, r=4, prec.tol=1e-12, prec.fail.stop=TRUE){
##   joint.meanSq = rem.joint.meanSq.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,
##     joint.mean=joint.mean,
##     delta=delta,n=n, r=r)
##   if (is.null(trans.prob))
##     trans.prob = process.prob.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,n=n)
##   prec.error.handler(joint.meanSq,trans.prob,prec.tol,prec.fail.stop,
##                      fnid="rem.cond.meanSq.one");
## }

## add.cond.meanSq.one <- function(t,lambda,mu,nu=0,X0=1,Xt,
##                                 trans.prob=NULL, joint.mean=NULL,
##                                 delta=0.001,n=1024, r=4,prec.tol=1e-12, prec.fail.stop=TRUE){
##   joint.meanSq = add.joint.meanSq.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,
##     joint.mean=joint.mean,
##     delta=delta,n=n,r=r)
##   if(is.null(trans.prob))
##     trans.prob = process.prob.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,n=n)
##   prec.error.handler(joint.meanSq,trans.prob,prec.tol,prec.fail.stop,
##                      fnid="rem.cond.meanSq.one");
## }


## add.cond.mean.one <- function(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024,r=4,
##                               trans.prob=NULL, joint.mean=NULL,
##                               prec.tol=1e-12, prec.fail.stop=TRUE){
##   if (is.null(joint.mean))
##     joint.mean <- add.joint.mean.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,
##                                      delta=delta,n=n,r=r);
##   if (is.null(trans.prob))
##     trans.prob = process.prob.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,n=n)
##   prec.error.handler(joint.mean,trans.prob,prec.tol,prec.fail.stop,
##                        fnid="add.cond.mean.one");
## }

## rem.cond.mean.one <- function(t,lambda,mu,nu=0,X0=1,Xt,
##                               trans.prob=NULL,
##                               delta=0.001,n=1024,r=4,
##                               prec.tol=1e-12, prec.fail.stop=TRUE){
##   joint.mean = rem.joint.mean.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,
##     delta=delta,n=n,r=r)
##   if(is.null(trans.prob))
##     trans.prob = process.prob.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,n=n)
##   prec.error.handler(joint.mean,trans.prob,prec.tol,prec.fail.stop,
##                      fnid="timeave.cond.mean.one");
## }


## timeave.cond.mean.one <- function(t,lambda,mu,nu=0,X0=1,Xt,
##                                   trans.prob=NULL,
##                                   delta=0.001,n=1024,r=4,
##                                   prec.tol=1e-12, prec.fail.stop=TRUE){
##   joint.mean = timeave.joint.mean.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,
##     delta=delta,n=n, r=r)
##   if(is.null(trans.prob))
##     trans.prob = process.prob.one(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,Xt=Xt,n=n)
##   prec.error.handler(joint.mean,trans.prob,prec.tol,prec.fail.stop,
##                      fnid="timeave.cond.mean.one");
## }

## hold.cond.mean.one <- timeave.cond.mean.one
























#
# Generating function: Ending state (s)
#
# lambda - birth rate
# mu - death rate
# nu - immigration rate (=lambda for TKF91)
# X0 - starting state
# careful! if lambda is too close to 0, things can go awry;
# if lambda==0, the correct result is programmed by hand.
process.generator <- function (s, time, lambda, mu, nu, X0) 
{
  if (lambda == mu)
    {return(process.generator.ratesEqual(s=s,time=time,lambda=lambda,mu=mu,nu=nu,X0=X0));}
  s <- matrix(s); #column
  time <- t(matrix(time));   #row
  timeMat <- matrix(rep(time,length(s)), nrow=length(s), byrow=TRUE)
  if (lambda < 1e-13){ ## ie "lambda==0". For tiny lambda, the inner term base of the exponent rounds to 1 and ruins the calculation.
    exponent <- (s%*%(1-exp(-mu*time)) +
                 exp(-mu*timeMat) - 1 )*nu/mu;
    factor1 <- exp(exponent);
  }
  else { #if lambda too small (say < 1e-15 ish) this will incorrectly evaluate to 1
    factor1 <- ((lambda - mu)/(lambda * (1 - s) %*% exp((lambda - mu) * time) + 
                               lambda * c(s) - mu))^(nu/lambda);
  }
  factor2 <- ((mu * (1 - s) %*% exp((lambda - mu) * time) +
               lambda * c(s) - mu)/(lambda * (1 - s) %*% exp((lambda - mu) * time) +
                                    lambda * c(s) - mu))^X0;
  return(factor2 * factor1);
}


### this works but slightly more functionality above.
## process.generator = function(s,t,lambda,mu,nu,X0) {
##  ((lambda-mu)/
##   (lambda*(1-s)*exp((lambda-mu)*t)+lambda*s-mu))^(nu/lambda)*
##     ((mu*(1-s)*exp((lambda-mu)*t)+lambda*s-mu)/
##      (lambda*(1-s)*exp((lambda-mu)*t)+lambda*s-mu))^X0
## }




#### vector functionality for t.  
process.prob.many <- function (t, lambda, mu, nu = 0, X0 = 1, n = 1024) 
{
  return(apply(as.matrix(t),1,
               function(z){
                 power.coef.many(process.generator, n = n, t = z, 
                                 lambda = lambda, mu = mu, nu = nu, X0 = X0);
               }));
#  return(power.coef.many(process.generator, n = n, t = t, lambda = lambda, 
#                         mu = mu, nu = nu, X0 = X0))  
}


## process.prob.many = function(t,lambda,mu,nu=0,X0=1,n=1024) {
##   return(power.coef.many(process.generator,
##                          n=n,
##                          t=t,
##                          lambda=lambda,
##                          mu=mu,
##                          nu=nu,
##                          X0=X0)
##          )
## }



add.cond.mean.many = function(t,lambda,mu,nu=0,X0=1,delta=0.001,n=1024, prec.tol=1e-12, prec.fail.stop=TRUE){
  joint.mean = add.joint.mean.many(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,delta=delta,n=n)
  trans.prob = process.prob.many(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,n=n)
  ##modify for "many" situation
  ##prec.error.handler(joint.mean,trans.prob,prec.tol,prec.fail.stop,
  ##                   fnid="rem.cond.mean.many");
  truncated.mean = joint.mean*as.numeric(abs(joint.mean)>10^-10)
  return(truncated.mean/trans.prob)
}

## rem.joint.mean.one = function(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024){
##   return(power.coef.one(num.deriv,
##                         ftn=rem.generator,
##                         var=1, ## differentiate rem.generator at r=1
##                         delta=delta,
##                         n=n,
##                         k=Xt,
##                         t=t,
##                         lambda=lambda,
##                         mu=mu,
##                         nu=nu,
##                         X0=X0)
##          )
## }


rem.cond.mean.many = function(t,lambda,mu,nu=0,X0=1,delta=0.001,n=1024, prec.tol=1e-12, prec.fail.stop=TRUE){
  joint.mean = rem.joint.mean.many(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,delta=delta,n=n)
  trans.prob = process.prob.many(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,n=n)
  ##modify for 'many' situation?
  ##prec.error.handler(joint.mean,trans.prob,prec.tol,prec.fail.stop,
  ##                   fnid="rem.cond.mean.many");
  
   truncated.mean = joint.mean*as.numeric(abs(joint.mean)>10^-10)
  
   return(truncated.mean/trans.prob)
}



timeave.cond.mean.many = function(t,lambda,mu,nu=0,X0=1,delta=0.001,n=1024, prec.tol=1e-12, prec.fail.stop=TRUE){
  joint.mean = timeave.joint.mean.many(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,delta=delta,n=n)
  trans.prob = process.prob.many(t=t,lambda=lambda,mu=mu,nu=nu,X0=X0,n=n)
#$#  prec.error.handler(joint.mean,trans.prob,prec.tol,prec.fail.stop,
##                     fnid="timeave.cond.mean.one");
  truncated.mean = joint.mean*as.numeric(abs(joint.mean)>10^-10)  
  return(truncated.mean/trans.prob)
}

#same as below but more functionality?  acts differently if lists are passed.
process.prob.one.fft <- function (t, lambda, mu, nu = 0, X0 = 1, Xt, n = 1024) 
{
  apply(as.matrix(t),1,
        function(z){
          max(power.coef.one(process.generator, n = n, k = Xt, t = z, 
                         lambda = lambda, mu = mu, nu = nu, X0 = X0),
              0);
        });
#  return(power.coef.one(process.generator, n = n, k = Xt, t = t, 
#                        lambda = lambda, mu = mu, nu = nu, X0 = X0))
}

#########See OPS.R for
########   process.prob.one <- function definitione






## ##tested all.cond.mean and all.cond.mean2
## all.cond.mean.PO.new <- function(data,lambda,mu,nu=0,
##                                     delta=0.001,n=1024, r=4,prec.tol=1e-12, prec.fail.stop=TRUE){
##   theArg <- CTMCPO2indepIntervals(data);  ##accepts either ctmcpo1 or ctmcpomany
##   #### SHOULD PRECOMPUTE transition probabilities AND PASS THEM IN 
##   Nplus <- sum(apply(theArg, 1, function(arg){
##     add.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,
##                       delta=delta,n=n, r=r, prec.tol=prec.tol, prec.fail.stop=prec.fail.stop)
##   }))
##   Nminus <- sum(apply(theArg, 1, function(arg){
##     rem.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n, prec.tol=prec.tol, r=r,prec.fail.stop=prec.fail.stop)
##   }))
##   Holdtime <- sum(apply(theArg, 1, function(arg){
##     hold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n, prec.tol=prec.tol, r=r,prec.fail.stop=prec.fail.stop)
##   }))
##   res <- c(Nplus,Nminus,Holdtime);
##   names(res) <- c("Nplus","Nminus","Holdtime");
##   return(res);
## }



## # '2' = second-order means, ie meansquares and crossproducts
## ## Well, apparently this is slower than the original "all.cond.mean2.PO",
## ## despite the fact that it does less actual slow computations. Apparently
## ## the increased calls to apply() are all that really matter.  Forgot the
## ## golden rule of never worrying about speed in R
## ## (Also, I tested 'mapply' vs. 'apply' and the mapply calls apparently seem to be faster
## ## (the time to do a cbind is insignificant).  I don't understand this.
## all.cond.mean2.PO.new <- function(data,lambda,mu,nu=0,
##                                   delta=0.001,n=1024, r=4,prec.tol=1e-12, prec.fail.stop=TRUE){
##   theArg <- CTMCPO2indepIntervals(data); ##accepts either ctmcpo1 or ctmcpomany
##   ##Get E(N+^2 | all data)

##   argList <- list(lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r,
##                    prec.tol=prec.tol,prec.fail.stop=prec.fail.stop);
##   trans.probs <- mapply(process.prob.one,
##                         t=theArg[,3],
##                         X0=theArg[,1],
##                         Xt=theArg[,2],
##                         MoreArgs=list(lambda=lambda, mu=mu,nu=nu),
##                         SIMPLIFY=TRUE) ## mapply or apply faster once we already have theArg?

##   ## ENplusSq.sumi <- sum(mapply(add.cond.meanSq.one,
##   ##                             t=theArg[,3],
##   ##                             X0=theArg[,1],
##   ##                             Xt=theArg[,2],
##   ##                             trans.prob=trans.probs,
##   ##                             MoreArgs=argList,
##   ##                             SIMPLIFY=TRUE))
  
##   theArg <- cbind(theArg, trans.probs) ##used for prod-means, e.g. addrem
  
##   ENplusi.joint <- apply(theArg, 1, function(arg){
##     add.joint.mean.one(t=arg[3],X0=arg[1],Xt=arg[2], 
##                        lambda=lambda,mu=mu,nu=nu,
##                        delta=delta,n=n,r=r)
##   })

##   ENminusi.joint <- apply(theArg,1,function(arg){
##     rem.joint.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],
##                        lambda=lambda,mu=mu,nu=nu,
##                        delta=delta,n=n,r=r)
##   })

##   EHoldtimei.joint <- apply(theArg,1,function(arg){
##     hold.joint.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],
##                         lambda=lambda,mu=mu,nu=nu,
##                         delta=delta,n=n,r=r)
##   })

  
##   ## ENplusi <- apply(theArg, 1, function(arg){
##   ##   add.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2], trans.prob=arg[4],
##   ##                     joint.mean=,
##   ##                     lambda=lambda,mu=mu,nu=nu,
##   ##                     delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   ## })
  
##   ENplusi <- mapply(add.cond.mean.one,
##                     t=theArg[,3],
##                     X0=theArg[,1],
##                     Xt=theArg[,2],
##                     trans.prob=trans.probs,
##                     joint.mean=ENplusi.joint,
##                     MoreArgs=argList,
##                     SIMPLIFY=TRUE)
##   ENplusSq.sumi <- sum(mapply(add.cond.meanSq.one,
##                     t=theArg[,3],
##                     X0=theArg[,1],
##                     Xt=theArg[,2],
##                     trans.prob=trans.probs,
##                     joint.mean=ENplusi.joint,
##                     MoreArgs=argList,
##                     SIMPLIFY=TRUE))
##   ## ENplusSq.sumi <- sum(apply(theArg, 1, function(arg){
##   ##   add.cond.meanSq.one(t=arg[3],X0=arg[1],Xt=arg[2],
##   ##                       trans.prob=arg[4],
##   ##                       lambda=lambda,mu=mu,nu=nu,                        
##   ##                       delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   ## }))
##   sqENplus.sumi <- sum(ENplusi^2)
##   ENplus <- sum(ENplusi)
##   NplusSq <- ENplusSq.sumi-sqENplus.sumi + ENplus^2; #final answer.
  
##   ##Get E(N-^2 | all data)

##   ENminusi <- mapply(rem.cond.mean.one,
##                      t=theArg[,3], X0=theArg[,1], Xt=theArg[,2],
##                      trans.prob=trans.probs, joint.mean=ENminusi.joint,
##                      MoreArgs=argList, SIMPLIFY=TRUE)
                        
##   ## ENminusi <- apply(theArg, 1, function(arg){
##   ##   rem.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,
##   ##                     delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   ## })

##   ENminusSq.sumi <- sum(mapply(rem.cond.meanSq.one,
##                            t=theArg[,3], X0=theArg[,1], Xt=theArg[,2],
##                            trans.prob=trans.probs, joint.mean=ENminusi.joint,
##                            MoreArgs=argList,
##                            SIMPLIFY=TRUE))
##   ## ENminusSq.sumi <- sum(apply(theArg, 1, function(arg){
##   ##   rem.cond.meanSq.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   ## }))
##   sqENminus.sumi <- sum(ENminusi^2)
##   ENminus <- sum(ENminusi)
##   NminusSq <- ENminusSq.sumi-sqENminus.sumi + ENminus^2; #final answer.

##   ##Get E(R^2 | all data)
##   ##   EHoldtimei <- apply(theArg, 1, function(arg){
##   ##   hold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   ## });
##   EHoldtimei <- mapply(hold.cond.mean.one,
##                        t=theArg[,3], X0=theArg[,1], Xt=theArg[,2],
##                        trans.prob=trans.probs,
##                        joint.mean=EHoldtimei.joint, ## not  necesary for holdtime to precompute..
##                        MoreArgs=argList, SIMPLIFY=TRUE)

##   EHoldtimeSq.sumi <- sum(mapply(hold.cond.meanSq.one,
##                                  t=theArg[,3], X0=theArg[,1], Xt=theArg[,2],
##                                  trans.prob=trans.probs,
##                                  MoreArgs=argList, SIMPLIFY=TRUE))
  
##   ## EHoldtimeSq.sumi <- sum(apply(theArg, 1, function(arg){
##   ##   hold.cond.meanSq.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   ## }))  
##   sqEHoldtime.sumi <- sum(EHoldtimei^2)
##   EHoldtime <- sum(EHoldtimei)
##   HoldtimeSq <- EHoldtimeSq.sumi-sqEHoldtime.sumi + EHoldtime^2; #final answer.
  
##   ##Get E( (N+)(N-) | all data )
##   ENplusNminus.sumi <- sum(apply(theArg, 1, function(arg){
##     addrem.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2], trans.prob=arg[4],
##                          lambda=lambda,mu=mu,nu=nu,
##                          delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   }))
##   ENplusENminus.sumi <- sum(ENplusi*ENminusi)
##   ##ENplus,ENminus already computed
##   NplusNminus <- ENplusNminus.sumi - ENplusENminus.sumi + ENplus*ENminus;

##   ##Get E( (N+)R | all data )
##   ENplusHoldtime.sumi <- sum(apply(theArg, 1, function(arg){
##     addhold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2], trans.prob=arg[4],
##                           lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r,
##                           prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   }))
##   ENplusEHoldtime.sumi <- sum(ENplusi*EHoldtimei)
##   ##ENplus Eholdtime already computed
##   NplusHoldtime <- ENplusHoldtime.sumi - ENplusEHoldtime.sumi + ENplus*EHoldtime;

##   ##Get E( (N-)R | all data )
##   ENminusHoldtime.sumi <- sum(apply(theArg, 1, function(arg){
##     remhold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],trans.prob=arg[4],
##                           lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   }))
##   ENminusEHoldtime.sumi <- sum(ENminusi*EHoldtimei)
##   ##ENminus, EHoldtime already computed
##   NminusHoldtime <- ENminusHoldtime.sumi - ENminusEHoldtime.sumi + ENminus*EHoldtime;

##   ##NminusSq<-HoldtimeSq<- NplusNminus<- NplusHoldtime<- NminusHoldtime<-
##   ##         ENplus<- ENminus<- EHoldtime <- NULL;
  
##   res <- c(NplusSq,NminusSq,HoldtimeSq, NplusNminus, NplusHoldtime, NminusHoldtime,
##            ENplus, ENminus, EHoldtime);
##   names(res) <- c("NplusSq","NminusSq","HoldtimeSq",
##                   "NplusNminus", "NplusHoldtime", "NminusHoldtime",
##                   "Nplus", "Nminus", "Holdtime")
##   return(res);
## }


##tested all.cond.mean and all.cond.mean2
all.cond.mean.PO <- function(data,lambda,mu,nu=0,delta=0.001,n=1024, r=4,prec.tol=1e-12, prec.fail.stop=TRUE){
  theArg <- CTMCPO2indepIntervals(data);  ##accepts either ctmcpo1 or ctmcpomany
  Nplus <- sum(apply(theArg, 1, function(arg){
    add.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,
                      delta=delta,n=n, r=r, prec.tol=prec.tol, prec.fail.stop=prec.fail.stop)
  }))
  Nminus <- sum(apply(theArg, 1, function(arg){
    rem.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n, prec.tol=prec.tol, r=r,prec.fail.stop=prec.fail.stop)
  }))
  Holdtime <- sum(apply(theArg, 1, function(arg){
    hold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n, prec.tol=prec.tol, r=r,prec.fail.stop=prec.fail.stop)
  }))
  res <- c(Nplus,Nminus,Holdtime);
  names(res) <- c("Nplus","Nminus","Holdtime");
  return(res);
}

# '2' = second-order means, ie meansquares and crossproducts
all.cond.mean2.PO <- function(data,lambda,mu,nu=0,delta=0.001,n=1024, r=4,prec.tol=1e-12, prec.fail.stop=TRUE){
  theArg <- CTMCPO2indepIntervals(data); ##accepts either ctmcpo1 or ctmcpomany
  ##Get E(N+^2 | all data)

  ENplusSq.sumi <- sum(apply(theArg, 1, function(arg){
    add.cond.meanSq.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,
                        delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
  }))
  
  ENplusi <- apply(theArg, 1, function(arg){
    add.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,
                      delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
  })
  sqENplus.sumi <- sum(ENplusi^2)
  ENplus <- sum(ENplusi)
  NplusSq <- ENplusSq.sumi-sqENplus.sumi + ENplus^2; #final answer.
  ## print(ENplusSq.sumi)  
  ## print(ENplusi)
  ## print(sqENplus.sumi)
  ## print(ENplus)
  
  ##Get E(N-^2 | all data)  
  ENminusSq.sumi <- sum(apply(theArg, 1, function(arg){
    rem.cond.meanSq.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
  }))
  ENminusi <- apply(theArg, 1, function(arg){
    rem.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,
                      delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
  })
  sqENminus.sumi <- sum(ENminusi^2)
  ENminus <- sum(ENminusi)
  NminusSq <- ENminusSq.sumi-sqENminus.sumi + ENminus^2; #final answer.

  ##Get E(R^2 | all data)  
  EHoldtimeSq.sumi <- sum(apply(theArg, 1, function(arg){
    hold.cond.meanSq.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
  }))
  EHoldtimei <- apply(theArg, 1, function(arg){
    hold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
  });
  sqEHoldtime.sumi <- sum(EHoldtimei^2)
  EHoldtime <- sum(EHoldtimei)
  HoldtimeSq <- EHoldtimeSq.sumi-sqEHoldtime.sumi + EHoldtime^2; #final answer.
  
  ##Get E( (N+)(N-) | all data )
  ENplusNminus.sumi <- sum(apply(theArg, 1, function(arg){
    addrem.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
  }))
  ENplusENminus.sumi <- sum(ENplusi*ENminusi)
  ##ENplus,ENminus already computed
  NplusNminus <- ENplusNminus.sumi - ENplusENminus.sumi + ENplus*ENminus;

  ##Get E( (N+)R | all data )
  ENplusHoldtime.sumi <- sum(apply(theArg, 1, function(arg){
    addhold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
  }))
  ENplusEHoldtime.sumi <- sum(ENplusi*EHoldtimei)
  ##ENplus Eholdtime already computed
  NplusHoldtime <- ENplusHoldtime.sumi - ENplusEHoldtime.sumi + ENplus*EHoldtime;

  ##Get E( (N-)R | all data )
  ENminusHoldtime.sumi <- sum(apply(theArg, 1, function(arg){
    remhold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
  }))
  ENminusEHoldtime.sumi <- sum(ENminusi*EHoldtimei)
  ##ENminus, EHoldtime already computed
  NminusHoldtime <- ENminusHoldtime.sumi - ENminusEHoldtime.sumi + ENminus*EHoldtime;

  ##NminusSq<-HoldtimeSq<- NplusNminus<- NplusHoldtime<- NminusHoldtime<-
  ##         ENplus<- ENminus<- EHoldtime <- NULL;
  
  res <- c(NplusSq,NminusSq,HoldtimeSq, NplusNminus, NplusHoldtime, NminusHoldtime,
           ENplus, ENminus, EHoldtime);
  names(res) <- c("NplusSq","NminusSq","HoldtimeSq",
                  "NplusNminus", "NplusHoldtime", "NminusHoldtime",
                  "Nplus", "Nminus", "Holdtime")
  return(res);
}


##backup -- the current version _has_ been tested on a handful of values against this backup version.

## ## # '2' = second-order means, ie meansquares and crossproducts
## all.cond.mean2.PO.backup <- function(data,lambda,mu,nu=0,delta=0.001,n=1024, r=4,prec.tol=1e-12, prec.fail.stop=TRUE){
##   theArg <- CTMCPO2indepIntervals(data); ##accepts either ctmcpo1 or ctmcpomany
##   ##Get E(N+^2 | all data)
##   ENplusSq.sumi <- sum(apply(theArg, 1, function(arg){
##     add.cond.meanSq.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,
##                         delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   }))
##   sqENplus.sumi <- sum(apply(theArg, 1, function(arg){
##     add.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,
##                       delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)^2
##   }))
##   ENplus <- sum(apply(theArg, 1, function(arg){
##     add.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,
##                       delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   }))
##   NplusSq <- ENplusSq.sumi-sqENplus.sumi + ENplus^2; #final answer.

##   ##Get E(N-^2 | all data)  
##   ENminusSq.sumi <- sum(apply(theArg, 1, function(arg){
##     rem.cond.meanSq.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   }))
##   sqENminus.sumi <- sum(apply(theArg, 1, function(arg){
##     rem.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,
##                       delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)^2
##   }))
##   ENminus <- sum(apply(theArg, 1, function(arg){
##     rem.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,
##                       delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   }))
##   NminusSq <- ENminusSq.sumi-sqENminus.sumi + ENminus^2; #final answer.

##   ##Get E(R^2 | all data)  
##   EHoldtimeSq.sumi <- sum(apply(theArg, 1, function(arg){
##     hold.cond.meanSq.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   }))
##   sqEHoldtime.sumi <- sum(apply(theArg, 1, function(arg){
##     hold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)^2
##   }))
##   EHoldtime <- sum(apply(theArg, 1, function(arg){
##     hold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   }))
##   HoldtimeSq <- EHoldtimeSq.sumi-sqEHoldtime.sumi + EHoldtime^2; #final answer.
  
##   ##Get E( (N+)(N-) | all data )
##   ENplusNminus.sumi <- sum(apply(theArg, 1, function(arg){
##     addrem.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   }))
##   ENplusENminus.sumi <- sum(apply(theArg, 1, function(arg){
##     add.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop) *
##       rem.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop);
      
##   }))
##   ##ENplus,ENminus already computed
##   NplusNminus <- ENplusNminus.sumi - ENplusENminus.sumi + ENplus*ENminus;

##   ##Get E( (N+)R | all data )
##   ENplusHoldtime.sumi <- sum(apply(theArg, 1, function(arg){
##     addhold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   }))
##   ENplusEHoldtime.sumi <- sum(apply(theArg, 1, function(arg){
##     add.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop) *
##       hold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop);
      
##   }))
##   ##ENplus Eholdtime already computed
##   NplusHoldtime <- ENplusHoldtime.sumi - ENplusEHoldtime.sumi + ENplus*EHoldtime;

##   ##Get E( (N-)R | all data )
##   ENminusHoldtime.sumi <- sum(apply(theArg, 1, function(arg){
##     remhold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop)
##   }))
##   ENminusEHoldtime.sumi <- sum(apply(theArg, 1, function(arg){
##     rem.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop) *
##       hold.cond.mean.one(t=arg[3],X0=arg[1],Xt=arg[2],lambda=lambda,mu=mu,nu=nu,delta=delta,n=n,r=r, prec.tol=prec.tol,prec.fail.stop=prec.fail.stop);
      
##   }))
##   ##ENminus, EHoldtime already computed
##   NminusHoldtime <- ENminusHoldtime.sumi - ENminusEHoldtime.sumi + ENminus*EHoldtime;
  
##   res <- c(NplusSq,NminusSq,HoldtimeSq, NplusNminus, NplusHoldtime, NminusHoldtime);
##   names(res) <- c("NplusSq","NminusSq","HoldtimeSq",
##                   "NplusNminus", "NplusHoldtime", "NminusHoldtime")
##   return(res);
## }



