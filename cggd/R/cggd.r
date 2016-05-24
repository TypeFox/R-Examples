# CGGD package
# Copyright 2007 Cun-Hui Zhang and Ofer Melnik
    
###############################################
# cggd:
#
# This is the main CGGD function. It iteratively builds a CGGD model in up to kmax iterations.  

cggd <- function(x, y, beta0 = rep(0,2), kmax = 300, 
                      TRR=FALSE, t0 = 1, TRACE=FALSE,
                      alpha1 = 0, alpha2 = 0, w = 1, 
                      tau = 1, tautil = -1, eps = -1,fctr=1e8)
{
   
  n <- length(y) # n - number of data points 

  p <- ceiling(length(x)/n) # p - number of potential variables

  if (length(beta0)!=p) beta0 <- rep(0,p) # Make sure initial betas have correct dimensionality

  if ((tau>1)||(tau<=0)) tau <- 1  # Verify that tau is correct

  # scale and normalize variables
  mu <- mean(y)
  y <- drop(y - mu)

  one <- rep(1, n)
  meanx <- drop(one %*% x)/n
  x <- scale(x, meanx, FALSE) # centers x
  normx <- sqrt(drop(one %*% (x^2)))
  x <- scale(x, FALSE, normx) # scales x


  if (p==1) { # if one variable just do least squares 
    olse <- sum(x*y)/sum(x^2)
    beta.tk <- c(beta0, olse)
    tk <-  0
    a.set.tk <- rep(TRUE,2)
    kmax <- 1      
    obj <-list( beta.tk=beta.tk, tk=tk, a.set.tk=a.set.tk, kmax=kmax, p=p)

    return(obj)
  }  

  # now p > 1
  eps <- eps * median( abs(t(x)%*%y)/apply(x^2,2,sum) ) # median of least squares solution for individual least squares

  # variable declaration, beta is dynamic 
  beta.tk <- matrix(beta0, p, kmax+1)
  tk <- rep(0, kmax)   
  a.set.tk <- matrix(FALSE, p, kmax+1) 

  # compute the active set of model 1 
  a.set <- rep(FALSE, p)
  g0 <- t(x)%*%(y - x%*%beta0)  # O: gradient at t0
  gmax <- max(abs(g0))
  a.set[ abs(g0) >= tau * gmax ]  <- TRUE # active set for model 1 
  a.set.tk[,1] <- a.set 

  if (TRR==FALSE)  {   tk[1] <- 0}
  else        {   tk[1] <- t0}


  k <- 2 
  max.tk <- 0

  success <- TRUE
  while (k <= kmax+1 && success) { # compute the trajectory for kmax -1 models
    # step k 
    if (TRR==FALSE)  {

      step.k<-cggd.mgd.one.step(x, y, a.set, beta.tk[,k-1], # [,k] is the kth column of beta.tk
               alpha1=alpha1, alpha2=alpha2, 
              w=w, tau=tau, tautil=tautil, eps=eps, TRACE)
    }
    else {
      if (k>1) t0 <- tk[k-1]
      step.k<-cggd.trr.one.step(x, y, a.set, beta.tk[,k-1],   
                                t0=t0, tau=tau, tautil=tautil, eps=eps)
    }



    # update a.set, a.set.tk, tk and k 

    if(TRACE) cat("k =", k," model size = ",sum(a.set), " gmax = ", step.k$gmax, "\n")

    success <- step.k$success        
    if(success)
    {
      beta.tk[,(k):(kmax+1)] <- matrix(step.k$beta[,1], p, kmax+2-k)       
      a.set <- step.k$a.set
      a.set.tk[,k] <- a.set
      tk[k] <- step.k$t1
      if (TRR==FALSE)  { max.tk <- max.tk + step.k$t1 }
      else        { max.tk <- step.k$t1 }
      k <- k+1
    }
  } # end while (k <= kmax) 
  # kmax may have changed due to early stopping

  k <- k-1
  beta.tk <- beta.tk[,1:k]  
  tk <- tk[1:k]
  a.set.tk <- a.set.tk[,1:k]


  obj <-list( beta.tk=beta.tk, tk=tk, a.set.tk=a.set.tk,  max.tk=max.tk, min.tk=tk[1], 
    kmax=length(tk), p=p, TRR=TRR, alpha1=alpha1, alpha2=alpha2,tau=tau,tautil=tautil,w=w,
    x=x,y=y,mu=mu, normx = normx, meanx = meanx)
  class(obj) <- "cggd"

  return(obj)
}
  
  
  
###############################################
# get.switch.time:
#
# This function finds the location of the next switch time
# It works in two stages. First it uses bump following to find consecutive local minimas
# of the switch time function, until it finds a negative minima.
# Then it uses root finding to narrow in on the switch time.


get.switch.time <- function(switch.func,e, boundary=0,terminate_value=1e15,tol = .Machine$double.eps^0.5, 
  fctr=1e8, TRACE=FALSE)
{

 # looking for negative values of the switch.function


 last_val <-0
 find_minimum<-TRUE
 more<-TRUE
 last_boundary <- 0
 push_boundary <- FALSE

 

 while(more)
 {

    if(push_boundary)
    {
      boundary <- boundary + 1e0
    }

    if(find_minimum) {
     res<-optim(par=boundary,fn=switch.func,method="L-BFGS-B",lower=boundary,control=list(factr=fctr,pgtol=0),upper=Inf, gr=NULL,hessian=FALSE,e)
    }
    else { # maximize over the next hump
     res<-optim(par=boundary,fn=switch.func,method="L-BFGS-B",lower=boundary,control=list(fnscale=-1,factr=fctr,pgtol=0),upper=Inf, gr=NULL,hessian=FALSE,e)       
    }

    if(TRACE)
    { cat("res par =", res$par," res value = ",res$value," counts =", res$counts, "fctr = ",fctr," convergence =", res$convergence, " message = ",res$message,"\n") }

    if(res$par >= terminate_value || res$value >= terminate_value ){ return(list(switch.time=0,success=FALSE)) }


    if(res$value <0) { more<-FALSE }
    else
    {
     if(round(last_val-res$value,digits=5)==0 || round(res$par-boundary,digits=5)==0) 
     { # stuck in place, so tighten search parameters
       fctr <- 1e-1
     }
     if(!res$convergence)
     { # got some error so push forward a little
       push_boundary <- TRUE   
     }
    }

    find_minimum <- !find_minimum
    if(!push_boundary)
    { last_boundary<-boundary }
    boundary <- res$par
    last_val <-res$value
 }

 switch.time<-res$par
 
 # Find the exact switch time
 # can further optimize by figuering out it is in the m.in or m.out and only doing those
 # can also check which variables in particular are negative and only do it for them
 if( res$value <0  && res$par !=last_boundary) {
     r<-uniroot(f=switch.func,int=c(last_boundary,boundary),lower=last_boundary,upper=boundary,maxiter=1000,tol=tol,e)

     switch.time <- r$root + r$estim.prec 
 }

 return(list(switch.time=switch.time,success=TRUE))
}
 

###############################################
# cggd.mgd.var.init:
#
# Initializes variables which are used by the cggd.mgd.gen.beta and cggd.switch.func functions.

cggd.mgd.var.init <- function(x,y,  a.set, beta_init, e,  
                    alpha1 = 0, alpha2 = 0, w = 1)
{

  e$q <-sum(a.set)
  e$n <- length(y)
  e$p <- length(beta_init)
  e$beta_init <- beta_init
  e$a.set <- a.set
  

  # the decomposition of P(t)GP(t)
  e$txa<-t(x[,a.set])
  e$txna<-t(x[,!a.set])
  e$svdG <-svd(e$txa %*% x[,a.set])

  Q1 <-  e$svdG$d  # diagonal, eigenvalues
  Q2 <-  Q1      # do twice for kernel calculations
  zeros <-  Q1 == 0 
  Q1[!zeros] <- Q1[!zeros]^alpha1
  if (alpha1==0) Q1[zeros]  <- 1
  Q2[!zeros] <- Q2[!zeros]^alpha2
  if (alpha2==0) Q2[zeros]  <- 1
  Q <- w*Q1 + (1-w)*Q2  # kernel calculation

  e$g0 <- t(x)%*%(y - x%*%beta_init) 

  # this is the kernel x g0
  
  tsvdGv <- t(e$svdG$v)
  e$Qg0<-(e$svdG$u)%*%(Q * c(tsvdGv%*%e$g0[a.set]))

  e$sQg0 <-c(tsvdGv%*%e$Qg0)

  # H = Q*G in eigen space
  e$H<-Q * (e$svdG$d)

  e$gen.beta<-cggd.mgd.gen.beta
}
    
 
###############################################
# cggd.mgd.gen.beta:
#
# This function calculates beta values at a particular t for the mgd model.

cggd.mgd.gen.beta<-function(t,e)
{
   beta <- matrix(e$beta_init,e$p,1)
    
   if (e$q==1) {
     Fstar <-  - t*e$H
     Fstar <- cggd.Fstar(Fstar)
     beta_idx <- beta[e$a.set] + t*Fstar*e$Qg0
     beta[e$a.set] <- beta_idx
   } # end if (q==1)
   if (e$q>1) {
   # H matrix for every eigen value x time points calculate 
     Fstar <-  - e$H*t
     Fstar <- cggd.Fstar(Fstar)

     ### on the active set   beta(t) = beta0 + t F*(-t H) Q g0, H = QG (in eigen space)
     beta_idx <- beta[e$a.set] + (e$svdG$u)%*%(Fstar*e$sQg0)*t # eq 6.65
     beta[e$a.set] <- beta_idx
   } # end if (q>1)

  return(list(beta=beta,beta_idx=beta_idx));
}

###############################################
# cggd.switch.func:
#
# Calculates which variables should be added/removed at a particular time.

cggd.switch.func <- function(t, e, do_in=TRUE, do_out=TRUE)
{
   r=e$gen.beta(t,e)

   C <- (matrix(e$y,e$n,1)-e$x%*%r$beta)
   g.in<-abs(e$txa%*%C)  
   g.in.max <- max(g.in)

   # check for enters
   if(e$q < e$p && do_in)
   {
     g.out <-abs(e$txna%*%C)
     g.out.max <- max(g.out)

     m.in <- e$tau*g.in.max - g.out.max 
   }
   else m.in <- 1e99

   # check for exits
   if (do_out && length(g.in)>1) 
   {
     i <- which.min(g.in)
     g.in.min <- g.in[i]
     if( abs( r$beta_idx[i] ) <e$eps ) {
       m.out <- g.in.min - e$tautil*g.in.max
     }
     else m.out <- 1e99 
   }
   else m.out <- 1e99

   return(list(m.in=m.in,m.out=m.out))     
}

   
cggd.switch.func.min <- function(t,e)
{
  a <- cggd.switch.func(t,e);
  m<-min(c(a$m.in,a$m.out))

  return(m)
}
    

###############################################
# cggd.mgd.one.step:
#
# Finds the next iteration of an mgd model.

cggd.mgd.one.step <- function(x, y, a.set, beta_init,    t1 = 1, 
                      alpha1 = 0, alpha2 = 0, w = 1, tau = 1, tautil = -1, eps = -1, TRACE=FALSE, fctr=1e8)
{
   
   g0 <- t(x)%*%(y - x%*%beta_init) 
   gmax <- max(abs(g0))
   a.set[ abs(g0) >= tau * gmax ]  <- TRUE
   
   # It cant be easier to exit than enter model.
   if (tautil >= tau) {
     tautil <- tau
     w <- 1
     alpha1 <-  -1
   }
      
   my.env <-environment()
   
   cggd.mgd.var.init(x=x,y=y, a.set=a.set, beta_init=beta_init, alpha1=alpha1, alpha2=alpha2, w=w, e=my.env )

    ##############################
    
    
   switch.time <- get.switch.time (cggd.switch.func.min,my.env,fctr=fctr)  
   
   info.init <- cggd.mgd.gen.beta(0,my.env)
   switch.info <- cggd.mgd.gen.beta(switch.time$switch.time,my.env)
   
   g.init <- t(x)%*%(y-x%*%info.init$beta) # the starting g (not sure why we need it)
   g <- t(x)%*%(y-x%*%switch.info$beta)  # the g at the switch time
   g.abs <- abs(g)
   gmax <- max(g.abs[a.set])
   a.set.init <- a.set
   a.set[g.abs >= tau*gmax] <- TRUE
   a.set[ (g.abs <= tautil*gmax) ]  <- FALSE  # & (abs(beta[,m+1]) <= eps) 
   return(list(success=switch.time$success,beta=switch.info$beta, a.set.init=a.set.init, a.set=a.set, 
    g.init=g.init, g=g, g.abs=g.abs, gmax=gmax, t1=switch.time$switch.time))
    
}

###############################################
# cggd.trr.var.init:
#
# Initializes variables which are used by the cggd.trr.gen.beta and cggd.switch.func functions.
 
cggd.trr.var.init <- function(x,y,  a.set, beta_init, e)
{

  e$q <-sum(a.set)
  e$n <- length(y)
  e$p <- length(beta_init)
  e$beta_init <- beta_init
  e$a.set <- a.set
  e$g0 <- t(x)%*%(y - x%*%beta_init) 

  # the decomposition of P(t)GP(t)
  e$txa<-t(x[,a.set])
  e$txna<-t(x[,!a.set])
  e$svdG <-svd(e$txa %*% x[,a.set])
  e$tsvdGv<-t(e$svdG$v)

  G.inverse <- e$svdG$d
  G.inverse[e$svdG$d<1e-05] <- 1e-06 
  G.inverse <- 1/G.inverse
  e$G.inverse.g0 <-  (e$svdG$u)%*%(G.inverse * c(e$tsvdGv%*%e$g0[a.set]))

  e$svGinvg0 <- e$tsvdGv%*%e$G.inverse.g0
  e$gen.beta<-cggd.trr.gen.beta
}

###############################################
# cggd.trr.gen.beta:
#
# This function calculates beta values at a particular t for the trr model.

cggd.trr.gen.beta<-function(t,e)
{


    beta <- matrix(e$beta_init,e$p,1)
    if (e$q==1) {
      Kstar <-  1 - (e$t0*e$svdG$d+1)/(t*e$svdG$d+1)
      beta_idx <- beta[e$a.set,] + Kstar * e$G.inverse.g0
      beta[e$a.set] <- beta_idx
    } # end if (q==1)
    if (e$q>1) {
      Kstar <- 1-(e$t0*e$svdG$d+1)/(t*e$svdG$d+1)
      Kstarg0 <- (e$svdG$u)%*%(Kstar * e$svGinvg0)   
      beta_idx <- beta[e$a.set,] + Kstarg0
      beta[e$a.set] <- beta_idx
    } # end if (q>1)

  return(list(beta=beta,beta_idx=beta_idx));
}

###############################################
# cggd.trr.one.step:
#
# Finds the next iteration of an trr model.
 
cggd.trr.one.step <- function(x, y, a.set, beta_init, t1 = 1, 
                     t0 = 1, tau = 1, tautil = -1, eps = -1, fctr=1e8)
{   
   g0 <- t(x)%*%(y - x%*%beta_init) 
   gmax <- max(abs(g0))
   a.set[ abs(g0) >= tau * gmax ]  <- TRUE
         
   my.env <-environment()
   
   cggd.trr.var.init(x=x,y=y, a.set=a.set, beta_init=beta_init, e=my.env )

   ##############################
    
   switch.time <- get.switch.time (switch.func=cggd.switch.func.min,e=my.env,boundary=t0,fctr=fctr)  
   
   info.init <- cggd.trr.gen.beta(t0,my.env)
   switch.info <- cggd.trr.gen.beta(switch.time$switch.time,my.env)
   
   g.init <- t(x)%*%(y-x%*%info.init$beta) # the starting g (not sure why we need it)
   g <- t(x)%*%(y-x%*%switch.info$beta)  # the g at the switch time
   g.abs <- abs(g)
   gmax <- max(g.abs[a.set])
   a.set.init <- a.set
   a.set[g.abs >= tau*gmax] <- TRUE
   a.set[ (g.abs <= tautil*gmax) ]  <- FALSE  # & (abs(beta[,m+1]) <= eps) 
   return(list(success=switch.time$success,beta=switch.info$beta, a.set.init=a.set.init, a.set=a.set, 
    g.init=g.init, g=g, g.abs=g.abs, gmax=gmax, t1=switch.time$switch.time))
}

 
 
###############################################

cggd.Fstar <- function(x, eps1=1e-6) {
  if (length(x)==1) 
    if (abs(x) > eps1) return( (exp(x)-1)/x )
    else return( 1 + x/2 + x^2/6 + x^3/24 + x^4/120)
  if (length(x) > 1) {
    Fstar <- x
    Fstar[ abs(x)>eps1 ] <- (exp(x[ abs(x)>eps1 ]) - 1)/x[ abs(x)>eps1 ]
    Fstar[ abs(x) <= eps1 ] <- 1 + x[ abs(x)<=eps1 ]/2 +  x[ abs(x)<=eps1 ]^2/6 + 
                               x[ abs(x)<=eps1 ]^3/24 +  x[ abs(x)<=eps1 ]^4/120
    return(Fstar)
  }
}



###############################################
# plot.cggd:
#
# Plots a graph showing the paths of the beta values of the cggd model at different interations.

plot.cggd <- function(x, steps=5, breaks = TRUE, first_k=1,last_k=Inf, xvar=c("step","t"),omit.zeros = TRUE, eps = 1e-10, ...)
{
  o<-x
  xvar <- match.arg(xvar)
  
  # generate the intermediate beta values
  d <- dim(o$beta.tk)  
  if(last_k>d[2]) { last_k=d[2] }
  
  
  intvl <- last_k-first_k +1
  num_points <- (intvl-1)*steps+intvl
  
  b <- matrix(0,d[1],num_points)
  s <- matrix(0,num_points,1)
  brk <- matrix(0,intvl,1);
  
  my.env <-environment()

  s_cum_tk <- 0
  
  b[,1]<-o$beta.tk[,1]
  s[1] <- o$tk[1] 
  
  brk[1] <- o$tk[1] 
  s_cum_tk <- o$tk[1]
  
  j <-2
  for (k in (first_k+1):last_k) 
  {
    
    if(o$TRR==FALSE)
    {
      cggd.mgd.var.init(x=o$x,y=o$y, a.set=o$a.set.tk[,k-1],beta_init=o$beta.tk[,k-1], alpha1=o$alpha1, alpha2=o$alpha2, w=o$w, e=my.env )      
    }
    else
    {
      t0 <- o$tk[k-1]
      cggd.trr.var.init(x=o$x,y=o$y, a.set=o$a.set.tk[,k-1],beta_init=o$beta.tk[,k-1], e=my.env )      
    }
    
    for (i in 1:steps)
    {
      if(o$TRR==FALSE) {
        tm <-  i*o$tk[k]/(steps+1)
      }
      else {
        tm <- o$tk[k-1] + i*(o$tk[k]-o$tk[k-1])/(steps+1)
      }
      r=gen.beta(tm,my.env)
      b[o$a.set[,k-1],j] <- r$beta[o$a.set[,k-1]]
      s[j] <- tm + s_cum_tk
      j <- j+1
    }

    b[o$a.set[,k],j]<-o$beta.tk[o$a.set[,k],k]
    s[j] <- o$tk[k] + s_cum_tk
    brk[k] <- o$tk[k] + s_cum_tk
    j <- j + 1

    if(o$TRR==FALSE) {
      s_cum_tk <- s_cum_tk + o$tk[k]
    }

  }
  
  b<-t(b)
  
  
  if(omit.zeros) {
      c1 <- drop(rep(1, nrow(b)) %*% abs(b))
      nonzeros <- c1 > eps
      cnums <- seq(nonzeros)[nonzeros]
      b <- b[, nonzeros]
   }
   else cnums <- seq(nrow(b))


  if(xvar=="step") { 
    s <- seq(from=first_k,to=last_k,by=1/(steps+1) ) 
    brk <- seq(from=first_k,to=last_k)    
  }
  
  matplot(s, b, xlab = xvar, ..., type = "b", 
          pch = "*", ylab = "Coefficients")
  #title(object$type,line=2.5)
  abline(h = 0, lty = 3)
  axis(4, at = b[nrow(b),  ], label = paste(cnums
                                        ), cex = 0.80000000000000004, adj = 0)
  if(breaks) {
    axis(3, at = brk, labels = paste(round(brk,digits=2)),cex=.8)
    abline(v = brk)
  }
}


#plot(o,xvar="step")
# diabetes.ntgd.2a <- cggd.lm(x.noise,y,m=1,kmax=11,tau=3/4)
 
#diabetes.ntgd.2a <- cggd.lm(x,y,kmax=11,tau=3/4)


###############################################
# predict.cggd:
#
# Generates the betas at a particular iteration(k)/t and or a prediction (y value)

predict.cggd <- function(object, newx, t, type = c("fit", "coefficients"), mode=c("k","t"), ...)
{
  o<-object
  # check parameters
  
  
  if( (mode=="t" &&( t < o$min.tk || t> max.tk) ) || (mode=="k" && ( t <1 ||t > length(o$tk))))
    { stop("t out of range.") }
  
  # find coefficients

  my.env <-environment()

  if(mode=="k")
  {
    it <- floor(t)
    if(it==t)
    {
      beta <- o$beta.tk[,t]
    }
    else
    {
      if(o$TRR==FALSE)
      {
        cggd.mgd.var.init(x=o$x,y=o$y, a.set=o$a.set.tk[,it],beta_init=o$beta.tk[,it], alpha1=o$alpha1, alpha2=o$alpha2, w=o$w, e=my.env )      
        tm <-  (t-it)*o$tk[it+1]
      }
      else 
      {
        t0 <- o$tk[it]
        cggd.trr.var.init(x=o$x,y=o$y, a.set=o$a.set.tk[,it],beta_init=o$beta.tk[,it], e=my.env )      
        tm <- o$tk[it] + (t-it)*(o$tk[it+1]-o$tk[it])
      }
      r <- gen.beta(tm,my.env)
      beta <- r$beta
    }
  }
  else  # mode == "t"
  {
    if(o$TRR==FALSE)
    {
      s <- 0
      for (k in 1:length(o$tk)) 
      {
        s <- s + o$tk[k]
        if(s==t)
        {
          beta <- o$beta.tk[,k]         
          break
        }
        if(s>t)
        {
          cggd.mgd.var.init(x=o$x,y=o$y, a.set=o$a.set.tk[,k-1],beta_init=o$beta.tk[,k-1], alpha1=o$alpha1, alpha2=o$alpha2, w=o$w, e=my.env )      
          r <- gen.beta(t-(s-o$a.set.tk[,k-1]), my.env)
          beta <- r$beta
          break
        }
      }
    }
    else
    {
      for (k in 1:length(o$tk)) 
      {
        if(o$tk[k]==t)
        {
          beta <- o$beta.tk[,k]         
          break
        }       
        if(o$tk[k]>t)
        {       
          t0 <- o$tk[k-1]
          cggd.trr.var.init(x=o$x,y=o$y, a.set=o$a.set.tk[,k-1],beta_init=o$beta.tk[,k-1], e=my.env )      
          r <- gen.beta(t,my.env)
          beta <- r$beta
          break
        }
      } # for
    } # TRR
  } # mode = "t"
  
    
  if(type=="coefficients") { return(beta) }
  else  { return(list(beta=beta, fit=scale(newx,o$meanx,o$normx)%*%beta+o$mu)) }
  
}
 
#predict(o,x,5,type="fit",mode="k")

###############################################
# cv.cggd:
#
# Performs a cross validation analysis of particular dataset. 

cv.cggd <-
function(x, y, nfolds = 6, kmax=40 , 
           trace = FALSE, plot.it = TRUE, se=TRUE, ...)
{
  
  ly <- length(y)
  all.folds <-  split(sample(1:ly), rep(1:nfolds, length = ly)) 
  residmat <- matrix(0, kmax+1, nfolds)
  
  min_kmax <- 1e50
  
  for(i in seq(nfolds)) {
    
    # browser()
    omit <- all.folds[[i]]    
    notomit <- rep(TRUE,ly)
    notomit[omit] <- FALSE
    o <- cggd(x[ notomit,  ], y[ notomit ], TRACE = trace, kmax=kmax, ...)

    fit <- matrix(0, length(omit), o$kmax)
    
    for(k in seq(o$kmax)) {
      #browser()
      a <-  predict(o, x[omit,  ,drop=FALSE], k, type ="fit", mode="k")
      fit[,k] <- a$fit
    }
    
    if(length(omit)==1)fit<-matrix(fit,nrow=1)
    
    residmat[1:o$kmax, i] <- apply((y[omit] - fit)^2, 2, mean)
    
    if(min_kmax>o$kmax) {min_kmax <- o$kmax}

    if(trace)
      cat("\n CV Fold", i, "\n\n")
  }
  cv <- apply(residmat[1:min_kmax,], 1, mean)
  cv.error <- sqrt(apply(residmat[1:min_kmax,], 1, var)/nfolds)
  object<-list(cv = cv, cv.error = cv.error)
  if(plot.it) plotCVcggd(object,se=se)
  
  return(object)
}

###############################################
# plotCVcggd :
#
# plots the results of a cross validation analysis, where each x corresponds to another model iteration.

plotCVcggd <-
function(cv.lars.object,se=TRUE){
  attach(cv.lars.object)
      plot(1:length(cv), cv, type = "b", ylim = range(cv, cv + cv.error, 
                                     cv - cv.error),xlab="step")
    if(se)
      error.bars(1:length(cv), cv + cv.error, cv - cv.error, 
                 width = 1/length(cv))
  detach(cv.lars.object)
  
invisible()
}
error.bars <-
function(x, upper, lower, width = 0.02, ...)
{
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

#cv.cggd(x2,y,kmax=40)
