checkin <- function(x,interval){
  nx = length(x)
  ni = length(interval)
  x[nx:(nx-ni*2+1)] = c(-interval,interval)
  return(x)
}
chain.name <- function (nc.chain, lq, lb, var.lab = NA){
  symb.q.p = symb.b.p = list
  if(lq==1) {
    symb.q = 'q'
    symb.q.p = lapply(0, function(i)bquote(q))
  }else{
    symb.q   = paste('Theta',0:(lq-1),sep = '')
    symb.q.p = lapply(0:(lq-1), function(i)bquote(theta[.(i)]))
  }

  if(lb==1) {
    symb.b   = 'B'
    symb.b.p = lapply(0, function(i)bquote(beta))
  }else{
    symb.b   = paste('Gamma',0:(lb-1),sep = '')
    symb.b.p = lapply(0:(lb-1), function(i)bquote(gamma[.(i)]))
  }

  symb     = c(symb.q   , symb.b)
  symb.p   = c(symb.q.p , symb.b.p)

  if((lb+lq) < nc.chain){
    symb.l   = paste('l', 1:(nc.chain - (lq+lb)) , sep = '')
    symb.l.p =    (lapply(1:(nc.chain - (lq+lb)) , function(i)bquote(lambda[.(i)])))
    symb     = c  (symb   , symb.l  )
    symb.p   = c  (symb.p , symb.l.p)
  }
  if(length(var.lab) > 1){
    symb  [1:(lb+lq)] = var.lab
    symb.p = symb
    las = 2
    marb = max(nchar(x = var.lab))/2
  }else{
    symb.p = as.expression(symb.p)
    las = 1
    marb = 5
  }
  return(list(symb = symb     ,
              symb.p = symb.p ,
              las = las       ,
              marb = marb)
  )
}

extract.vars = function (dw.object , reg = FALSE){
  nc.chain = ncol(dw.object$chain)
  lq       = dw.object$lq
  lb       = dw.object$lb
  para.q   = dw.object$para.q
  para.b   = dw.object$para.b
  names = chain.name (nc.chain, lq, lb, var.lab = NA)
  symb.p= names$symb.p
  symb  = names$symb

  if(reg){
    if(para.q & para.b){
      symb.p = names$symb.p[1:(lb+lq)]
      symb   = names$symb  [1:(lb+lq)]
    }else if (para.q && !para.b){
      symb.p = names$symb.p[1:lq]
      symb   = names$symb  [1:lq]
    }else if (!para.q && para.b){
      symb.p = names$symb.p[(lq+1):(lb+lq)]
      symb   = names$symb  [(lq+1):(lb+lq)]
    }else{
      symb.p = names$symb.p[1:2]
      symb   = names$symb  [1:2]
    }
  }
  return(list(symb=symb,symb.p=symb.p))
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}

# find mode of population
Mode <- function(x , adj = 1 , na.rm = TRUE , remove.outliers = TRUE) {
  if(length(x)<50000){
    if (remove.outliers){
      x   = remove_outliers(x,na.rm = na.rm)
    }
    x   = na.omit(x)
  }
  ux <- x #unique(x)
  if(length(ux)>20){
    den<-density(ux,adjust = adj)
    ux=den$x[which(den$y==max(den$y))]
  }else{
    ux = mean(ux,na.rm = 1)
    if(is.na(ux)) ux=0
  }
  return(ux[1])
}

# improper prior!
imp.d<-function(par, ...){
  #return(1/abs(par))
  #return(1/(1+par^2))
  #return(sign(par)/(1+par^2))
  return(par)
}

# laplace distribution
plaplace <- function(q, m=0, s=1){
  if(any(q<0))stop("q must contain non-negative values")
  if(any(s<=0))stop("s must be positive")
  u <- (q-m)/s
  t <- exp(-abs(u))/2
  ifelse(u<0,t,1-t)
}

dlaplace <- function(x, m=0, s=1, log=FALSE , normalized = TRUE){
  #if(any(s<=0))stop("s must be positive")
  s=abs(s)
  if(normalized){
    ss=s/2
  }else{
    ss=1
  }
  r=ss * exp(-abs(x-m)* s)
  if(log==TRUE) {r=log(r)}
  return (r)
}

rlaplace = function(n,mu,sigma){
  U = runif(n,0,1)
  #This will give negative value half of the time
  sign = ifelse(rbinom(n,1,.5)>.5,1,-1)
  y = mu + sign*(sigma)/sqrt(2)*log(1-U)
  y
}

ddw<-function(x, q=exp(-1),beta=1)
{
  if(max(q)> 1 | min(q) < 0){
    stop('q must be between 0 and 1', call. = FALSE)
  }
  if(min(beta) <= 0){
    stop('beta must be positive', call. = FALSE)
  }
  is.wholenumber <-function(x, tol = sqrt(.Machine$double.eps)) { abs(x - round(x)) < tol}
  for(i in 1:length(beta)){
    if(beta==1){
      res<-dgeom(x,1-q)
    } else{
      if(is.wholenumber(x)){
        res<-q^(x^(beta))-q^((x+1)^(beta))
      }else{
        res<- 1 #<<<< log(1)=0
        print("Warning message: non-integer value in ddw")
      }
    }
  }
  return(res)
}






pdw<-function(x,q=exp(-1),beta=1)
{
  if(x<0)
    res<-0
  else
    res<-1-q^((floor(x)+1)^(beta))
  return(res)
}

qdw<-function(p,q=exp(-1),beta=1)
{
  if(max(q)>1 | min(q) <0)
    stop('q must be between 0 and 1', call. = FALSE)
  if(min(beta) <= 0)
    stop('beta must be positive', call. = FALSE)
  if(p<=1-q)
    res<-0
  else
    res<-ceiling((log(1-p)/log(q))^(1/beta)-1)
  return(res)
}




rdw<-function(n,q=exp(-1),beta=1)
{
  if(max(q)>1 | min(q)<0){
    stop('q must be between 0 and 1', call. = FALSE)
  }
  if(min(beta) <=0){
    stop('beta must be positive', call. = FALSE)
  }
  y<-runif(n,0,1)
  x<-as.numeric(lapply(y,qdw,q=q,beta=beta))
  return(x)
}


digamma <- function (x, shape, scale = 1, log = FALSE)
{
  if (shape <= 0 | scale <= 0) {
    stop("Shape or scale parameter negative in dinvgamma().\n")
  }
  x=abs(x)
  alpha <- shape
  beta <- scale
  log.density <- alpha * log(beta+10^-10) - lgamma(alpha) - (alpha +
                                                               1) * log(x+10^-10) - (beta/x)
  if(log){
    return(log.density)
  }else{
    return(exp(log.density))
  }
}

dw.cov.matrix <- function(x,i){
  #per = round(i*exp(-sqrt(i)/1*runif(1,0,.05)))
  per = round(i*exp(-sqrt(i)/15*runif(1,0,.05)))
  #per=100
  if((is.matrix(x)==TRUE) & (per>5) ){
    res = cov(x[(i-per ):i,])
    res = make.positive.definite(res)
  }else{
    res = 1
  }
  return(res)
}


BF <- function (fun,  chain , est.stat = Mode , ...)
{
  options(warn = -1)
  mode.v = apply (chain , 2 , est.stat )
  fit = optim( par = mode.v , fn = fun ,
               gr = NULL, hessian = TRUE      ,
               control = list(fnscale = -1)        ,
               ...
  )
  options(warn = 0)
  mode = fit$par
  diag(fit$hessian)=diag(fit$hessian)+10^-8
  h = -solve(fit$hessian)
  p = length(mode)
  int = p/2 * log(2 * pi) + 0.5 * log(abs(det(h))) +
    fun(mode,...)

  stuff = list(mode = mode,
               var = h,
               int = int,
               converge = fit$convergence == 0
  )
  return(stuff)
}

make.positive.definite <-  function(m, tol)
{
  # Author:
  #   Copyright 2003-05 Korbinian Strimmer
  #   Rank, condition, and positive definiteness of a matrix
  #   GNU General Public License, Version 2

  # Method by Higham 1988

  if (!is.matrix(m)) {
    m = as.matrix(m)
  }

  d = dim(m)[1]
  if ( dim(m)[2] != d ) {
    stop("Input matrix is not square!")
  }

  es = eigen(m)
  esv = es$values

  if (missing(tol)) {
    tol = d*max(abs(esv))*.Machine$double.eps
  }
  delta =  2*tol
  # factor two is just to make sure the resulting
  # matrix passes all numerical tests of positive definiteness

  tau = pmax(0, delta - esv)
  dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)

  # print(max(DA))
  # print(esv[1]/delta)

  return( m + dm )
}

#generalized negative binomial density
dgnbinom <- function(y,mu,theta,log=FALSE){
  d = gamma(theta+y)/gamma(theta) * (1/factorial(y))  *
    (1/(mu+theta))^(y+theta) *
    (mu^y )* (theta^theta)
  if(log)d=log(d)
  return(d)
}
