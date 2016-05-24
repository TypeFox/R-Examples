cusum <- function(x,...)
  UseMethod("cusum")

cusum.default <- function(x, ...)
{
  if(length(levels(as.factor(x)))==2)
    out = bcusum(x,...)
  else
    out = ccusum(x,...)
  invisible(out)
}

Restart <- function(x,side='upper', k,h,...)
  UseMethod("Restart")

Restart.default <- function(x,side='upper', k,h,...)
{
  stopifnot(class(x)=="CCUSUM")
  Restart.CCUSUM(x,side=side,k=k,h=h,...)  
}

Restart.CCUSUM <- function(x,side='upper', k,h,...)
{
  if(missing(k)) k = x$k
  if(missing(h)) h = x$h
  side = match.arg(tolower(side),c("upper","lower"))
  n = length(x$x)
  if(side=='upper'){
    stopifnot(!is.na(x$ubrk))
    x0 = x$x[x$ubrk:n]
    sx = rev(x$sH[1:x$ubrk])
    nx = which(sx<=0)[1]-1
    stopifnot(nx>0)
    delta = x$K+sx[1]/nx
    mu = x$mu[x$ubrk:n]
    mu0 = mu+delta
  }else{
    stopifnot(!is.na(x$lbrk))
    x0 = x$x[x$lbrk:n]
    sx = rev(x$sL[1:x$lbrk])
    nx = which(sx<=0)[1]-1
    stopifnot(nx>0)
    delta = x$K+sx[1]/nx
    mu = x$mu[x$lbrk:n]
    mu0 = mu-delta
  }
  cusum(x=x0, mu=mu0,k=k, h=h,...)
}


ccusum <- function(x,mu=NULL, k=0.5,h=4, ...)
{
  name <- deparse(substitute(x))
  if(is.null(mu)) mu = mean(x,na.rm=TRUE)
  n = length(x)
  if(length(mu)==1) mu = rep(mu,n)    
  sele = is.na(x) | is.na(mu)
  x = x[!sele]; mu = mu[!sele]
  n = length(x)  
  se.x = sd(x)
  sH = rep(0,n); sL = sH;
  K = k * se.x
  H = h * se.x
  for(i in 2:n){
    sH[i] = max(0, x[i] - (mu[i] + K) + sH[i-1]);
    sL[i] = max(0, (mu[i] - K) - x[i] + sL[i-1]);
  }
    ##  find the first breaks
  up1 = which(sH>H)[1]
  if(!is.na(up1)) ubrk = up1
  else ubrk = NA
  dn1 = which(sL>H)[1]
  if(!is.na(dn1)) lbrk = dn1
  else lbrk = NA

  structure(list(x = x, sH=sH, sL=sL, n=n,mu=mu,
                 mean = mean(mu,na.rm=TRUE),
                 K=K, H=H, k=k, h=h,
                 stdev = se.x, ubrk=ubrk,lbrk=lbrk,
                 call = match.call(), data.name = name),
            class='CCUSUM')
}

bcusum <- function(x,prob=NULL, R0=1.0,Ra=2, ...)
{
  name <- deparse(substitute(x))
  if(is.null(prob)) prob = mean(x,na.rm=TRUE) #c0
  if(length(prob)==1){ # c0 = prob can be given or missing
    sele = is.na(x)
    x = as.character(as.factor(x[!sele]));
    xlevel = levels(as.factor(x[!sele]));
    stopifnot(length(xlevel)==2)
    c0 = prob
    cA = Ra * prob
    n = length(x)  
    sH = rep(0,n);
    for(i in 2:n){
      Wt1 = log(cA/c0)  # if y==1 log(cA/c0)
      Wt0 = log((1-cA)/(1-c0))  # if y==0, log((1-cA)/(1-c0))
      Wt = ifelse(x[i]==xlevel[2], Wt1, Wt0)
      sH[i] = max(0, Wt + sH[i-1]);
    }
  }else{
    stopifnot(length(x)==length(prob))
    sele = is.na(x) | is.na(prob)
    x = as.character(as.factor(x[!sele]));
    xlevel = levels(as.factor(x[!sele]));
    stopifnot(length(xlevel)==2)
    prob = prob[!sele]
    n = length(x)  
    sH = rep(0,n);
    for(i in 2:n){
      pt = prob[i]
      Wt1 = log((1-pt+R0*pt)/(1-pt+Ra*pt))  # if y==1 log(cA/c0)
      Wt0 = log((1-pt+R0*pt)*Ra/(1-pt+Ra*pt)/R0)  # if y==0, log((1-cA)/(1-c0))
      Wt = ifelse(x[i]==xlevel[2], Wt1, Wt0)
      sH[i] = max(0, Wt + sH[i-1]);
    }
  }
  ##  find the first breaks
  if(Ra>1){
    H = 4.5
    up1 = which(sH>H)[1]
    if(!is.na(up1)) ubrk = up1
    else ubrk = NA
  }else{
    H = -4
    sH = -sH
    up1 = which(sH<H)[1]
    if(!is.na(up1)) ubrk = up1
    else ubrk = NA
  }  
  structure(list(x = x, sH=sH,  n=n,prob=prob,
                 H=H,  ubrk=ubrk,
                 call = match.call(), data.name = name),
            class='BCUSUM')
}

plot.BCUSUM <- function(x,main = NULL, xlab = NULL,
                        ylab=NULL,ylim=NULL,col=NULL,...){
  if (is.null(col))
    col = ifelse(x$H>0, 2, 4)
  if (is.null(xlab)) 
    xlab <- paste("Cases")
  if (is.null(main)) 
    main <- deparse(x$call)
  if (is.null(ylab)) 
    ylab <- paste("Binary Cumulative Sum of", x$data.name)
  if (is.null(ylim)){
    ymin = max(-2*x$H, min(x$sH))
    ymax = min(2*x$H, max(x$sH))
    ylim <- range(min(-x$H,ymin), max(ymax,x$H))
  }
  x0 = c(1: x$n)
  plot(x$sH~x0,type='n',ylim=ylim,
       main=main, xlab=xlab,ylab=ylab,...)
  abline(h=0,col='gray')
  lines(x$sH~x0,col='red',...)
  abline(h=x$H,lty=2, col='gray')
  if(!is.na(x$ubrk)){
    abline(v=x$ubrk,col='red',lty=3)
    text(x$ubrk,x$H, paste(x$ubrk),col='red',pos=2)
  }

  invisible(x)
}

plot.CCUSUM <- function(x,main = NULL, xlab = NULL,
                        ylab=NULL,ylim=NULL,col=NULL,...){
  if (is.null(col)){
    col1 = 2
    col2 = 4
  }
  if (is.null(xlab)) 
    xlab <- paste("Cases")
  if (is.null(main)) 
    main <- deparse(x$call)
  if (is.null(ylab)) 
    ylab <- paste("Cumulative Sum of", x$data.name)
  if (is.null(ylim)){
    ymin = max(-2*x$H, min(c(x$sH,-x$sL)))
    ymax = min(2*x$H, max(c(x$sH,-x$sL)))
    ylim <- range(min(-x$H,ymin), max(ymax,x$H))
  }
  x0 = c(1: x$n)
  plot(x$sH~x0,type='n',ylim=ylim,
       main=main, xlab=xlab,ylab=ylab,...)
  abline(h=0,col='gray')
  lines(x$sH~x0,col=col1,...)
  lines(-x$sL~x0,col=col2,...)
  abline(h=-x$H,lty=2, col='gray')
  abline(h=x$H,lty=2, col='gray')
  if(!is.na(x$ubrk)){
    abline(v=x$ubrk,col=col1,lty=3)
    text(x$ubrk,x$H, paste(x$ubrk),col=col1,pos=2)
  }
  if(!is.na(x$lbrk)){
    abline(v=x$lbrk,col=col2,lty=3)
    text(x$lbrk,-x$H, paste(x$lbrk),col=col2,pos=2)
  }

  invisible(x)
}

print.CCUSUM <- function(x,...){
  cat("\nUse ccusum for continous CUSUM\n")
}

print.BCUSUM <- function(x,...){
  cat("\nUse bcusum for binary CUSUM\n")
}
