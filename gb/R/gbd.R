## Check validity of GLD parameters
.is.egld = function(lmd){
  out = TRUE
  if(is.null(lmd)) out=FALSE
  if(any(is.na(lmd))) out=FALSE
  if(length(lmd)!=4) out=FALSE
  if(lmd[2]<=0||lmd[3]<=0||lmd[4]<=0) out=FALSE
  return(out)
}
##########################################################################
###  EGLD basic functions
##########################################################################
pegld = function(x,lambda){
  if(!.is.egld(lambda)) stop("Invalid EGLD parameters.")
  pbeta((x-lambda[1])/lambda[2],lambda[3],lambda[4])
}

degld = function(x,lambda){
  if(!.is.egld(lambda)) stop("Invalid EGLD parameters.")
  xt = (x-lambda[1])/lambda[2]
  dbeta(xt,lambda[3],lambda[4])/lambda[2]
}

qegld = function(p,lambda){
  ## x is a probability
  if(!.is.egld(lambda)) stop("Invalid EGLD parameters.")
  qbeta(p,lambda[3],lambda[4])*lambda[2]+lambda[1]
}

regld = function(n, lambda){
  if(!.is.egld(lambda)) stop("Invalid EGLD parameters.")
  lambda[1]+lambda[2]*rbeta(n,lambda[3],lambda[4])
}

fit.egld = function(x,xmin=NULL,xmax=NULL)
{
  name <- deparse(substitute(x))
  if(!is.numeric(x)) stop("\n\n\t 'x' must be numeric values!\n\n\n");
  x = x[!is.na(x)]; n = length(x);
  if(n <6) stop("\n\n\t Sample size is too small!\n\n\n");
  if(n < 30) warning("\n Sample size might be too small to fit moment matching GLD!\n");
  s = sd(x); gpoints = seq(min(x)-s,max(x)+s,length=100);
  out = .bincount(x); os = out$y; nbins = length(out$y); xbin = out$x #
  if(nbins<3) stop("Data has less than 3 bins. Double-check the raw data.")

  if(is.null(xmin)){
    if(is.null(xmax)){
      type=0;xmin=min(x);xmax=max(x);
    }else{
      type=1;xmin=min(x);
    }
  }else{
    if(is.null(xmax)){
      type=2;xmax=max(x);
    }else{
      type=3;
    }
  }
  y = (x-xmin)/(xmax-xmin)
  s2 = var(y); mu=mean(y)
  b3 = (1./mu-1.-s2/mu/mu)*mu^3./s2;
  b4 = (1./mu-1.)*b3;
  
  p = c(xmin,xmax,b3,b4)
  out = .Fortran(.F_FitGBD, as.double(x), as.integer(length(x)),
    as.integer(type), llk=as.double(0),p=as.double(p));
  betas = out$p
  gpoints = seq(betas[1],betas[1]+betas[2],length=100);
  y = degld(gpoints,betas)

  DF = nbins - 5
  Fx = pegld(xbin[-1], betas); Fx = c(Fx[1], diff(Fx))
  es = n*Fx; es = c(es,n-sum(es))
  D2 = sum((os-es)^2/es)
  pv = ifelse(DF>0,pchisq(D2,DF,lower.tail=FALSE),NA)

  
  return(structure(list(x = gpoints,y = y,n = length(x),para=betas,
                        chisq = D2, p.value=pv, method="ml",
                        call = match.call(), data.name = name),
                   class = "egld"))
}

print.egld <- function (x, digits = NULL, ...) 
{
  tmp="";k=length(x$para)
  for(i in 1:(k-1))
    tmp = paste(tmp,formatC(x$para[i]),",",sep='')
  tmp = paste(tmp,formatC(x$para[k]),sep='')
  cat("\nCall:\n\t", deparse(x$call), "\n\nData: ", x$data.name, 
      " (", x$n, " obs.);", "\tMethod: = ", x$method,"\n",
      "Parameters: \t'para' = (",tmp, 
      ")\n\n", sep = "")
  cat("\nGoodness-of-fit: \n\tChi-Square = ",x$chisq, 
      ",\tp-value = ", format(x$p.value,6), "\n",sep = "")
  print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
  invisible(x)
}

plot.egld  <- 
function (x, main = NULL, xlab = NULL, ylab = "Density", type = "l", 
    zero.line = TRUE, ...) 
{
  tmp="";k=length(x$para)
  for(i in 1:(k-1))
    tmp = paste(tmp,formatC(x$para[i]),",",sep='')
  tmp = paste(tmp,formatC(x$para[k]),sep='')
  if (is.null(xlab)) 
    xlab <- paste("N =", x$n, "  Parameters =(",tmp,")")
  if (is.null(main)) 
    main <- deparse(x$call)
  sele = x$y>0
  plot.default(x$x[sele],x$y[sele], type = type, 
               main = main, xlab = xlab, ylab = ylab, ...)
  if (zero.line) 
    abline(h = 0, lwd = 0.1, col = "gray")
  invisible(NULL)
}

