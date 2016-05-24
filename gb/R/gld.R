## Check validity of GLD parameters
#.is.gld = function(lmd)
#  .Fortran(F_RIsGld, as.double(lmd),status=as.integer(0))$status>0

##########################################################################
###  GLD basic functions
##########################################################################
qrsgld = function(p,lambda){
  out = .Fortran(.F_RGldx,
    y=as.double(p), as.double(lambda),n=as.integer(length(p)))
  if(out$n<0) stop("Invalid RS-GLD parameters.")
  out$y
}

prsgld = function(x,lambda){
  out = .Fortran(.F_RGldFx,
    y=as.double(x), as.double(lambda),n=as.integer(length(x)))
  if(out$n<0) stop("Invalid RS-GLD parameters.")
  out$y
}

drsgld = function(x,lambda){
  out = .Fortran(.F_RGldfx,
    y=as.double(x), as.double(lambda),n=as.integer(length(x)))
  if(out$n<0) stop("Invalid RS-GLD parameters.")
  out$y
}

rrsgld = function(n, lambda){
  qrsgld(runif(n),lambda)
}

.lecount <- function(x,y) sum(y<=x)

.bincount <- function(x){
  x=x[!is.na(x)]; s = sd(x); mu = mean(x);
  n=length(x)
  x0 = seq(mu-5.7*s,mu+6*s,by=0.3*s)
  ps = apply(matrix(x0,ncol=1),1,FUN=.lecount,y=x)
  cts = diff(ps)
  bounds = NULL; counts=NULL;
  k = length(cts)
  for(i in 1:(k-1)){
    if(cts[i]>=5){
      bounds = c(bounds,x0[i])
      counts = c(counts,cts[i])
      ul = x0[i+1]
    }else{
      cts[i+1] = cts[i+1]+cts[i]
    }
  }
  m = n - sum(counts)
  if(m>5){
    bounds = c(bounds,ul)
    counts = c(counts,m)
  }else{
    counts[length(counts)] = counts[length(counts)]+m
  }
  list(x=bounds,y=counts);
}

.mom <- function(x){
  x = x[!is.na(x)]
  s = sd(x);
  m1 = mean(x);
  m2 = mean((x-m1)^2)
  m3 = mean((x-m1)^3)/s^3
  m4 = mean((x-m1)^4)/s^4
  c(m1,m2,m3,m4)
}

.mop <- function(x){
  pctls = quantile(x,c(.10,.25,.50,.75,.90))
  mp1 = pctls[3];
  mp2 = pctls[5]-pctls[1];
  mp3 = (pctls[3]-pctls[1])/(pctls[5]-pctls[3]);
  mp4 = (pctls[4]-pctls[2])/mp2;
  as.numeric(c(mp1,mp2,mp3,mp4))
}

#.lmom <- function(x){
#  .Fortran(.F_RLMoM, as.double(x), as.integer(length(x)),
#           para = as.double(rep(0,4)))$para
#}

.lmom = function(x){
# compute the sample first four L-moments
  x = sort(x);
  n = length(x);
  b0 = mean(x);
  w1 = (0:(n-1))/(n-1);
  b1 = mean(w1*x);
  s2 = 2:n;
  w2 = c(0,(s2-1)*(s2-2)/(n-1)/(n-2))
  b2 = mean(w2*x);
  s3 = 3:n;
  w3 = c(0,0,(s3-1)*(s3-2)*(s3-3)/(n-1)/(n-2)/(n-3));
  b3 = mean(w3*x);
  s4 = 4:n;
  w4 = c(0,0,0,(s4-1)*(s4-2)*(s4-3)*(s4-4)/(n-1)/(n-2)/(n-3)/(n-4));
  b4 = mean(w4*x);
  p10=1;p20=-1;p30=1;p40=-1;p50=1;
  p21=2; p31=-6;p41=12;p51=-20;
  p32=6;p42=-30;p52=90;
  p43 = 20;p53=-140;p54=70;
  l1=p10*b0;
  l2=p20*b0+p21*b1;
  l3=p30*b0+p31*b1+p32*b2;
  l4=p40*b0+p41*b1+p42*b2+p43*b3;
  return(c(l1,l2,l3,l4));
}

fit.gld <- function(x,method='LMoM'){
  name <- deparse(substitute(x))
  if(!is.numeric(x)) stop("\n\n\t 'x' must be numeric values!\n\n\n");
  x = x[!is.na(x)]; n = length(x);
  if(n < 6) stop("\n\n\t Sample size is too small!\n\n\n");
  if(n < 30) warning("\n Sample size might be too small to fit moment matching GLD!\n");
  s = sd(x); gpoints = seq(min(x)-s,max(x)+s,length=100);
  out = .bincount(x); os = out$y; nbins = length(out$y); xbin = out$x #
  if(nbins<3) stop("Data has less than 3 bins. Double-check the raw data.")

  method = match.arg(tolower(method), c("mom","mop","lmom"))
  out = switch(method,
    mom = .Fortran(.F_GLDMoM, b=as.double(.mom(x)), chisq = as.double(0.),
      as.integer(c(n,nbins)),as.double(os),as.double(xbin)),
    mop = .Fortran(.F_GLDMoP, b=as.double(.mop(x)), chisq = as.double(0.),
      as.integer(c(n,nbins)),as.double(os),as.double(xbin)),
    lmom = .Fortran(.F_GLDLMoM, b=as.double(.lmom(x)), chisq = as.double(0.),
      as.integer(c(n,nbins)),as.double(os),as.double(xbin))
    )

  DF = nbins - 5
  D2 = out$chisq
  pv = ifelse(DF>0,pchisq(D2,DF,lower.tail=FALSE),NA)

  y = drsgld(gpoints,out$b)
  sele = y>0
  return(structure(list(x = gpoints[sele],y = y[sele],n = n,para=out$b,
                        chisq = out$chisq, p.value=pv,method=method,
                        call = match.call(), data.name = name),
                   class = "gld"))
}

print.gld <- function (x, digits = NULL, ...) 
{
  tmp=NULL;k=length(x$para)
  for(i in 1:(k-1))
    tmp = paste(tmp,formatC(x$para[i]),",",sep='')
  tmp = paste(tmp,formatC(x$para[k]),sep='')
  cat("\nCall:\n\t", deparse(x$call), "\n\nData: ", x$data.name,
      " (", x$n, " obs.);", "\tMethod: = ", x$method,"\n",
      "Parameters: \t'para' = (",tmp, 
      ")", sep = "")
  cat("\nGoodness-of-fit: \n\tChi-Square = ",x$chisq, 
      ",\tp-value = ", format(x$p.value,6), "\n",sep = "")
  print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
  invisible(x)
}

plot.gld  <- 
function (x, main = NULL, xlab = NULL, ylab = "Density", type = "l", 
    zero.line = TRUE, ...) 
{
  tmp=NULL;k=length(x$para)
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

lines.gld  <- 
function (x, zero.line = TRUE, ...) 
{
  sele = x$y>0
  lines.default(x$x[sele],x$y[sele],  ...)
  invisible(NULL)
}
