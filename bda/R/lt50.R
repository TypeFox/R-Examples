.logit <- function(x){
  if(any(x<0))stop("Negative 'x' value")
  if(any(x>1))stop("Invalid probability (percentage/rate)")
  x[x==1] = 0.999999999999
  x[x==0] = 0.000000000001
  log(x/(1-x))
}

ld50.logitfit <- function(rate, dose, p = 0.5) {
  if(any(rate<0||rate>1))
    stop("Percentage can not be negative or greater than 1.")
  if(any(dose<0))
    stop("Dosage can not be negative.")
  rate0 = rate;  dose0 = dose
  dose[dose==0] = 0.00000000000001
  ldose = log(dose)
  ldose2 = ldose*ldose
  ldose3 = ldose2*ldose
  ldose5 = ldose2*ldose3
  lrate = .logit(rate)
  
  obj = lm(lrate~dose+ldose)
  b = coef(obj)
  x0 = seq(0.00001,max(dose0),length=100)
  lx0 = log(x0)
  lx02 = lx0*lx0
  lx03 = lx02*lx0
  lx05 = lx02*lx03
  y0 = 1/(1+exp(-b[1L]-b[2L]*x0-b[3L]*lx0))
  if(any(p>=1|p<=0))stop("Invalid 'p' in 'LC(p)'.")
  x.p <- .Fortran(.F_nrlogit, as.double(median(dose0)),
                  as.double(b),p=as.double(p),as.integer(length(p)))$p
  SE <- NA
  res <- structure(x.p, SE = SE, p = p,
                   dose=dose0,rate=rate0,
                   x=x0, y = y0)
  class(res) <- "glm.dose" 
  res
}

ld50.logit <- function(ndead, ntotal, dose, cf = 1:2, p = 0.5) {
  if(any(ndead>ntotal)||ntotal<1||any(ndead<0))stop("Wrong count(s)!")
  if(any(floor(ndead)!=ndead))stop("Counts must be integers")
  nalive = ntotal - ndead
  SF <- cbind(nalive,ndead)
  obj <- glm(SF ~ dose, family = binomial)
  if(obj$df.residual>obj$deviance){
    warning("Check for overdispersion problem!")
    print(paste("Residual deviance=",obj$deviance,
                ", Residual DF=",obj$df.residual))
  }
  eta <- family(obj)$linkfun(p)
  b <- coef(obj)[cf]
  x0 = seq(0,max(dose),length=100)
  y0 = 1/(1+exp(-b[1L]-b[2L]*x0))

  x.p <- (eta - b[1L])/b[2L]
  names(x.p) <- paste("p = ", format(p), ":", sep = "")
  pd <-  -cbind(1, x.p)/b[2L]
  SE <- sqrt(((pd %*% vcov(obj)[cf, cf]) * pd) %*% c(1, 1))
  res <- structure(x.p, SE = SE, p = p,
                   dose=dose,rate=nalive/ntotal,
                   x=x0, y = y0)
  class(res) <- "glm.dose" 
  res
}

print.glm.dose <- function(x, ...)
{
  M <- cbind(x, attr(x, "SE"))
  dimnames(M) <- list(names(x), c("Dose", "SE"))
  x <- M
  NextMethod("print")
}

plot.glm.dose <- function(x, lty=1,lwd=1,...)
{
  y = attr(x,"dose")
  z = attr(x,"rate")
  y0 = attr(x,"y")
  x0 = attr(x,"x")
  plot(y,z,ylab="Percent viability",ylim=c(0,1),
       xlab="Dose", xlim=c(0,max(y)))
  lines(x0,y0,new=TRUE,lty=lty,lwd=lwd)
  abline(h=0.5,lty=3)
}

lines.glm.dose <- function(x, lty=1,lwd=1,pch=pch,...)
{
  y = attr(x,"dose")
  z = attr(x,"rate")
  y0 = attr(x,"y")
  x0 = attr(x,"x")
  points(y,z,pch=pch)
  lines(x0,y0,lty=lty,lwd=lwd)
}
