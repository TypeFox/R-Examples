#.First.lib <- function(lib, pkg)
.onLoad <- function(lib, pkg)
{
  library.dynam("orsk", pkg, lib)
#  vers <- library(help=orsk)$info[[1]]
#  vers <- vers[grep("Version:",vers)]
#  vers <- rev(strsplit(vers," ")[[1]])[1]
#  packageStartupMessage("Loaded orsk",vers, appendLF = FALSE)
#  packageStartupMessage("Loaded orsk ", as.character(packageDescription("orsk")[["Version"]]),"\n")
}

orsk <-
  function(nctr, ntrt, a=NA, al=NA, au=NA, level=0.95, type="two-sided", method=c("grid","optim"), d=1e-4){
    call <- match.call()
    method <- match.arg(method)
    x <- nctr
    y <- ntrt
    if(x < 0 || y < 0) stop("x and y must be positive integers\n")
    if(type =="ci-only"){
    if(!is.na(a))
    warning("If type='ci=only', the point estimate of the odds ratio is not utilized\n")
    a <- NA
}
    if(type %in% c("lower", "two-sided","ci-only")){
      if(is.na(al))
      stop("missing al - lower bound of odds ratio confidence interval\n")
      else if(al < 0)
           stop("odds ratio or confidence intervals must be positive\n")
}
    if(type %in% c("upper", "two-sided", "ci-only")){
      if(is.na(au))
      stop("missing au - upper bound of odds ratio confidence interval\n")
      else if(au < 0)
           stop("odds ratio or confidence intervals must be positive\n")
}
    if(!is.na(a))
     if(a < 0) stop("odds ratio must be positive\n")
    if(level <= 0 || level >= 1) stop("confidence level must be between 0 and 1\n")
    if(method=="optim") d <- "NA   "
    else{
    if(d <= 0) stop("thereshold value must be positive\n")
    else if(d > 1) warning("thereshold value maybe too large\n")
}
    a1 <- a
    al1 <- al
    au1 <- au
    if(is.na(a)) a1 <- 0
    if(is.na(al)) al1 <- 0
    if(is.na(au)) au1 <- 0
    typef <- switch(type,
                 "lower" = 1,
                 "upper" = 2, 
                 "two-sided" = 3,
                 "ci-only" = 4)
    q <- qnorm(1-(1-level)/2)
    if(method=="grid"){
      z <- .Fortran("oddsratio",
                    x=as.integer(x),
                    y=as.integer(y),
                    a=as.double(a1),
                    al=as.double(al1),
                    au=as.double(au1),
                    q=as.double(q),
                    d=as.double(d),
                    t=as.integer(typef),
                    w=as.double(matrix(0, x-1, y-1)),
                    wind=as.integer(matrix(0, nrow=(x-1)*(y-1), ncol=2)),
                    package="orsk")
      w2 <- matrix(z$w, byrow=FALSE, nrow=x-1, ncol=y-1)
      windex <- matrix(z$wind, byrow=FALSE, nrow=(x-1)*(y-1), ncol=2)
      windex <- subset(windex, windex[,1] > 0)
      m <- dim(windex)[1]
      if(m == 0) 
        stop("Threshold value d=", d, " is too small\n")
### estimated odds ratio
      n01 <- windex[,1]
      n11 <- windex[,2]
      n10 <- y-n11
      n00 <- x-n01
      or.est <- n11*n00/(n10*n01)
      s <- sqrt(1/windex[,1] + 1/windex[,2] + 1/(x-windex[,1])+1/(y-windex[,2]))
      s <- exp(q*s)
      res.lo <- or.est/s
      res.up <- or.est*s
      rr <- (windex[,2]/y)/(windex[,1]/x)
                                        #standard error (SE) log relative risk
      slrr <- sqrt(1/windex[,2]-1/y+1/windex[,1]-1/x)
                                        #lower limit= the exponential of (log(Rel risk)-(1.96*SElogR))
                                        #upper limit= the exponential of (log(Rel risk)+(1.96*SElogR)) 
      rr.lo <- exp(log(rr) - q*slrr)
      rr.up <- exp(log(rr) + q*slrr)
      ctr.risk <- windex[,1]/nctr
      trt.risk <- windex[,2]/ntrt
      res <- cbind(ctr_yes=n01,ctr_no=n00, ctr_risk=ctr.risk, trt_yes=n11, trt_no=n10, trt_risk=trt.risk, OR=or.est, OR_lower=res.lo, OR_upper=res.up, RR=rr, RR_lower=rr.lo, RR_upper=rr.up, SS=w2[windex])
### sort R^2=w2
      k <- dim(res)[2]
      res <- res[order(res[,k]),]
      res <- as.data.frame(res)
    }
    else{
      fw <- function(p){
        n11 <- p[2]
        n01 <- p[1]
        n10 <- y - n11; n00 <- x - n01;
        v <- n11*n00/(n10*n01) ### compare to odds ratio
        f <- log(v) - log(a) ### compare to odds ratio
        s <- exp(q*sqrt(1/n11 + 1/n10 + 1/n01 + 1/n00))
        if(n11 <= 0 || n10 <= 0 || n01 <= 0 || n00 <= 0)
          print(c(n11, n10, n01, n00))  
        if(type!="upper")
        tmp1 <- log(v/s) - log(al)          ### compute confidence interval lower or upper bound
        if(type!="lower")
        tmp2 <- log(v*s) - log(au)         ### compute confidence interval lower or upper bound
        if(type=="lower")
          return(f^2 + tmp1^2)
        else if(type=="upper")
          return(f^2 + tmp2^2)
        else if(type=="two-sided")
          return(f^2 + tmp1^2 + tmp2^2)
        else if(type=="ci-only")
          return(tmp1^2 + tmp2^2)
      }

      tmp <- 1:max(2, ceiling(0.1*min(x-1, y-1))) # select number of random integer
      p0 <- cbind(sample(y-1)[tmp], sample(x-1)[tmp])
#      sink("/dev/null")
      ans <- multiStart(par=p0, fn=fw, lower=c(1,1), upper=c(x-1,y-1), action="optimize", control=list(trace=FALSE), quiet=TRUE, details=FALSE)
#      sink()
#      ans <- multiStartNew(par=p0, fn=fw, lower=c(1,1), upper=c(x-1,y-1), action="optimize", control=list(trace=FALSE), quiet=TRUE)
      pmat <- ans$par[ans$conv, ]
      pmat <- round(pmat)
      pmat <- pmat[!duplicated(pmat),]
      pmat <- matrix(pmat, ncol=2)
      pmat <- cbind(pmat, apply(pmat, 1, fw))
      ord1 <- order(pmat[, 3])
      ans <- pmat[ord1, ]
      ans <- matrix(ans, ncol=3)
      n11 <- ans[,2]; n01 <- ans[,1];
      n10 <- y - n11; n00 <- x - n01;
      or.est <- n11*n00/(n10*n01) ### compare to odds ratio
      s <- exp(q*sqrt (1/n11 + 1/n10 + 1/n01 + 1/n00)) 
      OR.lower <- or.est / s 
      OR.upper <- or.est * s
      rr <- (n11/y)/(n01/x)
      slrr <- sqrt(1/n11-1/y+1/n01-1/x)
                                        #lower limit= the exponential of (log(Rel risk)-(1.96*SElogR))
                                        #upper limit= the exponential of (log(Rel risk)+(1.96*SElogR)) 
      rr.lo <- exp(log(rr) - q*slrr)
      rr.up <- exp(log(rr) + q*slrr)
      ctr.risk <- n01/nctr
      trt.risk <- n11/ntrt
      res <- cbind(ctr_yes=n01,ctr_no=n00, ctr_risk=ctr.risk, trt_yes=n11, trt_no=n10, trt_risk=trt.risk, OR=or.est, OR_lower=OR.lower, OR_upper=OR.upper, RR=rr, RR_lower=rr.lo, RR_upper=rr.up, SS=ans[,3])
      res <- as.data.frame(res)
    }
#      colnames(res) <- c("ctr_yes","ctr_no","ctr_risk","trt_yes","trt_no","trt_risk", "OR","OR_lower","OR_upper", "RR", "RR_lower","RR_upper", "SS") 
#    if(sort){
    RET <- list(x=x, y=y, a=a, al=al, au=au, type=type, method=method, d=d, res=res) 
    RET$call <- call
    class(RET) <- "orsk"
    return(RET)
}

plot.orsk <- function(x, type=c("RR", "OR"), digits=2, factor=1, amount=NULL, ...){
type <- match.arg(type)
if(type=="OR"){
xla <- "Risk in the control group"
yla <- "Risk in the treatment group"
}
else {
xla <- "Relative risk"
}
maxdist <- pmax(abs(x$res$OR - x$a),
abs(x$res$OR_lower - x$al),
abs(x$res$OR_upper - x$au))
if(type=="OR"){
res1 <- x$res$ctr_risk[maxdist <= 10^-digits/2]
res2 <- x$res$trt_risk[maxdist <= 10^-digits/2]
if(!is.null(res1) && !is.null(res2))
plot(jitter(res1, factor=factor, amount=amount), jitter(res2, factor=factor, amount=amount), xlab=xla, ylab=yla)
else return(warnings("no results for the selected digits\n"))
}
else{
res1 <- x$res$RR[maxdist <= 10^-digits/2]
if(!is.null(res1))
dotPlot(res1, xlab=xla)
else return(warnings("no results for the selected digits\n"))
}
}
print.orsk <- function(x, ...) {
  if(class(x) != "orsk")
  stop("Not an object of class orsk\n")
  cat("\n")
  cat("\t Converting odds ratio to relative risk\n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Method: ", x$method, "\n")
  cat("Threshold value: ", x$d, "\n")
  cat("\n")
  invisible(x)
}

summary.orsk <- function(object, nlist=1:5, ...) {
  x <- object
  if(class(x) != "orsk")
  stop("Not an object of class orsk\n")
  cat("\n")
  cat("\t Converting odds ratio to relative risk\n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("type: ", x$type, "              method: ", x$method, "\n")
  cat("threshold value: ", x$d, "\n")
  cat(paste("The odds ratio utilized: ",x$a, ", confidence interval utilized: ", x$al, "-", x$au, "\n", sep=""))
  cat("\nThe following odds ratios and relative risks are for the scenarios created \nwith different numbers of events in control and treatment group that lead \nto comparable results for the above odds ratio and confidence interval\n")
  print(x$res[nlist,], digits=3)
  cat("\n")
  invisible(x)
}

