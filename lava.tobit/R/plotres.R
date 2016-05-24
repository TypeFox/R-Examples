##'Plot distribution of standardized residuals
##'
##'Plot empirical (KM) and model-specific cumulative distribution function of
##'standardized residuals
##'
##'
##'@param x Model, \code{lvmfit} object
##'@param var Character vector of (endogenous) variable names
##'@param ylab Label of x-axis
##'@param xlab Label of y-axis
##'@param main Title of plot
##'@param \dots Additional argument
##'@author Klaus K. Holst
##'@keywords models
##'@examples
##'
##'\dontrun{
##'
##'## Simulate data where (y01,y2)
##'## follows conditional bivariate normal distribution
##'## given covariate x. Instead of y01 we observe
##'## right censored version y2
##'n <- 200
##'m <- lvm(c(y01,y2) ~ x)
##'covariance(m) <- y01~y2
##'set.seed(1)
##'d <- sim(m,n)
##'d$cens1 <- rexp(n)
##'d$status1 <- with(d,y01<cens1)
##'d$y1 <- with(d, pmin(y01,cens1))
##'
##'## Estimate model parameters
##'d$S1 <- with(d, Surv(y1,status1))
##'m <- lvm(c(S1,y2)~x); covariance(m) <- S1~y2
##'e <- estimate(m,d,control=list(trace=1))
##'
##'## Plot cumulative distribution functions
##'par(mfrow=c(2,2)); plotres(e); plot(e)
##'}
##'
##' @export
plotres <- function(x,var=endogenous(x),
                    ylab="Cumulative Distribution Function",
                    xlab="Standardized residuals",
                    main,
                    ...) {
  require(survival)
  r <- residuals(x,std=TRUE)
  W <- Weight(x)
  
  for (v in var) {    
    if (v %in% colnames(W)) {
      S <- Surv(ifelse(W[,v]==-1,NA,r[,v]),
                ifelse(W[,v]==1,NA,r[,v]),
                type="interval2")      
    } else {
      S <- Surv(r[,v],rep(TRUE,length(r[,v])))
    }
    g <- survfit(S~1)
    mymain <- ifelse(!missing(main),main,v)    
    with(g,plot(1-surv~time,type="s",main=mymain,xlab=xlab,ylab=ylab))
    with(g,lines(1-upper~time,type="s",lty=2))
    with(g,lines(1-lower~time,type="s",lty=2))
    ro <- sort(r[,v]); 
    lines(ro,pnorm(ro),col="red",xlab=xlab,ylab=ylab)
    
   }
  invisible(x)
}
