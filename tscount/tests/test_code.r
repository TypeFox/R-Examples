#library(tscount)

begin <- Sys.time()

checkfit <- function(seed, n=50, model, param, xreg=NULL, distr="poisson", distrcoefs=NULL, link="identity", startestims=FALSE, extended=FALSE, interv=FALSE, ...){
  if(missing(seed)) seed <- 1945
  set.seed(seed)
  #Simulate a time series from the given model:
  timser <- tsglm.sim(n=n, param=param, model=model, xreg=xreg, link=link, distr=distr, distrcoefs=distrcoefs)$ts
  print(unlist(param))
  #Fit the given model to the time series:
  print(fit <- tsglm(ts=timser, model=model, xreg=xreg, link=link, distr=distr, ...))
  if(startestims){ #Try the different methods for start estimation:
    tsglm(ts=timser, model=model, xreg=xreg, link=link, distr=distr, start.control=list(method="iid"), final.control=NULL)
    tsglm(ts=timser, model=model, xreg=xreg, link=link, distr=distr, start.control=list(method="GLM"), final.control=NULL)
    tsglm(ts=timser, model=model, xreg=xreg, link=link, distr=distr, start.control=list(method="CSS"), final.control=NULL)
    tsglm(ts=timser, model=model, xreg=xreg, link=link, distr=distr, start.control=list(method="MM"), final.control=NULL)
  }
  if(extended){ #Apply the standard methods for the fitted model:
    summary(fit)
    residuals(fit, type="pearson")
    residuals(fit, type="anscombe")
    par(mfrow=c(3,2))
      plot(fit, ask=FALSE)
    fitted(fit)
    coef(fit)
    predict(fit, n.ahead=4)
    #logLik(fit) #already included in summary-method
    vcov(fit)
    #AIC(fit) #already included in summary-method
    #BIC(fit) #already included in summary-method
    se(fit, B=3)
    #pit(fit) #already included in plot-method
    #marcal(fit) #already included in plot-method
    scoring(fit)
  }
  if(interv){
    interv_test(fit, tau=floor(n/2), delta=0.8, external=FALSE, est_interv=TRUE, ...)
    interv_detect(fit, taus=floor(0.45*n):ceiling(0.55*n), delta=0.8, B=3, ...)
    interv_multiple(fit, taus=floor(0.45*n):ceiling(0.55*n), deltas=c(0,1), B=3, final.control_bootstrap=NULL, ...)
  }
  return(TRUE)
}

checkfit(model=list(past_obs=1, past_mean=1), param=list(intercept=2, past_obs=0.4, past_mean=0.3), startestims=TRUE, extended=TRUE)
checkfit(model=list(past_obs=1, past_mean=1), param=list(intercept=2, past_obs=0.4, past_mean=0.3, xreg=c(4)), xreg=matrix(rexp(1*50, rate=3), ncol=1), startestims=TRUE, extended=TRUE) #one covariate
checkfit(model=list(past_obs=1, past_mean=1), param=list(intercept=2, past_obs=0.4, past_mean=0.3, xreg=c(4,2)), xreg=matrix(rexp(2*50, rate=3), ncol=2), startestims=TRUE, extended=TRUE) #two covariates
checkfit(model=list(past_obs=1), param=list(intercept=2, past_obs=0.4))
checkfit(model=list(), param=list(intercept=2))
checkfit(model=list(past_obs=1:2, past_mean=1:2), param=list(intercept=2, past_obs=c(0.3,0.2), past_mean=c(0.2,0.1)))
checkfit(model=list(past_obs=1, past_mean=1), param=list(intercept=2, past_obs=0.4, past_mean=0.3), distr="nbinom", distrcoefs=c(size=2), extended=TRUE, start.control=list(use=20))
checkfit(model=list(past_obs=1, past_mean=1), param=list(intercept=2, past_obs=0.4, past_mean=0.3), distr="nbinom", distrcoefs=c(size=2), extended=TRUE, init.drop=TRUE)
checkfit(model=list(past_obs=1, past_mean=1), param=list(intercept=0.5, past_obs=0.4, past_mean=0.3), link="log", startestims=TRUE, extended=TRUE)
checkfit(model=list(past_obs=1, past_mean=1), param=list(intercept=0.5, past_obs=0.4, past_mean=0.3), link="log", distr="nbinom", distrcoefs=c(size=2), extended=TRUE, interv=TRUE)
checkfit(model=list(past_obs=1, past_mean=1, external=FALSE), param=list(intercept=2, past_obs=0.4, past_mean=0.3, xreg=2), xreg=cbind(lintrend=seq(0, 1, length.out=50)))
checkfit(model=list(past_obs=1, past_mean=1, external=TRUE), param=list(intercept=2, past_obs=0.4, past_mean=0.3, xreg=2), xreg=cbind(lintrend=seq(0, 1, length.out=50)))

#Functions for analytical mean, variance and autocorrelation:
ingarch.mean(intercept=2, past_obs=c(0.1,0.1), past_mean=0.1)
ingarch.acf(intercept=2, past_obs=0.1,  past_mean=0.1, type="acf", lag.max=5, plot=FALSE)
ingarch.var(intercept=2, past_obs=c(0.3,0.1))

Sys.time() - begin
