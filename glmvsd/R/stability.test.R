stability.test <- function(x, y, method=c("seq", "bs", "perturb"), 
                           penalty = c("LASSO", "SCAD", "MCP"), nrep = 50, 
                           remove = 0.2, tau = 0.5, nfolds = 5, 
					  family=c("gaussian","binomial")) {
  # model check
  method <- match.arg(method)
  penalty <- match.arg(penalty)
  family <- match.arg(family)
  y <- drop(y)
  y <- as.numeric(y)
  x <- as.matrix(x)
  p <- NCOL(x)
  n <- length(y)
  if (n != NROW(x)) 
    stop("x and y have different number of observations")
  if (nrep < 2) 
    stop("The number of repititions must be greater than 1.")
  if (remove<=0 || remove>=1)  
    stop("The proportion of data to be removed for sequential stability test must be in (0,1).")
  if (tau<=0 || tau>1) 
    stop("The perturbation size (tau) must be in (0,1].")
  # fit
  full <- modelfit(x, y, nfolds, penalty, family)
  fitted <- drop(as.matrix(cbind(rep(1,n),x)%*%full$coefit))
  sigmafit<-sqrt(sum((y-fitted)^2)/(n-sum(full$modelfit)))
  ####################sequential############
  model.sub<-matrix(NA, nrep, p)
  if(method == "seq"){
    m<-floor(n*(1-remove))
    for (i in 1:nrep){
      a <- sample(n,m,replace=F)
      model.sub[i, ] <- modelfit(x=x[a,],y=y[a],nfolds,penalty)$modelfit
    }
  }
  ####################Bootstrap####################
  if(method == "bs"){
    for (i in 1:nrep){
	    if(family=="gaussian") y.star <- rnorm(n,fitted,sigmafit)
	    if(family=="binomial") y.star <- rbinom(n, 1, exp(fitted)/(1+exp(fitted)))
      model.sub[i, ] <- modelfit(x=x,y=y.star,nfolds,penalty)$modelfit
    }
  }
  ####################Perturbation####################
  if(family=="gaussian" && method == "perturb"){
    for (i in 1:nrep){
      y.star<-rnorm(n,y,sqrt(tau)*sigmafit)
      model.sub[i, ] <- modelfit(x=x,y=y.star,nfolds,penalty)$modelfit
    }
  }
  # results collection
  out <- sum(abs(sweep(model.sub,MARGIN=2,full$modelfit,"-")))/nrep
  out
}
