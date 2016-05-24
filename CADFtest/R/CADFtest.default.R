CADFtest.default <- function(model, X=NULL, type=c("trend", "drift", "none"), 
                             data=list(), max.lag.y=1, min.lag.X=0, max.lag.X=0, dname=NULL,
                             criterion=c("none", "BIC", "AIC", "HQC", "MAIC"), ...)
{
# Author:       Claudio Lupi
# This version: July 22, 2009.
# This function computes Hansen's (1995) Covariate-Augmented Dickey-Fuller (CADF) test.
# The only required argument is y, the Tx1 time series to be tested (y can be a vector).
# If no time series of stationary covariates X is passed to the procedure, then an ordinary ADF test is performed.
# The test types are no-constant ("none"), constant ("drift"), constant plus trend ("trend", the default).
#
# max.lag.y >= 0
# min.lag.X <= 0
# max.lag.X >= 0

if (is.null(dname)){dname <- deparse(substitute(model))}

method <- "CADF test"
y <- model

if (is.null(X)) method <- "ADF test"

type      <- match.arg(type)

# modify argument type to be used by punitroot
switch(type,
		"trend" = urtype <- "ct",
		"drift" = urtype <- "c",
		"none"  =  urtype <- "nc")

criterion <- match.arg(criterion)

rho2 <- NULL # default value for rho^2
nX   <- 0    # default number of covariates. The exact number is computed below

if (is.ts(y)==FALSE) y <- ts(y)
trnd <- ts(1:length(y), start=start(y), frequency=frequency(y))

#############################################################################################################
#############################################################################################################
if (criterion=="none")  # no automatic model selection
{
test.results <- estmodel(y=y, X=X, trnd=trnd, type=type, 
			max.lag.y=max.lag.y, min.lag.X=min.lag.X, max.lag.X=max.lag.X, 
			dname=dname, criterion=criterion, obs.1=NULL, obs.T=NULL, ...)
}
#############################################################################################################
#############################################################################################################


#############################################################################################################
#############################################################################################################
if (criterion!="none")  # automatic model selection
{
  all.models <- expand.grid(max.lag.y:0, min.lag.X:0, max.lag.X:0)  # all possible models
  models.num <- dim(all.models)[1]                                  # number of models to be estimated 
  ICmatrix <- matrix(NA, models.num, 7)                             # matrix to store lag orders, & inf. crit.

  max.lag.y <- all.models[1, 1]
  min.lag.X <- all.models[1, 2]
  max.lag.X <- all.models[1, 3]
  interm.res <- estmodel(y=y, X=X, trnd=trnd, type=type, 
			  max.lag.y=max.lag.y, min.lag.X=min.lag.X, max.lag.X=max.lag.X, 
			  dname=dname, criterion=criterion, obs.1=NULL, obs.T=NULL, ...)
  ICmatrix[1, ] <- c(max.lag.y, min.lag.X, max.lag.X, interm.res$AIC, interm.res$BIC, interm.res$HQC, interm.res$MAIC)
  t.1 <- interm.res$est.model$index[1]
  t.T <- interm.res$est.model$index[length(interm.res$est.model$index)]

  for (modeln in 2:models.num)
  {
    max.lag.y <- all.models[modeln, 1]
    min.lag.X <- all.models[modeln, 2]
    max.lag.X <- all.models[modeln, 3]

	interm.res <- estmodel(y=y, X=X, trnd=trnd, type=type, 
				max.lag.y=max.lag.y, min.lag.X=min.lag.X, max.lag.X=max.lag.X, 
				dname=dname, criterion=criterion, obs.1=t.1, obs.T=t.T, ...)

	ICmatrix[modeln, ] <- c(max.lag.y, min.lag.X, max.lag.X, interm.res$AIC, interm.res$BIC, 
				interm.res$HQC, interm.res$MAIC)
  }
  
  if (criterion=="AIC")  selected.model <- which(ICmatrix[,4]==min(ICmatrix[,4]))
  if (criterion=="BIC")  selected.model <- which(ICmatrix[,5]==min(ICmatrix[,5]))
  if (criterion=="HQC")  selected.model <- which(ICmatrix[,6]==min(ICmatrix[,6]))
  if (criterion=="MAIC") selected.model <- which(ICmatrix[,7]==min(ICmatrix[,7]))
  
  if (length(selected.model) > 1) selected.model <- selected.model[length(selected.model)]

  max.lag.y <- ICmatrix[selected.model, 1]
  min.lag.X <- ICmatrix[selected.model, 2]
  max.lag.X <- ICmatrix[selected.model, 3]

  ################################## ESTIMATION & TEST WITH THE SELECTED MODEL ##############################

  test.results <- estmodel(y=y, X=X, trnd=trnd, type=type, 
			    max.lag.y=max.lag.y, min.lag.X=min.lag.X, max.lag.X=max.lag.X, 
			    dname=dname, criterion=criterion, obs.1=t.1, obs.T=t.T, ...)
}

class(test.results) <- c("CADFtest", "htest")
if (is.null(X)){names(test.results$statistic) <- paste("ADF(",max.lag.y,")",sep="")}
else{names(test.results$statistic) <- paste("CADF(",max.lag.y,",",max.lag.X,",",min.lag.X,")",sep="")}
test.results$estimate <- c("delta" = as.vector(test.results$est.model$coefficients[(2 - as.numeric(type=="none") +
			    as.numeric(type=="trend"))])) 
test.results$null.value <- c("delta" = 0) 
test.results$alternative <- "less" 
test.results$type <- type

return(test.results)
}

##############################################################################################################
##############################################################################################################
##############################################################################################################

estmodel <- function(y, X, trnd, type, max.lag.y, min.lag.X, max.lag.X, dname, criterion, obs.1, obs.T, ...)
{
  method <- "CADF test"
  if (is.null(X)) method <- "ADF test"
  rho2 <- NULL
  model <- "d(y) ~ "
  if (type=="trend") model <- paste(model, "trnd +", sep="")
  model <- paste(model, " L(y, 1)", sep="")

  if (max.lag.y > 0)
  {
    for (i in 1:max.lag.y) model <- paste(model, " + L(d(y), ",i,")", sep="")
  }

  if (is.null(X)==FALSE)
  {
    if (is.ts(X)==FALSE) X <- ts(X, start=start(y), frequency=frequency(y))
    nX <- 1; if (is.null(dim(X))==FALSE) nX <- dim(X)[2]  # number of covariates
    nX <- (max.lag.X - min.lag.X + 1)*nX              # number of X's (including the lags)
    if ((min.lag.X==0) & (max.lag.X==0)) model <- paste(model, " + L(X, 0)", sep="")
    if ((min.lag.X!=0) | (max.lag.X!=0))
    {
      for (i in min.lag.X:max.lag.X) model <- paste(model, " + L(X, ",i,")", sep="")
    }
  }

  if (type=="none") model <- paste(model, " -1", sep="")

  est.model      <- dynlm(formula=formula(model), start=obs.1, end=obs.T)
  summ.est.model <- summary(est.model)
  q              <- summ.est.model$df[1]
  TT             <- q + summ.est.model$df[2]

  sig2       <- sum(est.model$residuals^2)/TT
  lsig2      <- log(sig2)
  model.AIC  <- lsig2 + 2*q/TT
  model.BIC  <- lsig2 + q*log(TT)/TT
  model.HQC  <- lsig2 + 2*q*log(log(TT))/TT

  ytm1       <- est.model$model[, (2 + as.numeric(type=="trend"))]
  if (type=="drift") ytm1 <- ytm1 - mean(ytm1)
  if (type=="trend") 
    {
      dtrmod <- lsfit((1:TT), ytm1)
      ytm1   <- dtrmod$residuals
    }
  b0         <- est.model$coefficient[1 + as.numeric(type=="drift") + as.numeric(type=="trend")*2]
  sy2        <- sum(ytm1^2)
  tau        <- b0^2 * sy2 / sig2
  model.MAIC <- lsig2 + 2*(tau + q)/TT 

  t.value <- summ.est.model$coefficients[(2 - as.numeric(type=="none") + as.numeric(type=="trend")),3]

  if (is.null(X))
	{
		switch(type,
		"trend" = urtype <- "ct",
		"drift" = urtype <- "c",
		"none"  =  urtype <- "nc")
		p.value <- punitroot(t.value, N=TT, trend=urtype, statistic = "t") # MacKinnon p-values
	}

  if (is.null(X)==FALSE)
  {
    # Compute Hansen's p-value
    k <- length(est.model$coefficients)
    series  <- as.matrix(est.model$model)      # explanatory variables considered in the model (excluding constant)
    nseries <- dim(series)[2]                  # number of variables
    Xseries <- series[,(nseries-nX+1):nseries]     # the X's are the last nX columns of series
    if (nX==1) Xseries <- Xseries - mean(Xseries)  # demean the X's
    if (nX>1)  Xseries <- Xseries - apply(Xseries,2,mean)
    e <- as.matrix(est.model$residuals)
    if (nX==1) v <- Xseries * est.model$coefficients[k] + e
    if (nX>1)  v <- Xseries%*%est.model$coefficients[(k-nX+1):k] + e
    V <- cbind(e,v)
    mod <- lm(V~1)
    LRCM <- (kernHAC(mod, ...))*nrow(V)
    rho2 <- LRCM[1,2]^2/(LRCM[1,1]*LRCM[2,2])

    p.value <- CADFpvalues(t.value, rho2, type)
  }
    return(list(statistic=t.value,
           parameter=c("rho2" = rho2),
           method=method,
           p.value=as.vector(p.value),
           data.name=dname,
           max.lag.y=max.lag.y,
           min.lag.X=min.lag.X,
           max.lag.X=max.lag.X,
           AIC=model.AIC,
           BIC=model.BIC,
	   HQC=model.HQC,
	   MAIC=model.MAIC,
           est.model=est.model,
	   call=match.call(CADFtest)))
}
