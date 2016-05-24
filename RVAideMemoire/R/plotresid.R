# lme4 : getME, isLMM
# statmod : qresiduals

plotresid <- function(model,shapiro=FALSE) {
  res <- get.res(model)
  model.res <- res$residuals
  res.lab <- res$lab
  model.fit <- get.fit(model)
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  if (!inherits(model,"mlm")) {
    par(mfrow=c(1,2))
    plot(model.fit,model.res,xlab="Fitted values",ylab=res.lab,main=paste(res.lab,"vs fitted"))
    abline(h=0,col="grey",lty=3)
    panel.smooth(model.fit,model.res)
    qqPlot(model.res,lwd=1,grid=FALSE,xlab="Theoretical quantiles",ylab="Sample quantiles")
    if (shapiro) {
	shapiro.test(model.res)
    }
  } else {
    mqqnorm(model.res)
    if (shapiro) {
	mshapiro.test(model.res)
    }
  }
}

get.res <- function(x,...) {
  UseMethod("get.res")
}

get.res.default <- function(x,...) {stop("unknown model")}

get.res.lm <- function(x,...) {
  if (inherits(x,"mlm")) {get.res.mlm(x)} else
  if (inherits(x,"glm")) {get.res.glm(x)} else {
    list(residuals=rstudent(x),lab="Externally studentized residuals")
  }
}

get.res.mlm <- function(x,...) {
  list(residuals=resid(x),lab="")
}

get.res.glm <- function(x,...) {
  if (inherits(x,"negbin")) {get.res.negbin(x)} else {
    laws <- c("poisson","quasipoisson","binomial","quasibinomial")
    if (x$family[1] %in% laws) {
	list(residuals=statmod::qresiduals(x),lab="Quantile residuals")
    } else {
	list(residuals=rstudent(x),lab="Externally studentized residuals")
    }
  }
}

get.res.negbin <- function(x,...) {
  list(residuals=statmod::qresiduals(x),lab="Quantile residuals")
}

get.res.mer <- function(x,...) {
  stop(paste("for mixed models please update 'lmer' to version > 1.0 (actual: ",
    packageVersion("lme4"),")",sep=""))
}

get.res.glmmadmb <- function(x,...) {
  list(residuals=x$resid,lab="Residuals")
}

get.res.merMod <- function(x,...) {
  if (lme4::isLMM(x)) {
    list(residuals=residuals(x),lab="Residuals")
  } else {
    fam <- family(x)$family
    if (fam=="poisson") {
	y <- lme4::getME(x,"y")
	mu <- fitted(x)
	a <- ppois(y-1,mu)
	b <- ppois(y,mu)
	u <- runif(n=length(y),min=a,max=b)
	list(residuals=qnorm(u),lab="Quantile residuals")
    } else if (fam=="binomial") {
	p <- fitted(x)
	y <- lme4::getME(x,"y")
	mf <- model.frame(x)
	if ("(weights)" %in% colnames(mf)) { 
	  n <- mf$weights
	} else {
	  n <- rep(1,length(y))
	}
	y <- n*y
	a <- pbinom(y-1,n,p)
	b <- pbinom(y,n,p)
	u <- runif(n=length(y),min=a,max=b)
	list(residuals=qnorm(u),lab="Quantile residuals")
    } else if (grepl("Negative Binomial",fam)) {
	y <- lme4::getME(x,"y")
	size <- x@theta
	mu <- fitted(x)
	p <- size/(mu+size)
	a <- ifelse(y>0,pbeta(p,size,pmax(y,1)),0)
	b <- pbeta(p,size,y+1)
	u <- runif(n=length(y),min=a,max=b)
	list(residuals=qnorm(u),lab="Quantile residuals")
    } else {
	list(residuals=residuals(x),lab="Residuals")
    }
  }
}

get.res.lme <- get.res.nls <- get.res.gls <- function(x,...) {
  list(residuals=resid(x,type="pearson"),lab="Standardized residuals")
}

get.res.nlsList <- function(x,...) {
  list(residuals=resid(x,type="pooled"),lab="Standardized residuals")
}

get.res.survreg <- get.res.least.rect <- function(x,...) {
  list(residuals=residuals(x),lab="Residuals")
}


get.fit <- function(x,...) {
  UseMethod("get.fit")
}

get.fit.default <- function(x,...) {
  fitted(x)
}

get.fit.survreg <- function(x,...) {
  predict(x)
}
