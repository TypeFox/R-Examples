HOF.model <- function (
		occ, 
		grad, 
		M = max(occ), 
		y.name, 
		family = binomial, 
		lim = 100, 
		x.name, 
		...)  {
  if(max(occ) > M) stop('Maximum response value specified too low!')
  if(any(is.na(grad))) stop('No NA values in gradient allowed.')
  aic <- family()$aic
  dev.resids <- family()$dev.resids
  famname <- family()$family
  if (!(famname %in% c("binomial", "gaussian", "poisson"))) stop("Allowed families: binomial, gaussian, poisson")
  if(famname == "poisson" & sum(grep('.',as.character(occ), fixed=TRUE)) > 0)  stop('Occurrency data must be integer values for family poisson!')
  if(famname == "binomial" && !all(occ %in% c(0,1))) warning('Occurrency data should be 0 or 1 for family binomial!') 
  options(warn=-1)
#  if (getRversion() >= '2.15.1') globalVariables(c('II.res', 'IV.res', 'p', 'I.res','III.res', 'VII.res'), package='eHOF')
  
  div <- if(famname == "binomial") M else 1
  if(!exists('wt')) wt  <- if(famname == "binomial") M else 1
  if (length(wt) == 1) wt <- rep(wt, length(grad))
  x.orig <- grad
  x <- scale01(grad)
  nobs <- length(occ)

  mlHOF <- function (p, x, y, model, M,  ...) {
     mu <- HOF.fun(x, model, p, M)
     n <- wt
     if(famname == "binomial") {
          y <- y/M
          mu <- mu/M
      }
     dev <- sum(dev.resids(y, mu, wt))
     aic(y, n, mu, wt, dev)/2
  }

  trial.1 <- function(p, x, m, ...) { # PORT routine
    p <- svHOF(x, occ, M, model = m, mod = mod)
    temp <- nlminb(start=p, objective = mlHOF, lower = -lim, upper = lim, x = x, y = occ, M = M, model = m)
    temp$method <- 'nlminb'
    temp
    }
  trial.2 <- function(p, x, m, ...) { # Byrd et. al. (1995)
    p <- svHOF(x, occ, M, model = m, mod = mod)                                
    temp <- optim(par = p, mlHOF, x = x, method = "L-BFGS-B", lower = -lim, upper = lim, y = occ, M = M, model = m)
    temp$method <- 'Nelder-Mead'
    return(temp)
    }

  I.res <- NA; II.res <- NA; III.res <- NA;  IV.res <- NA; VII.res <- NA; p <- NA
## Optimization
for(m in eHOF.modelnames) {
  res <- paste(m, 'res',sep='.')
  mod <- switch(m,
    III  = II.res,
    V    = IV.res,
    VI   = IV.res,
    VII  = IV.res,            
    NA
    )
  try(assign(res, trial.1(p, x, m)), silent=TRUE)

  if(is.na(res)) try(
    { tmp <- trial.2(p, x, m)
    if(all(abs(tmp$par) <= lim)) assign(res, tmp) }, silent =TRUE)

# fitted values and deviance
  if(is.na(res)) {
  	assign(res, list(par=c(a=NA,b=NA,c=NA,d=NA), value=NA, counts=c(0,0),convergence=NULL, method='None', trial='No success', message=paste('No solution for model',m), fitted=rep(NA,length(x)),deviance=NA))
  	warning(paste('No solution for model', m)) } else {

  if(all(get(res)$par <= lim) & !any(is.na(get(res)$par))) {
      fv <- HOF.fun(x, m, get(res)$par, M)
      assign(res, c(get(res), deviance = sum(dev.resids(occ/div, fv/div, wt))))
      assign(res, c(get(res), fitted = if(is.nan(get(res)[['deviance']])) list(rep(NA, length(x))) else list(fv) ))
      } else {
  	assign(res, list(par=c(a=NA,b=NA,c=NA,d=NA), value=NA, counts=c(0,0),convergence=NULL, method='None', trial='No success', message=paste('No solution for model',m), fitted=rep(NA,length(x)), deviance=NA))
  	warning(paste('No solution for model',m)) }
   }
  }
## end of model loop

if(any(is.na(VI.res$fitted))) {
   VI.res$message <- 'Can not minimize bimodal model'
   VI.res$deviance <- NA
   }
if(!all(is.na(V.res$par))) 
  if(V.res$par[2] * V.res$par[4] < 0) {
    V.res$message <- 'response continuous, rejected'
    V.res$deviance <- NA
  }

  models <- list(I = I.res, II = II.res, III = III.res , IV = IV.res, V = V.res, VI = VI.res, VII = VII.res)
  out <- list(call = match.call(), x = x.orig, y = occ,  range = range(grad), M = M, family = famname, nobs = nobs, models = models)
  class(out) <- "HOF"
  options(warn=0)
  out
}
