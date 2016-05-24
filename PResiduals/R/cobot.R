ordinal.scores.logit <- function(y, X) {
  ## y is a numeric vector
  ## X is a vector or matrix with one or more columns.

  ## Ensure code works if called from somewhere else (not COBOT.scores()).
  ## Make X a matrix if it is a vector.  This makes later coding consistent.
  if(!is.matrix(X)) X = matrix(X, ncol=1)

  ## N: number of subjects; ny: number of y categories
  N = length(y)
  ny = length(table(y))

  ## na, nb: number of parameters in alpha and beta
  na = ny - 1
  nb = ncol(X)
  npar = na + nb
  ## Z is defined as in McCullagh (1980) JRSSB 42:109-142
  Z = outer(y, 1:ny, "<=")

  ## Fit proportional odds model and obtain the MLEs of parameters.
  mod = lrm(y ~ X, tol=1e-50, maxit=100)
  alpha = -mod$coeff[1:na]
  beta = -mod$coeff[-(1:na)]

  ## Scores are stored for individuals.
  dl.dtheta = matrix(NA, N, npar)
  ## Information matrices are stored as sums over all individuals.
  d2l.dalpha.dalpha = matrix(0,na,na)
  d2l.dalpha.dbeta = matrix(0,na,nb)
  d2l.dbeta.dbeta = matrix(0,nb,nb)
  d2l.dbeta.dalpha = matrix(0,nb,na)
  ## Predicted probabilities p0 and dp0.dtheta are stored for individuals.
  p0 = matrix(,N,ny)
  dp0.dtheta = array(0,c(N,ny,npar))
  ## Cumulative probabilities
  Gamma = matrix(0,N,na)
  dgamma.dtheta = array(0,c(N,na,npar))

  for (i in 1:N) {
    z = Z[i,]  ## z has length ny
    x = X[i,]

    ## gamma and phi are defined as in McCullagh (1980)
    gamma = 1 - 1/(1 + exp(alpha + sum(beta*x))) ## gamma has length na
    diffgamma = diff(c(gamma,1))
    invgamma = 1/gamma
    invgamma2 = invgamma^2
    invdiffgamma = 1/diffgamma
    invdiffgamma2 = invdiffgamma^2
    phi = log(gamma / diffgamma) ## phi has length na
    Gamma[i,] = gamma


#### Some intermediate derivatives
    ## g(phi) = log(1+exp(phi))
    dg.dphi = 1 - 1/(1 + exp(phi))
    ## l is the log likelihood (6.3) in McCullagh (1980)
    dl.dphi = z[-ny] - z[-1] * dg.dphi
    t.dl.dphi = t(dl.dphi)

    ## dphi.dgamma is a na*na matrix with rows indexed by phi
    ## and columns indexed by gamma
    dphi.dgamma = matrix(0,na,na)
    diag(dphi.dgamma) = invgamma + invdiffgamma
    if(na > 1)
      dphi.dgamma[cbind(1:(na-1), 2:na)] = -invdiffgamma[-na]

    dgamma.base = gamma * (1-gamma)
    dgamma.dalpha = diagn(dgamma.base)
    dgamma.dbeta = dgamma.base %o% x
    dgamma.dtheta[i,,] = cbind(dgamma.dalpha, dgamma.dbeta)

    d2gamma.base = gamma * (1-gamma) * (1-2*gamma)

    ##
    d2l.dphi.dphi = diagn(-z[-1] * dg.dphi * (1-dg.dphi))
    d2l.dphi.dalpha = d2l.dphi.dphi %*% dphi.dgamma %*% dgamma.dalpha
    d2l.dphi.dbeta = d2l.dphi.dphi %*% dphi.dgamma %*% dgamma.dbeta

    ##
    d2phi.dgamma.dalpha = array(0,c(na,na,na))
    d2phi.dgamma.dalpha[cbind(1:na,1:na,1:na)] = (-invgamma2 + invdiffgamma2) * dgamma.base
    if(na > 1) {
      d2phi.dgamma.dalpha[cbind(1:(na-1),1:(na-1),2:na)] = -invdiffgamma2[-na] * dgamma.base[-1]
      d2phi.dgamma.dalpha[cbind(1:(na-1),2:na,1:(na-1))] = -invdiffgamma2[-na] * dgamma.base[-na]
      d2phi.dgamma.dalpha[cbind(1:(na-1),2:na,2:na)] = invdiffgamma2[-na] * dgamma.base[-1]
    }

    ##
    d2phi.dgamma.dbeta = array(0,c(na,na,nb))
    rowdiff = matrix(0,na,na)
    diag(rowdiff) = 1
    if(na > 1) rowdiff[cbind(1:(na-1),2:na)] = -1
    d2phi.dgamma.dbeta.comp1 = diagn(-invdiffgamma2) %*% rowdiff %*% dgamma.dbeta
    d2phi.dgamma.dbeta.comp2 = diagn(-invgamma2) %*% dgamma.dbeta - d2phi.dgamma.dbeta.comp1
    for(j in 1:na) {
      d2phi.dgamma.dbeta[j,j,] = d2phi.dgamma.dbeta.comp2[j,]
      if(j < na)
        d2phi.dgamma.dbeta[j,j+1,] = d2phi.dgamma.dbeta.comp1[j,]
    }

    ##
    d2gamma.dalpha.dbeta = array(0,c(na,na,nb))
    for(j in 1:na)
      d2gamma.dalpha.dbeta[j,j,] = d2gamma.base[j] %o% x

    ##
    d2gamma.dbeta.dbeta = d2gamma.base %o% x %o% x


#### First derivatives of log-likelihood (score functions)
    dl.dalpha = dl.dphi %*% dphi.dgamma %*% dgamma.dalpha
    dl.dbeta = dl.dphi %*% dphi.dgamma %*% dgamma.dbeta
    dl.dtheta[i,] = c(dl.dalpha, dl.dbeta)


#### Second derivatives of log-likelihood
#### Since first derivative is a sum of terms each being a*b*c,
#### second derivative is a sum of terms each being (a'*b*c+a*b'*c+a*b*c').

#### d2l.dalpha.dalpha
    ## Obtain aprime.b.c
    ## Transpose first so that matrix multiplication is meaningful.
    ## Then transpose so that column is indexed by second alpha.
    aprime.b.c = t(crossprod(d2l.dphi.dalpha, dphi.dgamma %*% dgamma.dalpha))

    ## Obtain a.bprime.c
    ## run through the index of second alpha
    a.bprime.c = matrix(,na,na)
    for(j in 1:na)
      a.bprime.c[,j] = t.dl.dphi %*% d2phi.dgamma.dalpha[,,j] %*% dgamma.dalpha

    ## Obtain a.b.cprime
    ## cprime = d2gamma.dalpha.dalpha = 0 if indices of the two alphas differ.
    d2gamma.dalpha.dalpha = diagn(d2gamma.base)
    a.b.cprime = diagn(as.vector(dl.dphi %*% dphi.dgamma %*% d2gamma.dalpha.dalpha))

    ## summing over individuals
    d2l.dalpha.dalpha = aprime.b.c + a.bprime.c + a.b.cprime + d2l.dalpha.dalpha


#### d2l.dalpha.dbeta
    aprime.b.c = t(crossprod(d2l.dphi.dbeta, dphi.dgamma %*% dgamma.dalpha))
    a.bprime.c = a.b.cprime = matrix(,na,nb)
    for(j in 1:nb) {
      a.bprime.c[,j] = t.dl.dphi %*% d2phi.dgamma.dbeta[,,j] %*% dgamma.dalpha
      a.b.cprime[,j] = t.dl.dphi %*% dphi.dgamma %*% d2gamma.dalpha.dbeta[,,j]
    }
    d2l.dalpha.dbeta = aprime.b.c + a.bprime.c + a.b.cprime + d2l.dalpha.dbeta


#### d2l.dbeta.dalpha
#    dl.dbeta = dl.dphi %*% dphi.dgamma %*% dgamma.dbeta
#    aprime.b.c = t(crossprod(d2l.dphi.dalpha, dphi.dgamma %*% dgamma.dbeta))
#    a.bprime.c = a.b.cprime = matrix(,na,nb)
#    for(j in 1:nb) {
#      a.bprime.c[,j] = t.dl.dphi %*% d2phi.dgamma.dalpha[,,j] %*% dgamma.dbeta
#      a.b.cprime[,j] = t.dl.dphi %*% dphi.dgamma %*% d2gamma.dbeta.dalpha[,,j]
#    }
#    d2l.dbeta.dalpha = aprime.b.c + a.bprime.c + a.b.cprime + d2l.dbeta.dalpha


#### d2l.dbeta.dbeta
    aprime.b.c = t(crossprod(d2l.dphi.dbeta, dphi.dgamma %*% dgamma.dbeta))
    a.bprime.c = a.b.cprime = matrix(,nb,nb)
    for(j in 1:nb) {
      a.bprime.c[,j] = t.dl.dphi %*% d2phi.dgamma.dbeta[,,j] %*% dgamma.dbeta
      a.b.cprime[,j] = t.dl.dphi %*% dphi.dgamma %*% d2gamma.dbeta.dbeta[,,j]
    }
    d2l.dbeta.dbeta = aprime.b.c + a.bprime.c + a.b.cprime + d2l.dbeta.dbeta 


#### Derivatives of predicted probabilities
    p0[i,] = diff(c(0, gamma, 1))

    rowdiff = matrix(0,ny,na)
    diag(rowdiff) = 1
    rowdiff[cbind(2:ny,1:na)] = -1
    dp0.dalpha = rowdiff %*% dgamma.dalpha
    dp0.dbeta = rowdiff %*% dgamma.dbeta

    dp0.dtheta[i,,] = cbind(dp0.dalpha, dp0.dbeta)
  }

#### Final assembly
  d2l.dtheta.dtheta = rbind(
    cbind(d2l.dalpha.dalpha, d2l.dalpha.dbeta),
    cbind(t(d2l.dalpha.dbeta), d2l.dbeta.dbeta))

  ## sandwich variance estimate: ABA', where
  ## A = (-d2l.dtheta.dtheta/N)^(-1)
  ## B = B0/N
  ## One way of coding:
  ## A0 = solve(d2l.dtheta.dtheta)
  ## B0 = t(dl.dtheta) %*% dl.dtheta
  ## var.theta = A0 %*% B0 %*% t(A0)
  ## Suggested coding for better efficiency and accuracy
  SS = solve(d2l.dtheta.dtheta, t(dl.dtheta))
  var.theta = tcrossprod(SS, SS)

  ## The sum of scores should be zero at the MLE.
  ##   apply(dl.dtheta, 2, sum)
  ## Sandwich variance estimate should be similar to information matrix, I,
  ## which is the same as the lrm() output mod$var.
  ##   I = -solve(d2l.dtheta.dtheta)
  ##   print(I)
  ##   print(mod$var)
  ##   print(var.theta)

  ## dlow.dtheta and dhi.dtheta
  npar.z = dim(dl.dtheta)[2]
  dlow.dtheta = dhi.dtheta = matrix(, npar.z, N)
  for(i in 1:N) {
    if (y[i] == 1) {
      dlow.dtheta[,i] <- 0
    } else {
      dlow.dtheta[,i] <- dgamma.dtheta[i,y[i]-1,]
    }

    if (y[i] == ny) {
      dhi.dtheta[,i] <- 0
    } else {
      dhi.dtheta[,i] <- -dgamma.dtheta[i,y[i],]
    }
  }

  low.x = cbind(0, Gamma)[cbind(1:N, y)]
  hi.x = cbind(1-Gamma, 0)[cbind(1:N, y)]

  presid <- low.x - hi.x
  dpresid.dtheta <- dlow.dtheta - dhi.dtheta

  list(mod = mod,
       presid=presid,
       dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       var.theta = var.theta,
       p0 = p0,
       dp0.dtheta = dp0.dtheta,
       Gamma = Gamma,
       dgamma.dtheta = dgamma.dtheta,
       dlow.dtheta=dlow.dtheta,
       dhi.dtheta=dhi.dtheta,
       dpresid.dtheta=dpresid.dtheta)
}

#' @importFrom stats coef
ordinal.scores <- function(mf, mm, method) {
  ## mf is the model.frame of the data

  if (method[1]=='logit'){
    return(ordinal.scores.logit(y=as.numeric(model.response(mf)),X=mm))
  }

  ## Fit proportional odds model and obtain the MLEs of parameters.
  mod <- newpolr(mf, Hess=TRUE, method=method,control=list(reltol=1e-50,maxit=100))
  y <- model.response(mf)
  
  ## N: number of subjects; ny: number of y categories
  N = length(y)
  ny = length(mod$lev)

  ## na, nb: number of parameters in alpha and beta
  na = ny - 1

  x <- mm

  nb = ncol(x)

  npar = na + nb

  alpha = mod$zeta
  beta = -coef(mod)
  eta <- mod$lp

  ## Predicted probabilities p0 and dp0.dtheta are stored for individuals.
  p0 = mod$fitted.values
  dp0.dtheta = array(dim=c(N, ny, npar))
  Y <- matrix(nrow=N,ncol=ny)
  .Ycol <- col(Y)

  edcumpr <- cbind(mod$dcumpr, 0)
  e1dcumpr <- cbind(0,mod$dcumpr)

  for(i in seq_len(na)) {
    dp0.dtheta[,,i] <- (.Ycol == i) * edcumpr - (.Ycol == i + 1)*e1dcumpr
  }

  for(i in seq_len(nb)) {
    dp0.dtheta[,,na+i] <- mod$dfitted.values * x[,i]
  }
  
  ## Cumulative probabilities
  Gamma = mod$cumpr
  dcumpr <- mod$dcumpr
  
  ## Scores are stored for individuals.
  dl.dtheta = mod$grad

  ## dlow.dtheta and dhigh.dtheta
  y <- as.integer(y)
  dcumpr.x <- cbind(0, dcumpr, 0)
  dlow.dtheta <- t(cbind((col(dcumpr) == (y - 1L)) * dcumpr,
                         dcumpr.x[cbind(1:N,y)] * x))
  dhi.dtheta <- t(-cbind((col(dcumpr) == y) * dcumpr,
                         dcumpr.x[cbind(1:N,y + 1L)] * x))

  d2l.dtheta.dtheta <- mod$Hessian

  low.x = cbind(0, Gamma)[cbind(1:N, y)]
  hi.x = cbind(1-Gamma, 0)[cbind(1:N, y)]

  presid <- low.x - hi.x
  dpresid.dtheta <- dlow.dtheta - dhi.dtheta

  list(mod = mod,
       presid=presid,
       dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       p0 = p0,
       dp0.dtheta = dp0.dtheta,
       Gamma = Gamma,
       dcumpr=dcumpr,
       dlow.dtheta=dlow.dtheta,
       dhi.dtheta=dhi.dtheta,
       dpresid.dtheta=dpresid.dtheta)
}

#' Conditional ordinal by ordinal tests for association.
#'
#' \code{cobot} tests for independence between two ordered categorical
#' variables, \var{X} and \var{Y} conditional on other variables,
#' \var{Z}.  The basic approach involves fitting models of \var{X} on
#' \var{Z} and \var{Y} on \var{Z} and determining whether there is any
#' remaining information between \var{X} and \var{Y}.  This is done by
#' computing one of 3 test statistics.  \code{T1} compares empirical
#' distribution of \var{X} and \var{Y} with the joint fitted
#' distribution of \var{X} and \var{Y} under independence conditional
#' on \var{Z}. \code{T2} computes the correlation between ordinal
#' (probability-scale) residuals from both models and tests the null
#' of no residual correlation.  \code{T3} evaluates the
#' concordance--disconcordance of data drawn from the joint fitted
#' distribution of \var{X} and \var{Y} under conditional independence
#' with the empirical distribution. Details are given in \cite{Li C and
#' Shepherd BE, Test of association between two ordinal variables
#' while adjusting for covariates. Journal of the American Statistical
#' Association 2010, 105:612-620}.
#'
#' formula is specified as \code{\var{X} | \var{Y} ~ \var{Z}}.
#' This indicates that models of \code{\var{X} ~ \var{Z}} and
#' \code{\var{Y} ~ \var{Z}} will be fit.  The null hypothsis to be
#' tested is \eqn{H_0 : X}{H0 : X} independant of \var{Y} conditional
#' on \var{Z}.
#' 
#' Note that \code{T2} can be thought of as an adjusted rank
#' correlation.(\cite{Li C and Shepherd BE, A new residual for ordinal
#' outcomes. Biometrika 2012; 99:473-480})
#'
#' @param formula an object of class \code{\link{Formula}} (or one
#'   that can be coerced to that class): a symbolic description of the
#'   model to be fitted.  The details of model specification are given
#'   under \sQuote{Details}.
#' @param link The link family to be used for ordinal models of both
#' \var{X} and \var{Y}.  Defaults to \samp{logit}. Other options are
#' \samp{probit}, \samp{cloglog}, and \samp{cauchit}.
#' @param link.x The link function to be used for a model of the first
#' ordered variable. Defaults to value of \code{link}.
#' @param link.y The link function to be used for a model of the
#' second variable. Defaults to value of \code{link}.
#' @param data an optional data frame, list or environment (or object
#' coercible by \code{\link{as.data.frame}} to a data frame)
#' containing the variables in the model.  If not found in
#' \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{cobot} is called.
#' @param subset an optional vector specifying a subset of
#' observations to be used in the fitting process.
#' @param na.action how \code{NA}s are treated.
#' @param fisher  logical; if \code{TRUE}, Fisher transformation and delta method a
#' used to compute p value for the test statistic based on correlation of
#' residuals.
#' @param conf.int numeric specifying confidence interval coverage.
#' @return object of \samp{cobot} class.
#' @references Li C and Shepherd BE, Test of association between two
#' ordinal variables while adjusting for covariates. Journal of the
#' American Statistical Association 2010, 105:612-620.
#' @references Li C and Shepherd BE, A new residual for ordinal
#' outcomes. Biometrika 2012; 99:473-480
#' @import Formula
#' @importFrom stats na.fail model.matrix contrasts model.frame model.response cor
#' @export
#' @seealso \code{\link{Formula}}, \code{\link{as.data.frame}}
#' @include newPolr.R
#' @include diagn.R
#' @include GKGamma.R
#' @include pgumbel.R
#' @examples
#' data(PResidData)
#' cobot(x|y~z, data=PResidData)
cobot <- function(formula, link=c("logit", "probit", "cloglog", "cauchit"),
                  link.x=link,
                  link.y=link,
                  data, subset, na.action=na.fail,fisher=FALSE,conf.int=0.95) {
  F1 <- Formula(formula)
  Fx <- formula(F1, lhs=1)
  Fy <- formula(F1, lhs=2)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.action
  # We set xlev to a benign non-value in the call so that it won't get partially matched
  # to any variable in the formula. For instance a variable named 'x' could possibly get
  # bound to xlev, which is not what we want.
  mf$xlev <- integer(0) 
  mf[[1L]] <- as.name("model.frame")


  mx <- my <- mf

  # NOTE: we add the opposite variable to each model frame call so that
  # subsetting occurs correctly. Later we strip them off.
  mx[["formula"]] <- Fx
  yName <- paste0('(', all.vars(Fy[[2]])[1], ')')
  mx[[yName]] <- Fy[[2]]

  my[["formula"]] <- Fy
  xName <- paste0('(', all.vars(Fx[[2]])[1], ')')
  my[[xName]] <- Fx[[2]]

  mx <- eval(mx, parent.frame())
  mx[[paste0('(',yName,')')]] <- NULL
  
  my <- eval(my, parent.frame())
  my[[paste0('(',xName,')')]] <- NULL

  data.points <- nrow(mx)

  Terms <- attr(mx, "terms")
  zz <- model.matrix(Terms, mx, contrasts)
  zzint <- match("(Intercept)", colnames(zz), nomatch = 0L)
  if(zzint > 0L) {
    zz <- zz[, -zzint, drop = FALSE]
  }

  score.xz <- ordinal.scores(mx, zz, method=link.x)
  score.yz <- ordinal.scores(my, zz, method=link.y)

  npar.xz = dim(score.xz$dl.dtheta)[2]
  npar.yz = dim(score.yz$dl.dtheta)[2]


  xx = as.integer(model.response(mx)); yy = as.integer(model.response(my))
  nx = length(table(xx))
  ny = length(table(yy))
  N = length(yy)


#### Asymptotics for T3 = mean(Cprob) - mean(Dprob)
  ##  If gamma.x[0]=0 and gamma.x[nx]=1, then
  ##  low.x = gamma.x[x-1], hi.x = (1-gamma.x[x])
  low.x = cbind(0, score.xz$Gamma)[cbind(1:N, xx)]
  low.y = cbind(0, score.yz$Gamma)[cbind(1:N, yy)]
  hi.x = cbind(1-score.xz$Gamma, 0)[cbind(1:N, xx)]
  hi.y = cbind(1-score.yz$Gamma, 0)[cbind(1:N, yy)]
  Cprob = low.x*low.y + hi.x*hi.y
  Dprob = low.x*hi.y + hi.x*low.y
  mean.Cprob = mean(Cprob)
  mean.Dprob = mean(Dprob)
  T3 = mean.Cprob - mean.Dprob

  dCsum.dthetax = score.xz$dlow.dtheta %*% low.y + score.xz$dhi.dtheta %*% hi.y
  dCsum.dthetay = score.yz$dlow.dtheta %*% low.x + score.yz$dhi.dtheta %*% hi.x
  dDsum.dthetax = score.xz$dlow.dtheta %*% hi.y  + score.xz$dhi.dtheta %*% low.y
  dDsum.dthetay = score.yz$dlow.dtheta %*% hi.x  + score.yz$dhi.dtheta %*% low.x

  dT3sum.dtheta = c(dCsum.dthetax-dDsum.dthetax, dCsum.dthetay-dDsum.dthetay)

  ## Estimating equations for (theta, p3)
  ## theta is (theta.xz, theta.yz) and the equations are score functions.
  ## p3 is the "true" value of test statistic, and the equation is
  ## p3 - (Ci-Di)
  bigphi = cbind(score.xz$dl.dtheta, score.yz$dl.dtheta, T3-(Cprob-Dprob))

  ## sandwich variance estimate of var(thetahat)
  Ntheta = npar.xz + npar.yz + 1
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = score.xz$d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = score.yz$d2l.dtheta.dtheta
  A[Ntheta, -Ntheta] = -dT3sum.dtheta
  A[Ntheta, Ntheta] = N

  ## One way of coding:
  ##B = t(bigphi) %*% bigphi
  ##var.theta = solve(A) %*% B %*% solve(A)
  ## Suggested coding for better efficiency and accuracy:
  SS = solve(A, t(bigphi))
  var.theta = tcrossprod(SS, SS)
  varT3 = var.theta[Ntheta, Ntheta]

  pvalT3 = 2 * pnorm(-abs(T3)/sqrt(varT3))


#### Asymptotics for T4 = (mean(Cprob)-mean(Dprob))/(mean(Cprob)+mean(Dprob))
  T4 = (mean.Cprob - mean.Dprob)/(mean.Cprob + mean.Dprob)

  ## Estimating equations for (theta, P4)
  ## theta is (theta.xz, theta.yz) and the equations are score functions.
  ## P4 is a vector of (cc, dd, p4).  Their corresponding equations are:
  ## cc - Ci
  ## dd - Di
  ## p4 - (cc-dd)/(cc+dd)
  ## Then p4 is the "true" value of test statistic.
  bigphi = cbind(score.xz$dl.dtheta, score.yz$dl.dtheta,
    mean.Cprob - Cprob, mean.Dprob - Dprob, 0)

  ## sandwich variance estimate of var(thetahat)
  Ntheta = npar.xz + npar.yz + 3
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = score.xz$d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = score.yz$d2l.dtheta.dtheta
  A[Ntheta-3+(1:3), Ntheta-3+(1:3)] = diag(N, 3)
  A[Ntheta-2, 1:(npar.xz+npar.yz)] = -c(dCsum.dthetax, dCsum.dthetay)
  A[Ntheta-1, 1:(npar.xz+npar.yz)] = -c(dDsum.dthetax, dDsum.dthetay)

  revcpd = 1/(mean.Cprob + mean.Dprob)
  dT4.dcpd = (mean.Cprob-mean.Dprob)*(-revcpd^2)
  A[Ntheta, Ntheta-3+(1:2)] = -N * c(revcpd+dT4.dcpd, -revcpd+dT4.dcpd)

  ## One way of coding:
  ##B = t(bigphi) %*% bigphi
  ##var.theta = solve(A) %*% B %*% solve(A)
  ## Suggested coding for better efficiency and accuracy:
  SS = solve(A, t(bigphi))
  var.theta = tcrossprod(SS, SS)
  varT4 = var.theta[Ntheta, Ntheta]

  pvalT4 = 2 * pnorm(-abs(T4)/sqrt(varT4))


#### Asymptotics for T2 = cor(hi.x - low.x, hi.y - low.y)
  xresid = hi.x - low.x
  yresid = hi.y - low.y
  T2 = cor(xresid, yresid)

  xbyyresid = xresid * yresid
  xresid2 = xresid^2
  yresid2 = yresid^2
  mean.xresid = mean(xresid)
  mean.yresid = mean(yresid)
  mean.xbyyresid = mean(xbyyresid)
  ## T2 also equals numT2 / sqrt(varprod) = numT2 * revsvp
  numT2 = mean.xbyyresid - mean.xresid * mean.yresid
  var.xresid = mean(xresid2) - mean.xresid^2
  var.yresid = mean(yresid2) - mean.yresid^2
  varprod = var.xresid * var.yresid
  revsvp = 1/sqrt(varprod)
  dT2.dvarprod = numT2 * (-0.5) * revsvp^3

  ## Estimating equations for (theta, P5)
  ## theta is (theta.xz, theta.yz) and the equations are score functions.
  ## P5 is a vector (ex, ey, crossxy, ex2, ey2, p5).
  ## Their corresponding equations are:
  ## ex - (hi.x-low.x)[i]
  ## ey - (hi.y-low.y)[i]
  ## crossxy - ((hi.x-low.x)*(hi.y-low.y))[i]
  ## ex2 - (hi.x-low.x)[i]^2
  ## ey2 - (hi.y-low.y)[i]^2
  ## p5 - (crossxy-ex*ey)/sqrt((ex2-ex^2)*(ey2-ey^2))
  ## Then p5 is the "true" value of test statistic
  bigphi = cbind(score.xz$dl.dtheta, score.yz$dl.dtheta,
    mean.xresid - xresid, mean.yresid - yresid, mean.xbyyresid - xbyyresid,
    mean(xresid2) - xresid2, mean(yresid2) - yresid2, 0)

  ## sandwich variance estimate of var(thetahat)
  Ntheta = npar.xz + npar.yz + 6
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = score.xz$d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = score.yz$d2l.dtheta.dtheta
  A[Ntheta-6+(1:6), Ntheta-6+(1:6)] = diag(N, 6)

  dxresid.dthetax = score.xz$dhi.dtheta - score.xz$dlow.dtheta
  dyresid.dthetay = score.yz$dhi.dtheta - score.yz$dlow.dtheta
  bigpartial = rbind(c(dxresid.dthetax %*% rep(1, N), rep(0, npar.yz)),
    c(rep(0, npar.xz), dyresid.dthetay %*% rep(1, N)),
    c(dxresid.dthetax %*% yresid, dyresid.dthetay %*% xresid),
    c(dxresid.dthetax %*% (2*xresid), rep(0, npar.yz)),
    c(rep(0, npar.xz), dyresid.dthetay %*% (2*yresid)))
  A[Ntheta-6+(1:5), 1:(npar.xz+npar.yz)] = -bigpartial

  smallpartial = N *
    c(-mean.yresid * revsvp + dT2.dvarprod * (-2*mean.xresid*var.yresid),
      -mean.xresid * revsvp + dT2.dvarprod * (-2*mean.yresid*var.xresid),
      revsvp,
      dT2.dvarprod * var.yresid,
      dT2.dvarprod * var.xresid)
  A[Ntheta, Ntheta-6+(1:5)] = -smallpartial

  ## One way of coding:
  ##B = t(bigphi) %*% bigphi
  ##var.theta = solve(A) %*% B %*% solve(A)
  ## Suggested coding for better efficiency and accuracy:
  SS = solve(A, t(bigphi))
  var.theta = tcrossprod(SS, SS)
  varT2 = var.theta[Ntheta, Ntheta]

  if (fisher==TRUE){
    ####Fisher's transformation
    ## TS_f: transformed the test statistics
    ## var.TS_f: variance estimate for ransformed test statistics
    ## pvalT2: p-value based on transformed test statistics
    TS_f <- log( (1+T2)/(1-T2) )  
    var.TS_f <- varT2*(2/(1-T2^2))^2 
    pvalT2 <- 2 * pnorm( -abs(TS_f)/sqrt(var.TS_f))  
  } else {
    pvalT2 = 2 * pnorm(-abs(T2)/sqrt(varT2))
  }


#### Asymptotics for T1 = tau - tau0
  ## dtau0/dtheta
  ## P0 is the sum of product predicted probability matrix with dim(nx,ny)
  P0 = crossprod(score.xz$p0, score.yz$p0) / N

  cdtau0 = GKGamma(P0)
  C0 = cdtau0$scon
  D0 = cdtau0$sdis
  ## C0 = sum_{l>j,m>k} {P0[j,k] * P0[l,m]}
  ## D0 = sum_{l>j,m<k} {P0[j,k] * P0[l,m]}
  dC0.dP0 = matrix(,nx,ny)
  dD0.dP0 = matrix(,nx,ny)
  for(i in 1:nx)
    for(j in 1:ny) {
      dC0.dP0[i,j] = ifelse(i>1 & j>1, sum(P0[1:(i-1), 1:(j-1)]), 0) +
        ifelse(i<nx & j<ny, sum(P0[(i+1):nx, (j+1):ny]), 0)
      dD0.dP0[i,j] = ifelse(i>1 & j<ny, sum(P0[1:(i-1), (j+1):ny]), 0) +
        ifelse(i<nx & j>1, sum(P0[(i+1):nx, 1:(j-1)]), 0)
    }

  ## tau0 = (C0-D0)/(C0+D0)
  dtau0.dC0 = 2*D0/(C0+D0)^2
  dtau0.dD0 =-2*C0/(C0+D0)^2

  ## ## P0 is already a matrix
  dP0.dtheta.x = array(0, c(nx, ny, npar.xz))
  for(j in 1:ny) {
    aa = matrix(0, nx, npar.xz)
    for(i in 1:N)
      aa = aa + score.xz$dp0.dtheta[i,,] * score.yz$p0[i,j]
    dP0.dtheta.x[,j,] = aa/N
    ## simpler but mind-tickling version
    #dP0.dtheta.x[,j,] = (score.yz$p0[,j] %*% matrix(score.xz$dp0.dtheta,N))/N
  }

  dP0.dtheta.y = array(0, c(nx, ny, npar.yz))
  for(j in 1:nx) {
    aa = matrix(0, ny, npar.yz)
    for(i in 1:N)
      aa = aa + score.yz$dp0.dtheta[i,,] * score.xz$p0[i,j]
    dP0.dtheta.y[j,,] = aa/N
  }

  ## dC0.dtheta and dD0.dtheta
  dC0.dtheta.x = as.numeric(dC0.dP0) %*% matrix(dP0.dtheta.x, nx*ny)
  dD0.dtheta.x = as.numeric(dD0.dP0) %*% matrix(dP0.dtheta.x, nx*ny)
  dC0.dtheta.y = as.numeric(dC0.dP0) %*% matrix(dP0.dtheta.y, nx*ny)
  dD0.dtheta.y = as.numeric(dD0.dP0) %*% matrix(dP0.dtheta.y, nx*ny)

  ## dtau0/dtheta
  dtau0.dtheta.x = dtau0.dC0 * dC0.dtheta.x + dtau0.dD0 * dD0.dtheta.x
  dtau0.dtheta.y = dtau0.dC0 * dC0.dtheta.y + dtau0.dD0 * dD0.dtheta.y


  ## dtau/dPa
  ## tau = (C-D)/(C+D)
  Pa = table(xx, yy) / N
  cdtau = GKGamma(Pa)
  C = cdtau$scon
  D = cdtau$sdis
  dtau.dC = 2*D/(C+D)^2
  dtau.dD =-2*C/(C+D)^2

  ## Pa[nx,ny] is not a parameter, but = 1 - all other Pa parameters.
  ## Thus, d.Pa[nx,ny]/d.Pa[i,j] = -1.
  ## Also, d.sum(Pa[-nx,-ny]).dPa[i,j] = 1 when i<nx and j<ny, and 0 otherwise.
  ##
  ## In C = sum_{l>j,m>k} {Pa[j,k] * Pa[l,m]}, Pa[i,j] appears in
  ##   Pa[i,j] * XX (minus Pa[nx,ny] if i<nx & j<ny), and in
  ##   sum(Pa[-nx,-ny]) * Pa[nx,ny].
  ## So, dC.dPa[i,j] = XX (minus Pa[nx,ny] if i<nx & j<ny)
  ##                  + d.sum(Pa[-nx,-ny]).dPa[i,j] * Pa[nx,ny]
  ##                  - sum(Pa[-nx,-ny])
  ##                 = XX (with Pa[nx,ny] if present) - sum(Pa[-nx,-ny])
  ##
  ## D = sum_{l>j,m<k} {Pa[j,k] * Pa[l,m]} doesn't contain Pa[nx,ny]
  dC.dPa = matrix(,nx,ny)
  dD.dPa = matrix(,nx,ny)
  for(i in 1:nx)
    for(j in 1:ny) {
      dC.dPa[i,j] = ifelse(i>1 & j>1, sum(Pa[1:(i-1), 1:(j-1)]), 0) +
        ifelse(i<nx & j<ny, sum(Pa[(i+1):nx, (j+1):ny]), 0) - sum(Pa[-nx,-ny])
      dD.dPa[i,j] = ifelse(i>1 & j<ny, sum(Pa[1:(i-1), (j+1):ny]), 0) +
        ifelse(i<nx & j>1, sum(Pa[(i+1):nx, 1:(j-1)]), 0)
    }

  dtau.dPa = dtau.dC * dC.dPa + dtau.dD * dD.dPa
  dtau.dPa = dtau.dPa[-length(dtau.dPa)]  ## remove the last value


  ## Estimating equations for (theta, phi)
  ## theta is (theta.xz, theta.yz) and the equations are score functions.
  ## phi is (p_ij) for (X,Y), and the equations are
  ## I{subject in cell (ij)} - p_ij
  phi.Pa = matrix(0, N, nx*ny)
  phi.Pa[cbind(1:N, xx+(yy-1)*nx)] = 1
  phi.Pa = phi.Pa - rep(1,N) %o% as.numeric(Pa)
  phi.Pa = phi.Pa[,-(nx*ny)]

  bigphi = cbind(score.xz$dl.dtheta, score.yz$dl.dtheta, phi.Pa)

  ## sandwich variance estimate of var(thetahat, phihat)
  Ntheta = npar.xz + npar.yz + nx*ny-1
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = score.xz$d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = score.yz$d2l.dtheta.dtheta
  A[-(1:(npar.xz+npar.yz)), -(1:(npar.xz+npar.yz))] = -diag(N, nx*ny-1)

  ## One way of coding:
  ##B = t(bigphi) %*% bigphi
  ##var.theta = solve(A) %*% B %*% solve(A)
  ## Suggested coding for better efficiency and accuracy:
  ##SS = solve(A, t(bigphi))
  ##var.theta = SS %*% t(SS)
  ## Or better yet, no need to explicitly obtain var.theta.  See below.

  ## Test statistic T1 = tau - tau0
  T1 = cdtau$gamma - cdtau0$gamma
  ## dT.dtheta has length nx + ny + nx*ny-1
  dT1.dtheta = c(-dtau0.dtheta.x, -dtau0.dtheta.y, dtau.dPa)

  ## variance of T, using delta method
  ##varT = t(dT.dtheta) %*% var.theta %*% dT.dtheta
  SS = crossprod(dT1.dtheta, solve(A, t(bigphi)))
  varT1 = sum(SS^2)
  pvalT1 = 2 * pnorm(-abs(T1)/sqrt(varT1))

  ans <- structure(
          list(
            TS=list(
              T1=list(ts=T1, var=varT1, pval=pvalT1, label="Gamma(Obs) - Gamma(Exp)"), 
              T2=list(ts=T2, var=varT2, pval=pvalT2, label="Correlation of Residuals"),
              T3=list(ts=T3, var=varT3, pval=pvalT3, label="Covariance of Residuals")
            ),
            fisher=fisher, 
            conf.int=conf.int,
            data.points=data.points
          ),
          class="cobot")

  # Apply confidence intervals
  for (i in seq_len(length(ans$TS))){
    ts_ci <- getCI(ans$TS[[i]]$ts,ans$TS[[i]]$var,ans$fisher,conf.int)
    ans$TS[[i]]$lower <- ts_ci[1]
    ans$TS[[i]]$upper <- ts_ci[2]
  }

  ans
}
