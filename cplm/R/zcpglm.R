#######################################################
##    Zero-inflated compound Poisson GLM             ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################


zcpglm <- function(formula, link = "log", data, weights, offset, 
                  subset, na.action = NULL, contrasts = NULL, 
                  control = list(), optimizer = "nlminb") {

  call <- match.call()  
  if (missing(data)) 
    data <- environment(formula) 
  
  # get formula for the two parts: zero and tweedie 
  tmp <- strsplit(paste(deparse(formula), sep = "", collapse = ""), "\\|\\|")[[1]]
  if (length(tmp) == 1) {
    warning("formula for the zero part is missing: default to ~1")
    tmp <- c(tmp, "1")
  }
  ft <- as.formula(tmp[1])
  fz <- as.formula(paste(formula[[2]], "~", tmp[2]))
  # change environment so that model.frame will work (necessary for offset())!!!
  environment(ft) <- environment(fz) <- environment(formula)
  
  # generate model frame for each part
  call2 <- call 
  call2$formula <- ft # tweedie              
  frt <- cpglm.mf(call2, contrasts)
  call2$formula <- fz # zero
  frz <- cpglm.mf(call2, contrasts)
  fr <- list(Y = frt$Y, Xz = frz$X, Xt = frt$X, 
          offz = frz$off, offt = frt$off, wts = frt$wts)
  link.power <- make.link.power(link)
  control <- do.call("cplm.control", control)
  
  # estimate parameters 
  ans <- zcpglm.optim(fr, link.power, control, optimizer)  
  
  # return result
  ans@formula <- formula 
  ans@call <- call
  ans@na.action <- na.action
  ans@contrasts <- contrasts
  ans@model.frame <- list(zero = frz$mf, tweedie = frt$mf)
  return(ans)
}


# zero-inflated Tweedie using the Newton-Raphson method (EM is too slow)
zcpglm.optim <- function(fr, link.power = 0, 
                         control = list(), optimizer = "nlminb"){
  
  # dimensions
  nbz <- ncol(fr$Xz)
  nbt <- ncol(fr$Xt)
  nb <- nbz + nbt
  # link inverse functions
  tw <- tweedie(link.power = link.power)                 
  logit <- binomial()
  eps <- 0.001 #constant in numerical gradient/hessian
  yp <- fr$Y > 0 
  
  # statistics used in marginal loglikelihood and gradient
  llik_common <- function(parm){
    # extract parameters
    betaz <- parm[1:nbz]
    betat <- parm[(nbz + 1): nb]
    phi <- exp(parm[(nb + 1)]) #phi is on log scale
    p <- parm[(nb + 2)]
    etaz <- fr$Xz %*% betaz + fr$offz
    q <- as.numeric(logit$linkinv(etaz))
    etat <- fr$Xt %*% betat + fr$offt
    mu <- tw$linkinv(etat)
    ly <- as.numeric(log(dtweedie(y = fr$Y, mu = mu, 
                                  phi = phi, power = p)))
    ml <- q * as.integer(fr$Y == 0) + exp(log(1 - q) + ly) 
    ml[ml == 0] <- .Machine$double.eps
    llik <- rep(NA, length(ly))
    llik[!yp] <- log(ml[!yp])
    llik[yp] <- (log(1 - q) + ly)[yp]
    
    # dtweedie could result in zero!!!
    return(list(etaz = etaz, etat = etat, q = q, mu = mu,
                fy = exp(ly), ml = ml, llik = llik))
  }
  
  # marginal loglikelihood f(y) used in Newton-Raphson 
  llik_zcpglm <- function(parm){
    ll <- llik_common(parm)
    ##FIXME: the weight is probably wrong?
    - sum(fr$wts * ll$llik)
  }
  
  # gradient: analytical expressions for regression parameters
  grad_zcpglm <- function(parm){
    ll <- llik_common(parm)
    gr <- rep(NA, nb + 2)  
    gr[1:nbz] <- - colSums(fr$wts * as.numeric(logit$mu.eta(ll$etaz)) * 
        (as.integer(fr$Y == 0) - ll$fy) / ll$ml * fr$Xz)
    gr[(nbz + 1):nb] <- - colSums(fr$wts * (1 - ll$q) * ll$fy * tw$mu.eta(ll$etat)
       * (fr$Y - ll$mu) / (ll$ml * ll$mu^parm[nb + 2]) * fr$Xt) / exp(parm[nb + 1])
    llik_phip <- function(pm)
      llik_zcpglm(c(parm[1:nb], pm))
    gr[(nb + 1):(nb + 2)] <- grad(parm[-(1:nb)], llik_phip)
    gr
  }
 
  # generate starting values
  p <- 1.5   
  fitt <- glm.fit(fr$Xt, fr$Y, weights = fr$wts, 
            offset = fr$offt, family = tweedie(var.power = p, 
            link.power = link.power))
  phi <- sum((fitt$weights * fitt$residuals^2)) / fitt$df.residual
  pt0 <- exp(-fitt$fitted.values^(2 - p) / (phi * (2 - p)))
  yz <- 1 - as.integer(as.logical(fr$Y)) - pt0 * (fr$Y == 0)
  fitz <- glm.fit(fr$Xz, yz, weights = fr$wts, 
            offset = fr$offz, family = quasibinomial())
  parm <- as.numeric(c(fitz$coefficients, fitt$coefficients, log(phi), p))
  
  # run optimization 
  opt_ans <- cplm_optim(parm, llik_zcpglm, gr = grad_zcpglm, 
                      lower = c(rep(-Inf, nb + 1), control$bound.p[1]),
                      upper = c(rep(Inf, nb + 1), control$bound.p[2]),
                      control = control, optimizer = optimizer)
  if (opt_ans$convergence) warning(opt_ans$message)
  
  # compute hessian matrix for regression parameters  
  # (optim also computes those for phi and p)
  grad_zcpglm2 <- function(parm){
    ll <- llik_common(parm)
    gr <- rep(NA, nb)  
    gr[1:nbz] <- - colSums(fr$wts * as.numeric(logit$mu.eta(ll$etaz)) * 
        (as.integer(fr$Y == 0) - ll$fy) / ll$ml * fr$Xz)
    gr[(nbz + 1):nb] <- - colSums(fr$wts * (1 - ll$q) * ll$fy * tw$mu.eta(ll$etat)
       * (fr$Y - ll$mu) / (ll$ml * ll$mu^parm[nb + 2]) * fr$Xt) / exp(parm[nb + 1])
    gr
  }
  parm <- opt_ans$par
  hn <- matrix(0, nb, nb)
  for (i in 1:nb){
    parm[i] <- parm[i] - eps
    g1 <- grad_zcpglm2(parm)
    parm[i] <- parm[i] + 2 * eps
    g2 <- grad_zcpglm2(parm)
    hn[i,] <- (g2 - g1) / ( 2 * eps)
    parm[i] <- parm[i] - eps
  }
  
  # extract stats for output                                   
  nmz <- dimnames(fr$Xz)[[2]]
  nmt <- dimnames(fr$Xt)[[2]]
  betaz <- opt_ans$par[1:nbz]
  betat <- opt_ans$par[(nbz + 1):nb]
  names(betaz) <- nmz 
  names(betat) <- nmt
  q <- logit$linkinv(fr$Xz %*% betaz + fr$offz)
  mu <- tw$linkinv(fr$Xt %*% betat + fr$offt)
  Yhat <- as.numeric((1 - q) * mu)
  nmv <- c(paste("zero_", nmz, sep = ""), paste("tw_", nmt, sep = ""))
  vc <- solve(hn)
  dimnames(vc) <- list(nmv, nmv)
  
  # return results
  out <- new("zcpglm", 
             coefficients = list(zero = betaz, tweedie = betat),
             residuals = as.numeric(sqrt(fr$wts) * (fr$Y - Yhat)),
             fitted.values = Yhat, call = call("foo"),
             df.residual = as.integer(nrow(fr$Xt) - nb),             
             formula = ~ 1, control = control, contrasts = NULL,
             p = opt_ans$par[nb + 2], phi = exp(opt_ans$par[nb + 1]), 
             converged = as.logical(ifelse(opt_ans$convergence, 0, 1)),
             link.power= link.power,  na.action = NULL,
             model.frame = list(), llik = - opt_ans$value,
             offset = list(zero = fr$offz, tweedie = fr$offt), 
             prior.weights = fr$wts, y = fr$Y,
             inits = NULL, vcov = vc)
  return(out)  
}               

