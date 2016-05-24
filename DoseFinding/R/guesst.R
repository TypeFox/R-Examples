## function to calculate guesstimates for nonlinear parameters of the dose-response models

guesst <- function(d, p, model = c("emax", "exponential", "logistic", "quadratic",
                   "betaMod", "sigEmax"), less = TRUE,  local = FALSE, dMax, Maxd, scal){
  model <- match.arg(model)
  if(any(p <= 0) | any(p > 1)) stop("must have 0 < p <= 1")
  if(model == "emax"){
    if(!local){
      return(c(ed50 = mean(d * (1 - p)/p)))
    } else {
      if (any(p <= d/Maxd))
        stop("must have p > d/Maxd, for local version")
      val <- (d/p-d)/(1-d/(Maxd*p))
      return(c(ed50=mean(val)))
    }
  }
  if(model == "exponential"){
    if(any(p >= d/Maxd)) stop("must have p < d/Maxd")
    init <- d/log(1 + p)
    fooexp <- function(delta, d, p, Maxd){
      sum((exponential(d, 0, 1, delta)/
           exponential(Maxd, 0, 1, delta) - p)^2)           
    }
    val <- optimize(fooexp, c(0, 2*Maxd), d=d, p=p, Maxd=Maxd)$minimum
    return(c(delta = mean(val)))
  }
  if(model == "logistic"){
    if(length(d) == 1) {
      stop("logistic model needs at least two pairs (d,p)")
    }
    logit <- function(p) log(p/(1-p)) 
    if(length(d) == 2) {
      ed50  <- diff(rev(d)*logit(p))/diff(logit(p))
      delta <- diff(d)/diff(logit(p))
      return(c(ed50 = ed50, delta = delta))
    } else {
      m <- lm(logit(p)~d)
      par <- coef(m)
      names(par) <- NULL
      return(c(ed50 = -par[1]/par[2], delta = 1/par[2]))
    } 
    if(local){
      foolog <- function(par, d, p, Maxd) {
        e0 <- logistic(0,0,1,par[1],par[2])
        sum(((logistic(d,0,1,par[1],par[2]) - e0)/
             (logistic(Maxd,0,1,par[1],par[2])-e0)-p)^2)
      }
      res <- try(optim(par=res, fn=foolog, d=d, p=p, Maxd=Maxd))
      if(res$convergence > 0)
        stop("cannot find guesstimates for specified values")
      else res <- res$par
    }
    if(res[1] < 0)
      message("Message: specified values lead to negative ed50, which should be positive")
    return(res)
  }
  if(model == "quadratic"){
    aux <- sqrt(1 - p)
    if (less){
      return(c(delta = mean(-(1 - aux)/(2 * d))))
    } else {
      return(c(delta = mean(-(1 + aux)/(2 * d))))
    }
  }
  if(model == "betaMod"){
    if(scal <= dMax)
      stop("scal needs to be larger than dMax to calculate guesstimate")
    if(dMax > Maxd)
      stop("dose with maximum effect (dMax) needs to be smaller than maximum dose (Maxd)")
    k <- dMax/(scal-dMax)
    val <- d^k*(scal-d)/(dMax^k*(scal-dMax))
    beta <- log(p)/(log(val))
    return(c(delta1 = mean(k*beta), delta2 = mean(beta)))
  }
  if(model == "sigEmax"){
    if(length(d) == 1) {
      stop("sigmoid Emax model needs at least two pairs (d,p)")
    }
    if(length(d) == 2){
      num <- log((p[1]*(1-p[2]))/(p[2]*(1-p[1])))
      h <- num/log(d[1]/d[2])
      ed50 <- ((1-p[1])/p[1])^(1/h)*d[1]
      return(c(ed50=ed50, h=h))
    } else {
      y <- log((1-p)/p)
      x <- log(d)
      par <- coef(lm(y~x))
      names(par) <- NULL
      return(c(ed50 = exp(par[1]/-par[2]), delta = -par[2]))
    } 
    if(local) {
      fooSE <- function(par, d, p, Maxd) {
        sum((sigEmax(d,0,1,par[1],par[2])/
             sigEmax(Maxd,0,1,par[1],par[2])-p)^2)
      }
      res <- try(optim(par=res, fn=fooSE, d=d, p=p, Maxd=Maxd))
      if(res$convergence > 0)
        stop("cannot find guesstimates for specified values")
      else res <- res$par
    }
    if(res[1] < 0)
      message("Message: specified values lead to negative ed50, which should be positive")
    return(res)
  }
}
