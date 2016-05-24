## \int f(x)dg(x) with step functions f,g
## integral reduces to finite sum over all jump times
stepIntegrate <- function(f, g) {
    ##f[is.nan(f)] <- 0
    sum(f * (g - c(1, g[-length(g)])))
}

evalFun <- function(times, x, y) {
    f <- approxfun(x, y, method="constant", yleft=1, rule=2, f=0)
    f(times)
}

## get all jump points (i.e. events) up to time L
getTimes <- function(L, data) {
    obs <- c(0, data$time[(data$to != "cens") & (data$time <= L)], data$Y[(data$D == 1) & (data$Y <= L)], L)
    sort(unique(obs))
}

## example: Surv(V, Y, D) ~ Trt + strata(W)
parseFormula <- function(formula, data, one.sample=FALSE) {

  Terms <- terms(formula, "strata")
  mf <- model.frame(Terms, data)

  ## extract Surv object
  surv <- model.response(mf)
  if(!is.Surv(surv)) stop("Model response must be 'Surv' object")

  n <- nrow(surv)
  
  if(ncol(surv) == 2) {
      Start <- rep.int(0, n)
      Stop <- surv[,1]
      status <- surv[,2]
  } else {
      Start <- surv[,1]
      Stop <- surv[,2]
      status <- surv[,3]
  }
  
  ## extract stratum indicator
  st <- attr(Terms, "specials")$strata
  if(is.null(st)) {
      W <- rep.int(0, n)
      st.labels <- NULL
  } else {
      st.labels <- dimnames(attr(Terms, "factors"))[[1]][st]
      W <- as.integer(mf[[st.labels]])
  }

  ## extract group indicator
  if(!one.sample) {
      sdi <- setdiff(labels(Terms), st.labels)
      if(length(sdi) > 0) Trt <- as.factor(as.matrix(mf[sdi])) ## multiple covariates are combined into one factor
      else stop("No treatment groups specified!")
  } else Trt <- rep.int(0, n)
  
  data.frame(V=Start, Y=Stop, D=status, W=factor(W), Trt=factor(Trt))
}
