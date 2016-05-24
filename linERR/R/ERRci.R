ERRci <- function(object, prob=0.95)
{
  beta <- attr(object, "beta")
  max.exp <- attr(object, "max.exp") 
  covariates1 <- attr(object, "covariates1")
  covars1_2 <- ifelse(is.null(covariates1), NA, apply(covariates1, 2, as.list))
  if (all(is.na(covars1_2))) covars1_2 <- NULL
  covars1_2  <- lapply(covars1_2, as.numeric)
  data_2 <- attr(object, "data_2")
  rsets_2 <- attr(object, "rsets_2")
  doses_2 <- attr(object, "doses_2")
  ages_2 <- attr(object, "ages_2")
  lag <- attr(object, "Call")$lag
  failtimes <- as.list(names(rsets_2)[1:(length(rsets_2)-1)])
  failtimes <- lapply(failtimes, as.numeric)
  
  add.arguments <- function(f,n){
    # adds n arguments to a function f; returns that new function 
    t = paste("arg <- alist(",
              paste(sapply(1:n, function(i) paste("x",i, "=",sep="")), collapse=","),
              ")", sep="")
    formals(f) <- eval(parse(text=t))
    f
  }
  p.est <- function()
  {
    beta3  <- vector()
    for (i in 1:length(beta))
    {
      beta3[i] <- eval(parse(text = paste0("x", i)))[i]
    }
    beta_2   <- as.list(beta3)
    beta_2   <- lapply(beta_2, as.numeric)
    suppressWarnings(if (is.null(covariates1))
    {
      ncovs1 <- 0
    }else{
      ncovs1 <- dim(covariates1)[2]
    })
    res1 <- .Call("llhood", beta_2, data_2, rsets_2, ages_2, doses_2, covars1_2, as.integer(length(doses_2[[1]])),
                  as.integer(length(rsets_2)-1), as.integer(length(data_2)), failtimes, as.integer(ncovs1), 
                  as.integer(max.exp), as.numeric(lag))
    return(-sum(res1*is.finite(res1), na.rm=T))
  }
  p.est <- add.arguments(p.est, length(beta))
  
  up <- (attr(object, "aic")-2*length(beta))/2 + qchisq(prob, 1)/2
  g <- function(x1)
  {
    beta$x1 <<- x1
    if (length(beta) > 1)
    {
      return(-optim(par = unlist(beta), fn = p.est, x1 = x1)$value + up)
    }
    return(-p.est(unlist(beta)))
  }
  l1 <- uniroot(g, lower=beta$x1-0.001, upper=beta$x1+0.001, extendInt="yes")$root
  l2 <- ifelse(l1 < attr(object, "coef")[1], uniroot(g, lower=attr(object, "coef")[1], upper=attr(object, "coef")[1]+0.001, extendInt="downX")$root,
                                             uniroot(g, lower=attr(object, "coef")[1]-0.01, upper=attr(object, "coef")[1], extendInt="upX")$root)
  ci <- c(l1, l2) 
  if (l1 > l2) ci <- c(l2, l1)
  names(ci) <- c(paste0("lower ", (1-prob)*100/2, "%"), paste0("upper ", (1-(1-prob)/2)*100, "%"))
  return(ci)
}