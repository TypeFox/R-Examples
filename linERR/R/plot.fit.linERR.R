plot.fit.linERR <- function(x, ...)
{
  if (!exists(lower)) lower <- 0
  if (!exists(upper)) upper <- 1
  if (!exists(ci)) ci <- NULL
  
  beta <- attr(x, "beta")
  max.exp <- attr(x, "max.exp") 
  covariates1 <- attr(x, "covariates1")
  covars1_2 <- ifelse(is.null(covariates1), NA, apply(covariates1, 2, as.list))
  if (all(is.na(covars1_2))) covars1_2 <- NULL
  covars1_2  <- lapply(covars1_2, as.numeric)
  data_2 <- attr(x, "data_2")
  rsets_2 <- attr(x, "rsets_2")
  doses_2 <- attr(x, "doses_2")
  ages_2 <- attr(x, "ages_2")
  lag <- attr(x, "Call")$lag
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
  
  x <- seq(lower, upper, (upper-lower)/100)
  y <- vector()
  
  g <- function(x)
  {
    beta$x1 <<- x
    if (length(beta) > 1)
    {
      return(-optim(par = unlist(beta), fn = p.est, x1 = x)$value)
    }
    return(-p.est(unlist(beta)))
  }
  
  y <- vapply(x, g, 0.5)
  plot(x, y, type="l", ylab="Profile log-likelihood", xlab=expression(beta))
  if (!is.null(ci))
  {
    up <- (attr(x, "aic")-2*length(beta))/2
    abline(h=-up)
    abline(h=-up-qchisq(p=ci,df=1)/2,col="red")
    I <- which(y>=-up-qchisq(p=ci,df=1)/2)
    lines(x[I],rep(-up-qchisq(p=ci,df=1)/2,length(I)),
          lwd=5,col="red")
    abline(v=range(x[I]),lty=2,col="red")
  }
}