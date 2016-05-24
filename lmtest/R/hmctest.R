hmctest <- function(formula, point = 0.5, order.by = NULL, simulate.p = TRUE,
  nsim = 1000, plot = FALSE, data = list())
{
  dname <- paste(deparse(substitute(formula)))
  
  if(!inherits(formula, "formula")) {
    X <- if(is.matrix(formula$x))
           formula$x
         else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
           formula$y
         else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }  

  k <- ncol(X)
  n <- nrow(X)
  if(point < 1) point <- floor(point*n)
  if (point > n - k | point < k) stop("inadmissable breakpoint")

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    X <- as.matrix(X[order(z),])
    y <- y[order(z)]
  }


  resi <- lm.fit(X,y)$residuals
  hmc <- sum(resi[1:point]^2)/sum(resi^2)

  if(plot)
  {
    stats <- c(0,cumsum(resi^2))/sum(resi^2)
    stats <- ts(stats, start = 0, frequency = n)
    plot(stats, xlab="fraction", ylab="Harrison-McCabe statistics", xaxs="i",
      yaxs="i")
    abline(0,1)
  }

  names(hmc) <- "HMC"
  if (simulate.p)
  {
    stat <- rep(0, nsim)
    for (i in 1:nsim) {
      x <- rnorm(n)
      x <- (x - mean(x))/sqrt(var(x))
      stat[i] <- sum(x[1:point]^2)/sum(x^2)
    }
    PVAL <- mean(stat <= hmc)
  }
  else
    PVAL <- NA

  RVAL <- list(statistic = hmc,
      method = "Harrison-McCabe test",
      p.value= PVAL,
      data.name=dname)

  class(RVAL) <- "htest"
  return(RVAL)
}

