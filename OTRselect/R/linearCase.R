linearCase <- function(x, 
                       y,  
                       trt,  
                       propen,  
                       wgt,  
                       intercept) {

  p <- 1L:{ncol(x)+1L}
  n <- length(y)

  trtOpts <- sort(unique(trt))
  trtOpts <- trtOpts[-1L]
  k <- length(trtOpts)

  x1 <- cbind(1.0,x) * wgt
  y <- y * wgt

  xpro <- NULL
  for( i in 1L:k ) {
    xpro <- cbind(xpro, x1 * {{trt == trtOpts[i]} - propen[,i+1L]})
  }

  if( intercept ) {
    fit <- stats::lm(y ~ x + xpro, data = data.frame(x,xpro,y))$coef
    gamma <- fit[p]
    beta <- fit[-p]
  } else {
    fit <- stats::lm(y ~ -1 + x + xpro, data = data.frame(x,xpro,y))$coef
    gamma <- 0.0
    beta <- fit
  }

  if( any(is.na(fit)) ) {
    stop("NAs encountered in fit of baseline mean function.",
         call. = FALSE)
  }

  AL <- adaptiveLasso(x = x,
                      y = y,
                      gamma = gamma,
                      beta = beta,
                      propen = propen,
                      trt = trt,
                      wgt = wgt)


  return(AL)
}
