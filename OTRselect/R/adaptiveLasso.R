adaptiveLasso <- function(x, 
                          y, 
                          gamma, 
                          beta, 
                          propen, 
                          trt, 
                          wgt) {

  n <- length(y)

  trtOpts <- sort(unique(trt))
  base <- trtOpts[1L]
  trtOpts <- trtOpts[-1L]
  k <- length(trtOpts)

  x1 <- cbind(1.0,x) * wgt
  y <- y * wgt

  xpro <- NULL
  for( i in 1L:k ) {
    xpro <- cbind(xpro, x1 * {{trt == trtOpts[i]} - propen[,i+1L]})
  }

  penmat <- sweep(x = xpro,
                  MARGIN = 2L,
                  STATS = abs(beta),
                  FUN = "*")

  if( length(gamma) == 1L ) {
    YTemp <- y - gamma * wgt
  } else {
    YTemp <- y - {gamma[1L] + x %*% gamma[-1L]} * wgt
  }

  fit <- try(lars::lars(x = penmat,
                        y = YTemp,
                        normalize = FALSE,
                        intercept = FALSE),
             silent = TRUE)

  if( is(fit,"try-error") ) {
    stop("Unable to obtain lars fit.", call. = FALSE)
  }

  if( any( is.na(coef(fit)) ) ) {
    stop("NAs returned by lars.", call. = FALSE)
  }

  betaAL <- sweep(x = coef(fit),
                  MARGIN = 2L,
                  STATS = abs(beta),
                  FUN = "*")

  sumFit <- summary(fit)

  rss <- sumFit$Rss

  bic <- rss/rss[length(rss)] + sumFit$Df*log(n)/n

  betaAL <- betaAL[which.min(bic),]
  betaAL <- matrix(data = betaAL, ncol = k)
  rownames(betaAL) <- c("(Intercept)",colnames(x))

  temp <- cbind(1.0,x) %*% betaAL
  optimal <- max.col(temp, ties.method="first")
  useBase <- sapply(1:n,function(i){temp[i,optimal[i]] <= 0.0})
  index <- trtOpts[optimal]
  index[useBase] <- base

  return(list("beta" = betaAL,
              "optTx" = index))
}
