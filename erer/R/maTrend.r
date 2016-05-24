maTrend <- function(q, n = 300, nam.c, nam.d, simu.c = TRUE)
{
  # 1. Check inputs
  if (!inherits(q, "maBina")) {stop("Need an object from 'maBina()'.\n")}
  if (missing(nam.c)) {stop("Need a continous variable name'.\n")}
  if (identical(sort(unique(q$w$x[, nam.c])), c(0, 1))) {
    stop("nam.c must be a continuous variable.") }
  if (!missing(nam.d)) {
    if (!identical(sort(unique(q$w$x[, nam.d])), c(0, 1))) {
      stop("nam.d must be a binary variable.") } }
    
  xx <- q$w$x; b.est <- as.matrix(coef(q$w)); link <- q$w$family$link
  pfun <- switch(link, probit = pnorm, logit = plogis)
  dfun <- switch(link, probit = dnorm, logit = dlogis)
  result <- listn(q, nam.c)
  
  # 2. Compute probability for all and standard error
  # 2.1 Create a data matrix: simulation or original data
  if (simu.c) {
    mm <- matrix(colMeans(xx), ncol = ncol(xx), nrow = n, byrow = TRUE)
    colnames(mm) <- colnames(xx)
    ran <- range(xx[, nam.c])
    mm[, nam.c] <- seq(from = ran[1], to = ran[2], length.out = n)
    ndf <- n - ncol(xx) 
  } else {
    mm <- matrix(colMeans(xx), ncol = ncol(xx), nrow = nrow(xx), byrow = TRUE)
    colnames(mm) <- colnames(xx)
    mm[, nam.c] <- xx[, nam.c]
    mm <- mm[order(mm[, nam.c]), ]
    ndf <- q$w$df.residual  
  }

  # 2.2 Predicted prob and standard error
  pp <- pfun(mm %*% b.est)         
  dr <- diag(c(dfun(mm %*% b.est)), nrow = nrow(mm), ncol = nrow(mm)) %*% mm 
  va <- dr %*% vcov(q$w) %*% t(dr) 
  se <- sqrt(diag(va))
  t_value <- pp / se
  p_value <- 2 * (1 - pt(abs(t_value), ndf))
  trend <- data.frame(mm[, nam.c], pp, se, t_value, p_value) 
  colnames(trend) <- c(nam.c, 
    paste('all', c("pr", "se", "t_value", "p_value"), sep = "."))
  result <- listn(q, nam.c, mm, trend)
  
  # 4. Add probability by dummy variable (1 or 0) 
  if (!missing(nam.d)) {
    m1 <- mm; m1[, nam.d] <- 1
    m0 <- mm; m0[, nam.d] <- 0
    pp1 <- pfun(m1 %*% b.est)
    pp0 <- pfun(m0 %*% b.est)
    dr1 <- diag(c(dfun(m1 %*% b.est)), nrow = nrow(mm), ncol = nrow(mm)) %*% m1
    dr0 <- diag(c(dfun(m0 %*% b.est)), nrow = nrow(mm), ncol = nrow(mm)) %*% m0
    va1 <- dr1 %*% vcov(q$w) %*% t(dr1); se1 <- sqrt(diag(va1))
    va0 <- dr0 %*% vcov(q$w) %*% t(dr0); se0 <- sqrt(diag(va0))
    t_value1 <- pp1 / se1; p_value1 <- 2 * (1 - pt(abs(t_value1), ndf))
    t_value0 <- pp0 / se0; p_value0 <- 2 * (1 - pt(abs(t_value0), ndf))

    trend1 <- data.frame(mm[, nam.c], pp1, se1, t_value1, p_value1) 
    trend0 <- data.frame(mm[, nam.c], pp0, se0, t_value0, p_value0) 
    trend  <- data.frame(mm[, nam.c], pp, pp1, pp0)
    colnames(trend1) <- c(nam.c, paste(paste(nam.d, "_d1", sep = ""), 
      c("pr", "se", "t_value", "p_value"), sep = "."))
    colnames(trend0) <- c(nam.c, paste(paste(nam.d, "_d0", sep = ""), 
      c("pr", "se", "t_value", "p_value"), sep = "."))
    colnames(trend) <- c(nam.c, 
      paste(c("all", paste(nam.d, c("d1", "d0"), sep='_')), "pr", sep = "."))      
    result <- listn(q, nam.c, nam.d, mm, m1, m0, trend, trend1, trend0)
  }
 
  # 5. Output 
  class(result) <- "maTrend"; return(result)
}

# print and plot method for 'maTrend'
print.maTrend <- function(x, ...) {
  cat("\n===========================================")
  cat("\nCalculated Probability Matrix\n")
  cat("===========================================\n")
  print(dim(x$trend)); print(tail(x$trend))
}

plot.maTrend <- function(x, ...) { 
  pr <- x$trend; yvar <- toupper(as.character(x$q$w$formula)[2])
  if (is.null(x$nam.d)) {
    plot(pr[, 2] ~ pr[, 1], type="l", lty=1,
      ylim = range(pr[, 2]), 
      xlab = toupper(x$nam.c),
      ylab = paste("Probability (", yvar, " = 1)", sep=""), ... )
    grid()           
  } else {
    plot(pr[, 2] ~ pr[, 1], type="l", lty=1,
      ylim = range(pr[, 2:4]), 
      xlab = toupper(x$nam.c),
      ylab = paste("Probability (", yvar, " = 1)", sep=""), ... )
    lines(pr[, 3] ~ pr[, 1], type="l", lty=2)
    lines(pr[, 4] ~ pr[, 1], type="l", lty=3)
    abline(v=colMeans(x$q$w$x)[x$nam.c], lty=4, col = "red")
    grid() 
    legend(x = "topright", legend = colnames(pr)[-1], 
      lty = 1:3, bg = "white")
  }
}