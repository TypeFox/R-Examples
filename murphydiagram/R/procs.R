# Input checker needed below
input.check <- function(f, y, t, alpha){
  if ( (length(t) > 1) | (length(alpha) > 1) | (length(f) != length(y)) | any(c(!is.vector(f), !is.vector(y), !is.vector(t), !is.vector(alpha))) ) stop("invalid input") 
}

# Extremal score for quantiles
S.quantile <- function(f, y, t, alpha){
  input.check(f, y, t, alpha)
  ((y < f) - alpha) * ((t < f) - (t < y))
}

# Take the positive part of a vector
pos <- function(x) x*(x >= 0)

# Extremal score for expectiles
S.expectile <- function(f, y, t, alpha){
  input.check(f, y, t, alpha)
  c1 <- abs((y < f) - alpha) 
  c2 <- pos(y-t) - pos(f-t) - (y-f)*(t<f)
  out <- c1 * c2
  out
}
  
# Newey-West VCV matrix
# u: vector of data
# prewhite: logical, should prewhitening be done?
# k: truncation lag for autocorrelations. If set to NULL, will be chosen automatically.
# meth: either "qs" (Quadratic Spectral, Andrews 1991) or anything else (Bartlett Kernel, Newey/West)
vHAC <- function(u, prewhite = FALSE, k = NULL, meth = "qs"){
  
  if (!is.matrix(u)) u <- as.matrix(u)
  
  n <- nrow(u)
  nreg <- ncol(u)
  rho <- sigma <- matrix(0, nreg, 1)
  
  # do a VAR(1) prewhitening
  if (prewhite == TRUE){
    reg.x <- matrix(u[1:(n-1), ], n-1, nreg)
    reg.y <- matrix(u[2:n, ], n-1, nreg)
    aux <- lm(reg.y~reg.x-1)
    beta <- matrix(unname(aux$coefficients), nreg, nreg)
    v <- matrix(unname(aux$residuals), n-1, nreg)
  } else {
    v <- u
    beta <- matrix(0, nreg, nreg)
  }
  nv <- nrow(v)
  
  # choose nr of lags (if not provided)
  if (is.null(k)){
    
    for (i in 1:nreg){
      aux <- lm(v[2:nv, i]~v[1:(nv-1), i]-1)
      rho[i] <- unname(aux$coefficients)
      sigma[i] <- sum(unname(aux$residuals)^2) / nv
    }
    
    if (meth == "qs"){
	  # See Eq. (6.4) on page 835 of Andrews (1991) -> Note that his sigma^2 corresponds to our sigma
      top <- sum( (4*(rho^2) * (sigma^2)) / ((1-rho)^8) )
      bot <- sum( (sigma^2) / ((1-rho)^4) )
      k <- ceiling(1.3221*((top/bot)*n)^(0.2))    
    } else {
      top <- sum( (4*(rho^2) * (sigma^2)) / (((1-rho)^6)*((1+rho)^2)) )
      bot <- sum( (sigma^2) / ((1-rho)^4) )
      k <- min(c(ceiling(1.1447*((top/bot)*n)^(1/3)), round(0.5*n)))
    }
    
  }
  
  # compute HAC
  vcv <- (t(v) %*% v) / (n-1)
  
  if (k > 0){
    if (meth == "qs"){
      del <- ((6 * pi)/(5 * k)) * (1:(n-1))
      w <- 3 * (sin(del) / del - cos(del)) / (del^2)  
      if (prewhite == FALSE){
        mlag <- n - 1 
      } else {
        mlag <- nv - 1
      }
    } else {
      w <- 1 - (1:k)/(k+1)
      mlag <- k
    }
    for (i in 1:mlag){
      cov <- t(v[(i+1):nv, , drop = FALSE]) %*% (v[1:(nv-i), , drop = FALSE]) / (n-1)
      vcv <- vcv + w[i]*(cov + t(cov))
    }
  }
  
  d <- solve(diag(nreg) - t(beta))
  hac <- d %*% vcv %*% t(d)
  
  return(list(hac = hac, k = k))  
  
}

get_grid <- function(f1, f2, y, functional, alpha, literal = TRUE, lhdiff = 1e-10){
  if (literal == TRUE){
    # Choose range of t's (see Corollary 2a, 2b in paper)
    if (functional == "quantile"){
      tseq <- sort(unique(c(f1, f2, y)))
    } else if (functional == "expectile") {
	  # NOTE: Use y's for Murphy diagrams also in case alpha = 1/2, even if not strictly needed
      aux1 <- sort(unique(c(f1, f2, y)))
      aux2 <- sort(unique(c(f1, f2))) - lhdiff
      tseq <- sort(unique(c(aux1, aux2)))   
    }
  } else {
    # Simply choose equally spaced grid of fixed length
    aux <- c(f1, f2, y)
    tseq <- seq(from = min(aux) - 0.1*sd(aux), to = max(aux) + 0.1*sd(aux), length.out = 100)
  }
  tseq
}

murphydiagram <- function(f1, f2, y, functional = "expectile", alpha = 0.5, labels = c("Method 1", "Method 2"), 
						  colors = NULL, equally_spaced = FALSE){
  cex.gen <- 1.6
  # Define function for extremal score
  if (functional == "expectile"){
    g <- function(f, t) S.expectile(f, y, t, alpha)
  } else if (functional == "quantile"){
    g <- function(f, t) S.quantile(f, y, t, alpha)
  } else {
    stop("Please choose either expectile or quantile functional")
  }
  # Unless specified otherwise: Use colors as in paper
  if (is.null(colors)) colors <- c("#D55E00", "#56B4E9", "#000000")
  # Grid of theta values
  tseq <- get_grid(f1, f2, y, functional, alpha, literal =  1 - equally_spaced)
  # Data frame with score entries for all theta values
  df <- data.frame(tseq = tseq, s1 = numeric(length(tseq)), s2 = numeric(length(tseq)))
  for (j in 1:length(tseq)){
    aux1 <- g(f1, tseq[j])
    aux2 <- g(f2, tseq[j])
    df[j, 2:3] <- c(mean(aux1), mean(aux2))
  }
  # Plot: Scores for both methods
  if (all(y %in% c(0, 1))){
    xx <- c(-0.05, 1.05)
  } else {
    xx <- c(min(tseq) - 0.1, max(tseq) + 0.1)
  }
  matplot(x = tseq, y = df[,2:3], type = "l", lty = 1, lwd = 4, xlab = expression(paste("Parameter ", theta)), ylab = "", bty = "n", cex.lab = cex.gen, 
          cex.axis = cex.gen, xlim = xx, ylim = c(0, 1.2*max(df[,2:3])), col = colors)
  abline(h = 0, lty = 2)
  if (!is.null(labels)) legend("top", labels, col = colors, lwd = 4, bty = "n", horiz = TRUE, cex = 1.2)
  
}

murphydiagram_diff <- function(f1, f2, y, functional = "expectile", alpha = 0.5, equally_spaced = FALSE, lag_truncate = 0, conf_level = 0.95){
  cex.gen <- 1.6
  # Some variables
  nobs <- length(y)
  scl <- abs(qnorm(0.5*(1-conf_level)))
  # Define function
  if (functional == "expectile"){
    g <- function(f, t) S.expectile(f, y, t, alpha)
  } else if (functional == "quantile"){
    g <- function(f, t) S.quantile(f, y, t, alpha)
  } else {
    stop("Please choose either expectile or quantile functional")
  }
  # Choose range of t's
  tseq <- get_grid(f1, f2, y, functional, alpha, literal = 1 - equally_spaced)
  df <- data.frame(tseq = tseq, s1 = numeric(length(tseq)), s2 = numeric(length(tseq)), lb = numeric(length(tseq)), ub = numeric(length(tseq)))
  for (j in 1:length(tseq)){
    aux1 <- g(f1, tseq[j])
    aux2 <- g(f2, tseq[j])
    df[j, 2:3] <- c(mean(aux1), mean(aux2))
    aux.v <- try(vHAC(aux1 - aux2, k = lag_truncate, meth = "bartlett")$hac/nobs)
    # HAC estimator won't invert for some t in the tails -> Set variance to zero in these cases
    if (class(aux.v) == "try-error") aux.v <- 0
    df[j, 4:5] <- mean(aux1-aux2) + c(-1, 1)*scl*sqrt(aux.v)
  }
    
  # Plot: Score difference + confidence interval
  if (all(y %in% c(0, 1))){
    xx <- c(-0.05, 1.05)
  } else {
    xx <- c(min(tseq) - 0.1, max(tseq) + 0.1)
  }
  matplot(x = df$tseq, y = df[, 4:5], type = "n", ylab = "", xlab = expression(paste("Parameter ", theta)),
          bty = "n", cex.axis = cex.gen, cex.lab = cex.gen, xlim = xx, col = 1)
  polygon(c(df$t, rev(df$t)), c(df$ub, rev(df$lb)), col = "grey", border = NA)
  lines(x = df$t, y = (df$s1-df$s2), type = "l", col = 1, lwd = 2.5)
  abline(h = 0, lty = 2)
  
}

##############################################################################
#
# Helper functions and analytical expressions for simulation study
#
##############################################################################

# cdf and quantile function for unfocused forecaster (forecaster 3)
p3p <- function(q, tau) 0.5*(pnorm(q) + pnorm(q - tau))
p3m <- function(q, tau) 0.5*(pnorm(q) + pnorm(q + tau))
q3p <- function(p, tau) {
  if (p == 0) {
    return(-Inf)
  } else if (p == 1) {
    return(Inf)
  } else {
    optimize(f = function(x) abs(p - p3p(x, tau)), interval = c(0, tau) + qnorm(p))$minimum
  }
}
q3m <- function(p, tau) {
  if (p == 0) {
    return(-Inf)
  } else if (p == 1) {
    return(Inf)
  } else {
    optimize(f = function(x) abs(p - p3m(x, tau)), interval = c(0, -tau) + qnorm(p))$minimum
  }
}  

# analytic expected scores for mean functional
c_theta <- function(t) sqrt(2) * dnorm(t/sqrt(2)) + t * pnorm(t/sqrt(2))

aux_mean <- list()
aux_mean[[1]] <- function(t) 0.5*(c_theta(t) - t * pnorm(t) - dnorm(t)) # perfect (forecaster 1)
aux_mean[[2]] <- function(t) 0.5*(c_theta(t) - t * (t >= 0)) # climatological (forecaster 2)
aux_mean[[3]] <- function(t) { # unfocused (forecaster 3)
  tau <- 2
  c1 <- c_theta(t)
  c2 <- t*(pnorm(t-0.5*tau) + pnorm(t+0.5*tau)) + dnorm(t-0.5*tau) + dnorm(t+0.5*tau)
  return(0.5*(c1 - 0.5 * c2))
}
aux_mean[[4]] <- function(t) 0.5*(c_theta(t) - t * pnorm(t) + dnorm(t)) # sign-reversed (forecaster 4)

expected_score_mean <- function(theta, forecaster = "P"){
  forecaster_names <- c("P", "C", "U", "SR")
  if (! forecaster %in% forecaster_names) stop("Forecaster must be either P, C, U, or SR")
  ind <- which(forecaster_names == forecaster)
  aux_mean[[ind]](t = theta)
}


# analytic expected scores for quantile / prob functional
aux_quant <- list()

aux_quant[[1]] <- function(t, a) { # perfect (forecaster 1)
  bound <- t - qnorm(a)
  if (is.finite(bound)) {
    c1 <- -a * pnorm(t / sqrt(2))
    c2 <- a * pnorm(bound)
    c3 <- integrate(function(x) pnorm(t-x) * dnorm(x), bound, Inf)$value
    return(c1 + c2 + c3)
  } else {
    return(0)
  }
}
aux_quant[[2]] <- function(t, a) { # climatological (forecaster 2)
  x <- pnorm(t/sqrt(2))
  if (x < a) {
    return(x * (1 - a))
  } else {
    return(a * (1 - x))
  }
}
aux_quant[[3]] <- function(t, a) { # unfocused (forecaster 3)
  tau <- 2
  bound1 <- t - q3p(a, tau)
  bound2 <- t - q3m(a, tau)
  if (is.finite(bound1)) {
    c1 <- -a * pnorm(t / sqrt(2))
    c2 <- a * (pnorm(bound1) + pnorm(bound2))
    c31 <- integrate(function(x) pnorm(t-x)*dnorm(x), bound1, Inf)$value
    c32 <- integrate(function(x) pnorm(t-x)*dnorm(x), bound2, Inf)$value
    return(c1 + 0.5 * (c2 + c31 + c32))
  } else {
    return(0)
  }
}
aux_quant[[4]] <- function(t, a) { #  sign-reversed (forecaster 4)
  bound <- qnorm(a) - t
  if (is.finite(bound)) {
    c1 <- -a * pnorm(t / sqrt(2))
    c2 <- a * pnorm(-bound)
    c3 <- integrate(function(x) pnorm(t-x)*dnorm(x), -Inf, bound)$value
    return(c1 + c2 + c3)
  } else {
    return(0)
  }
}

expected_score_quantile <- function(theta, alpha, forecaster = "P") {
  forecaster_names <- c("P", "C", "U", "SR")
  if (! forecaster %in% forecaster_names) stop("Forecaster must be either P, C, U, or SR")
  if (any(alpha < 0 | alpha > 1)) stop("Please provide quantile levels between zero and one")
  ind <- which(forecaster_names == forecaster)
  E <- function(theta, alpha) aux_quant[[ind]](t = theta, a = alpha)
  if (length(alpha) >= 1 & length(theta) == 1) {
    out <- sapply(alpha, E, t = theta)
  } else if (length(theta) > 1 & length(alpha) == 1) {
    out <- sapply(theta, E, a = alpha)
  } else {
    out <- NA
  }
  return(out)
}

fluctuation_test <- function(loss1, loss2, mu = 0.5, dmv_fullsample = TRUE,
                             lag_truncate = 0, time_labels = NULL,
                             conf_level = 0.05){
  
  # Some input checks
  if (length(loss1) != length(loss2)){
    stop("Vectors of losses must have the same length")
  }
  if (all(abs(seq(from = 0.1, to = 0.9, by = 0.1) - mu) > 1e-12)){
    stop("mu must be in {0.1, 0.2, ..., 0.9}")
  }
  if (!lag_truncate %in% 0:5){
    stop("lag_truncate must be in {0, 1, .., 5}")
  }
  if (!is.null(time_labels) & length(time_labels) != length(loss1)){
    warning("Specified time labels are inconsistent - simple integers used instead")
    time_labels <- NULL
  }
  if (!conf_level %in% c(0.05, 0.1)){
    stop("significance_level must be either 0.05 or 0.1")
  }
  
  # Critical values (two-sided test, from Table 1 of Giacomini/Rossi 2010)
  # Second/third column correspond to 5/10 percent confidence level
  CV <- cbind(seq(from = 0.1, to = 0.9, by = 0.1), 
              c(3.393, 3.179, 3.012, 2.890, 2.779, 2.634, 2.560, 2.433, 2.248),
              c(3.170, 2.948, 2.766, 2.626, 2.500, 2.356, 2.252, 2.130, 1.950))
  
  # Wrapper for HAC variance
  vHAC2 <- function(ld, lag_truncate){
    vHAC(ld, k = lag_truncate, meth = "Bartlett")$hac
  }
  
  # Size for fonts etc
  cex_gen <- 1.6
  
  # Window length m
  P <- length(loss1) # Notation as in GR2010
  m <- round(mu*P)
  
  # Loss differences
  ld <- loss1 - loss2
  
  # Vector of statistics
  dm_num <- dm_den <- rep(0, P - m + 1)
  
  # Loop over time
  for (jj in m:P){
    ind <- which(m:P == jj)
    ld_tmp <- ld[(jj-m+1):jj]
    dm_num[ind] <- mean(ld_tmp)
    dm_den[ind] <- sqrt(vHAC2(ld_tmp, lag_truncate)/m)
  }
  
  # Construct statistics
  
  # Variant 1: Use full sample for HAC variance
  s2hat <- vHAC2(ld, lag_truncate)
  dm1 <- sqrt(m)*dm_num/sqrt(s2hat)
  
  # Variant 2: Use rolling samples
  dm2 <- dm_num/dm_den
  
  # Choose statistic
  if (dmv_fullsample){
    dm_final <- dm1 
  } else {
    dm_final <- dm2
  }
  
  # Select critical value
  if (conf_level == 0.05){
    CVs <- CV[abs(CV[, 1] - mu) < 1e-12, 2] * c(-1, 1)  
  } else if (conf_level == 0.1){
    CVs <- CV[abs(CV[, 1] - mu) < 1e-12, 3] * c(-1, 1)
  }
  
  # Plot
  plot(x = m:length(loss1), y = dm_final, ylim = c(-7, 7), bty = "n", ylab = "", 
       xlab = "Time (End of Rolling Window)",
       type = "l", col = "cornflowerblue", lwd = 2.5, axes = FALSE, cex.lab = cex_gen)
  axis(2, cex.axis = cex_gen)
  abline(h = CVs, lwd = 3.5)
  abline(h = 0, lwd = 1.8, lty = 2)
  
  if (is.null(time_labels)){
    time_labels <- 1:length(loss1)
  }
  inds <- floor(seq(from = m, to = length(time_labels), length.out = 5))
  axis(1, at = inds, labels = time_labels[inds], cex.axis = cex_gen)
  
  # Return two-element list which summarizes the results
  list(df = data.frame(time = time_labels[m:length(loss1)], dmstat = dm_final), CV = CVs)
  
}

# Wrapper function for extremal scores
extremal_score <- function(x, y, theta, functional = "expectile", alpha = 0.5){
  if (functional == "expectile"){
    S.expectile(f = x, y = y, t = theta, alpha = alpha)
  } else if (functional == "quantile"){
    S.quantile(f = x, y = y, t = theta, alpha = alpha)
  } else {
    stop("invalid input for functional")
  }
}

input.check2 <- function(x, y){
  if (length(y) != length(x) | !is.vector(x) | !is.vector(y)) stop("x and y must be vectors of the same length")
}

apl_score <- function(x, y, alpha = 0.5){
  input.check2(x, y)
  out <- (x - y)
  ind <- (out > 0)
  out[ind] <- (1-alpha) * out[ind]
  out[!ind] <- -alpha * out[!ind]
  out
}

ase_score <- function(x, y, alpha = 0.5){
  input.check2(x, y)
  out <- (x - y)
  ind <- (out > 0)
  out[ind] <- (1-alpha) * (out[ind]^2)
  out[!ind] <- alpha * (out[!ind]^2)
  out
}


