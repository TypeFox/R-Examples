CGaMRes <-
function(data, type.t = 1, K = 5, alpha = rep(0.001, K.aux),
                    beta = rep(0.0001, K.aux), c.r = rep(0, K.aux - 1), 
                    type.c = 4, epsilon = 1, iterations = 1000, 
                    burn.in = floor(iterations * 0.2), thinning = TRUE,
                    thpar = 3, printtime = TRUE) {
  tInit <- proc.time()
  times <- data[, 1]
  delta <- data[, 2]
  covar <- as.matrix(data[, -c(1, 2)])
  if(mean(covar) == "NaN") {
    stop("'data' must contain at least one explanatory variable.")
  }
  if (min(times) < 0) {
    stop ("Invalid argument: 'times' must be a nonnegative vector.")
  }
  if (min((delta  ==  0) + (delta  ==  1 )) == 0) {
    stop ("Invalid argument: 'delta' must have 0 - 1 entries.")
  }
  if (length(times) != length(delta)) {
    stop ("Invalid argument: 'times' and 'delta' must have same length.")
  }
  if (type.t == 2) {
    K.aux <- ceiling(max(times))
    if (K != K.aux) {
      warning (c("type.t = 2 requires K=", K.aux, ". K(", K, ") fixed at ", 
                 K.aux, "."))
      K <- K.aux
    }
  }
  if (type.t == 1 || type.t == 3) {
    if (class(try(K != 0, TRUE)) == "try-error") {
      K.aux <- 5
      warning ("'K' value not specified. 'K' fixed at ", K.aux, ".")
    } else {K.aux <- K}
    K <- K.aux
  }
  tol <- .Machine$double.eps ^ 0.5
  if (abs(type.t - round(type.t)) > tol || type.t < 1 || type.t > 3) {
    stop ("Invalid argument: 'type.t' must be an integer between 1 and 3.")
  }
  if (K <= 2 || abs(K - round(K)) > tol) {
    stop ("Invalid argument: 'K' must be an integer greater than 2.")
  }
  if (length(alpha) != K || length(beta) != K) {
    stop (c("Invalid argument: 'alpha', 'beta', must have length "), K, 
          (". Note that, if type.t = 2 is selected, 'K' will be redefined as 
           ceiling(max(times))."))
  }
  if (min(c(alpha, beta)) < 0) {
    stop ("Invalid argument: 'alpha' and 'beta' must have nonnegative entries.")
  } 
  if (abs(type.c - round(type.c)) > tol || type.c < 1 || type.c > 4) {
    stop ("Invalid argument: 'type.c' must be an integer between 1 and 4.")
  }
  if (type.c == 1 || type.c == 2) {
    if (length(c.r) != (K - 1)) {
      stop (c("Invalid argument: 'c.r' must have length, ", K - 1, ". Note that,
              if type.t = 2, 'K' is redefined as ceiling(max(times)) = ", K,
              "."))
    }
    for(k in 1:(K - 1)) {
      if (abs(c.r[k] - round(c.r[k])) > tol || min(c.r) < 0) {
        stop ("Invalid argument: 'c.r' entries must be nonnegative integers.")
      }
    }
  }
  if (type.c == 1 && sum(abs(c.r)) != 0 ) {
    c.r <- rep(0, K - 1)
    warning (c("'c.r' redefined as rep(0,", K - 1, ") because type.c = 1."))
  }
  if ((type.c == 3 || type.c == 4) && epsilon < 0) {
    stop ("Invalid argument: 'epsilon' must be nonnegative.")
  }
  if (iterations <= 0 || abs(iterations - round(iterations)) > tol 
      || iterations < 50) {
    stop ("Invalid argument: 'iterations' must be an integer greater than 50.")
  }
  if (burn.in <= 0 || abs(burn.in - round(burn.in)) > tol 
      || burn.in > iterations*0.9) {
    stop ("Invalid argument: 'burn.in' must be a postitive integer smaller than 
          iterations = ", iterations * 0.9, ".")
  }
  if (thinning != TRUE && thinning != FALSE) {
    stop ("Invalid argument: 'thinning' must be a logical value.")
  }
  if (thpar <= 0 || abs(thpar - round(thpar)) > tol 
      || thpar > 0.1 * iterations) {
    stop ("Invalid argument: 'thpar' must be a postitive integer smaller than 
          iterations * 0.10 = ", iterations * 0.1, ".")
  }
  if (printtime != TRUE && printtime != FALSE) {
    stop ("Invalid argument: 'printtime' must be a logical value.")
  }
  tao <- GaTao(times, delta, type.t, K)
  n <- CGaN(times, delta, type.t, K, covar)
  t.unc <- sort(times[delta==1])
  if (thinning == TRUE) {
    cols <- floor((iterations - burn.in) / thpar) + 1
  } else {
    cols <- iterations
  }
  a <- 0
  if (type.c == 3) {
    c.r <- rep(5, (K - 1))
  }
  if (type.c == 4) {
    a <- 1
  }  
  cat(paste("Iterating...", "\n"), sep = "")
  p <- length(covar[1, ])
  theta <- rep(1, p)
  X <- matrix(NA, nrow = (3 * K - 2 + a + p), ncol = cols)
  lambda.r <- rep(0.1, K)
  for(j in 1:iterations) {
    if (floor(j / 200) == ceiling(j / 200)) {
      cat(c("Iteration ", j, " of ", iterations, ".", "\n"), sep = "")
    }
    u.r <- GaUpdU(alpha, beta, c.r, lambda.r)
    m <- CGaM(times, delta, type.t, K, covar, theta)
    lambda.r <- UpdLambda(alpha, beta, c.r, u.r, n, m)
    theta <- CUpdTheta(theta, lambda.r, times, delta, type.t, K, covar)
    if (type.c == 3 || type.c == 4) {
      if (type.c == 4) {
        epsilon <- rgamma(1, shape = 0.01 + K, scale = 1 / (0.01 + sum(c.r)))
      }
      c.r <- GaUpdC(alpha, beta, c.r, lambda.r, u.r, epsilon)
    }
    if (thinning == FALSE) {
      for(k in 1:K) {
        X[k, j] <- lambda.r[k]
      }
      for(k in (K + 1):(2 * K - 1)) {
        X[k, j] <- u.r[k - K]
      }
      for(k in (2 * K):(3 * K - 2)) {
        X[k, j] <- c.r[k - 2 * K + 1]
      }
      if(type.c == 4){
        if (thinning == FALSE) {
          X[3 * K - 1, j] <- epsilon
        }
      }
      for(k in (3 * K - 1 + a):(3 * K - 2 + a + p)) {
        X[k, j] <- theta[k - 3 * K + 1 - a + 1]
      }      
    }
    if (thinning == TRUE && j >= burn.in) {
      if (floor((j - burn.in) / thpar) == (j - burn.in) / thpar) {
        for(k in 1:K) {
          X[k, (j - burn.in) / thpar + 1] <- lambda.r[k]
        }
        for(k in (K + 1):(2 * K - 1)) {
          X[k, (j - burn.in) / thpar + 1] <- u.r[k - K]
        }
        for(k in (2 * K):(3 * K - 2)) {
          X[k, (j - burn.in) / thpar + 1] <- c.r[k - 2 * K + 1]
        }
        if (type.c == 4) {
          X[3 * K - 1, (j - burn.in) / thpar + 1] <- epsilon
        }
        for(k in (3 * K - 1 + a):(3 * K - 2 + a + p)) {
          X[k, (j - burn.in) / thpar + 1] <- theta[k - 3 * K + 1 - a + 1]
        }
      }
    }
  }
  cat(c("Done.", "\n", "Generating survival and cumulative hazard function estimates... 0%", "\n"), sep = "")
  iterations <- dim(X)[2]
  S <- matrix(1, ncol = iterations + 1, nrow = 101)
  H <- matrix(0, ncol = iterations + 1, nrow = 101)
  for(i in 0:100) {
    S[i + 1, 1] <- H[i + 1, 1] <- max(tao) * i / 100    
  }
  prog <- 0
  prog.by <- 1/5
  for(it in 1:iterations) {
    if (it/iterations - prog >= prog.by) {
      prog <- prog + prog.by
      cat(c("Generating survival and cumulative hazard function estimates... ", prog * 100, "%", "\n"), sep = "")
    }
    for(i in 2:101) {
      aux <- 0
      aux2 <- 1
      k <- 1
      while(aux2 == 1) {
        if (tao[k + 1] < S[i, 1]) {
          aux <- aux + (tao[k + 1] - tao[k]) * X[k, it]
        }
        if (tao[k] < S[i, 1] && S[i, 1] <= tao[k + 1]) {
          aux <- aux + (S[i, 1] - tao[k]) * X[k, it]
          aux2 <- 0
        }
        k <- k + 1
      }
      H[i, it + 1] <- aux
      S[i, it + 1] <- exp(-aux)
    }
  }
  cat(paste("Done.", "\n"), sep = "")
  if (printtime) {
    cat(">>> Total processing time (sec.):\n")
    print(procTime <- proc.time() - tInit)
  }
  out <- list(times = times, delta = delta, covar = covar, type.t = type.t, 
              tao = tao, K = K, t.unc = t.unc, iterations = iterations, 
              summary = X, S = S, H = H, p = p)
  return(out)
}
