BeMRes <-
function(times, delta = rep(1, length(times)), alpha = rep(0.0001, K), 
                   beta = rep(0.0001, K), c.r = rep(0, K), type.c = 4, 
                   epsilon =  1, iterations = 1000, 
                   burn.in = floor(iterations * 0.2), thinning = TRUE, 
                   thpar = 5, printtime = TRUE) {  
  tInit <- proc.time()
  K <- max(times)
  tol = .Machine$double.eps ^ 0.5
  if (min(times) < 0) {
    stop ("Invalid argument: 'times' must be a nonnegative integer vector.")
  }
  for(i in 1:length(times)) {
    if (abs(times[i] - round(times[i])) > tol) {
      stop ("Invalid argument: 'times' must be a nonnegative integer vector.")
    } 
  }
  if (min((delta  ==  0) + (delta  ==  1 )) == 0) {
    stop ("Invalid argument: 'delta' must have 0 - 1 entries.")
  }
  if (length(times) != length(delta)) {
    stop ("Invalid argument: 'times' and 'delta' must have same length.")
  }
  if (length(alpha) != K || length(beta) != K) {
    stop (c("Invalid argument: 'alpha', 'beta', must have length "), K, ("."))
  }
  if (min(c(alpha, beta)) < 0) {
    stop ("Invalid argument: 'alpha' and 'beta' must have nonnegative entries.")
  } 
  if (abs(type.c - round(type.c)) > tol || type.c < 1 || type.c > 4) {
    stop ("Invalid argument: 'type.c' must be an integer between 1 and 4.")
  }
  if (type.c == 2) {
    if (length(c.r) != (K - 1)) {
      stop (c("Invalid argument: 'c.r' must have length, ", K - 1, "."))
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
      || burn.in > iterations * 0.9) {
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
  nm <- BeNM(times, delta)
  n <- nm$n
  m <- nm$m
  tao <- nm$tao
  t.unc <- nm$t.unc
  if (thinning == TRUE) {
    cols <- floor((iterations - burn.in) / thpar) + 1
  } else {
    cols <- iterations
  }  
  a <- 0
  if (type.c == 3 || type.c == 4) {
    c.r <- rep(5, (K - 1))
  }
  if (type.c == 4) {
    a <- 1
  }
  cat(c("Iterating...", "\n"), sep = "")
  X <- matrix(NA, nrow = (3 * K - 2 + a), ncol = cols)
  Pi.r <- rep(0.1, K)
  for(j in 1:iterations) {
    if (floor(j / 200) == ceiling(j / 200)) {
      cat(c("Iteration ", j, " of ", iterations, ".", "\n"), sep = "")
    }
    u.r <- BeUpdU(alpha, beta, c.r, Pi.r)
    Pi.r <- UpdPi(alpha, beta, c.r, u.r, n, m)
    if (type.c == 3 || type.c == 4) {
      if (type.c == 4) {
        epsilon <- rgamma(1, shape = 0.01 + K, scale = 1 / (0.01 + sum(c.r)))
      }
      c.r <- BeUpdC(alpha, beta, Pi.r, u.r, epsilon)
    }
    if (thinning == TRUE && j >= burn.in) {
      if (floor((j - burn.in) / thpar) == (j - burn.in) / thpar) {
        for(i in 1:K) {
          X[i, (j - burn.in) / thpar + 1] <- Pi.r[i]
        }
        for(i in (K + 1):(2 * K - 1)) {
          X[i, (j - burn.in) / thpar + 1] <- u.r[i - K]
        }
        for(i in (2 * K):(3 * K - 2)) {
          X[i, (j - burn.in) / thpar + 1] <- c.r[i - 2 * K + 1]
        }
        if (type.c == 4) {
          X[3 * K - 1, (j - burn.in) / thpar + 1] <- epsilon
        }
      }
    }
    if (thinning == FALSE) {
      for(i in 1:K) {
        X[i, j] <- Pi.r[i]
      }
      for(i in (K + 1):(2 * K - 1)) {
        X[i, j] <- u.r[i - K]
      }
      for(i in (2 * K):(3 * K - 2)) {
        X[i, j] <- c.r[i - 2 * K + 1]
      }
      if (type.c == 4) {
        X[3 * K - 1, j] <- epsilon
      }
    }
  }
  cat(c("Done.", "\n", "Generating survival function estimates... 0%", "\n"),
      sep = "")
  iterations <- dim(X)[2]
  S <- matrix(1, ncol = iterations + 1, nrow = K + 1)
  S[,1] <- 0:K
  prog <- 0
  prog.by <- 1/5
  for(it in 1:iterations) {
    if (it/iterations - prog >= prog.by) {
      prog <- prog + prog.by
      cat(c("Generating survival function estimates... ", prog * 100, "%",
            "\n"), sep = "")
    }
    for(k in 1:K) {
      S[k, it + 1] <- exp(-sum(X[1:k, it]))
    }
  }
  cat(c("Done.", "\n"), sep = "")
  if (printtime) {
    cat(">>> Total processing time (sec.):\n")
    print(procTime <- proc.time() - tInit)
  }
  out <- list(times = times, delta = delta, tao = tao, K = K, t.unc = t.unc,
              iterations = iterations, summary = X, S = S)
  out
  }
