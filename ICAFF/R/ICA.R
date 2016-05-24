ICA <-
function(cost, nvar, ncountries = 80, nimp = 10, maxiter = 100,
                lb = rep(-10, nvar), ub = rep(10, nvar),
                beta = 2, P_revolve = 0.3, zeta = 0.02, ...) {
  if (missing(cost) || !is.function(cost)) {
    stop("No cost function supplied")
  }
  if (missing(nvar)) {
    stop("Number of variables not specified.")
  }
  if(!is.numeric(lb) || !is.numeric(ub) || any(!is.finite(lb)) || any(!is.finite(ub)) ||
     any(lb >= ub)) {
    stop("The region is invalidly specified.")
  }
  eps <- .Machine$double.eps  ## constant
  ptm <- proc.time()
  lb <- rep(lb, length.out = nvar)  ## corrected
  ub <- rep(ub, length.out = nvar)  ## corrected
  ## ---------------------- initial population algorithm
  var <- rep(0, nvar)
  fit <- 0
  totalfit <- 0
  y <- list(var = var, fit = fit)
  colony <- rep(list(y), ncountries)

  f <- rep(0, ncountries)
  for(i in 1:ncountries) {
    colony[[i]]$var <- vari <- runif(nvar, min = lb, max = ub)
    colony[[i]]$fit <- f[i] <- cost(vari, ...)
  }

  k <- order(f)
  f[] <- 0
  colony <- colony[k]

  ## --------
  d1 <- round((ncountries - nimp)/nimp)
  d2 <- (ncountries - nimp) - (d1 * (nimp - 1))

  ## ---------
  imp <- lapply(colony[1:nimp],
                function(x) c(x, list(colonyi = list(),
                                      totalfit = totalfit)))
  colony3 <- colony[-(1:nimp)]
  
  ## ---------
  r <- 1:(ncountries - nimp)
  r <- sample(r)
  w <- 0
  for (i in 1:(nimp - 1)) {
    for (j in 1:d1) {
      w <- w + 1
      imp[[i]]$colonyi[[j]] <- colony3[[r[w]]]
    }
  }
  for (i in 1:d2) {
    w <- w + 1
    imp[[nimp]]$colonyi[[i]] <- colony3[[r[w]]]
  }
  ## =========================================================Totla fitness
  nimp <- length(imp)
  k <- rep(0, nimp)
  for (i in 1:(nimp - 1)) {
    for (j in 1:d1) {
      k[i] <- k[i] + imp[[i]]$colonyi[[j]]$fit
    }
    k[i] <- k[i]/d1
    imp[[i]]$totalfit <- imp[[i]]$fit + zeta * k[i]
  }
  for (j in 1:d2) {
    k[nimp] <- k[nimp] + imp[[nimp]]$colonyi[[j]]$fit
  }
  k[nimp] <- k[nimp]/d2
  imp[[nimp]]$totalfit <- imp[[nimp]]$fit + zeta * k[nimp]
  ## ====================================================

  p <- sapply(imp, `[[`, "fit")
  k <- which.min(p)
  gimp <- imp[[k]]
  
  ## ==================================================== main loop algorithm
  BEST <- matrix(rep(0, maxiter), ncol = 1)
  for (iter in 1:maxiter) {
    ## =======================================Assimilation
    nimp <- length(imp)
    for (i in 1:nimp) {
      ncolony <- length(imp[[i]]$colonyi)
      for (j in 1:ncolony) {
        d <- imp[[i]]$var - imp[[i]]$colonyi[[j]]$var
        d <- d * runif(nvar) * beta
        imp[[i]]$colonyi[[j]]$var <- imp[[i]]$colonyi[[j]]$var + d
        for (k in 1:nvar) {
          imp[[i]]$colonyi[[j]]$var[k] <- max(imp[[i]]$colonyi[[j]]$var[k], lb[k])
          imp[[i]]$colonyi[[j]]$var[k] <- min(imp[[i]]$colonyi[[j]]$var[k], ub[k])
        }
        imp[[i]]$colonyi[[j]]$fit <- cost(imp[[i]]$colonyi[[j]]$var)
      }
    }
    ## =======================================================Revolution
    nimp <- length(imp)
    for (i in 1:nimp) {
      ncolony <- length(imp[[i]]$colonyi)
      for (j in 1:ncolony) {
        if (runif(1, 0, 1) < P_revolve) {
          k <- sample(nvar, 1)
          d <- ub[k] - lb[k]
          d <- 0.1 * runif(1, -1:1) * d
          imp[[i]]$colonyi[[j]]$var[k] <- imp[[i]]$colonyi[[j]]$var[k] + d
          for (z in 1:nvar) {
            imp[[i]]$colonyi[[j]]$var[z] <- max(imp[[i]]$colonyi[[j]]$var[z], lb[z])
            imp[[i]]$colonyi[[j]]$var[z] <- min(imp[[i]]$colonyi[[j]]$var[z], ub[z])
          }
          imp[[i]]$colonyi[[j]]$fit <- cost(imp[[i]]$colonyi[[j]]$var)
        }
      }
    }
    ## ========================================================Exchange
    nimp <- length(imp)
    for (i in 1:nimp) {

      p <- sapply(imp[[i]]$colonyi, `[[`, "fit")
      index <- which.min(p)
      value1 <- p[index]

      if (value1 < imp[[i]]$fit) {
        bestcolony <- imp[[i]]$colonyi[[index]]
        imp[[i]]$colonyi[[index]]$var <- imp[[i]]$var
        imp[[i]]$colonyi[[index]]$fit <- imp[[i]]$fit
        imp[[i]]$var <- bestcolony$var
        imp[[i]]$fit <- bestcolony$fit
      }
    }
    ## ===============================================imperialistic competition
    nimp <- length(imp)
    if (nimp >= 2) {

      p <- sapply(imp, `[[`, "totalfit")
      index1 <- which.max(p)
      wimp <- imp[[index1]]

      p <- sapply(wimp$colonyi, `[[`, "fit")
      index2 <- which.max(p)
      wcolony <- wimp$colonyi[[index2]]

      l <- length(imp[[index1]]$colonyi)
      if (index2 == l) {
        length(imp[[index1]]$colonyi) <- length(imp[[index1]]$colonyi) - 1
      }
      if (index2 != l) {
        l <- l - 1
        for (j in index2:l) {
          imp[[index1]]$colonyi[[j]] <- imp[[index1]]$colonyi[[j + 1]]
        }
        length(imp[[index1]]$colonyi) <- length(imp[[index1]]$colonyi) - 1
      }
      l2 <- length(imp)

      p <- sapply(imp[1:l2], `[[`, "totalfit")
      
      p[index1] <- 0
      p <- max(p) - p - eps  # eps=2.2204e-16
      p <- p/sum(p)
      p <- cumsum(p)
      d <- runif(1, 0, 1)
      for (i in 1:l2) {
        if ((d < p[i]) || (d <- p[i])) {
          k <- i
          break
        }
      }
      ## ------------------
      n1 <- length(imp[[k]]$colonyi)
      length(imp[[k]]$colonyi) <- length(imp[[k]]$colonyi) + 1
      n1 <- n1 + 1
      imp[[k]]$colonyi[[n1]]$var <- wcolony$var
      imp[[k]]$colonyi[[n1]]$fit <- wcolony$fit
      n <- length(imp[[index1]]$colonyi)
      if (n == 0) {
        m <- nimp
        if (index1 == m) {
          length(imp) <- length(imp) - 1
        }
        if (index1 != m) {
          m <- m - 1
          for (j in index1:m) {
            imp[[j]] <- imp[[j + 1]]
          }
          length(imp) <- length(imp) - 1
        }
        l3 <- length(imp)

        p <- sapply(imp, `[[`, "totalfit")

        p <- max(p) - p + eps
        p <- p/sum(p)
        p <- cumsum(p)
        d <- runif(1, 0, 1)
        for (i in 1:l3) {
          if ((d < p[i]) || (d <- p[i])) {
            k <- i
            break
          }
        }
        ## ----------------
        n2 <- length(imp[[k]]$colonyi)
        length(imp[[k]]$colonyi) <- length(imp[[k]]$colonyi) + 1
        n2 <- n2 + 1
        imp[[k]]$colonyi[[n2]]$var <- wimp$var
        imp[[k]]$colonyi[[n2]]$fit <- wimp$fit
      }
      ## n == 0
    }
    ## nimp>=2 ################################ End of imperialistic competition
    nimp <- length(imp)

    p <- sapply(imp, `[[`, "fit")
    index <- which.min(p)
    value2 <- p[index]
    
    if (value2 < gimp$fit) {
      gimp <- imp[[index]]
    }
    BEST[iter, 1] <- gimp$fit
    nimp <- length(imp)
    if (nimp == 1) {
      for(i in (iter+1):maxiter){BEST[i,1]=BEST[iter,1]}
      break
    }
  }

##### End results algorithm nBEST
  obj <- list(call = match.call(),
              position = gimp$var,
              value = gimp$fit,
              nimp = nimp,
              trace = BEST,
              time = proc.time() - ptm)
  class(obj) <- "ICA"
  return(obj)
}
