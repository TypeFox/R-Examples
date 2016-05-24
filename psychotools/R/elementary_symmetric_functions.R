elementary_symmetric_functions <- function(par,
  order = 0L, log = TRUE, diff = FALSE, engine = NULL)
{
  ## check order, set poly and engine
  stopifnot(order %in% 0L:2L)
  poly <- is.list(par)
  if (is.null(engine)) engine <- if (order < 2L && log) "C" else "R"

  ## if necessary: create ncat and unlist par
  if (engine == "C" || poly) {
    ncat <- sapply(par, length)
    par <- unlist(par)
    
    ## sanity checks
    stopifnot(log)
    if(order == 2L) {
      if (poly) stop("Second order ESFs are not available for polytomous items.")
      else {
        warning("Second order ESFs are not available in C, changed to R.")
        engine <- "R"
      }
    }
  }

  if(engine == "R" & !poly) {
    ## Michelle Liou (1994). More on the Computation of Higher-Order
    ## Derivatives of the Elementary Symmetric Functions in the
    ## Rasch Model. Applied Psychological Measurement, 18(1), 53-62.

    ## Use difference algorithm for orders > 0?
    diff <- rep(diff, length.out = 2L)

    ## derivatives up to order
    order <- round(order)[1L]
    rval <- list()[1L:(order+1L)]
    names(rval) <- 0L:order

    ## transformations
    par <- as.vector(par)
    beta <- if(log) par else -log(par)
    eps <- exp(-beta)
    n <- length(eps)
    stopifnot(n > 1L)
  
    ## Order: 0  
    ## initialization: gamma_1(eps_1) = eps_1
    gamma <- c(eps[1L], numeric(n-1L))

    ## recursion: Equation 3
    for(i in 2L:n) gamma[1L:i] <- c(eps[i] + gamma[1L],
      eps[i] * gamma[(2L:i) - 1L] + gamma[2L:i])
  
    ## gamma_0 = 1
    gamma <- c(1, gamma)

    ## return value
    rval[[1L]] <- gamma
    if(order < 1L) return(rval)

    ## Order: 1
    if(diff[1L]) {
      ## initialization: gamma_1^(j) = 1
      gamma1 <- matrix(0, nrow = n+1, ncol = n)
      gamma1[2L,] <- 1
      ## recursion: Equation 4
      for(q in 3L:(n+1)) gamma1[q,] <- gamma[q-1L] - eps * gamma1[q-1L,]
    } else {
      ## re-call self, omitting i-th parameter
      gamma1 <- sapply(1L:n, function(i)
        c(0, elementary_symmetric_functions(eps[-i], order = 0L, log = FALSE, engine = "R")[[1L]]))
    }
    ## if input on log scale: include inner derivative
    if(log) gamma1 <- exp(t(t(log(gamma1)) - beta))
    ## return value
    rval[[2L]] <- gamma1
    if(order < 2L) return(rval)

    ## Order: 2
    if(diff[2L]) {
      ## initialization: gamma_2^(i,j) = 1
      gamma2 <- array(0, c(n+1L, n, n))
      gamma2[3L,,] <- 1
      ## auxiliary variables
      eps_plus_eps <- outer(eps, eps, "+")
      eps_times_eps <- exp(-outer(beta, beta, "+"))
      ## recursion: Jansen's Equation (Table 1, Forward, Second-Order)
      if(n > 2L) for(q in 4L:(n+1L)) gamma2[q,,] <- gamma[q-2L] -
        (eps_plus_eps * gamma2[q-1L,,] + eps_times_eps * gamma2[q-2L,,])
    } else {
      ## re-call self, omitting i-th and j-th parameter
      gamma2 <- array(0, c(n + 1L, n, n))
      for(i in 1L:(n-1L)) {
        for(j in (i+1L):n) {
          gamma2[, i, j] <- gamma2[, j, i] <- c(0, 0,
	    elementary_symmetric_functions(eps[-c(i, j)], order = 0L, log = FALSE, engine = "R")[[1L]])
        }
      }
      if(log) eps_times_eps <- exp(-outer(beta, beta, "+"))  
    }
    ## if input on log scale: include inner derivative
    if(log) for(q in 1L:(n+1L)) gamma2[q,,] <- eps_times_eps * gamma2[q,,]
    ## each diagonal is simply first derivative
    for(q in 2L:(n+1L)) diag(gamma2[q,,]) <- gamma1[q,]
    ## return value
    rval[[3L]] <- gamma2
    return(rval)
  }

  if(engine == "R" && poly) {
    ## ESF summation/difference algorithm for polytmous items.
    ## See: Fischer & Ponocny (FP), 1994 and 1995
    
    ## construct item parameter list and result vector
    m <- length(ncat)
    par <- split.default(par, rep.int(1L:m, ncat))
    eps <- lapply(par, function (x) c(1, exp(-x)))
    rval <- vector("list", (1L+order))
    names(rval) <- 0L:order

    ## zero derivative, summation algorithm
    ## setup result list, possible scores (r) with i items (i=2, ..., m)  and categories (ncat, including zero)
    r <- cumsum(ncat) + 1L
    rmax <- r[m] - 1L
    ncat <- ncat + 1L
    rxncat <- r*ncat
    rval[[1L]] <- vector("list", m)

    ## initialization: gamma_r(Item1) = eps_r, r = 0, ..., o_1
    rval[[1L]][[1L]] <- eps[[1L]]
  
    ## FP (1994), EQ (11); FP (1995) EQ (19.19)
    for (i in 2L:m) { 
      gamma <- rep.int(c(rval[[1L]][[i-1L]], rep.int(0, ncat[i])), ncat[i])[1L:rxncat[i]]
      dim(gamma) <- c(r[i], ncat[i])
      rval[[1L]][[i]] <- gamma %*% eps[[i]]
    }

    ## first derivative, preparations, setup result list: derivatives of first parameter for all items (m)
    if (order) {
    
      ## create result object
      rval[[2L]] <- lapply(1L:m, function (x) vector("list", x))

      if (!diff) {
        ## helper variables
        ## nOld: number of items in previous round
        ## r: possible scores with 2:m Items (cols) without item j (row)
        nOld <- 0L:(m-1L)
        r <- vapply(r, "-", numeric(m) , ncat) + 1L

        ## initialization, gamma^(1)_0(I1) = 1
        rval[[2L]][[1L]] <- 1
        
        ## FP (1994), EQ (12)/FP (1995) EQ (19.21)
        for (i in 2L:m) {
          for (j in 1L:nOld[i]) { #derivivatives of 'old' items
            gamma1 <- rep.int(c(rval[[2L]][[j]], rep.int(0, ncat[i])), ncat[i])[1L:(ncat[i] * r[j,i])]
            dim(gamma1) <- c(r[j,i], ncat[i])
            rval[[2L]][[j]] <- gamma1 %*% eps[[i]]
          }
          rval[[2L]][[i]] <- rval[[1L]][[nOld[i]]] #derivative of newly added item
        }
      } else{
        ## helper variables
        ## maximum score with all Items but without item j (+1 because zero is first element)
        r <- r[m] - ncat + 1L

        ## initialization, gamma^(i)_0(I1) = 1, i = 1, ..., m
        rval[[2L]][1L:m] <- 1
        
        ## FP (1994), EQ (14/13); FP (1995) EQ (19.22/19.23/19.24)    
        for (j in 1L:m) {
          for (i in 2L:r[j]) {
            rval[[2L]][[j]] <- c(rval[[2L]][[j]], rval[[1L]][[m]][i] - sum(rval[[2L]][[j]][(i-1L) : (i-min(i-1L, ncat[j]-1L))]
                                                                     * eps[[j]][2L:min(i, ncat[j])]))
          }
        }
      }
    
      ## include inner derivative (only for identified parameters) & construct result matrix
      eps_indent <- lapply(eps, function (x) x[-1L])
      rval[[2L]] <- lapply(rval[[2L]], function (x) if (length(x) < rmax) c(0, x, rep.int(0, rmax-length(x))) else c(0, x))
      rval[[2L]] <- mapply("%o%", rval[[2L]], eps_indent, SIMPLIFY=FALSE)
      rval[[2L]] <- do.call("cbind", rval[[2L]])

      ## shift derivatives to the correct position with helper function
      shift_val <- unlist(lapply(ncat - 2L, seq, from = 0L))
      shift_row <- function (x, y) if (y) c(rep.int(0, y), x[1L:(length(x)-y)]) else x
      for (j in seq_len(rmax)) rval[[2L]][, j] <- shift_row(rval[[2L]][, j], shift_val[j])
    }

    ## over and out
    rval[[1L]] <- rval[[1L]][[m]][, , drop=TRUE] # return only full zero gammas, rest was only for first order derivatives needed
    return(rval)
  }
  
  if(engine == "C") {
    rval <- .Call(esf, as.double(par), as.integer(ncat), as.integer(order), as.integer(diff))
    names(rval) <- 0L:order
    return(rval)
  }
}
