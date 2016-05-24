################################################################################
#       Function for sample size computation in Non Exchangeable Case          #
#       Decision Rules: At least r endpoints significant among m               #
#       Authors: Delorme P., Lafaye de Micheaux P., Liquet B. and Riou J.      #
#       Date: Le 19/09/2015   Version 1                                        #
################################################################################

## The function .indive.rm.ssc() is only for a weak control of the type-II q-gFWER, namely called with p = m.
## We will have to generalize it to strong control.
## We will have to deal with the unilateral and bilateral tests!!!


.indivne.rm.ssc <- function(method, asympt, r, m, p, nCovernE, delta, SigmaC, SigmaE, R, power, alpha, interval, q, maxpts, abseps, releps, nbcores, LB, orig.Hochberg) {
    if (method == "Bonferroni") {
        result <- uniroot.integer(function(nE) .Psirmsne(r = r, m = m, p = m, nE = as.integer(nE), nCovernE = nCovernE, delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R, alpha = alpha, q = q, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps, nbcores = nbcores, LB = LB)$pow - power, interval, pos.side = TRUE)$root     
    }
    if (method == "Hochberg") {
        result <- uniroot.integer(function(nE) .Psirmune(r = r, m = m, p = m, nE = as.integer(nE), nCovernE = nCovernE, delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R, alpha = alpha, q = q, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps, nbcores = nbcores, LB = LB, orig.Hochberg = orig.Hochberg)$pow - power, interval, pos.side = TRUE)$root 
    }
    if (method == "Holm") {
        result <- uniroot.integer(function(nE) .Psirmdne(r = r, m = m, p = m, nE = as.integer(nE), nCovernE = nCovernE, delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R, alpha = alpha, q = q, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps, nbcores = nbcores, LB = LB)$pow - power, interval, pos.side = TRUE)$root
    }
    return(result)  
}

################################################################################
# Computation of the power in Theorem 2.1
################################################################################

# Function for the computation of all possible combinations of size m taken among the elements of the vector x
.combn.mod <- function (x, m, FUN = NULL, simplify = TRUE, ...)  {
    stopifnot(length(m) == 1L)
      
    if (m < 0) {
        stop("m < 0")
    }
    if (m == 0) {
        return(if (simplify) vector(mode(x), 0L) else list())
    }
      
    n <- length(x)
    
    if (n < m) {
        stop("n < m")
    }
    
    m <- as.integer(m)
    e <- 0
    h <- m
    a <- 1L:m
      
    nofun <- is.null(FUN)
    
    if (!nofun && !is.function(FUN)) {
        stop("'FUN' must be a function or NULL.")
    }
      
    len.r <- length( r <- if (nofun) { x[a] } else { FUN(x[a], ...) })
      
    count <- as.integer(round(choose(n, m)))
      
    if (simplify) {
        dim.use <- if (nofun)
                       c(m, count)
                   else {
                       d <- dim(r)
                       if (length(d) > 1L)
                           c(d, count)
                       else if (len.r > 1L)
                           c(len.r, count)
                       else c(d, count)
                   }
    }
      
    if (simplify) {
        out <- matrix(r, nrow = len.r, ncol = count)
    } else {
        out <- vector("list", count)
        out[[1L]] <- r
    }
      
    i <- 2L
    nmmp1 <- n - m + 1L
      
    while (a[1L] != nmmp1) {
        if (e < n - h) {
            h <- 1L
            e <- a[m]
            j <- 1L
        } else {
            e <- a[m - h]
            h <- h + 1L
            j <- 1L:h
        }
          
        a[m - h + j] <- e + j
        r <- if (nofun)
                 x[a]
             else FUN(x[a], ...)
        if (simplify)
            out[, i] <- r
        else out[[i]] <- r
        i <- i + 1L
    }
    
    if (simplify)
        array(out, dim.use)
    else out
}

.Psirmsne <- function(r, m, p = m, nE, nCovernE, delta, SigmaC, SigmaE, R, alpha, q, asympt, maxpts, abseps, releps, nbcores, LB) {
  		
    nC <- as.integer(nCovernE * nE)
    df <- df.compute(nE, nC, SigmaE, SigmaC)
    
    # Computation of the critical value cjm for the Bonferroni method (Example 2.3 p. 3).
    # q-gFWER control
    if (asympt) {
        cjm <- rep(qnorm(1 - q * alpha / m, lower.tail = TRUE), m)
    } else {
        cjm <- rep(qt(1 - q * alpha / m, df, lower.tail = TRUE), m)
    }
    sighat <- sqrt(diag(SigmaE) / nE + diag(SigmaC) / nC)
    effectsize <- delta / sighat
    # Because P[T_j > c_{jm} | H_1] = P[T_j > c_{jm} - effectsize | H_0]
    thresh <- cjm - effectsize

    if ((p - r + 1) > r) { # Theorem 2.1, p. 3, first formula.
        pow <- 1
        error <- 0
        for (k in (p - r + 1):p) {
            coef <- choose(k - 1, p - r)
            # We are under H_1:
            if (asympt) {
                if (nbcores == 1) {
                    tmp <- apply(.combn.mod(1:p, k, FUN = function(Jset) {prob <- pmvnorm(lower = rep(-Inf, k), upper = thresh[Jset], mean = rep(0.0, k), corr = R[Jset,Jset], maxpts = maxpts, abseps = abseps / coef, releps = releps) ; if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.") ; return(list(prob, attr(prob, "error")))}, simplify = TRUE), MARGIN = 1, FUN = function(x) sum(as.numeric(x)))
                } else {
                    cl <- parallel::makeCluster(getOption("cl.cores", nbcores))
                    clusterExport(cl = cl, varlist = c("k", "thresh", "R", "maxpts", "abseps", "coef", "releps"), envir = environment())
                    if (LB) {
                        tmp <- apply(parallel::parSapplyLB(cl, unlist(apply(.combn.mod(1:p, k)[,sample(choose(p, k)), drop = FALSE], MARGIN = 2, FUN= list), recursive = FALSE), FUN = function(Jset) {prob <- pmvnorm(lower = rep(-Inf, k), upper = thresh[Jset], mean = rep(0.0, k), corr = R[Jset,Jset], maxpts = maxpts, abseps = abseps / coef, releps = releps) ; if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.") ; return(list(prob, attr(prob, "error")))}), MARGIN = 1, FUN = function(x) sum(as.numeric(x)))
                    } else {
                        tmp <- apply(parallel::parApply(cl, .combn.mod(1:p, k)[,sample(choose(p, k)), drop = FALSE], MARGIN = 2, FUN = function(Jset) {prob <- pmvnorm(lower = rep(-Inf, k), upper = thresh[Jset], mean = rep(0.0, k), corr = R[Jset,Jset], maxpts = maxpts, abseps = abseps / coef, releps = releps) ; if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.") ; return(c(prob, attr(prob, "error")))}), MARGIN = 1, FUN = function(x) sum(as.numeric(x)))
                    }
                    parallel::stopCluster(cl)
                }
            } else {
                if (nbcores == 1) {
                    tmp <- apply(.combn.mod(1:p, k, FUN = function(Jset) {prob <- pmvt(lower = rep(-Inf, k), upper = thresh[Jset], delta = rep(0.0, k), corr = R[Jset,Jset], df = df, type = "Kshirsagar", maxpts = maxpts, abseps = abseps / coef, releps = releps) ; if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.") ; return(list(prob, attr(prob, "error")))}, simplify = TRUE), MARGIN = 1, FUN = function(x) sum(as.numeric(x)))
                } else {
                    cl <- parallel::makeCluster(getOption("cl.cores", nbcores))
                    clusterExport(cl = cl, varlist = c("k", "thresh", "R", "df", "maxpts", "abseps", "coef", "releps"), envir = environment())
                    if (LB) {
                        tmp <- apply(parallel::parSapplyLB(cl, unlist(apply(.combn.mod(1:p, k)[,sample(choose(p, k)), drop = FALSE], MARGIN = 2, FUN= list), recursive = FALSE), FUN = function(Jset) {prob <- pmvt(lower = rep(-Inf, k), upper = thresh[Jset], delta = rep(0.0, k), corr = R[Jset,Jset], df = df, type = "Kshirsagar", maxpts = maxpts, abseps = abseps / coef, releps = releps) ; if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.") ; return(list(prob, attr(prob, "error")))}), MARGIN = 1, FUN = function(x) sum(as.numeric(x)))
                    } else {
                        tmp <- apply(parallel::parApply(cl, X = .combn.mod(1:p, k)[,sample(choose(p, k)), drop = FALSE], MARGIN = 2, FUN = function(Jset) {prob <- pmvt(lower = rep(-Inf, k), upper = thresh[Jset], delta = rep(0.0, k), corr = R[Jset,Jset], df = df, type = "Kshirsagar", maxpts = maxpts, abseps = abseps / coef, releps = releps) ; if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.") ; return(c(prob, attr(prob, "error")))}), MARGIN = 1, FUN = function(x) sum(as.numeric(x)))
                    }
                    parallel::stopCluster(cl)
                }
            }
            pow <- pow + (-1) ^ (k - p + r) * coef * tmp[1]
            error <- error + coef * tmp[2]
        }
    } else { # Theorem 2.1, p. 3, second formula.
        pow <- 0
        error <- 0
        for (k in r:p) {
            cork <- R[1:k, 1:k]
            coef <- choose(k - 1, k - r)
            # We are under H_1:
            if (asympt) {
                if (nbcores == 1) {
                    tmp <- apply(.combn.mod(1:p, k, FUN = function(Jset) {prob <- pmvnorm(lower = thresh[Jset], upper = rep(Inf, k), mean = rep(0.0, k), corr = R[Jset,Jset], maxpts = maxpts, abseps = abseps / coef, releps = releps) ; if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.") ; return(list(prob, attr(prob, "error")))}, simplify = TRUE), MARGIN = 1, FUN = function(x) sum(as.numeric(x)))
                } else {
                    cl <- parallel::makeCluster(getOption("cl.cores", nbcores))
                    clusterExport(cl = cl, varlist = c("k", "thresh", "R", "maxpts", "abseps", "coef", "releps"), envir = environment())
                    if (LB) {
                        tmp <- apply(parallel::parSapplyLB(cl, unlist(apply(.combn.mod(1:p, k)[,sample(choose(p, k)), drop = FALSE], MARGIN = 2, FUN= list), recursive = FALSE), FUN = function(Jset) {prob <- pmvnorm(lower = thresh[Jset], upper = rep(Inf, k), mean = rep(0.0, k), corr = R[Jset,Jset], maxpts = maxpts, abseps = abseps / coef, releps = releps) ; if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.") ; return(list(prob, attr(prob, "error")))}), MARGIN = 1, FUN = function(x) sum(as.numeric(x)))
                    } else {
                        tmp <- apply(parallel::parApply(cl, .combn.mod(1:p, k)[,sample(choose(p, k)), drop = FALSE], MARGIN = 2, FUN = function(Jset) {prob <- pmvnorm(lower = thresh[Jset], upper = rep(Inf, k), mean = rep(0.0, k), corr = R[Jset,Jset], maxpts = maxpts, abseps = abseps / coef, releps = releps) ; if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.") ; return(c(prob, attr(prob, "error")))}), MARGIN = 1, FUN = function(x) sum(as.numeric(x)))
                    }
                    parallel::stopCluster(cl)
                }
            } else {
                if (nbcores == 1) {
                    tmp <- apply(.combn.mod(1:p, k, FUN = function(Jset) {prob <- pmvt(lower = thresh[Jset], upper = rep(Inf, k), delta = rep(0.0, k), corr = R[Jset,Jset], df = df, type = "Kshirsagar", maxpts = maxpts, abseps = abseps / coef, releps = releps) ; if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.") ; return(list(prob, attr(prob, "error")))}, simplify = TRUE), MARGIN = 1, FUN = function(x) sum(as.numeric(x)))
                } else {
                    cl <- parallel::makeCluster(getOption("cl.cores", nbcores))
                    clusterExport(cl = cl, varlist = c("k", "thresh", "R", "df", "maxpts", "abseps", "coef", "releps"), envir = environment())
                    if (LB) {
                        tmp <- apply(parallel::parSapplyLB(cl, unlist(apply(.combn.mod(1:p, k)[,sample(choose(p, k)), drop = FALSE], MARGIN = 2, FUN= list), recursive = FALSE), FUN = function(Jset) {prob <- pmvt(lower = thresh[Jset], upper = rep(Inf, k), delta = rep(0.0, k), corr = R[Jset,Jset], df = df, type = "Kshirsagar", maxpts = maxpts, abseps = abseps / coef, releps = releps) ; if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.") ; return(list(prob, attr(prob, "error")))}), MARGIN = 1, FUN = function(x) sum(as.numeric(x)))
                    } else {
                        tmp <- apply(parallel::parApply(cl, .combn.mod(1:p, k)[,sample(choose(p, k)), drop = FALSE], MARGIN = 2, FUN = function(Jset) {prob <- pmvt(lower = thresh[Jset], upper = rep(Inf, k), delta = rep(0.0, k), corr = R[Jset,Jset], df = df, type = "Kshirsagar", maxpts = maxpts, abseps = abseps / coef, releps = releps) ; if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.") ; return(c(prob, attr(prob, "error")))}), MARGIN = 1, FUN = function(x) sum(as.numeric(x)))
                    }
                    parallel::stopCluster(cl)
                }
            }
            pow <- pow + (-1) ^ (k - r) * coef * tmp[1]
            error <- error + coef * tmp[2]
        }
    }
    return(list(pow = as.numeric(pow), error = as.numeric(error)))
}


################################################################################
# Computation of the power in Theorem 2.4
################################################################################

# \mathcal{J}(\underline{a}, p), section 2.2 p. 3 in our Stat. Med. paper (elements involved in \Sum' in P. Delorme's M. Sc. thesis  )
.ind.sumprime <- function(a, p) {
  
    # We remove the repetitions in a
    a <- unique(a)
    
    q <- length(a)
    a <- c(0, a)

    # Computation of the number of possibilities
    card <- .coeff.mnom(diff(c(a, p))) # Note: diff(c(a, p))) is \Delta^*_{\underline{a}}
    if (is.na(card) || (card == 0)) stop(".coeff.mnom() returned a non-integer value; either 'p - r + 1' or 'p' is too large. This package should not be used when the number of endpoints is too large.")
  
    # Will contain all the vectors \underline{j} in \mathcal{J}(\underline{a}, p)
    Matvecj <- matrix(NA, ncol = card, nrow = a[q + 1])

## Not completely (re)verified from here!!
    bloc.ind <- as.list(1:q)
    bloc.lengths <- vector("integer", q)

    # a[1] is a_0    a[h + 1] is a_h
    for (h in 1:q) {
        bloc.ind[[h]] <- if ((a[h + 1] - 1) < (a[h + 1 - 1] + 1)) a[h + 1 - 1] + 1 else (a[h + 1 - 1] + 1):(a[h + 1] - 1 + 1)
        bloc.lengths[h] <-  if ((a[h + 1] - 1) < (a[h + 1 - 1] + 1)) 1 else length(bloc.ind[[h]])
    }
    
    mat.temp <- .combn.mod(1:p, bloc.lengths[1]) # Note: bloc.lengths[1] is a_1 - 1
    Matvecj[bloc.ind[[1]], ] <- apply(mat.temp, MARGIN = 2, FUN = function(x) rep(x, card / ncol(mat.temp))) # Note: bloc.ind[[1]] is 1,...,a_1-1

    card <- card / ncol(mat.temp)
      
    if (length(bloc.ind) > 1) {
        for (bloc in 2:q) {
            tmp2 <- NULL
            for (j in 1:ncol(mat.temp)) {            
                tmp <- .combn.mod(setdiff(1:p, mat.temp[, j]), bloc.lengths[bloc])
                Matvecj[bloc.ind[[bloc]], 1:card + (j - 1) * card] <- apply(tmp, MARGIN = 2, FUN = function(x) rep(x, card / ncol(tmp)))
                tmp2 <- cbind(tmp2, tmp)
            }
            card <- card / ncol(tmp) 
            tmp2 <- rbind(matrix(apply(mat.temp, MARGIN = 2, FUN = function(x) rep(x, ncol(tmp))), ncol = ncol(tmp2)), tmp2)
            mat.temp <- tmp2
        }
    }
    return(Matvecj)
}


.Pane <- function(a.vec, v, effectsize, p, df, R, asympt, maxpts, abseps, releps) {
    MatJap <- .ind.sumprime(a.vec, p)
    # Indices of v associated to the T_k's in the formula of Corollary 2.5.
    indices <- rep(1:length(v), c(a.vec[1], diff(a.vec)))
    if (asympt) {
        res <- apply(MatJap, MARGIN = 2, FUN = function(jtuple) {
                      mask <- !duplicated(jtuple)
                      jtuple <- unique(jtuple)
                      lenjtuple <- length(jtuple)
                      if (lenjtuple == 1) {
                                        # Because P[T_j <= v_i | H_1] = P[T_j <= v_i - delta_j / sighat_j | H_0]
                          return(c(prob = pnorm(v[indices[mask]] - effectsize[jtuple]), error = 0))
                      } else {
                                        # Because P[T_j <= v_i | H_1] = P[T_j <= v_i - delta_j / sighat_j | H_0]
                          prob <- pmvnorm(lower = rep(-Inf, lenjtuple), upper = v[indices[mask]] - effectsize[jtuple], mean = rep(0.0, lenjtuple), corr = R[jtuple,jtuple], maxpts = maxpts, abseps = abseps, releps = releps)
                          if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.")
                          error <- attr(prob, "error")
                          return(c(prob = as.numeric(prob), error = error))
                      }
                  })
        error <- sum(res[2, ])
        res <- sum(res[1, ])
    } else {
        res <- apply(MatJap, MARGIN = 2, FUN = function(jtuple) {
                      mask <- !duplicated(jtuple)
                      jtuple <- unique(jtuple)
                      lenjtuple <- length(jtuple)
                      if (lenjtuple == 1) {
                                        # Because P[T_j <= v_i | H_1] = P[T_j <= v_i - delta_j / sighat_j | H_0]
                          return(c(prob = pt(v[indices[mask]] - effectsize[jtuple], df = df), error = 0))
                      } else {
                                        # Because P[T_j <= v_i | H_1] = P[T_j <= v_i - delta_j / sighat_j | H_0]
                          prob <- pmvt(lower = rep(-Inf, lenjtuple), upper = v[indices[mask]] - effectsize[jtuple], delta = rep(0.0, lenjtuple), corr = R[jtuple,jtuple], df = df, type = "Kshirsagar", maxpts = maxpts, abseps = abseps, releps = releps)
                          if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.")
                          error <- attr(prob, "error")
                          return(c(prob = as.numeric(prob), error = error))
                      }
                  })
        error <- sum(res[2, ])
        res <- sum(res[1, ])
    }
    return(c(prob = as.numeric(res), error = error))
}

.sumastarne <- function(a.vec, v, effectsize, p, r, df, R, asympt, maxpts, abseps, releps) {
    lena <- length(a.vec)
    a.vec <- c(0, a.vec) # So that a_0 = 0    
    prod <- 1
    for (h in 1:(p - r  + 1)) prod <- prod * choose(a.vec[h + 1] - a.vec[h] - 1, a.vec[h + 1] - h)
    a.vec <- a.vec[-1]
    if (prod == 0) {
        res <- 0
        error <- 0
    } else { # Portion of equ. (1) p. 3
        res <- .Pane(a.vec, v, effectsize, p, df, R, asympt, maxpts, abseps / prod, releps)
        error <- res["error"] * prod
        res <- ((-1) ^ sum(a.vec)) * res["prob"] * prod
    }
    return(c(prob = as.numeric(res), error = error))
}


.Psirmune <- function(r, m, p = m, nE, nCovernE, delta, SigmaC, SigmaE, R, alpha, q, asympt, maxpts, abseps, releps, nbcores, LB, orig.Hochberg) {

    nC <- as.integer(nCovernE * nE)
    df <- df.compute(nE, nC, SigmaE, SigmaC)

    if (orig.Hochberg && (q == 1)) D1qm <- 1 else D1qm <- .D1(q, m)
    
    if (q == 1) {
        # FWER control (see for example the slides for my talk at McGill, 2012, and also Romano (2006) equ. (13) and first line of the following paragraph).
        alphais <- 1 / (m:1)
        alphaprime <- alpha * alphais / D1qm
        if (asympt) {
            u <- qnorm(1 - rev(alphaprime), lower.tail = TRUE)
        } else {
            u <- qt(1 - rev(alphaprime), df, lower.tail = TRUE)
        }
    } else { 
        # q-gFWER control; see Theorem 3.1 p. 6 in Romano (2006). It is described in our Example 2.6 p. 4.
        if (q < m){
            alphais <- c(rep(q / m, q), q / ((m - 1):q))
            alphaprime <- alpha * alphais / D1qm
            if (asympt) {
                u <- qnorm(1 - rev(alphaprime))
            } else {
                u <- qt(1 - rev(alphaprime), df)
            }
        } else {
            alphais <- rep(q / m, q)
            alphaprime <- alpha * alphais / D1qm
            if (asympt) {
                u <- qnorm(1 - rev(alphaprime))
            } else{
                u <- qt(1 - rev(alphaprime), df)
            }
        }
    }
    v <- u[(m - p + 1):(m - r  + 1)]
    sighat <- sqrt(diag(SigmaE) / nE + diag(SigmaC) / nC)
    effectsize <- delta / sighat
    
    # Equ. (1) p. 3.
    w <- 1:(p - r + 1)
    if (nbcores == 1) {
        res <- apply(.ind.sumastar(w, p)$Mat, FUN = .sumastarne, MARGIN = 2, v = v, effectsize = effectsize, p = p, r = r, df = df, R = R, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps)
    } else {
        cl <- parallel::makeCluster(getOption("cl.cores", nbcores))
        if (LB) {
            res <- parallel::parSapplyLB(cl, unlist(apply(.ind.sumastar(w, p)$Mat[, sample(.nbsumastar(w, p)), drop = FALSE], MARGIN = 2, FUN= list), recursive = FALSE), FUN = .sumastarne, v = v, effectsize = effectsize, p = p, r = r, df = df, R = R, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps)
        } else {
            res <- parallel::parApply(cl, .ind.sumastar(w, p)$Mat[, sample(.nbsumastar(w, p)), drop = FALSE], FUN = .sumastarne, MARGIN = 2, v = v, effectsize = effectsize, p = p, r = r, df = df, R = R, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps)
        }
        parallel::stopCluster(cl)
    }
    error <- sum(res[2, ])
    res <-  1 - (-1) ^ ((p - r + 1) * (p - r + 2) / 2) * sum(res[1, ])

    return(list(pow = as.numeric(res), error = as.numeric(error)))
}


################################################################################
# Computation of the power in Theorem 2.8
################################################################################


.Patildene <- function(a.vec, d, effectsize, p, df, R, asympt, maxpts, abseps, releps) {
    r <- length(a.vec)
    MatJap <- .ind.sumprime(a.vec, p)
    # Indices of v associated to the T_k's in the formula of Corollary 2.9.
    indices <- rep(p:(p - r  + 1), c(a.vec[1], diff(a.vec)))
    if (asympt) {
        res <- apply(MatJap, MARGIN = 2, FUN = function(jtuple) {
                      mask <- !duplicated(jtuple)
                      jtuple <- unique(jtuple)
                      lenjtuple <- length(jtuple)
                      if (lenjtuple == 1) {
                                        # Because P[T_j > d_i | H_1] = P[T_j > d_i - delta_j / sighat_j | H_0]
                          return(c(prob = 1 - pnorm(d[indices[mask]] - effectsize[jtuple]), error = 0))
                      } else {
                                        # Because P[T_j > d_i | H_1] = P[T_j > d_i - delta_j / sighat_j | H_0]
                          prob <- pmvnorm(lower = d[indices[mask]] - effectsize[jtuple], upper = rep(Inf, lenjtuple), mean = rep(0.0, lenjtuple), corr = R[jtuple,jtuple], maxpts = maxpts, abseps = abseps, releps = releps)
                          if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.")
                          error <- attr(prob, "error")
                          return(c(prob = as.numeric(prob), error = error))
                      }
                  })
        error <- sum(res[2, ])
        res <- sum(res[1, ])

    } else {
        res <- apply(MatJap, MARGIN = 2, FUN = function(jtuple) {
                      mask <- !duplicated(jtuple)
                      jtuple <- unique(jtuple)
                      lenjtuple <- length(jtuple)
                      if (lenjtuple == 1) {
                                        # Because P[T_j > d_i | H_1] = P[T_j > d_i - delta_j / sighat_j | H_0]
                          return(c(prob = 1 - pt(d[indices[mask]] - effectsize[jtuple], df = df), error = 0))
                      } else {
                                        # Because P[T_j > d_i | H_1] = P[T_j > d_i - delta_j / sighat_j | H_0]
                          prob <- pmvt(lower = d[indices[mask]] - effectsize[jtuple], upper = rep(Inf, lenjtuple), delta = rep(0.0, lenjtuple), corr = R[jtuple,jtuple], df = df, type = "Kshirsagar", maxpts = maxpts, abseps = abseps, releps = releps)
                          if (attr(prob, "msg") != "Normal Completion") warning("You should consider increasing the value of 'maxpts'.")
                          error <- attr(prob, "error")
                          return(c(prob = as.numeric(prob), error = error))
                      }
                  })
        error <- sum(res[2, ])
        res <- sum(res[1, ])
    }
    return(c(prob = as.numeric(res), error = error))
}


.sumastartildene <- function(a.vec, d, effectsize, p, r, df, R, asympt, maxpts, abseps, releps) {
    lena <- length(a.vec)
    a.vec <- c(0, a.vec) # So that a_0 = 0    
    prod <- 1
    for (h in 1:r) prod <- prod * choose(a.vec[h + 1] - a.vec[h] - 1, a.vec[h + 1] - h)
    a.vec <- a.vec[-1]
    if (prod == 0) {
        res <- 0
        error <- 0
    } else { # Portion of equ. (2) p. 4
        res <- .Patildene(a.vec, d, effectsize, p, df, R, asympt, maxpts, abseps / prod, releps)
        error <- res["error"] * prod
        res <- ((-1) ^ sum(a.vec)) * res["prob"] * prod
    }
    return(c(prob = as.numeric(res), error = error))
}




.Psirmdne <- function(r, m, p = m, nE, nCovernE, delta, SigmaC, SigmaE, R, alpha, q, asympt, maxpts, abseps, releps, nbcores, LB) {

    nC <- as.integer(nCovernE * nE)
    df <- df.compute(nE, nC, SigmaE, SigmaC)
    
    if (q == 1) {
        # FWER control (see for example the slides for my talk at McGill, 2012, and also Romano (2006) equ. (13) and first line of the following paragraph).
        alphais <- 1 / (m:1)
        alphadbleprime <- alpha * alphais
        if (asympt) {
            h <- qnorm(1 - rev(alphadbleprime), lower.tail = TRUE)
        } else {
            h <- qt(1 - rev(alphadbleprime), df, lower.tail = TRUE)
        }
    } else { 
        # q-gFWER control; see Theorem 3.1 p. 6 in Romano (2006). It is described in our Example 2.10 p. 4.
        if (q < m) {
            alphais <- c(rep(q / m, q), q / ((m - 1):q))
            alphadbleprime <- alpha * alphais
            if (asympt) {
                h <- qnorm(1 - rev(alphadbleprime))
            } else {
                h <- qt(1 - rev(alphadbleprime), df)
            }
        } else {
            alphais <- rep(q / m, q)
            alphadbleprime <- alpha * alphais
            if (asympt) {
                h <- qnorm(1 - rev(alphadbleprime))
            } else{
                h <- qt(1 - rev(alphadbleprime), df)
            }
        }
    }
    d <- h[(m - p + 1):m]
    sighat <- sqrt(diag(SigmaE) / nE + diag(SigmaC) / nC)
    effectsize <- delta / sighat

    # Equ. (2) p. 4.
    tvec <- 1:r
    if (nbcores == 1) {
        res <- apply(.ind.sumastar(tvec, p)$Mat, FUN = .sumastartildene, MARGIN = 2, d = d, effectsize = effectsize, p = p, r = r, df = df, R = R, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps)
    } else {
        cl <- parallel::makeCluster(getOption("cl.cores", nbcores))
        if (LB) {
            res <- parallel::parSapplyLB(cl, unlist(apply(.ind.sumastar(tvec, p)$Mat[, sample(.nbsumastar(tvec, p)), drop = FALSE], MARGIN = 2, FUN= list), recursive = FALSE), FUN = .sumastartildene, d = d, effectsize = effectsize, p = p, r = r, df = df, R = R, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps)
        } else {
            res <- parallel::parApply(cl, .ind.sumastar(tvec, p)$Mat[, sample(.nbsumastar(tvec, p)), drop = FALSE], FUN = .sumastartildene, MARGIN = 2, d = d, effectsize = effectsize, p = p, r = r, df = df, R = R, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps)
        }
        parallel::stopCluster(cl)
    }
    error <- sum(res[2, ])
    res <-  (-1) ^ (r * (r + 1) / 2) * sum(res[1, ])

    return(list(pow = as.numeric(res), error = as.numeric(error)))
}
