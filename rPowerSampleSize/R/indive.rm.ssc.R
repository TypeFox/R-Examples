################################################################################
#       Function for sample size computation in Exchangeable Case              #
#       Decision Rules: At least r endpoints significant among m               #
#       Authors: Delorme P., Lafaye de Micheaux P., Liquet B. and Riou J.      #
#       Date: Le 17/09/2015   Version 1                                        #
################################################################################

## The function .indive.rm.ssc() is only for a weak control of the type-II q-gFWER, namely called with p = m.
## We will have to generalize it to strong control.
## We will have to deal with the unilateral and bilateral tests!!!




.indive.rm.ssc <- function(method, asympt, r, m, p, nCovernE, delta, SigmaC, SigmaE, R, power, alpha, interval, q, maxpts, abseps, releps, nbcores, LB, orig.Hochberg) {
    if (method == "Bonferroni") {
        result <- uniroot.integer(function(nE) .Psirmse(r = r, m = m, p = m, nE = as.integer(nE), nCovernE = nCovernE, delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R = R, alpha = alpha, q = q, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps)$pow - power, interval = interval, pos.side = TRUE)$root 
    }
    if (method == "Hochberg") {
        result <- uniroot.integer(function(nE) .Psirmue(r = r, m = m, p = m, nE = as.integer(nE), nCovernE = nCovernE, delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R = R, alpha = alpha, q = q, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps, nbcores = nbcores, LB = LB, orig.Hochberg = orig.Hochberg)$pow - power, interval, pos.side = TRUE)$root 
  }
    if (method == "Holm") {
        result <- uniroot.integer(function(nE) .Psirmde(r = r, m = m, p = m, nE = as.integer(nE), nCovernE = nCovernE, delta = delta, SigmaC = SigmaC, SigmaE = SigmaE, R = R, alpha = alpha, q =  q, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps, nbcores = nbcores, LB = LB)$pow - power, interval = interval, pos.side = TRUE)$root 
    }
    return(result)
}

################################################################################
# Computation of the power in Corollary 2.2
################################################################################

.Psirmse <- function(r, m, p, nE, nCovernE, delta, SigmaC, SigmaE, R, alpha, q, asympt, maxpts, abseps, releps) {

    nC <- as.integer(nCovernE * nE)
    df <- df.compute(nE, nC, SigmaE, SigmaC)
    
    # Computation of the critical value cm for the Bonferroni method (Example 2.3 p. 3).
    # q-gFWER control
    if (asympt) {
        cm <- qnorm(1 - q * alpha / m, lower.tail = TRUE)
    } else {
        cm <- qt(1 - q * alpha / m, df, lower.tail = TRUE)
    }
    sighat <- sqrt(diag(SigmaE) / nE + diag(SigmaC) / nC)
    effectsize <- delta / sighat
    # Because P[T_j > c_m | H_1] = P[T_j > c_m - effectsize | H_0]
    thresh <- cm - effectsize
    
    if ((p - r + 1) > r) { # Corollary 2.2, p. 3, first formula.
        pow <- 1
        error <- 0
        for (k in (p - r + 1):p) {
            cork <- R[1:k, 1:k]
            coef <- choose(p, k) * choose(k - 1, p - r)
            # We are under H_1:
            if (asympt) {
                prob <- pmvnorm(lower = rep(-Inf, k), upper = thresh[1:k], mean = rep(0.0, k), corr = cork, maxpts = maxpts, abseps = abseps / coef, releps = releps)
                if (attr(prob, "msg") != "Normal Completion") warning(paste("You should consider increasing the value of 'maxpts'. Note: coef = ", coef))
            } else {
                prob <- pmvt(lower = rep(-Inf, k), upper = thresh[1:k], delta = rep(0.0, k), corr = cork, df = df, type = "Kshirsagar", maxpts = maxpts, abseps = abseps / coef, releps = releps)
                if (attr(prob, "msg") != "Normal Completion") warning(paste("You should consider increasing the value of 'maxpts'. Note: coef = ", coef))
            }
            error <- error + coef * attr(prob, "error")
            pow <- pow + (-1) ^ (k - p + r) * coef * prob
        }
    } else { # Corollary 2.2, p. 3, second formula.
        pow <- 0
        error <- 0
        for (k in r:p) {
            cork <- R[1:k, 1:k]
            coef <- choose(k - 1, k - r) * choose(p, k)
            # We are under H_1:
            if (asympt) {
                prob <- pmvnorm(lower = thresh[1:k], upper = rep(Inf, k), mean = rep(0.0, k), corr = cork, maxpts = maxpts, abseps = abseps / coef, releps = releps)
                if (attr(prob, "msg") != "Normal Completion") warning(paste("You should consider increasing the value of 'maxpts'. Note: coef = ", coef))
            } else {
                prob <- pmvt(lower = thresh[1:k], upper = rep(Inf, k), delta = rep(0.0, k), corr = cork, df = df, type = "Kshirsagar", maxpts = maxpts, abseps = abseps / coef, releps = releps)
                if (attr(prob, "msg") != "Normal Completion") warning(paste("You should consider increasing the value of 'maxpts'. Note: coef = ", coef))
            }
            error <- error + coef * attr(prob, "error")
            pow <- pow + (-1) ^ (k - r) * coef * prob
        }  
    }
    return(list(pow = as.numeric(pow), error = error))
}

################################################################################
# Computation of the power in Theorem 2.4 using Corollary 2.5
################################################################################

.nbsumastar <- function(l.vec, p, q = length(l.vec), k = 0) {
    # Uses recursivity. 
    if (all.equal(l.vec, 1:q)) return(choose(p + q, q) - choose(p + q, q - 1)) # Catalan's triangle numbers.
    if (q == 0) {
        k <- k + 1
        return(k = k)
    }
    for (a in l.vec[q]:p) {
        k <- .nbsumastar(l.vec, a, q - 1, k)
    }
    return(k)
}

.ind.sumastar <- function(l.vec, p, q = length(l.vec), Mat = matrix(NA, nrow = q, ncol = .nbsumastar(l.vec, p)), k = 1, indices = NULL, oldas = rep(NA, q)) {
    # Uses recursivity. 
    if (q == 0) {
        Mat[, k] <- oldas
        k <- k + 1
        return(list(Mat = Mat, k = k, indices = indices, oldas = oldas))
    }
  
    for (a in l.vec[q]:p) {
        indices <- c(indices, q)
        oldas[q] <- a
        tmp <- .ind.sumastar(l.vec, a, q - 1, Mat, k, indices, oldas)
        Mat <- tmp$Mat
        k <- tmp$k
        indices <- tmp$indices
        oldas <- tmp$oldas
    }
    return(list(Mat = Mat, k = k, indices = indices, oldas = oldas))
}

.coeff.mnom <- function(x) { # Multinomial coefficient.
    if (length(x) == 1) {
        res <- choose(sum(x), x)
    } else {
        res <- exp(lgamma(sum(x) + 1) - sum(lgamma(x + 1)))
    }
    # We use round() because .coeff.mnom(c(1,2)) does not give the same answer as as.integer(.coeff.mnom(c(1,2))) !
    return(as.integer(round(res)))
}

.Pae <- function(a.vec, v, effectsize, p, df, R, asympt, maxpts, abseps, releps) {
    # Indices of v associated to the T_k's in the formula of Corollary 2.5.
    indices <- rep(1:length(a.vec), c(a.vec[1], diff(a.vec)))
    lasta <- a.vec[length(a.vec)] # This is a_{p-r+1}. This is also length(indices)
    # Correlation matrix.
    cor <- R[1:lasta, 1:lasta]
    coef <- .coeff.mnom(diff(c(0, a.vec, p)))
    if (asympt) {
        if (lasta == 1) { # i.e., p = r.
                                        # Because P[T_j <= v_i | H_1] = P[T_j <= v_i - delta_j / sighat_j | H_0]
            prob <- pnorm(v - effectsize[1])
            error <- 0
        } else {
                                        # Because P[T_j <= v_i | H_1] = P[T_j <= v_i - delta_j / sighat_j | H_0]
            prob <- pmvnorm(lower = rep(-Inf, lasta), upper = v[indices] - effectsize[1:lasta], mean = rep(0.0, lasta), corr = cor, maxpts = maxpts, abseps = abseps / coef, releps = releps)
            error <- attr(prob, "error")
            if (attr(prob, "msg") != "Normal Completion") warning(paste("You should consider increasing the value of 'maxpts'. Note: coef = ", coef))
        }
    } else {
        if (lasta == 1) {
                                        # Because P[T_j <= v_i | H_1] = P[T_j <= v_i - delta_j / sighat_j | H_0]
            prob <- pt(v - effectsize[1], df = df)
            error <- 0
        } else {
                                        # Because P[T_j <= v_i | H_1] = P[T_j <= v_i - delta_j / sighat_j | H_0]
            prob <- pmvt(lower = rep(-Inf, lasta), upper = v[indices] - effectsize[1:lasta], delta = rep(0.0, lasta), corr = cor, df = df, type = "Kshirsagar", maxpts = maxpts, abseps = abseps / coef, releps = releps)
            error <- attr(prob, "error")
            if (attr(prob, "msg") != "Normal Completion") warning(paste("You should consider increasing the value of 'maxpts'. Note: coef = ", coef))
        }
    }
    # Equ. in Corollary 2.5 p. 4.
    res <- coef * c(prob, error)
    return(c(prob = res[1], error = res[2]))
}

.sumastare <- function(a.vec, v, effectsize, p, r, df, R, asympt, maxpts, abseps, releps) {
    lena <- length(a.vec)
    a.vec <- c(0, a.vec) # So that a_0 = 0    
    prod <- 1
    for (h in 1:(p - r  + 1)) prod <- prod * choose(a.vec[h + 1] - a.vec[h] - 1, a.vec[h + 1] - h)
    a.vec <- a.vec[-1]
    if (prod == 0) {
        res <- c(prob = 0, error = 0)
    } else { # Portion of equ. (1) p. 3
        res <- ((-1) ^ sum(a.vec)) * .Pae(a.vec, v, effectsize, p, df, R, asympt, maxpts, abseps, releps) * prod
        res <- c(prob = res["prob"], error = abs(res["error"]))
    }
    return(res)
}


.Psirmue <- function(r, m, p, nE, nCovernE, delta, SigmaC, SigmaE, R, alpha, q, asympt, maxpts, abseps, releps, nbcores, LB, orig.Hochberg) {
    
    nC <- as.integer(nCovernE * nE)
    df <- df.compute(nE, nC, SigmaE, SigmaC)

    if (orig.Hochberg && (q == 1)) D1qm <- 1 else D1qm <- .D1(q, m)

    if (q == 1) {
        # FWER control (see for example the slides for my talk at McGill, 2012, and also Romano (2006) equ. (13) and first line of the following paragraph).
        if (asympt) {
            u <- qnorm(1 - (alpha / (D1qm * (1:m))), lower.tail = TRUE)
        } else {
            u <- qt(1 - (alpha / ( D1qm * (1:m))), df, lower.tail = TRUE)
        }
    } else { 
        # q-gFWER control; see Theorem 3.1 p. 6 in Romano (2006). It is described in our Example 2.6 p. 4.
        if (q < m){
            alphaprime <- alpha * c(rep(q / m, q), q / ((m - 1):q)) / D1qm
            if (asympt) {
                u <- qnorm(1 - rev(alphaprime))
            } else {
                u <- qt(1 - rev(alphaprime), df)
            }
        } else {
            alphaprime <- alpha * rep(q / m, q) / D1qm
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
        res <- apply(.ind.sumastar(w, p)$Mat, FUN = .sumastare, MARGIN = 2, v = v, effectsize = effectsize, p = p, r = r, df, R = R, asympt = asympt, maxpts, abseps, releps)
    } else {
        cl <- parallel::makeCluster(getOption("cl.cores", nbcores))
        if (LB) {
            res <- parallel::parSapplyLB(cl, unlist(apply(.ind.sumastar(w, p)$Mat[, sample(.nbsumastar(w, p)), drop = FALSE], MARGIN = 2, FUN= list), recursive = FALSE), FUN = .sumastare, v = v, effectsize = effectsize, p = p, r = r, df = df, R = R, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps)
        } else {
            res <- parallel::parApply(cl, .ind.sumastar(w, p)$Mat[, sample(.nbsumastar(w, p)), drop = FALSE], FUN = .sumastare, MARGIN = 2, v = v, effectsize = effectsize, p = p, r = r, df = df, R = R, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps)
        }
        parallel::stopCluster(cl)
    }
    error <- sum(res[2, ])
    res <-  1 - (-1) ^ ((p - r + 1) * (p - r + 2) / 2) * sum(res[1, ])

    return(list(pow = as.numeric(res), error = error))
}

################################################################################
# Computation of the power in Theorem 2.8 using Corollary 2.9
################################################################################

.Patildee <- function(a.vec, d, effectsize, p, r, df, R, asympt, maxpts, abseps, releps) {
    # Indices of v associated to the T_k's in the formula of Corollary 2.9.
    indices <- rep(p:(p - r  + 1), c(a.vec[1], diff(a.vec)))
    lasta <- a.vec[length(a.vec)] # This is a_r. This is also length(indices).
    # Correlation matrix.
    cor <- R[1:lasta, 1:lasta]
    coef <- .coeff.mnom(diff(c(0, a.vec, p)))
    if (asympt) {
        if (lasta == 1) { # i.e., r = 1.
                                        # Because P[T_j > d_i | H_1] = P[T_j > d_i - delta_j / sighat_j | H_0]
            prob <- 1 - pnorm(d[p] - effectsize[1])
            error <- 0
        } else {
                                        # Because P[T_j > d_i | H_1] = P[T_j > d_i - delta_j / sighat_j | H_0]
            prob <- pmvnorm(lower = d[indices] - effectsize[1:lasta], upper = rep(Inf, lasta), mean = rep(0, lasta), corr = cor, maxpts = maxpts, abseps = abseps / coef, releps = releps)
            error <- attr(prob, "error")
            if (attr(prob, "msg") != "Normal Completion") warning(paste("You should consider increasing the value of 'maxpts'. Note: coef = ", coef))
        }
    } else {
        if (lasta == 1) { # i.e., r = 1.
                                        # Because P[T_j > d_i | H_1] = P[T_j > d_i - delta_j / sighat_j | H_0]
            prob <- 1 - pt(d[p] - effectsize[1], df = df)
            error <- 0
        } else {
                                        # Because P[T_j > d_i | H_1] = P[T_j > d_i - delta_j / sighat_j | H_0]
            prob <- pmvt(lower = d[indices] - effectsize[1:lasta], upper = rep(Inf, lasta), delta = rep(0, lasta), corr = cor, df = df, type = "Kshirsagar", maxpts = maxpts, abseps = abseps / coef, releps = releps)
            error <- attr(prob, "error")
            if (attr(prob, "msg") != "Normal Completion") warning(paste("You should consider increasing the value of 'maxpts'. Note: coef = ", coef))
        }
    }
    # Equ. in Corollary 2.5 p. 4.
    res <- coef * c(prob, error)
    return(c(prob = res[1], error = res[2]))
}

.sumastartildee <- function(a.vec, d, effectsize, p, r, df, R, n, asympt, maxpts, abseps, releps) {
    lena <- length(a.vec)
    a.vec <- c(0, a.vec) # So that a_0 = 0
    prod <- 1
    for (h in 1:r) prod <- prod * choose(a.vec[h + 1] - a.vec[h] - 1, a.vec[h + 1] - h)
    a.vec <- a.vec[-1]
    if (prod == 0) {
        res <- c(prob = 0, error = 0)
    } else { # Portion of equ. (2) p. 4
        res <- ((-1) ^ sum(a.vec)) * .Patildee(a.vec, d, effectsize, p, r, df, R, asympt, maxpts, abseps, releps) * prod
        res <- c(prob = res["prob"], error = abs(res["error"]))
    }
    return(res)
}


.Psirmde <- function(r, m, p, nE, nCovernE, delta, SigmaC, SigmaE, R, alpha, q, asympt, maxpts, abseps, releps, nbcores, LB) {
    
    nC <- as.integer(nCovernE * nE)
    df <- df.compute(nE, nC, SigmaE, SigmaC)

    if (q == 1) {
        # FWER control (see for example the slides for my talk at McGill, 2012, and also Romano (2006) equ. (13) and first line of the following paragraph).
        if (asympt) {
            h <- qnorm(1 - (alpha / (1:m)), lower.tail = TRUE)
        } else {
            h <- qt(1 - (alpha / (1:m)), df, lower.tail = TRUE)
        }
    } else { 
        # q-gFWER control; see Theorem 3.1 p. 6 in Romano (2006). It is described in our Example 2.10 p. 4.
        if (q < m){
            alphadbleprime <- alpha * c(rep(q / m, q), q / ((m - 1):q))
            if (asympt) {
                h <- qnorm(1 - rev(alphadbleprime))
            } else {
                h <- qt(1 - rev(alphadbleprime), df)
            }
        } else {
            alphadbleprime <- alpha * rep(q / m, q)
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
        res <- apply(.ind.sumastar(tvec, p)$Mat, FUN = .sumastartildee, MARGIN = 2, d = d, effectsize = effectsize, p = p, r = r, df, R = R, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps)
    } else {
        cl <- parallel::makeCluster(getOption("cl.cores", nbcores))
        if (LB) {
            res <- parallel::parSapplyLB(cl, unlist(apply(.ind.sumastar(tvec, p)$Mat[, sample(.nbsumastar(tvec, p)), drop = FALSE], MARGIN = 2, FUN= list), recursive = FALSE), FUN = .sumastartildee, d = d, effectsize = effectsize, p = p, r = r, df = df, R = R, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps)
        } else {        
            res <- parallel::parApply(cl, .ind.sumastar(tvec, p)$Mat[, sample(.nbsumastar(tvec, p)), drop = FALSE], FUN = .sumastartildee, MARGIN = 2, d = d, effectsize = effectsize, p = p, r = r, df = df, R = R, asympt = asympt, maxpts = maxpts, abseps = abseps, releps = releps)
        }
        parallel::stopCluster(cl)
    }
    error <- sum(res[2, ])
    res <-  (-1) ^ (r * (r + 1) / 2) * sum(res[1, ])

    return(list(pow = as.numeric(res), error = error))
}
