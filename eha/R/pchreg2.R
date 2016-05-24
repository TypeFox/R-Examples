pchreg2 <- function(X, Y, cuts, offset, strata, init, control, center){
    ## Test of stratification of the pch-model ("piecewise constant")
    ## strata \in {1, ..., ns}

    nn <- NROW(Y)
    if (!is.matrix(X)) X <- matrix(X, ncol = 1)
    if (NROW(X) != nn) stop("[pchreg2]: error in X")
    if (missing(strata) || is.null(strata)){
        strata <- rep(1, nn)
        ns <- 1
    }else{
        strata <- factor(strata)
        name.s <- levels(strata)
        strata <- as.integer(strata)
        ns <- max(strata)
    }

    ns <- length(unique(strata))
    ncov <- NCOL(X) # Assume here that ncov >= 0!
    
    n.ivl <- length(cuts) + 1
    ## Some sane checks .....
    ##cuts <- c(0, cuts, Inf)
    if (control$trace){
        cat("[pchreg2]    ns = ", ns, "\n")
        cat("[pchreg2] n.ivl = ", n.ivl, "\n")
    }
    if (length(cuts)){
        split <- SurvSplit(Y, cuts)
        
        strat <- strata[split$idx]
        offset <- offset[split$idx]
        X <- X[split$idx, ,drop = FALSE]
        T <- split$Y$exit - split$Y$enter
        d <- split$Y$event
        ivl <- split$ivl
    }else{
        strat <- strata
        T <- Y[, 2] - Y[, 1]
        d <- Y[, 3]
        ivl <- rep(1, NROW(Y))
    }
    
    Dtot <- sum(d)

    loglik0 <- function(){
        alpha <- matrix(0, nrow = ns, ncol = n.ivl)
        res <- -Dtot
        for (i in seq_len(n.ivl)){
            for (j in 1:ns){
                indx <- (strat == j) & (ivl == i)
                D <- sum(d[indx])
                if (D > 0){
                    sumT <- sum(T[indx])
                    alpha[j, i] <- D / sumT 
                    res <- res + D * (log(D) - log(sumT))
                }
            }
        }
        ##cat("res = ", res, "\n")
        list(loglik = res, hazards = alpha)
    }
    
    loglik <- function(beta){
        ##cat("beta = ", beta, "\n") 
        zb <- offset + X %*% beta
        ##cat("zb = ", zb, "\n")
        tezb <- T * exp(zb)
        ##cat("tezb = ", tezb, "\n")
        ##res <- sum(d * zb) - Dtot
        res <- drop(d %*% zb) - Dtot
        ##cat("res.first = ", res, "\n")
        ##cat("res1 = ", res, "\n")
        for (i in seq_len(n.ivl)){
            for (j in 1:ns){
                indx <- (strat == j) & (ivl == i)
                D <- sum(d[indx])
                ##cat("i = ", i, " j = ", j, " D = ", D, "\n")
                if (D > 0){
                    res <- res + D * (log(D) - log(sum(tezb[indx])))
                }
            }
        }
        ##cat("res = ", res, "\n")
        res
    }

    dloglik <- function(beta){
        zb <- offset + X %*% beta
        tezb <- T * exp(zb)

        res <- drop(d %*% X) # 10 feb 2014
        ##res <- numeric(ncov)
        ##for (j in seq_len(ncov)){
          ##  res[j] <- sum(d * X[, j]) 
        ##}
        ##cat("res1[1] = ", res[1], "\n")
        for (i in seq_len(n.ivl)){
            for (j in 1:ns){
                indx <- (strat == j) & (ivl == i)
                D <- sum(d[indx])
                ##cat("i = ", i, " j = ", j, " D = ", D, "\n")
                if (D > 0){
                    Stezb <- sum(tezb[indx]) 
                    for (j in seq_len(ncov)){ # Do smarter!
                        res[j] <- res[j] -
                            D * sum(X[indx, j] * tezb[indx]) /
                             Stezb
                    }
                }
            } # end for (j in ...
        } # end for (i in ...
        ##cat("res[1] = ", res[1], "\n")
        res
    } # end dloglik

    getAlpha <- function(beta){
        ## Note: f = (alpha * exp(zb)^d * exp(-t * alpha * exp(zb))
        tezb <- T * exp(X %*% beta)
        alpha <- matrix(0, nrow = ns, ncol = n.ivl)
        for (i in seq_len(n.ivl)){
            for (j in 1:ns){
                indx <- (strat == j) & (ivl == i)
                D <- sum(d[indx])
                alpha[j, i] <- D / sum(tezb[indx])
            }
        }
    alpha
    } # end getAlpha

    fit <- loglik0()
    
    if (ncov){
        means <- colMeans(X)
        for (i in seq_len(ncov)){
            X[, i] <- X[, i] - means[i]
        }
        beta <- init
        res <- optim(beta, loglik, dloglik,
                     method = "BFGS", hessian = TRUE, 
                     control = list(fnscale = -1, reltol = 1e-10))
        beta <- res$par
        fit$gradient <- dloglik(beta)
        ##cat("score = ", deriv, " at solution\n")
        fit$loglik <- c(fit$loglik, res$value)
        fit$coefficients <- beta
        fit$var <- solve(-res$hessian)
        fit$hazards <- getAlpha(beta)
        fit$hazards <- fit$hazards * exp(-drop(means %*% fit$coefficients))
        colnames(fit$var) <- rownames(fit$var) <- colnames(X)
    }else{# No covariates
        fit$loglik <- rep(fit$loglik, 2)
        fit$coefficients <- NULL
    }
    
    fit$df <- ncov
    fit$cuts <- cuts
    fit$fail <- FALSE # Optimism...
    cn <- character(n.ivl)
    if (n.ivl > 1){
        cn[1] <- paste("(.., ", cuts[1], "]", sep = "")
        cn[n.ivl] <- paste("(", cuts[n.ivl - 1], ", ...]", sep = "")
        if (n.ivl > 2){
            for (i in 2:(n.ivl - 1)){
                cn[i] <- paste("(", cuts[i-1], ", ", cuts[i], "]", sep = "")
            }
        }
        colnames(fit$hazards) <- cn
    }
    
    if (ns > 1) rownames(fit$hazards) <- name.s
    fit$n.strata <- ns

    fit
}
