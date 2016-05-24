# Fitting multivariate generalized linear model.  
# Author: Hua Zhou and Yiwen Zhang 
##============================================================## 

##============================================================## 
## set class
##============================================================##
setClass("MGLMreg", representation(call = "function", data = "list", coefficients = "list", 
    SE = "list", test = "matrix", Hessian = "matrix", logL = "numeric", BIC = "numeric", 
    AIC = "numeric", iter = "numeric", LRT = "numeric", distribution = "character", 
    fitted = "matrix"))

##============================================================## 
## Regression function 
##============================================================##

MGLMreg <- function(formula, data, dist, init, weight, epsilon = 1e-08, maxiters = 150, 
    display = FALSE, LRT = FALSE, parallel = FALSE, cores, regBeta = FALSE) {
    
    ##----------------------------------------## 
    ## Creating the environment
    ##----------------------------------------## 
    ow <- getOption("warn")
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weight"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame(n = 1))
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    weight <- mf$`(weight)`
    X <- model.matrix(mt, mf, contrasts)
    if (is.null(colnames(Y))) 
        colnames(Y) <- paste("Col", 1:ncol(Y), sep = "_")
    
    ##----------------------------------------## 
    ## Give warnings about zero rows
    ##----------------------------------------## 
    if (dist != "NegMN") {
        if (any(rowSums(Y) == 0)) {
            rmv <- rowSums(Y) == 0
            warning(paste(sum(rmv), " rows are removed because the row sums are 0", 
                "\n", sep = ""))
            Y <- Y[!rmv, ]
            X <- X[!rmv, ]
        }
        if (any(colSums(Y) == 0)) {
            rmv <- colSums(Y) == 0
            warning(paste(sum(rmv), " columns are removed because the column sums are 0", 
                "\n", sep = ""))
            Y <- Y[, !rmv]
        }
    }
    
    ## ----------------------------------------## 
    ## Check dimensions
    N <- nrow(Y)
    d <- ncol(Y)
    p <- ncol(X)
    
    ## ----------------------------------------## 
    ## Check weight length and values
    if (nrow(X) != N) 
        stop("Unequal numbers of observations in X and Y")
    if (is.null(weight)) 
        weight <- rep(1, N)
    if (!is.null(weight) && length(weight) != N) 
        stop("Length of weights doesn't match with the sample size")
    if (!is.null(weight) && any(weight < 0)) 
        stop("Negative weights are not allowed")
    if (missing(weight)) 
        weight <- rep(1, N)
    
    ## ----------------------------------------## 
    ## Check distribution
    if (!is.element(dist, c("MN", "DM", "NegMN", "GDM"))) 
        stop(paste("Dist '", dist, "' not valid\n", sep = ""))
    
    ## ----------------------------------------## 
    ## Input missing weight and cores
    if (missing(cores)) {
        if (parallel) 
            cores <- floor(detectCores()/2) else cores <- NULL
    }
    
    ## ----------------------------------------## 
    ## Fit the model
    if (dist != "GDM") {
        if (parallel) {
            cl <- makeCluster(getOption("cl.cores", cores))
            clusterExport(cl, "glm.private")
            sys <- Sys.info()
            sys <- sys[1]
        } else {
            cl <- NULL
            sys <- NULL
        }
        if (dist == "MN") {
            ## ----------------------------------------## 
            ## MN
            if (N < p * (d - 1)) 
                warning(paste("Sample size is smaller than the number of parameters.", 
                  "\n", sep = ""))
            
            if (missing(init)) {
                init <- matrix(0, p, (d - 1))
                for (i in 1:(d - 1)) {
                  init[, i] <- glm.fit(X, Y[, i], family = poisson(link = "log"), 
                    weights = weight)$coefficients
                }
            } else if (any(dim(init) != c(p, (d - 1)))) 
                stop("Dimension of the initial values is not compatible with the data")
            
            est <- eval(call("DMD.MN.reg", Y = Y, X = X, weight = weight, init = init, 
                epsilon = epsilon, maxiters = maxiters, display = display, parallel = parallel, 
                cores = cores, cl = cl, sys = sys))
        } else if (dist == "DM") {
            ## ----------------------------------------## 
            ## DM
            if (N < p * d) 
                warning(paste("Sample size is smaller than the number of parameters", 
                  "\n", sep = ""))
            
            if (missing(init)) {
                # alphahat <- as.vector( DMD.DM.fit(data=Y, weight=weight,
                # epsilon=epsilon)$estimate) alphahat <- alphahat/sum(alphahat) init <-
                # rbind(alphahat,matrix(0,(p-1),d))
                init <- matrix(0.1, p, d)
                options(warn = -1)
                for (j in 1:d) {
                  fit <- glm.fit(x = X, y = Y[, j]/rowSums(Y), family = binomial(link = "logit"))
                  init[, j] <- fit$coefficients
                }
                options(warn = ow)
            } else if (any(dim(init) != c(p, d))) {
                stop("Dimension of the initial values is not compatible with the data")
            }
            est <- eval(call("DMD.DM.reg", Y = Y, X = X, weight = weight, init = init, 
                epsilon = epsilon, maxiters = maxiters, display = display, parallel = parallel, 
                cores = cores, cl = cl, sys = sys))
        } else if (dist == "NegMN") {
            ## ----------------------------------------## 
            ## NegMN
            if (N < p * (d + 1)) 
                warning(paste("Sample size is smaller than the number of parameters", 
                  "\n", sep = ""))
            
            if (regBeta) {
                if (missing(init)) {
                  alpha_init <- matrix(0, p, d)
                  for (i in 1:d) {
                    alpha_init[, i] <- glm.fit(X, Y[, i], family = poisson(link = "log"), 
                      weights = weight)$coefficients
                  }
                  beta_init <- glm.fit(X, (rowSums(Y) + 10), family = quasipoisson(link = "log"), 
                    weights = weight)$coefficients
                  init <- cbind(alpha_init, beta_init)
                } else {
                  if (any(dim(init) != c(p, (d + 1)))) 
                    stop("Dimension of the initial values is not compatible with the data.") else {
                    alpha_init <- init[, 1:d]
                    beta_init <- init[, (d + 1)]
                  }
                }
                est <- eval(call("DMD.NegMN.reg", Y = Y, init = init, X = X, weight = weight, 
                  epsilon = epsilon, maxiters = maxiters, display = display, parallel = parallel, 
                  cores = cores, cl = cl, sys = sys))
            } else {
                if (missing(init)) {
                  init <- matrix(0, p, d)
                  for (i in 1:d) {
                    init[, i] <- glm.fit(X, Y[, i], family = poisson(link = "log"), 
                      weights = weight)$coefficients
                  }
                } else {
                  if (any(dim(init) != c(p, d))) 
                    stop("Dimension of the initial values is not compatible with the data.")
                }
                est <- eval(call("DMD.NegMN.Alpha.reg", Y = Y, init = init, X = X, 
                  weight = weight, epsilon = epsilon, maxiters = maxiters, display = display, 
                  parallel = parallel, cores = cores, cl = cl, sys = sys))
            }
        }
        if (parallel) 
            stopCluster(cl)
        ## ----------------------------------------## 
        ## GDM
    } else if (dist == "GDM") {
        if (d == 2) 
            stop("When d=2, GDM model is equivilant to DM model, please use dist='DM'.")
        if (parallel) {
            cl <- makeCluster(getOption("cl.cores", cores))
            clusterExport(cl, "DMD.DM.reg")
            clusterExport(cl, "ddirm")
            sys <- Sys.info()
            sys <- sys[1]
        } else {
            cl <- NULL
            sys <- NULL
        }
        
        if (missing(init)) {
            # init <- as.vector( DMD.GDM.fit(data=Y, weight=weight,
            # epsilon=epsilon)$estimate) init <- init/sum(init) init <-
            # rbind(rep(init[1:(d-1)],2), matrix(0, (p-1),2*(d-1)))
            init <- NULL
            
        } else if (any(dim(init) != c(p, 2 * (d - 1)))) {
            stop("Dimension of the initial values is not compatible with the data")
        } else if (N < p * (d - 1) * 2) 
            warning(paste("Sample size is smaller than the number of parameters", 
                "\n", sep = ""))
        
        est <- eval(call("DMD.GDM.reg", Y = Y, X = X, weight = weight, init = init, 
            epsilon = epsilon, maxiters = maxiters, display = display, parallel = parallel, 
            cores = cores, cl = cl, sys = sys))
    }
    
    ## ----------------------------------------## 
    ## Clean up the results
    ## ----------------------------------------## 
    if (dist == "MN") {
        colnames(est$coefficients) <- colnames(Y)[1:(ncol(Y) - 1)]
        colnames(est$SE) <- colnames(Y)[1:(ncol(Y) - 1)]
    } else if (dist == "DM") {
        colnames(est$coefficients) <- colnames(Y)
        colnames(est$SE) <- colnames(Y)
    } else if (dist == "GDM") {
        colnames(est$coefficients) <- c(paste("alpha_", colnames(Y)[1:(d - 1)], sep = ""), 
            paste("beta_", colnames(Y)[1:(d - 1)], sep = ""))
        colnames(est$SE) <- c(paste("alpha", colnames(Y)[1:(d - 1)], sep = ""), paste("beta_", 
            colnames(Y)[1:(d - 1)], sep = ""))
        colnames(est$gradient) <- c(paste("alpha_", colnames(Y)[1:(d - 1)], sep = ""), 
            paste("beta_", colnames(Y)[1:(d - 1)], sep = ""))
    } else if (dist == "NegMN") {
        if (regBeta) {
            colnames(est$coefficients) <- c(colnames(Y), "phi")
            colnames(est$SE) <- c(colnames(Y), "phi")
        } else {
            colnames(est$coefficients$alpha) <- colnames(Y)
            colnames(est$SE$SE.alpha) <- colnames(Y)
        }
    }
    
    if (dist == "NegMN" & !regBeta) {
        rownames(est$coefficients[[1]]) <- colnames(X)
        rownames(est$SE[[1]]) <- colnames(X)
    } else {
        rownames(est$coefficients) <- colnames(X)
        rownames(est$SE) <- colnames(X)
    }
    wald.value <- c(est$wald.value)
    wald.p <- c(est$wald.p)
    est$test <- cbind(wald.value, wald.p)
    rownames(est$test) <- colnames(X)
    colnames(est$test) <- c("wald value", "Pr(>wald)")
    
    
    ## ----------------------------------------## 
    ## LRT test the hypothesis beta_j=0
    ## ----------------------------------------## 
    logL <- est$logL
    if (LRT) {
        options(warn = -1)
        LRT.value <- rep(NA, p)
        LRT.p <- rep(NA, p)
        
        for (t in 1:p) {
            subX <- X[, -t]
            if (dist == "MN") {
                subest <- eval(call("DMD.MN.reg", Y = Y, X = X[, -t], weight = weight, 
                  init = init[-t, ], epsilon = epsilon, maxiters = maxiters, display = FALSE, 
                  parallel = FALSE, cores = cores, cl = cl, sys = sys))
            } else if (dist == "DM") {
                subest <- eval(call("DMD.DM.reg", Y = Y, X = X[, -t], weight = weight, 
                  init = init[-t, ], epsilon = epsilon, maxiters = maxiters, display = FALSE, 
                  parallel = FALSE, cores = cores, cl = cl, sys = sys))
            } else if (dist == "NegMN" & regBeta) {
                subest <- eval(call("DMD.NegMN.reg", Y = Y, X = X[, -t], weight = weight, 
                  init = init[-t, ], epsilon = epsilon, maxiters = maxiters, display = FALSE, 
                  parallel = FALSE, cores = cores, cl = cl, sys = sys))
            } else if (dist == "NegMN" & (!regBeta)) {
                subest <- eval(call("DMD.NegMN.Alpha.reg", Y = Y, X = X[, -t], weight = weight, 
                  init = init[-t, ], epsilon = epsilon, maxiters = maxiters, display = FALSE, 
                  parallel = FALSE, cores = cores, cl = cl, sys = sys))
            } else if (dist == "GDM") {
                subest <- eval(call("DMD.GDM.reg", Y = Y, X = subX, weight = weight, 
                  init = init[-p, ], epsilon = epsilon, maxiters = maxiters, display = FALSE, 
                  parallel = FALSE, cores = cores, cl = cl, sys = sys))
            }
            sublogL <- subest$logL
            LRT.value[t] <- 2 * (logL - sublogL)
        }
        ## ----------------------------------------## 
        ## Calculate the degrees of freedom
        ## ----------------------------------------## 
        Dof <- ifelse(dist == "MN", d - 1, ifelse(dist == "DM", d, ifelse(dist == 
            "GDM", 2 * (d - 1), ifelse(regBeta, d + 1, d))))
        LRT.p <- pchisq(LRT.value, Dof, lower.tail = FALSE)
        test <- cbind(LRT.value, LRT.p)
        colnames(test) <- c("LRT value", "Pr(>LRT)")
        est$test <- cbind(est$test, test)
        options(warn = ow)
    }
    
    ## ----------------------------------------## 
    ## More things to output
    ## ----------------------------------------## 
    est$call <- match.call()
    est$data <- list(Y = Y, X = X)
    est$distribution <- ifelse(dist == "MN", "Multinomial", ifelse(dist == "DM", 
        "Dirichlet Multinomial", ifelse(dist == "GDM", "Generalized Dirichlet Multinomial", 
            "Negative Multinomial")))
    
    ## ----------------------------------------## 
    ## Set class
    ## ----------------------------------------## 
    class(est) <- "MGLMreg"
    return(est)
}

## ============================================================## 
## Function 1 fit Multinomial Reg The dimension of Y is n by d (d>=3)
## ============================================================##

DMD.MN.reg <- function(Y, init, X, weight, epsilon, maxiters, display, parallel, 
    cores, cl, sys) {
    
    ## ----------------------------------------## 
    ## Function to calculate log-likelihood
    ## ----------------------------------------## 
    dmultn_L <- function(B1) {
        B <- cbind(B1, 0)
        alpha <- exp(X %*% B)
        return(sum(Y * (X %*% B) * weight) - sum(Y * log(rowSums(alpha)) * weight) + 
            sum(lgamma(m + 1) * weight) - sum(lgamma(Y + 1) * weight))
    }
    
    ##----------------------------------------## 
    ## Function to calculate gradients
    ##----------------------------------------## 
    dl.MN_fctn <- function(i, beta) {
        x <- (weight * X)[i, ]
        y <- Y[i, ]
        prob <- exp(colSums(x * beta))/(sum(exp(colSums(x * beta))) + 1)
        dl <- y[-ncol(Y)] - sum(y) * prob
        return(matrix(c(outer(x, dl)), (length(y) - 1) * length(x), 1))
    }
    
    ## ----------------------------------------## 
    ## Function to calculate hessian matrix
    ## ----------------------------------------## 
    H.MN_fctn <- function(i, beta) {
        x <- (weight * X)[i, ]
        y <- Y[i, ]
        prob <- exp(colSums(x * beta))/(sum(exp(colSums(x * beta))) + 1)
        if (d > 2) 
            dprob <- diag(prob) else if (d == 2) 
            dprob <- prob
        d2l <- sum(y) * (outer(prob, prob) - dprob)
        o.x <- outer(x, x)
        H <- d2l %x% o.x
        return(H)
    }
    
    ## ----------------------------------------## 
    ## Keep some original values
    ## ----------------------------------------## 
    fitted <- matrix(NA, nrow(Y), ncol(Y))
    emptyRow <- rowSums(Y) == 0
    Xf <- X
    ow <- getOption("warn")
    
    weight <- weight[rowSums(Y) != 0]
    X <- as.matrix(X[rowSums(Y) != 0, ])
    Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0])
    d <- ncol(Y)
    m <- rowSums(Y)
    p <- ncol(X)
    N <- nrow(Y)
    beta <- init
    lliter <- rep(0, maxiters)
    lliter[1] <- dmultn_L(init)
    niter <- 1
    options(warn = -1)
    
    ## ----------------------------------------## 
    ## Begin the main loop
    ## ----------------------------------------## 
    while (((niter <= 2) || ((ll2 - ll1)/(abs(ll1) + 1) > epsilon)) & (niter < maxiters)) {
        
        niter <- niter + 1
        ll1 <- lliter[niter - 1]
        
        ## ----------------------------------------## 
        ## Newton update
        ## ----------------------------------------## 
        P <- exp(X %*% beta)
        P <- P/(rowSums(P) + 1)
        dl <- Y[, -d] - P * m
        if (d > 2) {
            dl <- colSums(kr(dl, X, weight))
        } else if (d == 2) {
            dl <- colSums(dl * weight * X)
        }
        # dl <- Reduce('+', lapply(1:N, dl.MN_fctn, beta) )
        H <- Reduce("+", lapply(1:N, H.MN_fctn, beta))
        update <- NULL
        try(update <- solve(H, dl), silent = TRUE)
        if (is.null(update)) {
            ll.Newton <- NA
        } else if (is.numeric(update)) {
            beta_Newton <- beta - matrix(update, p, d - 1)
            ll.Newton <- dmultn_L(beta_Newton)
            ## ----------------------------------------## 
            ## Half stepping
            ## ----------------------------------------## 
            if (is.nan(ll.Newton) || ll.Newton >= 0) {
                ll.Newton <- NA
            } else if (!is.na(ll.Newton) & ll1 >= ll.Newton) {
                for (step in 1:20) {
                  beta_N <- beta - matrix(update/(2^step), p, d - 1)
                  llnew <- dmultn_L(beta_N)
                  if (is.nan(llnew) || llnew > 0) {
                    next
                  } else if (llnew > ll.Newton) {
                    ll.Newton <- llnew
                    beta_Newton <- beta_N
                  }
                  if (is.na(llnew) | llnew > ll1) {
                    break
                  }
                }
            }
            if (is.nan(ll.Newton) | ll.Newton >= 0) 
                ll.Newton <- NA
        } else {
            ll.Newton <- NA
        }
        
        if (is.na(ll.Newton) || ll.Newton < ll1) {
            ## ----------------------------------------## MM update
            denominator <- 1 + rowSums(exp(X %*% beta))
            a <- m/denominator
            beta_MM <- matrix(0, p, (d - 1))
            if (!parallel) {
                for (i in 1:(d - 1)) {
                  beta_MM[, i] <- glm.fit(X, Y[, i]/a, weights = weight * a, family = poisson(link = "log"))$coefficients
                }
            } else {
                y.list <- split(Y/a, rep(1:d, each = nrow(Y)))
                y.list[[d]] <- NULL
                if (sys == "Windows") {
                  fit.list <- clusterMap(cl, glm.private, y.list, .scheduling = "dynamic", 
                    MoreArgs = list(weights = weight * a, X = X, family = poisson(link = "log")))
                } else if (sys != "Windows") {
                  fit.list <- mcmapply(cl, glm.private, y.list, MoreArgs = list(weights = weight * 
                    a, X = X, family = poisson(link = "log")), mc.cores = cores, 
                    mc.preschedule = FALSE)
                }
                beta_MM <- do.call("cbind", fit.list)
            }
            ll.MM <- dmultn_L(beta_MM)
            ## ----------------------------------------## Choose the update
            if (is.na(ll.Newton) | (ll.MM < 0 & ll.MM > ll1)) {
                if (display) 
                  print(paste("Iteration ", niter, " MM update", ll.MM, sep = ""))
                beta <- beta_MM
                lliter[niter] <- ll.MM
                ll2 <- ll.MM
            }
        } else {
            if (display) 
                print(paste("Iteration ", niter, " Newton's update", ll.Newton, sep = ""))
            beta <- beta_Newton
            lliter[niter] <- ll.Newton
            ll2 <- ll.Newton
        }
    }
    ## ----------------------------------------## 
    ## End of the main loop
    ## ----------------------------------------## 
    options(warn = ow)
    
    ## ----------------------------------------## 
    ## Compute output statistics
    ## ----------------------------------------## 
    BIC <- -2 * lliter[niter] + log(N) * p * (d - 1)
    AIC <- 2 * p * (d - 1) - 2 * lliter[niter]
    P <- exp(X %*% beta)
    P <- cbind(P, 1)
    P <- P/rowSums(P)
    if (d > 2) {
        H <- kr(P[, 1:(d - 1)], X)
    } else if (d == 2) {
        H <- P[, 1] * X
    }
    H <- t(H) %*% (weight * m * H)
    for (i in 1:(d - 1)) {
        id <- (i - 1) * p + (1:p)
        H[id, id] = H[id, id] - t(X) %*% (weight * m * P[, i] * X)
    }
    
    SE <- matrix(NA, p, (d - 1))
    wald <- rep(NA, p)
    wald.p <- rep(NA, p)
    
    
    ## ----------------------------------------## 
    ## Check the estimate
    ## ----------------------------------------## 
    dl <- Reduce("+", lapply(1:N, dl.MN_fctn, beta))
    if (mean(dl^2) > 0.001) {
        warning(paste("The algorithm doesn't converge within", niter, "iterations. Please check gradient.", 
            "\n", sep = " "))
    }
    eig <- eigen(H)$values
    if (any(is.complex(eig)) || any(eig > 0)) {
        warning("The estimate is a saddle point.")
    } else if (any(eig == 0)) {
        warning(paste("The Hessian matrix is almost singular.", "\n", sep = ""))
    } else {
        ## ----------------------------------------## 
        ## If Hessian is negative definite, estimate SE and Wald
        ## ----------------------------------------## 
        Hinv <- chol2inv(chol(-H))
        SE <- matrix(sqrt(diag(Hinv)), p, (d - 1))
        wald <- rep(0, p)
        wald.p <- rep(0, p)
        for (j in 1:p) {
            id <- c(0:(d - 2)) * p + j
            wald[j] <- beta[j, ] %*% solve(Hinv[id, id]) %*% beta[j, ]
            wald.p[j] <- pchisq(wald[j], (d - 1), lower.tail = FALSE)
        }
    }
    
    
    ## ----------------------------------------## 
    ## Fitted value
    ## ----------------------------------------## 
    fitted <- exp(Xf %*% beta)/(rowSums(exp(Xf %*% beta)) + 1)
    fitted <- cbind(fitted, (1 - rowSums(fitted)))
    fitted[emptyRow, ] <- 0
    
    list(coefficients = beta, SE = SE, Hessian = H, wald.value = wald, wald.p = wald.p, 
        DoF = p * (d - 1), logL = lliter[niter], BIC = BIC, AIC = AIC, iter = niter, 
        gradients = dl, fitted = fitted, cl = cl)
    
}

## ============================================================## 
## Function 2: Dirichlet Multinomial regression
## ============================================================##
objfun <- function(alpha, x, y, d, p) {
    alpha <- matrix(alpha, p, d)
    alpha <- exp(x %*% alpha)
    m <- rowSums(y)
    
    logl <- rowSums(lgamma(y + alpha) - lgamma(alpha)) + lgamma(rowSums(alpha)) - 
        lgamma(rowSums(alpha) + m) + lgamma(m + 1) - rowSums(lgamma(y + 1))
    return(-sum(logl))
}

objfun.grad <- function(alpha, x, y, d, p) {
    alpha <- matrix(alpha, p, d)
    Beta <- exp(x %*% alpha)
    m <- rowSums(y)
    tmpvector <- digamma(rowSums(Beta) + m) - digamma(rowSums(Beta))
    tmpvector[is.nan(tmpvector)] <- 0
    tmpmatrix <- digamma(Beta + y) - digamma(Beta)
    tmpvector2 <- trigamma(rowSums(Beta)) - trigamma(m + rowSums(Beta))
    tmpmatrix2 <- trigamma(Beta) - trigamma(Beta + y)
    tmpmatrix2 <- Beta * tmpmatrix - Beta^2 * tmpmatrix2
    dalpha <- Beta * (tmpmatrix - tmpvector)
    
    expr <- paste("rbind(", paste(rep("dalpha", p), collapse = ","), ")", sep = "")
    dalpha <- eval(parse(text = expr))
    dalpha <- matrix(c(dalpha), nrow(x), ncol = p * d)
    expr2 <- paste("cbind(", paste(rep("x", d), collapse = ","), ")", sep = "")
    x <- eval(parse(text = expr2))
    dl <- colSums(dalpha * x)
    return(-dl)
}

objfun.hessian <- function(alpha, x, y, d, p) {
    alpha <- matrix(alpha, p, d)
    Beta <- exp(x %*% alpha)
    m <- rowSums(y)
    tmpvector <- digamma(rowSums(Beta) + m) - digamma(rowSums(Beta))
    tmpvector[is.nan(tmpvector)] <- 0
    tmpmatrix <- digamma(Beta + y) - digamma(Beta)
    tmpvector2 <- trigamma(rowSums(Beta)) - trigamma(m + rowSums(Beta))
    tmpmatrix2 <- trigamma(Beta) - trigamma(Beta + y)
    tmpmatrix2 <- Beta * tmpmatrix - Beta^2 * tmpmatrix2
    
    ## kr
    expr <- paste("rbind(", paste(rep("Beta", p), collapse = ","), ")", sep = "")
    Beta1 <- eval(parse(text = expr))
    Beta1 <- matrix(c(Beta1), nrow(x), ncol = p * d)
    expr2 <- paste("cbind(", paste(rep("x", d), collapse = ","), ")", sep = "")
    x1 <- eval(parse(text = expr2))
    Hessian <- Beta1 * x1
    Hessian <- t(Hessian) %*% (tmpvector2 * Hessian)
    
    for (i in 1:d) {
        idx <- (i - 1) * p + (1:p)
        Hessian[idx, idx] <- Hessian[idx, idx] - t(x) %*% ((tmpvector * Beta[, i] - 
            tmpmatrix2[, i]) * x)
    }
    return(-Hessian)
}


DMD.DM.reg <- function(Y, init, X, weight, epsilon, maxiters, display, parallel, 
    cores, cl, sys) {
    
    ## ----------------------------------------## 
    ## Keep some original values
    ## ----------------------------------------## 
    ow <- getOption("warn")
    fitted <- matrix(NA, nrow(Y), ncol(Y))
    emptyRow <- rowSums(Y) == 0
    Xf <- X
    
    weight <- weight[rowSums(Y) != 0]
    X <- as.matrix(X[rowSums(Y) != 0, ])
    Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0])
    d <- ncol(Y)
    p <- ncol(X)
    m <- rowSums(Y)
    N <- nrow(Y)
    beta <- init
    Beta <- exp(X %*% init)
    lliter <- rep(0, maxiters)
    lliter[1] <- sum(ddirm(Y, Beta) * weight, na.rm = TRUE)
    ll2 <- lliter[1]
    options(warn = -1)
    niter <- 1
    
    ## ----------------------------------------##
    ## Begin the main loop
    ## ----------------------------------------##
    while ((niter <= 2 || abs(ll2 - ll1)/(abs(ll1) + 1) > epsilon) & (niter < maxiters)) {
        niter <- niter + 1
        ll1 <- lliter[niter - 1]
        tmpvector <- digamma(rowSums(Beta) + m) - digamma(rowSums(Beta))
        tmpvector[is.nan(tmpvector)] <- 0
        tmpmatrix <- digamma(Beta + Y) - digamma(Beta)
        ## ----------------------------------------## 
        ## Newton update tmpvector2 <-
        ## trigamma(rowSums(Beta)) - trigamma(m+rowSums(Beta)) tmpmatrix2 <-
        ## trigamma(Beta) - trigamma(Beta+Y) tmpmatrix2 <- Beta*tmpmatrix -
        ## Beta^2*tmpmatrix2 Hessian <- kr(Beta, X) Hessian <- t(Hessian)%*%(
        ## (tmpvector2*weight)*Hessian ) for(i in 1:d){ idx <- (i-1)*p+(1:p) Hessian[idx,
        ## idx] <- Hessian[idx, idx] - t(X) %*% ( weight*((tmpvector*Beta[,i])-
        ## tmpmatrix2[,i])*X) } dalpha <- Beta*(tmpmatrix - tmpvector) dl <- colSums(
        ## kr(dalpha, X, weight) )
        Hessian <- -objfun.hessian(c(beta), X, Y, d, p)
        dl <- -objfun.grad(c(beta), X, Y, d, p)
        if (all(!is.na(dl)) & mean(dl^2) < 1e-04) 
            break
        temp.try <- NULL
        try(temp.try <- solve(Hessian, dl), silent = TRUE)
        if (is.null(temp.try) | any(is.nan(temp.try))) {
            ll.Newton <- NA
        } else if (is.numeric(temp.try)) {
            beta_Newton <- beta - matrix(temp.try, p, d)
            ll.Newton <- sum(ddirm(Y, exp(X %*% beta_Newton)) * weight, na.rm = TRUE)
            ## ----------------------------------------## Half stepping
            if (is.nan(ll.Newton) || ll.Newton >= 0) {
                ll.Newton <- NA
            } else if (!is.na(ll.Newton) & ll1 >= ll.Newton) {
                for (st in 1:20) {
                  beta_N <- beta - matrix(temp.try * (0.5^st), p, d)
                  llnew <- sum(ddirm(Y, exp(X %*% beta_N)) * weight)
                  if (is.na(llnew) | is.nan(llnew) | llnew > 0) {
                    next
                  } else if (llnew > ll.Newton) {
                    ll.Newton <- llnew
                    beta_Newton <- beta_N
                  }
                  if (!is.na(llnew) & llnew > ll1) {
                    break
                  }
                }
            }
        } else {
            ll.Newton <- NA
        }
        if (is.na(ll.Newton) || ll.Newton < ll1) {
            ## ----------------------------------------## 
            ## MM update
            ## ----------------------------------------## 
            beta_MM <- beta
            weight.fit <- weight * tmpvector
            wnz <- weight.fit != 0 & weight.fit != Inf
            weight.fit <- weight.fit[wnz]
            X1 <- X[wnz, ]
            Y_new <- Beta * tmpmatrix
            Y_new <- Y_new[wnz, ]/weight.fit
            if (!parallel) {
                for (j in 1:d) {
                  ## Surrogate 1 Poisson Regression
                  wy <- Y_new[, j]
                  wy[is.na(wy)] <- Y_new[is.na(wy), j]
                  ff <- glm.fit(X1, Y_new[, j], weights = weight.fit, family = poisson(link = "log"), 
                    control = list(epsilon = 1e-08))
                  # a <- b <- NA ff <- nlminb(rep(0.1, p), ll.obj, gradient=ll.grad,
                  # #hessian=ll.hessian, y=wy, X=X1, w=weight.fit) b <- -ff$objective a <-
                  # -ll.obj(beta[,j], wy, X1, weight.fit) b <- -ll.obj(ff$coefficients, wy, X1,
                  # weight.fit)
                  if (ff$converged) 
                    beta_MM[, j] <- ff$coefficients  #par
                  # print(paste( ff$converged, sep=' *** '))
                }
            } else {
                weight.fit[weight.fit == 0] <- 1
                y.list <- split(Y_new/weight.fit, rep(1:d, each = nrow(Y)))
                if (sys[1] == "Windows") {
                  fit.list <- clusterMap(cl, glm.private, y.list, .scheduling = "dynamic", 
                    MoreArgs = list(weights = weight.fit, X = X, family = poisson(link = "log")))
                } else if (sys[1] != "Windows") {
                  fit.list <- mcmapply(cl, glm.private, y.list, MoreArgs = list(weights = weight.fit, 
                    X = X, family = poisson(link = "log")), mc.cores = cores, mc.preschedule = FALSE)
                }
                beta_MM <- do.call("cbind", fit.list)
            }
            ll.MM <- sum(ddirm(Y, exp(X %*% beta_MM)) * weight, na.rm = TRUE)
            # print(paste(ll1, ll.MM, ll1-ll.MM, sep=' '))
            ## ----------------------------------------## 
            ## Choose the update
            ## ----------------------------------------## 
            if (is.na(ll.Newton) | (ll.MM < 0 & ll.MM > ll1)) {
                if (display) 
                  print(paste("Iteration ", niter, " MM update, log-likelihood ", 
                    ll.MM, sep = ""))
                beta <- beta_MM
                Beta <- exp(X %*% beta_MM)
                lliter[niter] <- ll.MM
                ll2 <- ll.MM
            }
        } else {
            if (display) 
                print(paste("Iteration ", niter, " Newton's update, log-likelihood", 
                  ll.Newton, sep = ""))
            beta <- beta_Newton
            Beta <- exp(X %*% beta_Newton)
            lliter[niter] <- ll.Newton
            ll2 <- ll.Newton
        }
    }
    ## ----------------------------------------## 
    ## End of the main loop
    ## ----------------------------------------## 
    options(warn = ow)
    
    ## ----------------------------------------## 
    ## Compute some output statistics
    ## ----------------------------------------## 
    BIC <- -2 * ll2 + log(N) * p * d
    AIC <- -2 * ll2 + 2 * p * d
    tmpvector <- digamma(rowSums(Beta) + m) - digamma(rowSums(Beta))
    tmpmatrix <- digamma(Beta + Y) - digamma(Beta)
    tmpvector2 <- trigamma(rowSums(Beta)) - trigamma(m + rowSums(Beta))
    tmpmatrix2 <- trigamma(Beta) - trigamma(Beta + Y)
    tmpmatrix2 <- Beta * tmpmatrix - Beta^2 * tmpmatrix2
    SE <- matrix(NA, p, d)
    wald <- rep(NA, p)
    wald.p <- rep(NA, p)
    H <- matrix(NA, p * d, p * d)
    
    ## ---------------------------------------------------------------## 
    ## Check diverge
    ## ---------------------------------------------------------------## 
    if (any(Beta == Inf, is.nan(Beta), is.nan(tmpvector2), is.nan(tmpvector), is.nan(tmpmatrix2))) {
        warning(paste("Out of range of trigamma. \n                  No standard error or tests results reported.\n                  Regression parameters diverge.  Recommend multinomial logit model.", 
            "\n", sep = ""))
    } else {
        ## ----------------------------------------## 
        ## Check gradients
        ## ----------------------------------------## 
        dalpha <- Beta * (tmpmatrix - tmpvector)
        dl <- colSums(kr(dalpha, X, weight))
        if (mean(dl^2) > 1e-04) {
            warning(paste("The algorithm doesn't converge within", niter, "iterations. The norm of the gradient is", 
                sqrt(sum(dl^2)), ". Please interpret hessian matrix and MLE with caution.", 
                "\n", sep = " "))
        }
        ## -----------------------------------## 
        ## Check whether H is negative definite, H
        ## -----------------------------------## 
        ## <- sapply(1:N, function(i, a, b) return(a[i, ]%x%b[i, ]),Beta, X)
        H <- kr(Beta, X)
        H <- t(H) %*% (H * (tmpvector2))
        for (i in 1:d) {
            id <- (i - 1) * p + c(1:p)
            H[id, id] <- H[id, id] + t(X) %*% ((-tmpvector * Beta[, i] + tmpmatrix2[, 
                i]) * X)
        }
        eig <- eigen(H)$values
        if (any(is.complex(eig)) || any(eig > 0) || any(abs(eig) < 1e-04)) {
            warning("The estimate is a saddle point. The hessian matrix is not negative definite.")
        } else if (all(eig < 0)) {
            ## --------------------------------------## 
            ## Calculate SE and wald
            ## --------------------------------------## 
            Hinv <- chol2inv(chol(-H))
            SE <- matrix(sqrt(diag(Hinv)), p, d)
            wald <- rep(NA, p)
            wald.p <- rep(NA, p)
            for (j in 1:p) {
                id <- c(0:(d - 1)) * p + j
                invH <- NULL
                try(invH <- solve(Hinv[id, id], beta[j, ]), silent = TRUE)
                if (is.numeric(invH)) {
                  wald[j] <- sum(beta[j, ] * invH)
                  wald.p[j] <- pchisq(wald[j], d, lower.tail = FALSE)
                }
            }
        }
    }
    
    fitted <- exp(Xf %*% beta)/rowSums(exp(Xf %*% beta))
    fitted[emptyRow, ] <- 0
    
    list(coefficients = beta, SE = SE, Hessian = H, wald.value = wald, wald.p = wald.p, 
        DoF = p * d, logL = ll2, BIC = BIC, AIC = AIC, iter = niter, gradients = dl, 
        fitted = fitted, logLiter = lliter)
}

## ============================================================## 
## Function 3 fit GDM Reg 
## ============================================================##

DMD.GDM.reg <- function(Y, init, X, weight, epsilon, maxiters, display, parallel, 
    cores, cl, sys) {
    
    ## ----------------------------------------## 
    ## Keep some original values
    ## ----------------------------------------## 
    fitted <- matrix(NA, nrow(Y), ncol(Y))
    emptyRow <- rowSums(Y) == 0
    Xf <- X
    ow <- getOption("warn")
    
    weight <- weight[rowSums(Y) != 0]
    X <- as.matrix(X[rowSums(Y) != 0, ])
    Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0])
    d <- ncol(Y)
    p <- ncol(X)
    m <- rowSums(Y)
    N <- nrow(Y)
    
    colOrder <- order(colSums(Y), decreasing = TRUE)
    Y <- Y[, colOrder]
    outOrder <- order(colOrder[-d])
    Ys <- t(apply(apply(apply(Y, 1, rev), 2, cumsum), 2, rev))
    
    if (is.null(init)) {
        init <- matrix(0, p, 2 * (d - 1))
        options(warn = -1)
        for (j in 1:(d - 1)) {
            den <- Ys[, j]
            den[den == 0] <- m[den == 0]
            fit <- glm.fit(X, Y[, j]/den, family = binomial(link = "logit"))
            init[, c(j, j + d - 1)] <- fit$coefficients
        }
        options(warn = ow)
    }
    
    alpha <- beta <- alpha_se <- beta_se <- matrix(0, p, (d - 1))
    gradients <- matrix(0, 2 * p, (d - 1))
    niter <- rep(0, (d - 1))
    options(warn = -1)
    if (!parallel) {
        for (j in 1:(d - 1)) {
            YGDM <- cbind(Y[, j], Ys[, (j + 1)])
            reg.result <- eval(call("DMD.DM.reg", Y = YGDM, init = init[, c(j, j + 
                d - 1)], X = X, weight = weight, epsilon = epsilon, maxiters = maxiters, 
                display = display, parallel = FALSE, cores = NULL, cl = NULL, sys = NULL))
            alpha[, j] <- reg.result$coefficients[, 1]
            beta[, j] <- reg.result$coefficients[, 2]
            alpha_se[, j] <- reg.result$SE[, 1]
            beta_se[, j] <- reg.result$SE[, 2]
            niter[j] <- reg.result$iter
            gradients[, j] <- reg.result$gradients
        }
    } else {
        Y.list <- list()
        init.list <- list()
        for (i in 1:(d - 1)) {
            Y.list[[i]] <- cbind(Y[, i], Ys[, (i + 1)])
            init.list[[i]] <- init[, c(i, i + d - 1)]
        }
        if (sys[1] == "Windows") {
            fit.list <- clusterMap(cl, DMD.DM.reg, Y.list, init.list, .scheduling = "dynamic", 
                MoreArgs = list(X = X, weight = weight, epsilon = epsilon, maxiters = maxiters, 
                  display = display, parallel = FALSE))
        } else if (sys[1] != "Windows") {
            fit.list <- mcmapply(cl, DMD.DM.reg, Y.list, init.list, MoreArgs = list(X = X, 
                weight = weight, epsilon = epsilon, maxiters = maxiters, display = display, 
                parallel = FALSE), mc.cores = cores, mc.preschedule = FALSE)
        }
        
        for (j in 1:(d - 1)) {
            alpha[, j] <- fit.list[[j]]$coefficients[, 1]
            beta[, j] <- fit.list[[j]]$coefficients[, 2]
            alpha_se[, j] <- fit.list[[j]]$SE[, 1]
            beta_se[, j] <- fit.list[[j]]$SE[, 2]
            niter[j] <- fit.list[[j]]$iter
            gradients[, j] <- fit.list[[j]]$gradient
        }
    }
    options(warn = ow)
    ## ----------------------------------------## 
    ## End of the main loop
    ## ----------------------------------------## 
    
    ll2 <- sum(weight * dgdirm(Y, exp(X %*% alpha), exp(X %*% beta)), na.rm = TRUE)
    
    ## ----------------------------------------## 
    ## Compute output statistics
    ## ----------------------------------------## 
    BIC <- -2 * ll2 + log(N) * p * (d - 1) * 2
    AIC <- -2 * ll2 + 2 * p * (d - 1) * 2
    gradients <- cbind(gradients[1:p, ], gradients[(p + 1):(2 * p), ])
    wald <- rep(NA, p)
    wald.p <- rep(NA, p)
    SE <- matrix(Inf, p, 2 * (d - 1))
    H <- matrix(NA, 2 * (d - 1), 2 * (d - 1))
    
    ## ----------------------------------------## 
    ## Check the gradients
    ## ----------------------------------------##     
    if (mean(gradients^2) > 1e-04) {
        warning(paste("The algorithm doesn't converge within", sum(niter), "iterations. The norm of the gradient is ", 
            sum(gradients^2), " Please interpret hessian matrix and MLE with caution.", 
            sep = " "))
    }
    
    ## ----------------------------------------## 
    ## Check whether H is negative definite
    ## ----------------------------------------## 
    B1 <- exp(X %*% alpha)
    B2 <- exp(X %*% beta)
    a1 <- B1 * (digamma(B1 + Y[, -d]) - digamma(B1)) - B1^2 * (-trigamma(B1 + Y[, 
        -d]) + trigamma(B1))
    a2 <- B1 * (digamma(B1 + B2 + Ys[, -d]) - digamma(B1 + B2)) - B1^2 * (-trigamma(B1 + 
        B2 + Ys[, -d]) + trigamma(B1 + B2))
    b <- B1 * B2 * (-trigamma(B1 + B2 + Ys[, -d]) + trigamma(B1 + B2))
    d1 <- B2 * (digamma(B2 + Ys[, -1]) - digamma(B2)) - B2^2 * (-trigamma(B2 + Ys[, 
        -1]) + trigamma(B2))
    d2 <- B2 * (digamma(B1 + B2 + Ys[, -d]) - digamma(B1 + B2)) - B2^2 * (-trigamma(B1 + 
        B2 + Ys[, -d]) + trigamma(B1 + B2))
    
    if (any(is.nan(a1), is.nan(a2), is.nan(b), is.nan(d1), is.nan(d2))) {
        warning(paste("Out of range of trigamma. \n\t\t\tNo standard error or tests results reported.\n\t\t\tRegression parameters diverge. Recommend multinomial logit model.", 
            "\n", sep = ""))
        SE <- matrix(NA, p, 2 * (d - 1))
        wald <- rep(NA, p)
        wald.p <- rep(NA, p)
        H <- matrix(NA, p * 2 * (d - 1), p * 2 * (d - 1))
    } else {
        Ha <- matrix(0, p * (d - 1), p * (d - 1))
        for (i in 1:(d - 1)) {
            id <- (i - 1) * p + c(1:p)
            Ha[id, id] <- t(X) %*% (weight * (a1[, i] - a2[, i]) * X)
        }
        Hb <- matrix(0, p * (d - 1), p * (d - 1))
        for (i in 1:(d - 1)) {
            id <- (i - 1) * p + c(1:p)
            Hb[id, id] <- t(X) %*% (weight * b[, i] * X)
        }
        Hd <- matrix(0, p * (d - 1), p * (d - 1))
        for (i in 1:(d - 1)) {
            id <- (i - 1) * p + c(1:p)
            Hd[id, id] <- t(X) %*% (weight * (d1[, i] - d2[, i]) * X)
        }
        H <- rbind(cbind(Ha, Hb), cbind(Hb, Hd))
        eig <- eigen(H)$values
        
        if (any(is.complex(eig)) || any(eig > 0)) {
            warning(paste("The estimate is a saddle point.\n\t\t\tNo standard error estimate or test results reported.", 
                sep = ""))
        } else if (any(eig == 0)) {
            warning(paste("The hessian matrix is almost singular.\n\t\t\tNo standard error estimate or test results reported.", 
                sep = ""))
        } else if (all(eig < 0)) {
            ## ----------------------------------------------## 
            ## Calculate the Wald statistic
            ## ----------------------------------------------## 
            ## Hinv <- chol2inv(chol(-H) )
            Hinv <- solve(-H)
            SE <- matrix(sqrt(diag(Hinv)), p, 2 * (d - 1))
            for (j in 1:p) {
                id <- c(0:(2 * (d - 1) - 1)) * p + j
                invH <- NULL
                try(invH <- solve(Hinv[id, id], c(alpha[j, ], beta[j, ])), silent = TRUE)
                if (is.numeric(invH)) {
                  wald[j] <- sum(c(alpha[j, ], beta[j, ]) * invH)
                  wald.p[j] <- pchisq(wald[j], 2 * (d - 1), lower.tail = FALSE)
                }
            }
        }
    }
    
    fitted[, 1] <- exp(Xf %*% alpha[, 1])/(exp(Xf %*% alpha[, 1]) + exp(Xf %*% beta[, 
        1]))
    fitted[, 2] <- (1 - fitted[, 1]) * exp(Xf %*% alpha[, 2])/(exp(Xf %*% alpha[, 
        2]) + exp(Xf %*% beta[, 2]))
    if (d > 3) {
        for (f in 3:(d - 1)) {
            fitted[, f] <- (1 - rowSums(fitted[, 1:(f - 1)])) * exp(Xf %*% alpha[, 
                f])/(exp(Xf %*% alpha[, f]) + exp(Xf %*% beta[, f]))
        }
    } else if (d == 3) {
        fitted[, d] <- (1 - rowSums(fitted[, 1:(d - 1)]))
    }
    fitted[emptyRow, ] <- 0
    
    list(coefficients = cbind(alpha[, outOrder], beta[, outOrder]), SE = SE[, c(outOrder, 
        outOrder + d - 1)], Hessian = H[c(outOrder, outOrder + d - 1), c(outOrder, 
        outOrder + d - 1)], wald.value = wald, wald.p = wald.p, logL = ll2, BIC = BIC, 
        AIC = AIC, iter = sum(niter), DoF = 2 * p * (d - 1), gradients = gradients, 
        fitted = fitted[, order(colOrder)])
    
}


## ============================================================## 
## Function 4 fit NegMN Reg 
## ============================================================##
DMD.NegMN.reg <- function(Y, init, X, weight, epsilon, maxiters, display, parallel, 
    cores, cl, sys) {
    ow <- getOption("warn")
    fitted <- matrix(NA, nrow(Y), ncol(Y))
    emptyRow <- rowSums(Y) == 0
    Ys <- t(apply(apply(apply(Y, 1, rev), 2, cumsum), 2, rev))
    d <- ncol(Y)
    p <- ncol(X)
    m <- rowSums(Y)
    N <- nrow(Y)
    alpha <- init[, 1:d]
    beta <- init[, (d + 1)]
    Alpha <- exp(X %*% alpha)
    rowsum_Alpha <- rowSums(Alpha) + 1
    Beta <- c(exp(X %*% beta))
    lliter <- rep(NA, maxiters)
    lliter[1] <- sum(weight * dneg(Y, Alpha, Beta), na.rm = TRUE)
    niter <- 1
    options(warn = -1)
    while (((niter <= 2) || ((ll2 - ll1)/(abs(ll1) + 1) > epsilon)) & (niter < maxiters)) {
        niter <- niter + 1
        ll1 <- lliter[niter - 1]
        Alpha <- exp(X %*% alpha)
        rowsum_Alpha <- rowSums(Alpha) + 1
        Beta <- c(exp(X %*% beta))
        Prob <- cbind(Alpha/rowsum_Alpha, 1/rowsum_Alpha)
        tmpBeta <- digamma(Beta + m) - digamma(Beta)
        tmpBeta[is.nan(tmpBeta)] <- 0
        w_beta <- log(rowsum_Alpha)
        dlbeta <- colSums((tmpBeta - w_beta) * Beta * X)
        hbeta_w <- (trigamma(Beta + m) - trigamma(Beta) + tmpBeta - w_beta) * Beta
        hbeta <- sapply(1:N, function(i, A, B, w) return(w[i] * A[i, ] %x% B[i, ]), 
            X, X, hbeta_w)
        hbeta <- matrix(rowSums(hbeta), p, p)
        if (all(eigen(hbeta)$value < 0)) {
            beta_MM <- beta - solve(hbeta, dlbeta)
            Beta_MM <- c(exp(X %*% beta_MM))
            lltemp <- sum(weight * dneg(Y, Alpha, Beta_MM), na.rm = TRUE)
            if (lltemp < ll1) {
                Y_reg <- weight * Beta * tmpBeta/w_beta
                beta_MM <- glm.fit(X, Y_reg, weights = w_beta, family = poisson(link = "log"), 
                  start = beta)$coefficients
            }
        } else {
            Y_reg <- weight * Beta * tmpBeta/w_beta
            beta_MM <- glm.fit(X, Y_reg, weights = w_beta, family = poisson(link = "log"), 
                start = beta)$coefficients
        }
        Beta_MM <- c(exp(X %*% beta_MM))
        alpha_MM <- matrix(0, p, d)
        w_alpha <- (Beta_MM + m)/(rowSums(Alpha) + 1)
        if (!parallel) {
            for (j in 1:d) {
                alpha_MM[, j] <- glm.fit(X, weight * Y[, j]/w_alpha, weights = w_alpha, 
                  family = poisson(link = "log"), start = alpha[, j])$coefficients
            }
        } else {
            y.list <- split(weight * Y/w_alpha, rep(1:d, each = nrow(Y)))
            start.list <- split(alpha, rep(1:d, each = nrow(alpha)))
            if (sys[1] == "Windows") {
                fit.list <- clusterMap(cl, glm.private, y.list, start.list, .scheduling = "dynamic", 
                  MoreArgs = list(X = X, family = poisson(link = "log"), weights = w_alpha))
            } else if (sys[1] != "Windows") {
                fit.list <- mcmapply(cl, glm.private, y.list, start.list, MoreArgs = list(X = X, 
                  family = poisson(link = "log"), weights = w_alpha), mc.cores = cores, 
                  mc.preschedule = FALSE)
            }
            alpha_MM <- do.call("cbind", fit.list)
        }
        Alpha_MM <- exp(X %*% alpha_MM)
        Prob_MM <- cbind(Alpha_MM/(1 + rowSums(Alpha_MM)), 1/(1 + rowSums(Alpha_MM)))
        Beta_MM <- c(exp(X %*% beta_MM))
        ll.MM <- sum(weight * dneg(Y, Alpha_MM, Beta_MM), na.rm = TRUE)
        deta <- matrix(0, N, (d + 1))
        deta[, 1:d] <- Y - Alpha * Beta - Prob[, 1:d] * (m - Beta * (rowsum_Alpha - 
            1))
        deta[, (d + 1)] <- Beta * (digamma(Beta + m) - digamma(Beta) - log(rowsum_Alpha))
        score <- rowSums(sapply(1:N, function(i, A, B, weight) return(weight[i] * 
            A[i, ] %x% B[i, ]), deta, X, weight))
        hessian <- matrix(0, p * (1 + d), p * (1 + d))
        upleft <- sapply(1:N, function(i, A, B) return(A[i, ] %x% B[i, ]), cbind(Prob[, 
            1:d], -Beta/(Beta + m)), X)
        hessian <- t(t(upleft) * (weight * (Beta + m))) %*% t(upleft)
        for (j in 1:d) {
            idx <- (j - 1) * p + (1:p)
            hessian[idx, idx] <- hessian[idx, idx] - t(X) %*% (X * (weight * (Beta + 
                m) * Prob[, j]))
        }
        tmpvector2 <- as.vector(Beta * (tmpBeta + Beta * (trigamma(Beta + m) - trigamma(Beta)) - 
            log(rowsum_Alpha) - Beta/(Beta + m)))
        idx <- d * p + c(1:p)
        hessian[idx, idx] = hessian[idx, idx] + t(X) %*% (X * (weight * tmpvector2))
        temp.try <- NULL
        try(temp.try <- solve(hessian, score), silent = TRUE)
        if (is.null(temp.try)) {
            ll.Newton <- NA
        } else if (is.numeric(temp.try)) {
            B_Newton <- cbind(alpha, beta) - matrix(temp.try, p, (d + 1))
            B <- exp(X %*% B_Newton)
            Alpha_Newton <- B[, 1:d]
            Beta_Newton <- B[, (d + 1)]
            ll.Newton <- sum(weight * dneg(Y, Alpha_Newton, Beta_Newton), na.rm = TRUE)
            if (!is.na(ll.Newton) & ll.MM > ll.Newton) {
                for (step in 1:20) {
                  B_N <- cbind(alpha, beta) - matrix(temp.try/(2^step), p, d + 1)
                  B <- exp(X %*% B_N)
                  Alpha_N <- B[, 1:d]
                  Beta_N <- B[, (d + 1)]
                  llnew <- sum(weight * dneg(Y, Alpha_N, Beta_N), na.rm = T)
                  if (is.na(llnew) | is.nan(llnew) | llnew > 0) {
                    next
                  } else if (llnew > ll.Newton) {
                    ll.Newton <- llnew
                    B_Newton <- B_N
                  }
                  if (!is.na(llnew) & llnew > ll1) {
                    break
                  }
                }
            }
            if (is.nan(ll.Newton)) 
                ll.Newton <- NA else if (ll.Newton >= 0) 
                ll.Newton <- NA
        } else {
            ll.Newton <- NA
        }
        if (is.na(ll.Newton) | ll.MM > ll.Newton) {
            if (display) 
                print(paste("Iteration ", niter, " MM update", sep = ""))
            alpha <- alpha_MM
            beta <- beta_MM
            ll2 <- ll.MM
        } else {
            if (display) 
                print(paste("Iteration ", niter, " Newton's update", sep = ""))
            alpha <- B_Newton[, 1:d]
            beta <- B_Newton[, (d + 1)]
            ll2 <- ll.Newton
        }
        lliter[niter] <- ll2
    }
    options(warn = ow)
    BIC <- -2 * ll2 + log(N) * p * (d + 1)
    AIC <- -2 * ll2 + 2 * p * (d + 1)
    A <- exp(X %*% cbind(alpha, beta))
    Alpha <- A[, 1:d]
    Beta <- A[, (d + 1)]
    Prob <- cbind(Alpha, 1)
    Prob <- Prob/(rowSums(Prob))
    tmpv2 <- Beta + m
    tmpv1 <- digamma(tmpv2) - digamma(Beta)
    SE <- matrix(NA, p, d + 1)
    wald <- rep(NA, p)
    wald.p <- rep(NA, p)
    H <- matrix(NA, p * (d + 1), p * (d + 1))
    if (any(A == Inf, is.nan(A))) {
        warning("Out of range of trigamma.  No SE or tests results reported.\n\t\t\t\tRegression parameters diverge.Recommend multinomial logit model")
    } else {
        deta <- matrix(0, N, (d + 1))
        deta[, 1:d] <- Y - Alpha * Beta - Prob[, 1:d] * (m - Beta * (rowsum_Alpha - 
            1))
        deta[, (d + 1)] <- Beta * (digamma(Beta + m) - digamma(Beta) + log(Prob[, 
            d + 1]))
        score <- rowSums(sapply(1:N, function(i, A, B, weight) return(weight[i] * 
            A[i, ] %x% B[i, ]), deta, X, weight))
        H <- sapply(1:N, function(i, a, x) return(a[i, ] %x% x[i, ]), cbind(Prob[, 
            1:d], -Beta/tmpv2), X)
        H <- H %*% (weight * tmpv2 * t(H))
        for (i in 1:d) {
            id <- (i - 1) * p + c(1:p)
            H[id, id] <- H[id, id] - t(X) %*% (weight * tmpv2 * Prob[, i] * X)
        }
        id <- d * p + c(1:p)
        tmpv3 <- Beta * (tmpv1 + Beta * (trigamma(tmpv2) - trigamma(Beta)) + log(Prob[, 
            (d + 1)]) - Beta/tmpv2)
        H[id, id] <- H[id, id] + t(X) %*% (weight * tmpv3 * X)
        if (mean(score^2) > 1e-04) {
            warning(paste("The algorithm doesn't converge within", niter, "iterations. The norm of the gradient is ", 
                sum(score^2), " Please interpret hessian matrix and MLE with caution.", 
                sep = " "))
        }
        eig <- eigen(H)$values
        if (any(eig > 0)) {
            warning("The estimate is a saddle point.")
        } else if (any(eig == 0)) {
            warning("The hessian matrix is almost singular.")
        } else if (all(eig < 0)) {
            Hinv <- chol2inv(chol(-H))
            SE <- matrix(sqrt(diag(Hinv)), p, (d + 1))
            wald <- rep(0, p)
            wald.p <- rep(0, p)
            for (j in 1:p) {
                id <- c(0:d) * p + j
                wald[j] <- c(alpha[j, ], beta[j]) %*% chol2inv(chol(Hinv[id, id])) %*% 
                  c(alpha[j, ], beta[j])
                wald.p[j] <- pchisq(wald[j], (d + 1), lower.tail = FALSE)
            }
        }
    }
    fitted <- c(exp(X %*% beta)) * exp(X %*% alpha)/rowSums(exp(X %*% alpha))
    list(coefficients = cbind(alpha, beta), SE = SE, Hessian = H, BIC = BIC, AIC = AIC, 
        wald.value = wald, wald.p = wald.p, logL = lliter[niter], iter = (niter), 
        gradient = score, fitted = fitted)
}

# DMD.NegMN.reg <- function(Y, init, X, weight, epsilon, maxiters, display,
# parallel, cores, cl, sys){ ##----------------------------------------## ## Keep
# some original values ow <- getOption('warn') fitted <- matrix(NA, nrow(Y),
# ncol(Y)) emptyRow <- rowSums(Y)==0 Ys <- t( apply(apply(
# apply(Y,1,rev),2,cumsum),2,rev) ) d <- ncol(Y) p <- ncol(X) m <- rowSums(Y) N
# <- nrow(Y) alpha <- init[, 1:d] beta <- init[, (d+1)] Alpha <- exp(X%*%alpha)
# rowsum_Alpha <- rowSums(Alpha)+1 Beta <- c(exp(X%*%beta)) lliter <- rep(NA,
# maxiters) lliter[1] <- sum(weight*dneg(Y, Alpha,Beta), na.rm=TRUE) ll2 <-
# lliter[1] niter <- 1 options(warn=-1)
# ##---------------------------------------## ## Begin the main loop while(
# ((niter <=2)|| ((ll2-ll1)/(abs(ll1)+1) > epsilon))&(niter<maxiters) ){ niter <-
# niter+1 ll1 <- lliter[niter-1] Alpha <- exp(X%*%alpha) rowsum_Alpha <-
# rowSums(Alpha)+1 Beta <- c(exp(X%*%beta)) Prob <-
# cbind(Alpha/rowsum_Alpha,1/rowsum_Alpha) tmpBeta <-
# digamma(Beta+m)-digamma(Beta) tmpBeta[is.nan(tmpBeta)] <- 0
# ##----------------------------------------## ## Newton Update deta <- matrix(0,
# N, (d+1)) deta[, 1:d] <- Y - Alpha*Beta - Prob[, 1:d]*(m-Beta*(rowsum_Alpha-1)
# ) deta[, (d+1)] <- Beta*(digamma(Beta+m)-digamma(Beta)-log(rowsum_Alpha)) score
# <- colSums( kr(deta, X, weight) ) hessian <- matrix(0, p*(1+d), p*(1+d)) upleft
# <- kr(cbind(Prob[,1:d], -Beta/(Beta+m)),X) hessian <-
# t(upleft*(weight*(Beta+m)))%*%upleft for(j in 1:d){ idx <- (j-1)*p + (1:p)
# hessian[idx, idx] <- hessian[idx, idx]- t(X)%*%(X*(weight*(Beta+m)*Prob[,j])) }
# tmpvector2 <- as.vector(Beta*(tmpBeta+Beta*(trigamma(Beta+m)-trigamma(Beta))-
# log(rowsum_Alpha)- Beta/(Beta+m) )) idx <- d*p + c(1:p) hessian[idx, idx]=
# hessian[idx, idx]+t(X)%*%(X*(weight*tmpvector2)) temp.try <- NULL try(temp.try
# <- solve(hessian,score), silent=TRUE) if(is.null(temp.try)){ ll.Newton <- NA
# }else if(is.numeric(temp.try)){ B_Newton <- cbind(alpha, beta)-matrix(temp.try,
# p, (d+1)) B <- exp(X%*%B_Newton) Alpha_Newton <- B[, 1:d] Beta_Newton <- B[,
# (d+1)] ll.Newton <- sum(weight*dneg(Y,Alpha_Newton, Beta_Newton), na.rm=TRUE)
# ## ----------------------------------------## ## Half stepping
# if(is.nan(ll.Newton) || ll.Newton>0){ ll.Newton <- NA }else
# if(!is.na(ll.Newton)&ll1 >= ll.Newton){ for(step in 1:40){ B_N <- cbind(alpha,
# beta) - matrix(temp.try*(0.5^step), p, d+1) B <- exp(X%*%B_N) Alpha_N <- B[,
# 1:d] Beta_N <- B[, (d+1)] llnew <- sum(weight*dneg(Y, Alpha_N, Beta_N),na.rm=T)
# if(is.nan(llnew) | is.na(llnew) | llnew>0){ next }else if(llnew > ll.Newton){
# ll.Newton <- llnew B_Newton <- B_N } if(is.na(llnew) |llnew>ll1){ break } } }
# if(is.nan(ll.Newton)) ll.Newton <- NA else if(ll.Newton>=0) ll.Newton <- NA
# }else{ ll.Newton <- NA } if(is.na(ll.Newton) || ll.Newton<ll1){ ##
# ----------------------------------------## ## MM Update w_beta <-
# log(rowsum_Alpha) dlbeta <- colSums((tmpBeta - w_beta)*Beta*X) hbeta_w <-
# (trigamma(Beta+m)-trigamma(Beta)+tmpBeta-w_beta)*Beta hbeta <- kr(X, X,
# hbeta_w) hbeta <- matrix(colSums(hbeta), p, p) if( all(eigen(hbeta)$value<0) ){
# beta_MM <- beta - solve(hbeta, dlbeta) Beta_MM <- c(exp(X%*%beta_MM)) lltemp <-
# sum(weight*dneg(Y, Alpha, Beta_MM), na.rm=TRUE) if(lltemp < ll1){ Y_reg <-
# weight*Beta*tmpBeta/w_beta beta_MM <- glm.fit(X, Y_reg, weights=w_beta,
# family=poisson(link='log'), start=beta)$coefficients } }else{ Y_reg <-
# weight*Beta*tmpBeta/w_beta beta_MM <- glm.fit(X, Y_reg, weights=w_beta,
# family=poisson(link='log'), start=beta)$coefficients } Beta_MM <-
# c(exp(X%*%beta_MM)) alpha_MM <- matrix(0,p,d) w_alpha
# <-(Beta_MM+m)/(rowSums(Alpha)+1) if(!parallel){ for(j in 1:d){ alpha_MM[,j] <-
# glm.fit(X, weight*Y[,j]/w_alpha, weights=w_alpha, family=poisson(link='log'),
# start=alpha[,j])$coefficients } }else{ y.list <- split(weight*Y/w_alpha,
# rep(1:d, each=nrow(Y))) start.list <- split(alpha, rep(1:d, each=nrow(alpha)))
# if(sys[1]=='Windows'){ fit.list <- clusterMap(cl, glm.private, y.list,
# start.list, .scheduling='dynamic', MoreArgs=list(X=X, family=poisson(link
# ='log'), weights=w_alpha)) }else if(sys[1]!='Windows'){ fit.list <-
# mcmapply(cl, glm.private, y.list, start.list, MoreArgs=list(X=X,
# family=poisson(link ='log'), weights=w_alpha), mc.cores = cores, mc.preschedule
# = FALSE) } alpha_MM <- do.call('cbind', fit.list) } Alpha_MM <-
# exp(X%*%alpha_MM) Prob_MM
# <-cbind(Alpha_MM/(1+rowSums(Alpha_MM)),1/(1+rowSums(Alpha_MM))) Beta_MM <-
# c(exp(X%*%beta_MM)) ll.MM <- sum(weight*dneg(Y, Alpha_MM, Beta_MM), na.rm=TRUE)
# ## ----------------------------------------## ## Choose update if(
# is.na(ll.Newton)|ll.MM<0&ll.MM>ll1 ){ if(display) print(paste('Iteration ',
# niter, ' MM update', ll.MM, sep='')) alpha <- alpha_MM beta <- beta_MM
# lliter[niter] <- ll.MM ll2 <- ll.MM } }else{ if(display) print(paste('Iteration
# ', niter, ' Newton's update', ll.Newton, sep='')) alpha <- B_Newton[,1:d] beta
# <- B_Newton[, (d+1)] lliter[niter] <- ll.Newton ll2 <- ll.Newton } } ##
# ----------------------------------------## ## End of the main loop
# options(warn=ow) ##----------------------------------------## ## Compute output
# statistics BIC <- -2*ll2 + log(N)*p*(d+1) AIC <- -2*ll2 + 2*p*(d+1) A <-
# exp(X%*%cbind(alpha, beta)) Alpha <- A[, 1:d] Beta <- A[,(d+1)] Prob <-
# cbind(Alpha, 1) Prob <- Prob/(rowSums(Prob)) tmpv2 <- Beta+m tmpv1 <-
# digamma(tmpv2) - digamma(Beta) SE <- matrix(NA, p,d+1) wald <- rep(NA, p)
# wald.p <- rep(NA, p) H <- matrix(NA, p*(d+1), p*(d+1) )
# ##----------------------------------------## ## Check diverge if(any(A==Inf,
# is.nan(A))){ warning(paste('Out of range of trigamma.  No SE or tests results
# reported.  Regression parameters diverge.Recommend multinomial logit
# model','\n', sep='')) }else{ ##----------------------------------------## ##
# Calculate dl deta <- matrix(0, N, (d+1)) deta[, 1:d] <- Y - Alpha*Beta - Prob[,
# 1:d]*(m-Beta*(rowsum_Alpha-1) ) deta[, (d+1)] <-
# Beta*(digamma(Beta+m)-digamma(Beta)+log(Prob[, d+1])) score <- colSums(kr(deta,
# X, weight)) ##----------------------------------------## ## Calculate H H <-
# kr(cbind(Prob[, 1:d],-Beta/tmpv2), X) H <- t(H)%*%(weight*tmpv2*H) for(i in
# 1:d){ id <- (i-1)*p + c(1:p) H[id, id] <- H[id, id]-t(X)%*%(weight*tmpv2*Prob[,
# i]*X) } id <- d*p+c(1:p) tmpv3 <- Beta*( tmpv1+
# Beta*(trigamma(tmpv2)-trigamma(Beta))+ log(Prob[, (d+1)])-Beta/tmpv2) H[id, id]
# <- H[id, id]+t(X)%*%(weight*tmpv3*X)
# ##----------------------------------------## ## Check gradients
# if(mean(score^2)>1e-4){ warning(paste('The algorithm doesn't converge within',
# niter, 'iterations. The norm of the gradient is ', sum(score^2), ' Please
# interpret hessian matrix and MLE with caution.', '\n', sep=' ') ) }
# ##-----------------------------------## ## Check whether H is negative
# definite, eig <- eigen(H)$values if(any(is.complex(eig)) || any(eig>0)){
# warning(paste('The estimate is a saddle point.', '\n', sep='')) }else
# if(any(eig==0)){ warning(paste('The hessian matrix is almost singular.','\n',
# sep='')) }else if(all(eig<0)){ ##--------------------------------------## ##
# Calculate SE and wald Hinv <-chol2inv(chol(-H) ) SE <-matrix( sqrt(diag(Hinv )
# ), p, (d+1)) wald <- rep(0, p) wald.p<- rep(0, p) for(j in 1:p){ id<- c(0:
# d)*p+j wald[j] <- c(alpha[j,] ,beta[j])%*%chol2inv(chol(Hinv[id,id]))%*%
# c(alpha[j,] ,beta[j]) wald.p[j] <- pchisq(wald[j], (d+1), lower.tail=FALSE) } }
# } fitted <- c(exp(X%*%beta)) * exp(X%*%alpha)/rowSums(exp(X%*%alpha))
# list(coefficients=cbind(alpha, beta), SE=SE, Hessian=H, BIC=BIC, AIC=AIC,
# wald.value=wald, wald.p=wald.p, logL=lliter[niter], iter=(niter),
# gradients=score, fitted=fitted) }


## ============================================================## 
## Function 5 fit NegMN Reg, but do not link the over-dispersion parameter
## ============================================================##

DMD.NegMN.Alpha.reg <- function(Y, init, X, weight, epsilon, maxiters, display, parallel, 
    cores, cl, sys) {
    
    ## ----------------------------------------## 
    ## Keep some original values
    ## ----------------------------------------## 
    ow <- getOption("warn")
    fitted <- matrix(NA, nrow(Y), ncol(Y))
    emptyRow <- rowSums(Y) == 0
    
    Ys <- t(apply(apply(apply(Y, 1, rev), 2, cumsum), 2, rev))
    d <- ncol(Y)
    p <- ncol(X)
    m <- rowSums(Y)
    N <- nrow(Y)
    
    alpha <- init
    Alpha <- exp(X %*% alpha)
    rowsum_Alpha <- rowSums(Alpha) + 1
    beta <- DMD.NegMN.fit(Y, weight = weight)
    beta <- beta$estimate[(d + 1)]
    if (beta < 0) 
        beta <- 1
    Beta <- rep(beta, N)
    lliter <- rep(NA, maxiters)
    lliter[1] <- sum(weight * dneg(Y, Alpha, Beta), na.rm = TRUE)
    ll2 <- lliter[1]
    betaiter <- rep(NA, maxiters)
    betaiter[1] <- beta
    niter <- 1
    div <- FALSE
    options(warn = -1)
    ## ----------------------------------------## 
    ## Begin the main loop
    ## ----------------------------------------## 
    while (((niter <= 2) || ((ll2 - ll1)/(abs(ll1) + 1) > epsilon)) & (niter < maxiters)) {
        
        niter <- niter + 1
        ll1 <- lliter[niter - 1]
        Alpha <- exp(X %*% alpha)
        rowsum_Alpha <- rowSums(Alpha) + 1
        Prob <- cbind(Alpha/rowsum_Alpha, 1/rowsum_Alpha)
        tmpBeta <- digamma(beta + m) - digamma(beta)
        tmpBeta[is.nan(tmpBeta)] <- 0
        
        dlbeta <- sum(tmpBeta - log(rowsum_Alpha))
        hbeta <- sum(trigamma(beta + m) - trigamma(beta))
        ## ----------------------------------------## 
        ## Newton Update
        deta <- Y - (beta + m) * Prob[, 1:d]
        score <- colSums(kr(deta, X, weight))
        score <- c(score, dlbeta)
        upleft <- kr(Prob[, 1:d], X)
        upright <- -colSums(upleft)
        upleft <- t(upleft) %*% (upleft * (weight * (beta + m)))
        for (j in 1:d) {
            idx <- (j - 1) * p + 1:p
            upleft[idx, idx] <- upleft[idx, idx] - t(X) %*% (X * (weight * (beta + 
                m) * Prob[, j]))
        }
        hessian <- rbind(cbind(upleft, upright), c(upright, hbeta))
        temp.try <- NULL
        try(temp.try <- solve(hessian, score), silent = TRUE)
        if (is.null(temp.try)) {
            ll.Newton <- NA
        } else if (is.numeric(temp.try)) {
            beta_Newton <- beta - temp.try[length(temp.try)]
            B_Newton <- alpha - matrix(temp.try[1:p * d], p, d)
            Alpha_Newton <- exp(X %*% B_Newton)
            ll.Newton <- sum(weight * dneg(Y, Alpha_Newton, rep(beta_Newton, N)), 
                na.rm = TRUE)
        }
        ## ----------------------------------------## 
        ## Half stepping
        if (is.nan(ll.Newton) || ll.Newton >= 0) {
            ll.Newton <- NA
        } else if (!is.na(ll.Newton) & ll1 >= ll.Newton) {
            for (step in 1:20) {
                temp <- temp.try/(2^step)
                B_N <- alpha - matrix(temp[1:(p * d)], p, d)
                Alpha_N <- exp(X %*% B_N)
                beta_N <- beta - temp[p * d + 1]
                llnew <- sum(weight * dneg(Y, Alpha_N, rep(beta_N, N)), na.rm = T)
                if (is.na(llnew) | is.nan(llnew) | llnew > 0) {
                  next
                } else if (llnew > ll.Newton) {
                  ll.Newton <- llnew
                  B_Newton <- B_N
                  beta_Newton <- beta_N
                }
                if (!is.na(llnew) & llnew > ll1) {
                  break
                }
                if (is.nan(ll.Newton) | ll.Newton >= 0) 
                  ll.Newton <- NA
            }
            if (is.nan(ll.Newton)) 
                ll.Newton <- NA else if (ll.Newton >= 0) 
                ll.Newton <- NA
        }
        
        if (is.na(ll.Newton) || ll.Newton < ll1) {
            ## ----------------------------------------## 
            ## MM Update
            dlbeta <- sum(tmpBeta - log(rowsum_Alpha))
            hbeta <- sum(trigamma(beta + m) - trigamma(beta))
            beta_MM <- beta
            if (hbeta != 0) {
                temp_MM <- beta - dlbeta/hbeta
                ## check ll increase
                if (!is.na(temp_MM) && !is.nan(temp_MM) && temp_MM > 0) 
                  beta_MM <- temp_MM
            }
            alpha_MM <- matrix(0, p, d)
            w_alpha <- (beta_MM + m)/rowsum_Alpha
            if (sum(w_alpha^2) == 0) 
                break else w_alpha[w_alpha == 0] <- 1
            if (any(is.na(w_alpha))) 
                stop("The algorithm diverged. Please try other model.")
            if (!parallel) {
                for (j in 1:d) {
                  alpha_MM[, j] <- glm.fit(X, weight * Y[, j]/w_alpha, weights = w_alpha, 
                    family = poisson(link = "log"), start = alpha[, j])$coefficients
                }
            } else {
                y.list <- split(weight * Y/w_alpha, rep(1:d, each = nrow(Y)))
                start.list <- split(alpha, rep(1:d, each = nrow(alpha)))
                if (sys[1] == "Windows") {
                  fit.list <- clusterMap(cl, glm.private, y.list, start.list, .scheduling = "dynamic", 
                    MoreArgs = list(X = X, family = poisson(link = "log"), weights = w_alpha))
                } else if (sys[1] != "Windows") {
                  fit.list <- mcmapply(cl, glm.private, y.list, start.list, MoreArgs = list(X = X, 
                    family = poisson(link = "log"), weights = w_alpha), mc.cores = cores, 
                    mc.preschedule = FALSE)
                }
                alpha_MM <- do.call("cbind", fit.list)
            }
            Alpha_MM <- exp(X %*% alpha_MM)
            Prob_MM <- cbind(Alpha_MM/(1 + rowSums(Alpha_MM)), 1/(1 + rowSums(Alpha_MM)))
            ll.MM <- sum(weight * dneg(Y, Alpha_MM, rep(beta_MM, N)), na.rm = TRUE)
            
            ## ----------------------------------------## 
            ## Choose update
            ## ----------------------------------------## 
            if (is.na(ll.Newton) | ll.MM > ll1) {
                if (beta_MM <= 0) {
                  warning(paste("The estimate of overdispersion parameter is ,", 
                    beta, ". It is smaller or eqal to zero.  Please consider other models.", 
                    sep = " "))
                  div <- TRUE
                  break
                }
                beta <- beta_MM
                alpha <- alpha_MM
                lliter[niter] <- ll.MM
                ll2 <- ll.MM
                if (display) 
                  print(paste("Iteration ", niter, " MM update, log-likelihood ", 
                    ll.MM, sep = ""))
            }
        } else {
            if (beta_Newton <= 0) {
                warning(paste("The estimate of overdispersion parameter is ,", beta, 
                  ". It is smaller or eqal to zero.  Please consider other models.", 
                  sep = " "))
                div <- TRUE
                break
            }
            alpha <- B_Newton
            beta <- beta_Newton
            lliter[niter] <- ll.Newton
            ll2 <- ll.Newton
            if (display) 
                print(paste("Iteration ", niter, " Newton's update, log-likelihood ", 
                  ll.Newton, sep = ""))
        }
        if (div) 
            break
    }
    ## ----------------------------------------## 
    ## End of the main loop
    ## ----------------------------------------## 
    options(warn = ow)
    
    ## ----------------------------------------## 
    ## Compute output statistics
    ## ----------------------------------------## 
    BIC <- -2 * ll2 + log(N) * (p * d + 1)
    AIC <- -2 * ll2 + 2 * (p * d + 1)
    A <- exp(X %*% alpha)
    Alpha <- A[, 1:d]
    Prob <- cbind(Alpha, 1)
    Prob <- Prob/(rowSums(Prob))
    tmpv1 <- digamma(beta + m) - digamma(beta)
    SE <- matrix(NA, p, d)
    SE_beta <- NA
    wald <- rep(NA, p)
    wald.p <- rep(NA, p)
    H <- matrix(NA, p * (d + 1), p * (d + 1))
    
    ## ----------------------------------------## 
    ## Check diverge
    ## ----------------------------------------## 
    if (any(A == Inf, is.nan(A), is.nan(beta), is.nan(beta))) {
        warning(paste("Out of range of trigamma.  No SE or tests results reported.\n\t\t\t\tRegression parameters diverge.Recommend multinomial logit model", 
            "\n", sep = ""))
    } else {
        ## ----------------------------------------## 
        ## Calculate dl
        ## ----------------------------------------## 
        dlbeta <- sum(tmpv1 - log(rowSums(Alpha) + 1))
        deta <- Y - (beta + m) * Prob[, 1:d]
        # score <- rowSums(sapply(1:N, function(i, A, B, weight)
        # return(weight[i]*A[i,]%x%B[i,]), deta, X, weight))
        score <- colSums(kr(deta, X, weight))
        score <- c(score, dlbeta)
        ## ----------------------------------------## 
        ## Calculate H upleft <- sapply(1:N,
        ## function(i,A,B) return(A[i,]%x%B[i,]), Prob[,1:d], X)
        upleft <- kr(Prob[, 1:d], X)
        upright <- -colSums(kr(Prob[, 1:d], X))
        upleft <- t(upleft * (weight * (beta + m))) %*% upleft
        for (j in 1:d) {
            idx <- (j - 1) * p + 1:p
            upleft[idx, idx] <- upleft[idx, idx] - t(X) %*% (X * (weight * (beta + 
                m) * Prob[, j]))
        }
        H <- rbind(cbind(upleft, upright), c(upright, hbeta))
        
        ## ----------------------------------------## 
        ## Check gradients
        ## ----------------------------------------## 
        if (mean(score[1:(p * d)]^2) > 0.01) {
            warning(paste("The algorithm doesn't converge within", niter, "iterations. The norm of the gradient is ", 
                sum(score^2), " Please interpret hessian matrix and MLE with caution.", 
                "\n", sep = " "))
        }
        ## -----------------------------------## 
        ## Check whether H is negative definite,
        ## -----------------------------------## 
        eig <- eigen(H)$values
        if (any(is.complex(eig)) || any(eig > 0)) {
            warning(paste("The estimate is a saddle point.", "\n", sep = ""))
        } else if (any(eig == 0)) {
            warning(paste("The hessian matrix is almost singular.", "\n", sep = ""))
        } else if (all(eig < 0)) {
            ## --------------------------------------## 
            ## Calculate SE and wald
            ## --------------------------------------## 
            Hinv <- chol2inv(chol(-H))
            SE <- matrix(sqrt(diag(Hinv))[1:p * d], p, d)
            SE_beta <- Hinv[length(Hinv)]
            wald <- rep(0, p)
            wald.p <- rep(0, p)
            for (j in 1:p) {
                id <- c(0:(d - 1)) * p + j
                wald[j] <- alpha[j, ] %*% chol2inv(chol(Hinv[id, id])) %*% alpha[j, 
                  ]
                wald.p[j] <- pchisq(wald[j], d, lower.tail = FALSE)
            }
        }
    }
    
    fitted <- beta * exp(X %*% alpha)/rowSums(exp(X %*% alpha))
    
    list(coefficients = list(alpha = alpha, phi = beta), SE = list(SE.alpha = SE, 
        SE.beta = SE_beta), Hessian = H, BIC = BIC, AIC = AIC, wald.value = wald, 
        wald.p = wald.p, logL = lliter[niter], iter = (niter), gradients = score, 
        fitted = fitted)
}

## ============================================================## 
## re-organize the glm.fit function 
## ============================================================##
glm.private <- function(Y, start = NULL, weights, X, family) {
    fit <- glm.fit(x = X, y = Y, weights = weights, family = poisson(link = "log"), 
        start = start)
    return(fit$coefficients)
} 
