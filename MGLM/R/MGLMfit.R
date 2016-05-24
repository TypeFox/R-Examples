# Fit multivariate response distribution 
# Author: Hua Zhou and Yiwen Zhang
#============================================================##

##============================================================## 
## Set class
##============================================================##
setClass("MGLMfit", representation(estimate = "vector", SE = "vector", vcov = "matrix", 
    logL = "numeric", BIC = "numeric", AIC = "numeric", LRTpvalue = "numeric", iter = "numeric", 
    distribution = "character", gradient = "numeric", fitted = "vector"))

##============================================================## 
## Distribution fitting function 
##============================================================##
MGLMfit <- function(data, dist, init, weight, epsilon = 1e-08, maxiters = 150, display = FALSE) {
    
    N <- nrow(data)
    d <- ncol(data)
    
    ## ----------------------------------------## 
    ## Check weights
    ## ----------------------------------------##
    if (!missing(weight) && length(weight) != N) 
        stop("Length of weights doesn't match with the sample size")
    if (!missing(weight) && any(weight < 0)) 
        stop("Negative weights are not allowed")
    if (!is.element(dist, c("DM", "GDM", "NegMN"))) 
        stop(paste("Dist '", dist, "' not valid\n", sep = ""))
    
    ##----------------------------------------## 
    ## Give warnings about zero rows
    ##----------------------------------------##
    if (dist != "NegMN") {
        if (any(rowSums(data) == 0)) {
            rmv <- sum(rowSums(data) == 0)
            warning(paste(rmv, " rows are removed because the row sums are 0", sep = ""))
        }
        if (any(colSums(data) == 0)) {
            rmv <- sum(colSums(data) == 0)
            warning(paste(rmv, " columns are removed because the column sums are 0", 
                sep = ""))
        }
    }
    
    ##----------------------------------------## 
    ## Fit distribution
    ##----------------------------------------##
    if (dist == "DM") {
        if (!missing(init) && length(init) != d) 
            stop("Dimension of the initial values is not compatible with the data")
        est <- DMD.DM.fit(data = data, init = init, weight = weight, epsilon = epsilon, 
            maxiters = maxiters, display = display)
        
    } else if (dist == "GDM") {
        if (!missing(init) && length(init) != 2 * (d - 1)) 
            stop("Dimension of the initial values is not compatible with the data")
        est <- DMD.GDM.fit(data = data, init = init, weight = weight, epsilon = epsilon, 
            maxiters = maxiters, display = display)
        
    } else if (dist == "NegMN") {
        if (!missing(init) && length(init) != (d + 2)) 
            stop("Dimension of the initial values is not compatible with the data")
        est <- DMD.NegMN.fit(data = data, init = init, weight = weight, epsilon = epsilon, 
            maxiters = maxiters, display = display)
    }
    
    ##----------------------------------------## 
    ## Clean up the results
    ##----------------------------------------## 
    if (dist == "DM") {
        if (!is.null(colnames(data))) {
            names(est$estimate) <- paste("alpha", colnames(data), sep = "_")
        } else {
            names(est$estimate) <- paste("alpha", 1:d, sep = "_")
        }
    } else if (dist == "GDM") {
        if (!is.null(colnames(data))) {
            names(est$estimate) <- c(paste("alpha", colnames(data)[1:(d - 1)], sep = "_"), 
                paste("beta", colnames(data)[1:(d - 1)], sep = "_"))
        } else {
            names(est$estimate) <- c(paste("alpha", 1:(d - 1), sep = "_"), paste("beta", 
                1:(d - 1), sep = "_"))
        }
    } else if (dist == "NegMN") {
        if (!is.null(colnames(data))) {
            names(est$estimate) <- c(paste("p", colnames(data), sep = "_"), "phi")
        } else {
            names(est$estimate) <- c(paste("p", 1:d, sep = "_"), "phi")
        }
    }
    
    est$distribution <- ifelse(dist == "DM", "Dirichlet Multinomial", ifelse(dist == 
        "GDM", "Generalized Dirichlet Multinomial", "Negative Multinomial"))
    
    class(est) <- "MGLMfit"
    return(est)
}

##============================================================## 
## Function 1 fit Direchelet Multinomial
##============================================================##

DMD.DM.fit <- function(data, init, weight, epsilon = 1e-08, maxiters = 150, display = FALSE) {
    
    ##----------------------------------------## 
    ## Remove the zero rows and columns
    ##----------------------------------------## 
    if (!missing(weight)) 
        weight <- weight[rowSums(data) != 0]
    
    data <- data[rowSums(data) != 0, colSums(data) != 0]
    N <- nrow(data)
    d <- ncol(data)
    m <- rowSums(data)
    max <- max(m)
    k <- c(0:(max - 1))
    
    ##----------------------------------------## 
    ## Give default value to the missing variables
    ##----------------------------------------## 
    if (missing(weight)) {
        weight <- rep(1, N)
    }
    if (missing(init)) {
        
        rho <- sum((colSums(weight * (data/m)^2))/(colSums(weight * data/m)))
        if (rho == d) {
            init <- rep(1e-06, d)
        } else {
            init <- as.vector((1/N) * (colSums(weight * data/m)) * (d - rho)/(rho - 
                1))
        }
    }
    
    ##----------------------------------------## 
    ## Get prepared for the loop
    ##----------------------------------------## 
    Sj <- function(xj, k, weight) Sj <- colSums(weight * outer(xj, k, ">"))
    s <- apply(data, 2, Sj, k = k, weight = weight)
    r <- Sj(m, k, weight = weight)
    alpha_hat <- init
    niter <- 1
    log_term <- sum(lgamma(m + 1)) - sum(lgamma(data + 1))
    LL <- -sum(r * log(sum(alpha_hat) + k)) + sum(s * log(outer(k, alpha_hat, "+"))) + 
        log_term
    alpha_hat <- init
    DM.LL_iter <- rep(0, maxiters)
    DM.LL_iter[1] <- LL
    
    ##----------------------------------------## 
    ## The main loop
    ##----------------------------------------## 
    while (((niter <= 2) || ((DM.LL2 - DM.LL1)/(abs(DM.LL1) + 1) > epsilon)) & (niter < 
        maxiters)) {
        
        niter <- niter + 1
        DM.LL1 <- -sum(r * log(sum(alpha_hat) + k)) + sum(s * log(outer(k, alpha_hat, 
            "+"))) + log_term
        
        ##----------------------------------------## 
        ## MM update
        numerator <- colSums(s/(outer(rep(1, length(k)), alpha_hat) + outer(k, rep(1, 
            d))))
        denominator <- sum(r/(sum(alpha_hat) + k))
        alpha_MM <- alpha_hat * numerator/denominator
        
        ##----------------------------------------## 
        ## Newton update
        a <- sum(r/(sum(alpha_hat) + k)^2)
        b <- colSums(s/(outer(rep(1, length(k)), alpha_hat) + outer(k, rep(1, d)))^2)
        dl <- (colSums(s/(outer(rep(1, max), alpha_hat) + outer(k, rep(1, d)))) - 
            sum(r/(sum(alpha_hat) + k)))
        alpha_Newton <- alpha_hat + (1/b) * dl + (a * sum((1/b) * dl)/(1 - a * sum(1/b))) * 
            (1/b)
        
        ##----------------------------------------## 
        ## Choose the update
        if (any(is.na(alpha_Newton)) | any(alpha_Newton < 0)) {
            if (is.nan(alpha_MM) || is.na(alpha_MM) || any(alpha_MM == Inf)) {
                stop("DM model is not suitable for this dataset. \n\t\t\t\tPlease use anoter model or privide initial value.")
            } else {
                alpha_hat <- alpha_MM
                if (display) 
                  print(paste("Iteration", niter, "MM update", sep = " "))
            }
        } else {
            LL.MM <- -sum(r * log(sum(alpha_MM) + k)) + sum(s * log(outer(k, alpha_MM, 
                "+"))) + log_term
            LL.Newton <- -sum(r * log(sum(alpha_Newton) + k)) + sum(s * log(outer(k, 
                alpha_Newton, "+"))) + log_term
            if (LL.MM > LL.Newton) {
                alpha_hat <- alpha_MM
                if (display) 
                  print(paste("Iteration", niter, "MM update", sep = " "))
            } else {
                alpha_hat <- alpha_Newton
                if (display) 
                  print(paste("Iteration", niter, "Newton update", sep = " "))
            }
        }
        DM.LL2 <- -sum(r * log(sum(alpha_hat) + k)) + sum(s * log(outer(k, alpha_hat, 
            "+"))) + log_term
        DM.LL_iter[niter] <- DM.LL2
    }
    ##----------------------------------------## 
    ## End of the main loop
    ##----------------------------------------## 
    
    ##----------------------------------------## 
    ## Check the gradients
    ##----------------------------------------## 
    a <- sum(r/(sum(alpha_hat) + k)^2)
    b <- colSums(s/(outer(rep(1, max), alpha_hat) + outer(k, rep(1, d)))^2)
    dl <- (colSums(s/(outer(rep(1, max), alpha_hat) + outer(k, rep(1, d)))) - sum(r/(sum(alpha_hat) + 
        k)))
    if (mean(dl^2) > 1e-04) 
        warning(paste("The algorithm doesn't converge within", niter, sep = ""))
    
    ##----------------------------------------## 
    ## Compute output statistics 1)SE 
	## 2) LRT test against the MN model 3) BIC
    ##----------------------------------------## 
    invI <- diag(1/b) + a/(1 - a * sum(1/b)) * outer(1/b, 1/b, "*")
    SE <- sqrt(diag(invI))
    p_MN <- t(matrix(apply(data/m, 2, mean), d, N))
    logL_MN <- sum(weight * data * log(p_MN)) + sum(lgamma(m + 1)) - sum(lgamma(data + 
        1))
    LRT <- 2 * (DM.LL2 - logL_MN)
    p_value <- pchisq(q = LRT, df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
    BIC <- -2 * DM.LL2 + log(N) * d
    AIC <- -2 * DM.LL2 + 2 * d
    DoF <- d
    fitted <- alpha_hat/sum(alpha_hat)
    
    list(estimate = alpha_hat, SE = SE, vcov = invI, gradient = dl, logL = DM.LL2, 
        iter = niter, BIC = BIC, AIC = AIC, LRT = LRT, LRTpvalue = p_value, fitted = fitted, 
        DoF = DoF)
}

##============================================================## 
## Function 2 fit Generalized DM 
##============================================================##

DMD.GDM.fit <- function(data, init, weight, epsilon = 1e-08, maxiters = 150, display = FALSE) {
    
    ow <- getOption("warn")
    ##----------------------------------------## 
    ## Remove the zero rows and columns
    ##----------------------------------------## 
    if (!missing(weight)) 
        weight <- weight[rowSums(data) != 0]
    
    data <- data[rowSums(data) != 0, colSums(data) != 0]
    N <- nrow(data)
    d <- ncol(data)
    m <- rowSums(data)
    max <- max(m)
    k <- c(0:(max - 1))
    Y <- t(apply(apply(apply(data, 1, rev), 2, cumsum), 2, rev))
    Y2 <- as.matrix(Y[, (2:d)])
    
    ##----------------------------------------## 
    ## Give default value to the missing variables
    ##----------------------------------------## 
    if (missing(weight)) 
        weight <- rep(1, N)
    if (missing(init)) {
        y <- as.matrix(Y[apply(Y, 1, min) != 0, 1:(d - 1)])
        y2 <- as.matrix(Y[apply(Y, 1, min) != 0, 2:d])
        x <- as.matrix(Y[apply(Y, 1, min) != 0, 1:(d - 1)])
        rho <- (colSums((x/y)^2)/colSums(x/y)) + colSums((y2/y)^2)/colSums((y2/y))
        init_alpha <- rep(min((1/N) * (colSums(x/y)) * (2 - rho)/(rho - 1)), (d - 
            1))
        init_alpha[rho == 2] <- 1e-06
        init_beta <- init_alpha
    } else {
        init_alpha <- init[, 1:(d - 1)]
        init_beta <- init[, d:(2 * d - 2)]
    }
    ##----------------------------------------## 
    ## Get prepared for the loop
    ##----------------------------------------## 
    alpha_hat <- rep(0, (d - 1))
    beta_hat <- rep(0, (d - 1))
    gradients <- matrix(0, 2, d - 1)
    LL <- rep(0, (d - 1))
    SE <- matrix(0, 2, (d - 1))
    niter <- rep(0, (d - 1))
    
    ##----------------------------------------## 
    ## Distribution fitting
    ##----------------------------------------## 
    options(warn = -1)
    for (i in 1:(d - 1)) {
        fit <- NULL
        data.fit <- cbind(data[, i], Y2[, i])
        init.fit <- c(init_alpha[i], init_beta[i])
        try(fit <- DMD.DM.fit(data = data.fit, weight = weight, init = init.fit, 
            epsilon = epsilon, maxiters = maxiters, display), silent = TRUE)
        if (is.null(fit)) {
            stop("GDM model is not suitable for this dataset. \n\t\t\t\tPlease use anoter model or privide initial value.")
        }
        alpha_hat[i] <- fit$estimate[1]
        beta_hat[i] <- fit$estimate[2]
        gradients[, i] <- fit$gradient
        LL[i] <- fit$logL
        SE[, i] <- fit$SE
        niter[i] <- fit$iter
    }
    options(warn = ow)
    
    ##----------------------------------------## 
    ##Check the gradients
    ##----------------------------------------## 
    if (mean(gradients^2) > 1e-04) 
        warning(paste("The algorithm doesn't converge within", niter, sep = ""))
    
    ##----------------------------------------## 
    ## Compute output statistics
    ##----------------------------------------## 
    p_MN <- matrix(apply(data/m, 2, mean), N, d, byrow = TRUE)
    logL_MN <- sum(weight * data * log(p_MN)) + sum(lgamma(m + 1)) - sum(lgamma(data + 
        1))
    LRT <- 2 * (sum(LL) - logL_MN)
    p_value <- pchisq(q = LRT, df = (d - 1), ncp = 0, lower.tail = FALSE, log.p = FALSE)
    BIC <- -2 * sum(LL) + log(N) * (d - 1) * 2
    AIC <- -2 * sum(LL) + 2 * (d - 1) * 2
    fitted <- alpha_hat/(alpha_hat + beta_hat)
    for (f in 2:(d - 1)) {
        fitted[f] <- (1 - sum(fitted[1:(f - 1)])) * fitted[f]
    }
    fitted[d] <- (1 - sum(fitted[1:(d - 1)]))
    DoF <- 2 * (d - 1)
    
    list(estimate = c(alpha_hat, beta_hat), SE = c(SE), gradient = gradients, logL = sum(LL), 
        BIC = BIC, AIC = AIC, iter = sum(niter), LRT = LRT, LRTpvalue = p_value, 
        fitted = fitted, DoF = DoF)
}

##============================================================## 
## Function 3 fit Negative Multinomial
##============================================================##
DMD.NegMN.fit <- function(data, init, weight, epsilon = 1e-08, maxiters = 150, display = FALSE) {
    
    N <- nrow(data)
    d <- ncol(data)
    m <- rowSums(data)
    mm <- mean(m)
    max <- max(m)
    k <- c(0:(max - 1))
    cat_count <- colSums(data)
    data <- as.matrix(data)
    ## ----------------------------------------## 
    ## Give default value to the missing variables
    ## ----------------------------------------## 
    if (missing(weight)) 
        weight <- rep(1, N)
    Sj <- function(xj, k, weight) Sj <- colSums(weight * outer(xj, k, ">"))
    Xj <- colSums(data * weight)
    r <- Sj(m, k, weight = weight)
    mbar <- mean(m * weight)
    s2 <- sd(m)^2 * (N - 1)/N
    if (missing(init)) {
        if ((s2 - mbar) > 0) {
            init_beta <- mbar^2/(s2 - mbar)
            init_pi <- rep(0, (d + 1))
            init_pi[d + 1] <- mbar/s2
            init_pi <- c((init_pi[d + 1] * Xj)/(init_beta * N), mbar/s2)
        } else {
            warning("Warning: this is underdispersion, better use other models.\n\t\t\t\tHere we assign beta_init=1")
            init_beta <- 1
            init_pi <- rep(1/(d + 1), d + 1)
        }
    } else {
        init_beta <- init[(d + 2)]
        init_pi <- init[1:(d + 1)]
    }
    
    ##----------------------------------------## 
    ## Get prepared for the main loop
    ##----------------------------------------## 
    beta_hat <- init_beta
    pi_hat <- init_pi[1:d]
    pi2_hat <- init_pi[(d + 1)]
    LL1 <- sum(lgamma(beta_hat + m)) + sum(data %*% log(pi_hat)) + N * beta_hat * 
        log(pi2_hat) - sum(lgamma(data + 1))
    LL_iter <- rep(NA, maxiters)
    LL_iter[1] <- LL1
    niter <- 1
    
    ##----------------------------------------## 
    ## The main loop
    ##----------------------------------------## 
    while ((niter <= 2) || (((LL2 - LL1)/(abs(LL1) + 1) > epsilon) & (niter < maxiters))) {
        
        LL1 <- LL_iter[niter]
        niter <- niter + 1
        
        ##----------------------------------------## 
        ## MM Update
        beta_MM <- -(sum(r * beta_hat/(beta_hat + k))/(N * log(pi2_hat)))
        pi2_MM <- (N * beta_MM)/(sum(Xj) + N * beta_MM)
        pi_MM <- Xj/(sum(Xj) + N * beta_MM)
        LL_MM <- sum(r * log(beta_MM + k)) + sum(Xj * log(pi_MM)) + N * beta_MM * 
            log(pi2_MM) - sum(lgamma(data + 1))
        
        ##----------------------------------------## 
        ## Newton's Update
        score <- c(Xj/pi_hat - sum(weight) * beta_hat/pi2_hat, sum(r/(beta_hat + 
            c(0:(max - 1)))) + sum(weight) * log(pi2_hat))
        
        diaginv <- c(pi_hat^2/Xj, 1/(sum(r/(beta_hat + c(0:(max - 1)))^2) - sum(weight)/beta_hat))
        rank1v <- rep(sqrt(beta_hat * sum(weight))/pi2_hat, (d + 1))
        rank1v[d + 1] <- sqrt(sum(weight)/beta_hat)
        newton_iterate <- c(pi_hat, beta_hat) + diaginv * score - sum(diaginv * rank1v * 
            score)/(1 + sum(diaginv * rank1v^2)) * (diaginv * rank1v)
        pi_Newton <- newton_iterate[1:d]
        pi2_Newton <- 1 - sum(pi_Newton)
        beta_Newton <- newton_iterate[length(newton_iterate)]
        
        ##----------------------------------------## 
        ## Choose update
        check <- all(beta_Newton > 0) && all(pi_Newton > 0) && all(pi_Newton < 1) && 
            all(pi2_Newton < 1) && all(pi2_Newton > 0) && (!is.nan(beta_Newton)) && 
            (!is.nan(pi_Newton))
        if (check) {
            LL_Newton <- sum(r * log(beta_Newton + k)) + sum(Xj * log(pi_Newton)) + 
                N * beta_Newton * log(pi2_Newton) - sum(lgamma(data + 1))
            beta_hat <- beta_MM * (LL_MM >= LL_Newton) + beta_Newton * (1 - (LL_MM >= 
                LL_Newton))
            pi_hat <- pi_MM * (LL_MM >= LL_Newton) + pi_Newton * (1 - (LL_MM >= LL_Newton))
            pi2_hat <- pi2_MM * (LL_MM >= LL_Newton) + pi2_Newton * (1 - (LL_MM >= 
                LL_Newton))
            LL2 <- LL_MM * (LL_MM >= LL_Newton) + LL_Newton * (1 - (LL_MM >= LL_Newton))
        } else {
            beta_hat <- beta_MM
            pi_hat <- pi_MM
            pi2_hat <- pi2_MM
            LL2 <- LL_MM
        }
        
        ##----------------------------------------## 
        ## Display the method used, if asked
        if (display) {
            if (!check) 
                print(paste("iteration", niter, ", MM Algorithm, logL=", LL_MM, sep = "")) else if (LL_MM >= LL_Newton) 
                print(paste("iteration", niter, ", MM Algorithm, logL=", LL_MM, sep = "")) else print(paste("iteration", niter, ", Newton's Method, logL=", LL_Newton, 
                sep = ""))
        }
        LL_iter[niter] <- LL2
        
    }
    ##----------------------------------------## 
    ## End of the main iteration
    ##----------------------------------------## 
    
    ##----------------------------------------## 
    ## Check gradients
    ##----------------------------------------## 
    score <- c(Xj/pi_hat - sum(weight) * beta_hat/pi2_hat, sum(r/(beta_hat + c(0:(max - 
        1)))) + sum(weight) * log(pi2_hat))
    
    if (mean(score^2) > 1e-04) 
        warning(paste("The algorithm doesn't converge within ", niter, " iterations", 
            sep = ""))
    
    ##----------------------------------------## 
    ## Compute output statistics
    ##----------------------------------------## 
    BIC <- -2 * LL2 + log(N) * (d + 1)
    AIC <- -2 * LL2 + 2 * (d + 1)
    dinv <- c(pi_hat^2/Xj, 1/(sum(r/(beta_hat + k)^2) - sum(weight)/beta_hat))
    r1v <- rep(sqrt(beta_hat * sum(weight))/pi2_hat, d + 1)
    r1v[length(r1v)] <- sqrt(sum(weight)/beta_hat)
    SE <- sqrt(dinv - (dinv * r1v)^2/(1 + sum(dinv * r1v^2)))
    fitted <- beta_hat * pi_hat/sum(c(pi_hat, pi2_hat))
    
    list(estimate = c(pi_hat, beta_hat), DoF = (d + 1), gradient = score, SE = SE, 
        logL = LL2, lliter = LL_iter, BIC = BIC, AIC = AIC, itern = niter, LRT = NA, 
        LRTpvalue = NA, fitted = fitted)
} 
