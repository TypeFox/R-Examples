rsggm <- function(x, gamma, lambda, nlambda = 10, delta = 0.2, penalty.offdiag = FALSE, method = "glasso", maxit = 100, tol.plogL = 1e-08, msg = TRUE, Omega.init, mu.init) {
    
    if (!is.matrix(x)) 
        stop("\"x\" must be a matrix")
    if (!is.numeric(gamma)) 
        stop("\"gamma\" must be numeric")
    if (any(gamma < 0)) 
        stop("\"gamma\" must be positive")
    if (delta <= 0 || delta >= 1) 
        stop("\"delta\" must be in (0,1)")
    p <- ncol(x)
    N <- nrow(x)
    ngamma <- length(gamma)
    if (!missing(lambda)) {
        if (!is.matrix(lambda)) 
            stop("\"lambda\" must be a matrix")
        nlambda <- nrow(lambda)
        if (ngamma != ncol(lambda)) 
            stop("The number of column of lambda must be equal to length of gamma")
        for (i in 1:ngamma) lambda[, i] <- sort(lambda[, i], decreasing = TRUE)
    }
    
    
    Omega_rlasso_array <- array(0, dim = c(p, p, nlambda, ngamma))
    Omega_rlasso_list <- NULL
    mu_rlasso_list <- NULL
    w_rlasso_list <- NULL
    mu_rlasso_array <- array(0, dim = c(p, nlambda, ngamma))
    w_rlasso_array <- array(0, dim = c(N, nlambda, ngamma))
    w_initial_mat <- set.initial(x, gamma, 1e-05)
    nedges <- matrix(NA, nlambda, ngamma)
    ones <- rep(1, p)
    
    colnamesx <- colnames(x)
    if (is.null(colnamesx)) 
        colnamesx <- paste("V", 1:p, sep = "")
    
    
    if (missing(lambda)) {
        lambda <- matrix(NA, nlambda, ngamma)
        for (i_gamma in 1:ngamma) {
            w_gd <- w_initial_mat[, i_gamma]
            gamma_element <- gamma[i_gamma]
            # set initial values
            w_gd <- w_initial_mat[, i_gamma]
            mu <- c(crossprod(x, w_gd))
            # update Snew
            S_new <- crossprod(x, w_gd * x) - 2 * mu %*% (w_gd %*% x) + crossprod(t(mu))
            S_new2 <- (1 + gamma_element) * S_new
            S_new2adj <- S_new2
            diag(S_new2adj) <- 0
            lambdamax <- max(abs(S_new2adj))
            lambdamin <- lambdamax * delta
            lambda[, i_gamma] <- exp(seq(log(lambdamax), log(lambdamin), length.out = nlambda))
            
        }
    }
    
    
    
    
    # conduct the gamma-lasso
    for (i_gamma in 1:ngamma) {
        w_gd <- w_initial_mat[, i_gamma]
        Omega_rlasso_list_lam <- NULL
        mu_rlasso_list_lam <- NULL
        w_rlasso_list_lam <- NULL
        for (i_lambda in 1:nlambda) {
            lambda_element <- lambda[i_lambda, i_gamma]
            gamma_element <- gamma[i_gamma]
            if (msg) {
                cat(paste(i_lambda + (i_gamma - 1) * nlambda, "/", ngamma * nlambda, ": ", sep = ""))
                cat(paste("lambda: ", sprintf("%5e", lambda_element), ", ", sep = ""))
                cat(paste("gamma: ", sprintf("%5e", gamma_element), "\n", sep = ""))
                cat(paste("     computing", sep = ""))
            }
            
            # set initial values
            w_gd <- w_initial_mat[, i_gamma]
            
            # mu
            if (missing(mu.init)) 
                mu <- c(crossprod(x, w_gd))
            if (!missing(mu.init)) 
                mu <- mu.init
            # update Snew
            S_new <- crossprod(x, w_gd * x) - 2 * mu %*% (w_gd %*% x) + crossprod(t(mu))
            S_new2 <- (1 + gamma_element) * S_new
            # Sigma
            if (missing(Omega.init)) {
                Omega <- diag(1/(diag(S_new) * (1 + gamma_element)))
            } else {
                Omega <- Omega.init
            }
            
            
            
            plogL_weight0_current <- -Inf
            plogL_weight0 <- NULL
            
            # iteration via the weighted glasso
            for (i_weight in 1:maxit) {
                
                # update the weights and Snew2
                if (i_weight != 1) {
                  xx <- x - matrix(mu, N, p, byrow = T)
                  w_gd <- (xx * (xx %*% Omega)) %*% ones
                  w_gd <- exp(-gamma_element/2 * w_gd)
                  w_gd <- w_gd/sum(w_gd)
                  w_gd <- c(w_gd)
                  mu <- c(crossprod(x, w_gd))
                  S_new <- crossprod(x, w_gd * x) - 2 * mu %*% (w_gd %*% x) + crossprod(t(mu))
                  S_new2 <- (1 + gamma_element) * S_new
                }
                if (sum((1 - w_gd) < 1e-05) == 1) {
                  stop("gamma is too large!")
                }
                
                lambda_mat <- matrix(lambda_element, p, p)
                if (!penalty.offdiag) 
                  diag(lambda_mat) <- 0
                
                
                if (method == "QUIC") {
                  result_rlasso_QUIC <- QUIC(S_new2, (1 + gamma_element) * lambda_mat, msg = 0, tol = 1e-05)
                  Omega <- result_rlasso_QUIC$X
                  plogL_weight0[i_weight] <- result_rlasso_QUIC$regloglik * 2/p
                }
                
                if (method == "glasso") {
                  result_rlasso_glasso <- glasso(S_new2, (1 + gamma_element) * lambda_mat)
                  Omega <- result_rlasso_glasso$wi
                  plogL_weight0[i_weight] <- result_rlasso_glasso$loglik * 2/p
                }
                
                tol.plogL_weight0 <- plogL_weight0[i_weight] - plogL_weight0_current
                plogL_weight0_current <- plogL_weight0[i_weight]
                if (abs(tol.plogL_weight0) < tol.plogL) 
                  break
                
                if (msg) 
                  cat(".")
                
            }
            
            colnames(Omega) <- colnamesx
            rownames(Omega) <- colnamesx
            names(mu) <- colnamesx
            
            Omega_rlasso_list_lam <- c(Omega_rlasso_list_lam, list(as(Omega, "sparseMatrix")))
            names(Omega_rlasso_list_lam)[i_lambda] <- paste("lambda", i_lambda, sep = "")
            
            mu_rlasso_list_lam <- c(mu_rlasso_list_lam, list(mu))
            names(mu_rlasso_list_lam)[i_lambda] <- paste("lambda", i_lambda, sep = "")
            
            w_rlasso_list_lam <- c(w_rlasso_list_lam, list(w_gd))
            names(w_rlasso_list_lam)[i_lambda] <- paste("lambda", i_lambda, sep = "")
            
            tol_hantei <- 1e-05
            n_nonzero <- sum(abs(Omega) > tol_hantei)
            nedges[i_lambda, i_gamma] <- trunc((n_nonzero - p)/2)
            
            if (msg) {
                cat(paste("\n     number of edges: ", nedges[i_lambda, i_gamma], "\n\n", sep = ""))
            }
        }
        Omega_rlasso_list <- c(Omega_rlasso_list, list(Omega_rlasso_list_lam))
        mu_rlasso_list <- c(mu_rlasso_list, list(mu_rlasso_list_lam))
        w_rlasso_list <- c(w_rlasso_list, list(w_rlasso_list_lam))
        names(Omega_rlasso_list)[i_gamma] <- paste("gamma", i_gamma, sep = "")
        names(mu_rlasso_list)[i_gamma] <- paste("gamma", i_gamma, sep = "")
        names(w_rlasso_list)[i_gamma] <- paste("gamma", i_gamma, sep = "")
        
        
    }
    
    ans <- list(Omega = Omega_rlasso_list, mu = mu_rlasso_list, weight = w_rlasso_list, nedges = nedges, lambda = lambda, gamma = gamma)
    class(ans) <- "rsggm"
    ans$call <- match.call()
    ans
    
}


set.initial <- function(x, gamma, tol = 1e-05) {
    ngamma <- length(gamma)
    ans <- matrix(NA, nrow(x), ngamma)
    ones <- rep(1, ncol(x))
    N <- nrow(x)
    p <- ncol(x)
    for (i_gamma in 1:ngamma) {
        gamma_element <- gamma[i_gamma]
        w0 <- rep(1/N, N)
        for (n_iter in 1:100) {
            
            mu0 <- c(crossprod(x, w0))
            S_new <- crossprod(x, w0 * x) - 2 * mu0 %*% (w0 %*% x) + crossprod(t(mu0))
            
            Omega00_diag <- 1/((1 + gamma_element) * diag(S_new))
            w0old <- w0
            xx <- x - matrix(mu0, N, p, byrow = T)
            w0 <- (xx^2 * matrix(Omega00_diag, N, p, byrow = T)) %*% ones
            w0 <- exp(-gamma_element/2 * w0)
            w0 <- w0/sum(w0)
            w0 <- c(w0)
            tol0 <- sum((w0old - w0)^2)
            if (sum((1 - w0) < 1e-05) == 1) {
                stop("gamma is too large!")
            }
            if (tol0 < tol) 
                break
        }
        ans[, i_gamma] <- w0
    }
    ans
}


print.rsggm <- function(x, digits = max(3, getOption("digits") - 3), num.result = 20, ...) {
    lambda <- x$lambda
    gamma <- x$gamma
    gamma <- c(gamma)
    nedges <- x$nedges
    colnames(lambda) <- gamma
    colnames(nedges) <- gamma
    cat("\nCall:", paste(deparse(x$call)), "\n")
    cat("\ngamma:\n")
    print(gamma, digits = digits)
    cat("\nlambda:\n")
    print(lambda, digits = digits)
    cat("\nnumber of edges:\n")
    print(nedges, digits = digits)
    cat("\n")
    invisible(x)
}



rsggm.generator <- function(N = 200, p = 20) {
    
    if (N < 50) 
        stop("\"N\" must be greater than or equal to 50")
    if (p < 3) 
        stop("\"p\" must be greater than or equal to 3")
    # generate random network
    Net <- matrix(0, p, p)
    Net0 <- rbinom(p * (p - 1)/2, size = 1, prob = 0.01)
    Net[upper.tri(Net)] <- Net[lower.tri(Net)] <- Net0
    
    # generate covariance matrix
    E0 <- matrix(runif(p^2, min = -0.5, max = 0.5), p, p)
    E <- matrix(NA, p, p)
    E[E0 > 0] <- E0[E0 > 0] + 0.25
    E[E0 <= 0] <- E0[E0 <= 0] - 0.25
    E[Net == 0] <- 0
    Ebar <- (E + t(E))/2
    L <- Ebar + (0.1 - eigen(Ebar)$values[p]) * diag(p)
    Linvhalf <- diag(diag(solve(L))^(1/2))
    Omega_true <- Linvhalf %*% L %*% Linvhalf
    
    # generate dataset
    x <- mvrnorm(N - 20, rep(0, p), solve(Omega_true))
    # generate outliers
    x <- rbind(x, mvrnorm(20, rep(5, p), diag(p)))
    x
}


out <- function(object, gamma, lambda) {
    if (sum(object$gamma == gamma) == 0) {
        tmp <- NULL
        for (i in 1:length(object$gamma)) {
            if (i != length(object$gamma)) 
                tmp <- paste(tmp, object$gamma[i], ", ", sep = "") else tmp <- paste(tmp, object$gamma[i], sep = "")
        }
        stop(paste("\"gamma\" must be one of ", tmp, sep = ""))
    }
    i_gamma <- which(object$gamma == gamma)
    lambda0 <- object$lambda[, i_gamma]
    if (min(lambda0) >= lambda || max(lambda0) <= lambda) 
        stop(paste("\"lambda\" must be in (", min(lambda0), ", ", max(lambda0), ")", sep = ""))
    i_lambda <- which.min((lambda0 - lambda)[lambda0 - lambda > 0])
    t <- (lambda0[i_lambda] - lambda)/(lambda0[i_lambda] - lambda0[i_lambda + 1])
    if (t <= 0 || t >= 1) 
        stop("Error in calculating the value of t")
    Omega <- (1 - t) * object$Omega[[i_gamma]][[i_lambda]] + t * object$Omega[[i_gamma]][[i_lambda + 1]]
    mu <- (1 - t) * object$mu[[i_gamma]][[i_lambda]] + t * object$mu[[i_gamma]][[i_lambda + 1]]
    weight <- (1 - t) * object$weight[[i_gamma]][[i_lambda]] + t * object$weight[[i_gamma]][[i_lambda + 1]]
    
    tol_hantei <- 1e-05
    n_nonzero <- sum(abs(Omega) > tol_hantei)
    nedges <- trunc((n_nonzero - length(object$mu[[1]][[1]]))/2)
    ans <- list(Omega = Omega, mu = mu, weight = weight, nedges = nedges, lambda = lambda, gamma = gamma)
    ans
}