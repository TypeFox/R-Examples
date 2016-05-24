# Initial version of lasso estimation for SDEs
# fixes a small bug in passing arguments to the
# inner optim function

lasso <- function (yuima, lambda0, start, delta = 1, ...)
{
    call <- match.call()
    if (missing(yuima)) 
    yuima.stop("yuima object 'yuima' is missing.")
    pars <- yuima@model@parameter@all
    npars <- length(pars)
    if (missing(lambda0)) {
        lambda0 <- rep(1, npars)
        names(lambda0) <- pars
        lambda0 <- as.list(lambda0)
    }
    if(!is.list(lambda0)){
        lambda0 <- as.numeric(lambda0)
        lambda0 <- as.numeric(matrix(lambda0, npars,1))
        names(lambda0) <- pars
        lambda0 <- as.list(lambda0)
    }
    if (missing(start)) 
    yuima.stop("Starting values for the parameters are missing.")
    fail <- lapply(lambda0, function(x) as.numeric(NA))
    cat("\nLooking for MLE estimates...\n")
    fit <- try(qmle(yuima, start = start, ...), silent = TRUE)
    if (class(fit) == "try-error"){
      tmp <- list(mle = fail, sd.mle = NA, lasso = fail, sd.lasso = NA)
      class(tmp) <- "yuima.lasso"
      return(tmp)
    }
    theta.mle <- coef(fit)
    SIGMA <- try(sqrt(diag(vcov(fit))), silent = TRUE)
    if (class(SIGMA) == "try-error"){
     tmp <- list(mle = theta.mle, sd.mle = NA, lasso = fail, sd.lasso = NA)
     class(tmp) <- "yuima.lasso"
     return(tmp)
    }

    H <- try(solve(vcov(fit)), silent = TRUE)
    if (class(H) == "try-error"){
      tmp <- list(mle = theta.mle, sd.mle = SIGMA, lasso = fail, sd.lasso = NA)
      class(tmp) <- "yuima.lasso"
      return(tmp)
    }

    lambda <- unlist(lambda0[names(theta.mle)])/abs(theta.mle)^delta
    
    #lambda1 <- unlist(lambda0[names(theta.mle)])/abs(theta.mle)
    idx <- which(lambda > 10000)
    lambda[idx] <- 10000
    f2 <- function(theta) as.numeric(t(theta - theta.mle) %*% 
    H %*% (theta - theta.mle) + lambda %*% abs(theta))
    cat("\nPerforming LASSO estimation...\n")
    args <- list(...)
    args$joint <- NULL
    args$par <- theta.mle
    args$fn <- f2
    args$hessian <- TRUE
    args$delta <- NULL
    args$control <- list(maxit = 30000, temp = 2000, REPORT = 500)
    #  fit2 <- try(optim(theta.mle, f2, hessian = TRUE, args, control = list(maxit = 30000, 
    #     temp = 2000, REPORT = 500)), silent = FALSE)
    fit2 <- try( do.call(optim, args = args), silent = TRUE)
    
    
    if (class(fit2) == "try-error"){
      tmp <- list(mle = theta.mle, sd.mle = SIGMA, lasso = fail, sd.lasso = NA)
      class(tmp) <- "yuima.lasso"
      return(tmp)
    }
    theta.lasso <- fit2$par
    SIGMA1 <- try(sqrt(diag(solve(fit2$hessian))), silent = TRUE)
    if (class(SIGMA1) == "try-error"){
      tmp <- list(mle = theta.mle, sd.mle = SIGMA, lasso = theta.lasso,
    sd.lasso = NA)
      class(tmp) <- "yuima.lasso"
      return(tmp)
    }

    tmp <- list(mle = theta.mle, sd.mle = SIGMA, lasso = theta.lasso,
    sd.lasso = SIGMA1, call = call, lambda0 = lambda0)
    class(tmp) <- "yuima.lasso"
    return(tmp)

}




print.yuima.lasso <- function(x,...){
    cat("Adaptive Lasso estimation\n")
    cat("\nCall:\n")
    print(x$cal)
    if(!is.null(x$mle) & !is.null(x$sd.mle)){
        qmle.tab <- rbind(x$mle, x$sd.mle)
        rownames(qmle.tab) <- c("Estimate", "Std. Error")
        cat("\nQMLE estimates\n")
        print(t(qmle.tab))
    }
    
    if(!is.null(x$lasso) & !is.null(x$sd.lasso)){
        lasso.tab <- rbind(x$lasso, x$sd.lasso)
        rownames(lasso.tab) <- c("Estimate", "Std. Error")
        cat("\nLASSO estimates\n")
        print(t(lasso.tab))
    }
}

