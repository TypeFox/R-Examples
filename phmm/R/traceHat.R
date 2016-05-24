#' Conditional Akaike Information Criterion for PHMM
#' 
#' Function calculating the conditional Akaike information criterion (Vaida and
#' Blanchard 2005) for PHMM fitted model objects, according to the formula
#' \eqn{-2*log-likelihood + k*rho}{-2*log-likelihood + k*rho}, where
#' \eqn{rho}{rho} represents the "effective degrees of freedom" in the sense of
#' Hodges and Sargent (2001). The function uses the log-likelihood conditional
#' on the estimated random effects; and trace of the "hat matrix", using the
#' generalized linear mixed model formulation of PHMM, to estimate
#' \eqn{rho}{rho}. The default k = 2, conforms with the usual AIC.
#' 
#' 
#' @aliases cAIC cAIC.phmm cAIC.coxph
#' @param object A fitted PHMM model object of class \code{phmm}.
#' @param method Passed to \code{\link{traceHat}}. Options include "direct",
#' "pseudoPois", or "HaLee". The methods "direct" and "HaLee" are algebraically
#' equivalent.
#' @param ... Optionally more fitted model objects.
#' @param k numeric, the penalty per parameter to be used; the default k = 2
#' conforms with the classical AIC.
#' @return Returns a numeric value of the cAIC corresonding to the PHMM fit.
#' @seealso \code{\link{phmm}}, \code{\link[stats]{AIC}}
#' @references Vaida, F, and Blanchard, S. 2005. Conditional Akaike information
#' for mixed-effects models. Biometrika, 92(2), 351-.
#' 
#' Donohue, MC, Overholser, R, Xu, R, and Vaida, F (January 01, 2011).
#' Conditional Akaike information under generalized linear and proportional
#' hazards mixed models. \emph{Biometrika}, 98, 3, 685-700.
#' 
#' Breslow, NE, Clayton, DG. (1993). Approximate Inference in Generalized
#' Linear Mixed Models. Journal of the American Statistical Association, Vol.
#' 88, No. 421, pp. 9-25.
#' 
#' Whitehead, J. (1980). Fitting Cox\'s Regression Model to Survival Data using
#' GLIM. Journal of the Royal Statistical Society. Series C, Applied
#' statistics, 29(3), 268-.
#' 
#' Hodges, JS, and Sargent, DJ. 2001. Counting degrees of freedom in
#' hierarchical and other richly-parameterised models. Biometrika, 88(2), 367-.
#' @keywords survival
#' @examples
#' 
#' \dontrun{
#' n <- 50      # total sample size
#' nclust <- 5  # number of clusters
#' clusters <- rep(1:nclust,each=n/nclust)
#' beta0 <- c(1,2)
#' set.seed(13)
#' #generate phmm data set
#' Z <- cbind(Z1=sample(0:1,n,replace=TRUE),
#'            Z2=sample(0:1,n,replace=TRUE),
#'            Z3=sample(0:1,n,replace=TRUE))
#' b <- cbind(rep(rnorm(nclust),each=n/nclust),rep(rnorm(nclust),each=n/nclust))
#' Wb <- matrix(0,n,2)
#' for( j in 1:2) Wb[,j] <- Z[,j]*b[,j]
#' Wb <- apply(Wb,1,sum)
#' T <- -log(runif(n,0,1))*exp(-Z[,c('Z1','Z2')]%*%beta0-Wb)
#' C <- runif(n,0,1)
#' time <- ifelse(T<C,T,C)
#' event <- ifelse(T<=C,1,0)
#' mean(event)
#' phmmd <- data.frame(Z)
#' phmmd$cluster <- clusters
#' phmmd$time <- time
#' phmmd$event <- event
#' 
#' fit.phmm <- phmm(Surv(time, event) ~ Z1 + Z2 + (-1 + Z1 + Z2 | cluster), 
#'    phmmd, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
#'    NINIT = 10, MAXSTEP = 100, CONVERG=90)
#' 
#' # Same data can be fit with lmer,
#' # though the correlation structures are different.
#' poisphmmd <- pseudoPoisPHMM(fit.phmm)
#' 
#' library(lme4)
#' fit.lmer <- lmer(m~-1+as.factor(time)+z1+z2+
#'   (-1+w1+w2|cluster)+offset(log(N)), 
#'   as.data.frame(as(poisphmmd, "matrix")), family=poisson)
#' 
#' fixef(fit.lmer)[c("z1","z2")]
#' fit.phmm$coef
#' 
#' VarCorr(fit.lmer)$cluster
#' fit.phmm$Sigma
#' 
#' logLik(fit.lmer)
#' fit.phmm$loglik
#' 
#' traceHat(fit.phmm)
#' 
#' summary(fit.lmer)@AICtab
#' AIC(fit.phmm)
#' cAIC(fit.phmm)
#' }
#'
cAIC <- function(object, ..., k = 2) UseMethod("cAIC")

cAIC.phmm <- function(object, method = "direct", ..., k = 2){
  as.numeric(-2*object$loglik["Conditional"])+k*traceHat(object, method = method)
}

cAIC.coxph <- function(object, ..., k = 2){
	-2*object$loglik[2] + k*length(object$coef)
}

AIC.coxph <- cAIC.coxph

#' Trace of the "hat" matrix from PHMM-MCEM fit
#' 
#' Compute trace of the ``hat'' matrix from PHMM-MCEM fit using a direct
#' approximation method (Donohue, et al, submitted), an approximation via
#' hierarchical likelihoods (Ha et al, 2007), or an approximation via a
#' generalized linear mixed-effects model (GLMM) (Donohue, et al, submitted).
#' 
#' 
#' @param x an object of class \code{phmm},
#' @param method acceptable values are \code{"direct"}, \code{"pseudoPois"}, or
#' \code{"HaLee"},
#' @return The trace of the "hat" matrix which can be used as a measure of
#' complexity of the model.
#' @seealso \code{\link{phmm}}, \code{\link{AIC.phmm}}
#' @references Breslow, NE, Clayton, DG. (1993). Approximate Inference in
#' Generalized Linear Mixed Models. Journal of the American Statistical
#' Association, Vol. 88, No. 421, pp. 9-25.
#' 
#' Donohue, M, Xu, R, Vaida, F, Haut R. Model Selection for Clustered Data:
#' Conditional Akaike Information under GLMM and PHMM. Submitted.
#' 
#' Ha, ID, Lee, Y, MacKenzie, G. (2007). Model Selection for multi-component
#' frailty models. Statistics in Medicine, Vol. 26, pp. 4790-4807.
#' 
#' Whitehead, J. (1980). Fitting Cox\'s Regression Model to Survival Data using
#' GLIM. Journal of the Royal Statistical Society. Series C, Applied
#' statistics, 29(3), 268-.
#' @keywords survival
#' @examples
#' \dontrun{
#' n <- 50      # total sample size
#' nclust <- 5  # number of clusters
#' clusters <- rep(1:nclust,each=n/nclust)
#' beta0 <- c(1,2)
#' set.seed(13)
#' #generate phmm data set
#' Z <- cbind(Z1=sample(0:1,n,replace=TRUE),
#'            Z2=sample(0:1,n,replace=TRUE),
#'            Z3=sample(0:1,n,replace=TRUE))
#' b <- cbind(rep(rnorm(nclust),each=n/nclust),rep(rnorm(nclust),each=n/nclust))
#' Wb <- matrix(0,n,2)
#' for( j in 1:2) Wb[,j] <- Z[,j]*b[,j]
#' Wb <- apply(Wb,1,sum)
#' T <- -log(runif(n,0,1))*exp(-Z[,c('Z1','Z2')]%*%beta0-Wb)
#' C <- runif(n,0,1)
#' time <- ifelse(T<C,T,C)
#' event <- ifelse(T<=C,1,0)
#' mean(event)
#' phmmd <- data.frame(Z)
#' phmmd$cluster <- clusters
#' phmmd$time <- time
#' phmmd$event <- event
#' 
#' fit.phmm <- phmm(Surv(time, event) ~ Z1 + Z2 + (-1 + Z1 + Z2 | cluster), 
#'    phmmd, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
#'    NINIT = 10, MAXSTEP = 100, CONVERG=90)
#' 
#' # Same data can be fit with lmer,
#' # though the correlation structures are different.
#' poisphmmd <- pseudoPoisPHMM(fit.phmm)
#' 
#' library(lme4)
#' fit.lmer <- lmer(m~-1+as.factor(time)+z1+z2+
#'   (-1+w1+w2|cluster)+offset(log(N)), 
#'   as.data.frame(as(poisphmmd, "matrix")), family=poisson)
#' 
#' fixef(fit.lmer)[c("z1","z2")]
#' fit.phmm$coef
#' 
#' VarCorr(fit.lmer)$cluster
#' fit.phmm$Sigma
#' 
#' logLik(fit.lmer)
#' fit.phmm$loglik
#' 
#' traceHat(fit.phmm)
#' }
traceHat <- function(x, method = "direct"){ 
  if(! method %in% c("direct", "pseudoPois", "HaLee", "approx")) 
    stop("Undefined traceHat method.")

    time <- x$Y[, 1]; delta <- x$Y[, 2]; z <- x$Z
    w <- x$W; b <- x$bhat.long; Sigma <- Matrix(x$Sigma)
    eventtimes <- unique(sort(time[delta == 1]))
    neventtimes <- length(eventtimes)
    cluster <- as.numeric(as.character(x$cluster))
    nclust <- length(unique(cluster))

  if(method%in%c("direct", "HaLee", "approx")){
    # use Ha, Lee, MacKenzie 2007 framework
    # to compute traceHat directly
    # z -- covs for fixed effects
    # w -- covs for random effects
    # time -- min(T, C)
    # Sigma -- varcov Matrix for random effects
    # fitted -- linear part

    fitted <- x$linear.predictors  
    xx <- cBind(ID = 1:length(time), time = time, delta = delta)
    Lambda = cumsum(x$lambda)[!duplicated(sort(time), fromLast = TRUE)]
    lambda = Lambda-c(0, Lambda[1:(length(Lambda)-1)])
    Lambda <- cBind(Lambda = Lambda, 
            lambda = lambda, 
                    time = unique(sort(time)))
    xx <- merge(xx, Lambda, by = "time", all.x = TRUE)
    Lambda <- xx[order(xx$ID), "Lambda"]
    lambda <- xx[order(xx$ID), "lambda"]
  
    xxx <- xx[xx$delta == 1, c("time", "lambda")]
    ulambda <- xxx[!duplicated(xxx$time), "lambda"]
  
    time <- time[order(cluster)]
    delta <- delta[order(cluster)]
    if(ncol(z) == 1){ z <- Matrix(z[order(cluster)])
      }else{ z <- Matrix(z[order(cluster), ]) }
    if(ncol(w) == 1){ 
      w <- Matrix(w[order(cluster)])
      b <- Matrix(b[order(cluster)])
      bhat <- b[!duplicated(sort(cluster))]
      }else{ 
        w <- Matrix(w[order(cluster), ]) 
        b <- Matrix(b[order(cluster), ])
        bhat <- b[!duplicated(sort(cluster)), ]
        }
    
    cluster <- sort(cluster)
    nclust <- length(unique(cluster))
    cluster <- rep(1:nclust, table(cluster))
  
    X <- z
    Z <- bdiag(lapply(unique(cluster), function(x){
         block <- w[cluster == x, ]
         if(ncol(w) > 1) return(rBind(block)) else return(cBind(block))
       }))
    D <- bdiag(rep(list(Sigma), nclust))
  
    W3 <- diag(as.vector(exp(fitted)))
    B <- diag(Lambda)
    W1 <- W3 %*% B
    
    M <- outer(time, eventtimes, FUN = function(x, y) ifelse(x >= y, 1, 0))
  
    dk <- rep(0, neventtimes)
    for(k in 1:neventtimes)
      dk[k] <- sum(delta[time == eventtimes[k]])
    C <- diag(dk/(ulambda^2))
    W2 <- (W3 %*% M) %*% solve(C) %*% t(W3 %*% M)
    W <- W1-W2
    
    if(nrow(Sigma) == 1 & Sigma[1, 1] == 0){ 
      return(sum(diag(solve(crossprod(X, W) %*% X %*% crossprod(X, W) %*% X))))
    }else if(method == "direct"){
      U <- cBind(as.matrix(X), as.matrix(Z))
      A <- Matrix(bdiag(matrix(0, x$nfixed, x$nfixed), solve(D)))
      UW <- crossprod(U, W)
      return(sum(diag(UW %*% U %*% 
                solve(UW %*% U + A))))
    }else if(method %in% c("HaLee", "approx")){
      U <- solve(D)
      XW <- crossprod(X, W); ZW <- crossprod(Z, W)
      J <- if(method == "HaLee") 0 else U %*% as.vector(t(bhat)) %*% t(as.vector(t(bhat))) %*% U
      return(sum(diag(solve(
      rBind(cBind(XW %*% X, XW %*% Z), 
            cBind(ZW %*% X, ZW %*% Z + U))) %*% 
      rBind(cBind(XW %*% X, XW %*% Z), 
            cBind(ZW %*% X, ZW %*% Z + J)))))
    }
  }else if(method == "pseudoPois"){
    xx <- pseudoPoisPHMM(x)
    
    xx <- cBind( ID = 1:nrow(xx), xx)
    xx <- xx[order(xx[, "cluster"], xx[, "ID"]), ]
    z <- xx[, c(paste("t", 1:neventtimes, sep = ''), 
                paste("z", 1:x$nfixed, sep = ''))]
    w <- Matrix(as.matrix(xx[, paste("w", 1:x$nrandom, sep = '')]))
    cluster <- xx[, "cluster"]
    fitted <- xx[, "linear.predictors"] + 
      xx[, paste("t", 1:neventtimes, sep = '')] %*% log(x$lambda[x$lambda != 0]) + 
      log(xx[, "N"])
    nclust <- length(unique(cluster))
    cluster <- rep(1:nclust, table(cluster))
    Sigma = Matrix(x$Sigma)

    X <- as.matrix(z)
    Z <- bdiag(lapply(unique(cluster), function(x) w[cluster == x, ]))
    D <- bdiag(rep(list(Sigma), nclust))
    W <- Diagonal(x = as.numeric(exp(fitted)))
    U <- solve(D)
    ZWZ <- crossprod(Z, W) %*% Z
    ZWX <- crossprod(Z, W) %*% X
    XWX <- crossprod(X, W) %*% X
    XWZ <- crossprod(X, W) %*% Z

    return(nclust * x$nrandom + x$nfixed - sum(diag(solve(ZWZ + U - ZWX %*% solve(XWX) %*% XWZ) %*% U)))
    }else return(NULL) # end pseudoPois
}
