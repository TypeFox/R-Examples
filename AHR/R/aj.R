## check if multi-state model is a competing risks model (no transient states)
isCmprsk <- function(tra) all(tra[1,-1]) && all(!tra[-1,])

## estimate asymptotic covariance function of log(1-CIF)
ajLogCovCIF <- function(times, fit, target=2, n) {

    n.times <- length(times)
    tt <- fit$time
    delta.na <- fit$delta.na

    Y <- fit$n.risk ## = number of patients in initial state for each event time
    dA1 <- delta.na[1, target,] ## Nelson-Aalen increments for transition 0 -> target
    dA <- colSums(matrix(delta.na[1, -c(1, target),], ncol=length(tt)))  ## sum of Nelson-Aalen increments for all other transitions

    F1u <- fit$est[1, target,]
    F1s <- evalFun(times, tt, F1u)
    Fu <- colSums(matrix(fit$est[1, -c(1, target),], ncol=length(tt)))
    
    a1 <- ((1 - Fu)^2 * dA1 + F1u^2 * dA) / Y
    a2 <- ((1 - Fu) * dA1 + F1u * dA) / Y
    a3 <- (dA1 + dA) / Y

    rlm <- do.call(rbind, lapply(1:n.times, function(s1) {
        J <- tt <= times[s1]
        if(any(J)) {
            b1 <- sum(a1[J])
            b2 <- sum(a2[J])
            b3 <- sum(a3[J])
            K <- s1:n.times
            c(rep.int(0, s1-1), (b1 - (F1s[s1] + F1s[K]) * b2 + F1s[s1] * F1s[K] * b3) / ((1 - F1s[s1])*(1-F1s[K])))
        } else rep.int(0, n.times)
    }))

    res <- rlm + t(rlm)
    diag(res) <- diag(rlm)
    res
}

#' Aalen-Johansen estimator (empirical transition matrix)
#'
#' Wrapper for 'etm' function from the 'etm' package to be used with \code{\link{ahrAJ}}
#' 
#' @title aj
#' @param times a vector of evaluation times
#' @param data data frame (see \code{\link{etm}} function documentation)
#' @param param list of parameters (target, states, transitions, censoring, s, t, cov) (see \link{etm} documentation)
#' @return a list containing
#'  \item{times}{the argument \code{times} passed to the function}
#'  \item{S}{vector of 1 minus transition probabilities at \code{times} (one for each element of \code{times})}
#'  \item{V}{vector of variances at 'times' (only if \code{param$cov} is TRUE)}
#'  \item{logCOV}{matrix containing estimated values of the log-covariance function evaluated for all pairs of elements of the vector \code{times} (only if \code{param$cov} is TRUE and model is a competing risks model)}
#' @export
#' @seealso \code{\link[etm]{etm}}
#' @details For a description of the parameters in the list \code{param} see the documentation of the \code{\link{etm}} function in package \code{etm}.
#' @examples
#' ## competing risks
#' T <- rexp(100)
#' C <- rexp(100)
#' r <- rbinom(100, 2, 0.5)
#' r[(r == 0) | (T > C)] <- "cens"
#' data <- data.frame(id=1:100, time=pmin(T,C), from=rep(0, 100), to=r)
#' data <- data[order(data$time),]
#' tra <- matrix(FALSE, nrow=3, ncol=3)
#' tra[1, 2:3] <- TRUE
#' # estimate cumulative incidence function for event type 1
#' fit <- aj(sort(data$time), data, list(target="0 1", states=c("0", "1", "2"), transitions=tra,
#'    censoring="cens", s=0, t="last", covariance=TRUE))
aj <- function(times, data, param) {

    states <- param$states
    tra <- param$transitions
    
    ## empirical transition matrix / Aalen-Johansen estimator
    fit <- etm::etm(data, states, tra, param$censoring, s=param$s, t=param$t, covariance=param$cov)
    
    ## convert target string to state numbers
    tmp <- strsplit(param$target, split=NULL)[[1]]
    from.to <- c(which(tmp[1] == states), which(tmp[3] == states))

    F <- pmax(0, 1 - fit$est[from.to[1],from.to[2],])    
    S <- evalFun(times, fit$time, F)
    
    ## cum.n.risk <- rowSums(as.matrix(fit$n.risk[, tra[, from.to[2]][-from.to[2]]]))
    ## cumulative number-at-risk (sum over all transient states, assuming only one absorbing state)
    ## cum.n.risk <- rowSums(as.matrix(fit$n.risk))   
    ## number-at-risk of transitioning to target state at each element of 'times'
    ## f <- approxfun(fit$time, cum.n.risk, method="constant", rule=2, f=0)
    ## nar <- f(times)
    
    obj <- list(times=times, S=S, n.atrisk=NULL)
    
    if(param$cov) {
        obj$V <- evalFun(times, fit$time, fit$cov[param$target, param$target, ])
        if(isCmprsk(tra)) obj$logCOV <- ajLogCovCIF(times, fit, from.to[2], length(unique(data$id)))
    }

    obj    
}
