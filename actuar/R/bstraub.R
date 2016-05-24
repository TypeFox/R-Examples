### ===== actuar: An R Package for Actuarial Science =====
###
### Buhlmann-Straub credibility model calculations.
###
### Computation of the between variance estimators has been moved to
### external functions bvar.unbiased() and bvar.iterative() to share
### with hache().
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Sebastien Auclair, Louis-Philippe Pouliot

bstraub <- function(ratios, weights, method = c("unbiased", "iterative"),
                    tol = sqrt(.Machine$double.eps), maxit = 100, echo = FALSE)
{
    ## If weights are not specified, use equal weights as in
    ## Buhlmann's model.
    if (missing(weights))
    {
        if (any(is.na(ratios)))
            stop("missing ratios not allowed when weights are not supplied")
        weights <- array(1, dim(ratios))
    }

    ## Check other bad arguments.
    if (ncol(ratios) < 2)
        stop("there must be at least one node with more than one period of experience")
    if (nrow(ratios) < 2)
        stop("there must be more than one node")
    if (!identical(which(is.na(ratios)), which(is.na(weights))))
        stop("missing values are not in the same positions in 'weights' and in 'ratios'")
    if (all(!weights, na.rm = TRUE))
        stop("no available data to fit model")

    ## Individual weighted averages. It could happen that a contract
    ## has no observations, for example when applying the model on
    ## claim amounts. In such a situation, we will put the total
    ## weight of the contract and the weighted average both equal to
    ## zero. That way, the premium will be equal to the credibility
    ## weighted average, as it should, but the contract will have no
    ## contribution in the calculations.
    weights.s <- rowSums(weights, na.rm = TRUE)
    ratios.w <- ifelse(weights.s > 0, rowSums(weights * ratios, na.rm = TRUE) / weights.s, 0)

    ## Size of the portfolio.
    ncontracts <- sum(weights.s > 0)
    ntotal <- sum(!is.na(weights))

    ## Collective weighted average.
    weights.ss <- sum(weights.s)

    ## Estimation of s^2
    s2 <-  sum(weights * (ratios - ratios.w)^2, na.rm = TRUE) / (ntotal - ncontracts)

    ## First estimation of a. Always compute the unbiased estimator.
    a <- bvar.unbiased(ratios.w, weights.s, s2, ncontracts)

    ## Iterative estimation of a. Compute only if
    ## 1. asked to in argument;
    ## 2. weights are not all equal (Buhlmann model).
    ## 3. the unbiased estimator is > 0;
    method <- match.arg(method)

    if (method == "iterative" &&
        diff(range(weights, na.rm = TRUE)) > .Machine$double.eps^0.5)
    {
        a <-
            if (a > 0)
                bvar.iterative(ratios.w, weights.s, s2, ncontracts, start = a,
                               tol = tol, maxit = maxit, echo = echo)
            else
                0
    }

    ## Final credibility factors and estimator of the collective mean.
    if (a > 0)
    {
        cred <- 1 / (1 + s2/(weights.s * a))
        ratios.zw <- drop(crossprod(cred, ratios.w)) / sum(cred)
    }
    else
    {
        cred <- numeric(length(weights.s))
        ratios.zw <- drop(crossprod(weights.s, ratios.w)) / sum(weights.s)
    }

    structure(list(means = list(ratios.zw, ratios.w),
                   weights = list(if (a > 0) sum(cred) else weights.ss, weights.s),
                   unbiased = if (method == "unbiased") c(a, s2),
                   iterative = if (method == "iterative") c(a, s2),
                   cred = cred,
                   nodes = list(nrow(weights))),
              class = "bstraub",
              model = "Buhlmann-Straub")
}

predict.bstraub <- function(object, levels = NULL, newdata, ...)
    structure(object$means[[1]] + object$cred * (object$means[[2]] - object$means[[1]]), ...)

bvar.unbiased <- function(x, w, within, n)
{
    w.s <- sum(w)
    x.w <- drop(crossprod(w, x)) / w.s

    w.s * (drop(crossprod(w, (x - x.w)^2)) - (n - 1) * within) / (w.s^2 - sum(w^2))
}

### codetools does not like the way 'a1' is defined in function
### 'bvar.iterative' below. Avoid false positive in R CMD check.
if(getRversion() >= "2.15.1")  utils::globalVariables(c("a1"))

bvar.iterative <- function(x, w, within, n, start,
                           tol = sqrt(.Machine$double.eps), maxit = 100,
                           echo = FALSE)
{
    if (echo)
    {
        cat("Iteration\tBetween variance estimator\n")
        exp <- expression(cat(" ", count, "\t\t ", a1 <- a, fill = TRUE))
    }
    else
        exp <- expression(a1 <-  a)

    a <- start
    count <- 0

    repeat
    {
        eval(exp)

        if (maxit < (count <- count + 1))
        {
            warning("maximum number of iterations reached before obtaining convergence")
            break
        }

        cred <- 1 / (1 + within/(w * a))
        x.z <- drop(crossprod(cred, x)) / sum(cred)
        a <- drop(crossprod(cred, (x - x.z)^2)) / (n - 1)

        if (abs((a - a1)/a1) < tol)
            break
    }
    a
}
