### ===== actuar: An R Package for Actuarial Science =====
###
### Auxiliary function to fit regression credibility model using the
### original Hachemeister model.
###
### AUTHORS: Tommy Ouellet, Vincent Goulet <vincent.goulet@act.ulaval.ca>

### codetools does not like the way 'coll1' is defined in function
### 'hache.origin' below. Avoid false positive in R CMD check.
if(getRversion() >= "2.15.1")  utils::globalVariables(c("coll1"))

hache.origin <- function(ratios, weights, xreg, tol, maxit, echo)
{
    ## Frequently used values
    weights.s <- rowSums(weights, na.rm = TRUE) # contract total weights
    has.data <- which(weights.s > 0)	# contracts with data
    ncontracts <- nrow(ratios)   	# number of contracts
    eff.ncontracts <- length(has.data)  # effective number of contracts
    p <- ncol(xreg)               	# rank (>= 2) of design matrix
    n <- nrow(xreg)                     # number of observations

    ## Fit linear model to each contract. For contracts without data,
    ## fit some sort of empty model to ease use in predict.hache().
    f <- function(i)
    {
    	z <-
            if (i %in% has.data)            # contract with data
            {
                y <- ratios[i, ]
                not.na <- !is.na(y)
                lm.wfit(xreg[not.na, , drop = FALSE], y[not.na], weights[i, not.na])
            }
            else                            # contract without data
                lm.fit(xreg, rep.int(0, n))
        z[c("coefficients", "residuals", "weights", "rank", "qr")]
    }
    fits <- lapply(seq_len(ncontracts), f)

    ## Individual regression coefficients
    ind <- sapply(fits, coef)
    ind[is.na(ind)] <- 0

    ## Individual variance estimators. The contribution of contracts
    ## without data is 0.
    S <- function(z)                    # from stats:::summary.lm
    {
        nQr <- NROW(z$qr$qr)
        r1 <- z$rank
        r <- z$residuals
        w <- z$weights
        sum(w * r^2) / (nQr - r1)
    }
    sigma2 <- sapply(fits[has.data], S)
    sigma2[is.nan(sigma2)] <- 0

    ## Initialization of a few containers: p x p x ncontracts arrays
    ## for the weight and credibility matrices; p x p matrices for the
    ## between variance-covariance matrix and total credibility
    ## matrix.
    cred <- W <- array(0, c(p, p, ncontracts))
    A <- cred.s <- matrix(0, p, p)

    ## Weight matrices: we use directly (X'WX)^{-1}. This is quite
    ## different from hache.barycenter().
    V <- function(z)                # from stats:::summary.lm
    {
        r1 <- z$rank
        if (r1 == 1L)
            diag(as.double(chol2inv(z$qr$qr[1L, 1L, drop = FALSE])), p)
        else
            chol2inv(z$qr$qr[1L:r1, 1L:r1, drop = FALSE])
    }
    W[, , has.data] <- sapply(fits[has.data], V)

    ## Starting credibility matrices and collective regression
    ## coefficients.
    cred[, , has.data] <- diag(p)         # identity matrices
    coll <- rowSums(ind) / eff.ncontracts # coherent with above

    ## === ESTIMATION OF WITHIN VARIANCE ===
    s2 <- mean(sigma2)

    ## === ESTIMATION OF THE BETWEEN VARIANCE-COVARIANCE MATRIX ===
    ##
    ## This is an iterative procedure similar to the Bischel-Straub
    ## estimator. Following Goovaerts & Hoogstad, stopping criterion
    ## is based in the collective regression coefficients estimates.
    ##
    ## If printing of iterations was asked for, start by printing a
    ## header and the starting values.
    if (echo)
    {
        cat("Iteration\tCollective regression coefficients\n")
        exp <- expression(cat(" ", count, "\t\t ", coll1 <- coll,
                              fill = TRUE))
    }
    else
        exp <- expression(coll1 <-  coll)

    ## Iterative procedure
    count <- 0
    repeat
    {
        eval(exp)

        ## Stop after 'maxit' iterations
        if (maxit < (count <- count + 1))
        {
            warning("maximum number of iterations reached before obtaining convergence")
            break
        }

        ## Calculation of the between variance-covariance matrix.
        A[] <- rowSums(sapply(has.data, function(i)
                              cred[, , i] %*% tcrossprod(ind[, i] - coll))) /
                                  (eff.ncontracts - 1)

        ## Symmetrize A
        A <- (A + t(A))/2

        ## New credibility matrices
        cred[, , has.data] <- sapply(has.data, function(i)
                                     A %*% solve(A + s2 * W[, , i]))

        ## New collective regression coefficients
        cred.s <- apply(cred[, , has.data], c(1L, 2L), sum)
        coll <- solve(cred.s,
                      rowSums(sapply(has.data, function(i)
                                     cred[, , i] %*% ind[, i])))

        ## Test for convergence
        if (max(abs((coll - coll1)/coll1)) < tol)
            break
    }

    ## Final calculation of the between variance-covariance matrix and
    ## credibility matrices.
    A[] <- rowSums(sapply(has.data, function(i)
                          cred[, , i] %*% tcrossprod(ind[, i] - coll))) /
                              (eff.ncontracts - 1)
    A <- (A + t(A))/2
    cred[, , has.data] <- sapply(has.data, function(i)
                                 A %*% solve(A + s2 * W[, , i]))

    ## Credibility adjusted coefficients. The coefficients of the
    ## models are replaced with these values. That way, prediction
    ## will be trivial using predict.lm().
    for (i in seq_len(ncontracts))
        fits[[i]]$coefficients <- coll + drop(cred[, , i] %*% (ind[, i] - coll))

    ## Add names to the collective coefficients vector.
    names(coll) <- rownames(ind)

    ## Results
    list(means = list(coll, ind),
         weights = list(cred.s, W),
         unbiased = NULL,
         iterative = list(A, s2),
         cred = cred,
         nodes = list(ncontracts),
         adj.models = fits)
}
