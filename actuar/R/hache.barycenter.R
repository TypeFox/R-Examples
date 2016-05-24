### ===== actuar: An R Package for Actuarial Science =====
###
### Auxiliary function to fit regression credibility model by
### positioning the intercept at the barycenter of time.
###
### AUTHORS: Xavier Milhaud, Vincent Goulet <vincent.goulet@act.ulaval.ca>

hache.barycenter <- function(ratios, weights, xreg, method,
                             tol, maxit, echo)
{
    ## Frequently used values
    weights.s <- rowSums(weights, na.rm = TRUE) # contract total weights
    has.data <- which(weights.s > 0)	# contracts with data
    ncontracts <- nrow(ratios)   	# number of contracts
    eff.ncontracts <- length(has.data)  # effective number of contracts
    p <- ncol(xreg)               	# rank (>= 2) of design matrix
    n <- nrow(xreg)                     # number of observations

    ## Putting the intercept at the barycenter of time amounts to use
    ## a "weighted orthogonal" design matrix in the regression (that
    ## is, X'WX = I for some diagonal weight matrix W). In theory,
    ## there would be one orthogonal design matrix per contract. In
    ## practice, we orthogonalize with a "collective" barycenter. We
    ## use average weights per period across contracts since these
    ## will be closest to the their individual analogues.
    ##
    ## We orthogonalize the original design matrix using QR
    ## decomposition. We also keep matrix R as a transition matrix
    ## between the original base and the orthogonal base.
    w <- colSums(weights, na.rm = TRUE)/sum(weights.s)
    Xqr <- qr(xreg * sqrt(w))   # QR decomposition
    R <- qr.R(Xqr)              # transition matrix
    x <- qr.Q(Xqr) / sqrt(w)    # weighted orthogonal matrix

    ## Fit linear model to each contract. For contracts without data,
    ## fit some sort of empty model to ease use in predict.hache().
    f <- function(i)
    {
    	z <-
            if (i %in% has.data)            # contract with data
            {
                y <- ratios[i, ]
                not.na <- !is.na(y)
                lm.wfit(x[not.na, , drop = FALSE], y[not.na], weights[i, not.na])
            }
            else                            # contract without data
                lm.fit(x, rep.int(0, n))
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
        rank <- z$rank
        r <- z$residuals
        w <- z$weights
        sum(w * r^2) / (nQr - rank)
    }
    sigma2 <- sapply(fits[has.data], S)
    sigma2[is.nan(sigma2)] <- 0

    ## Initialization of a few containers: p x p x ncontracts arrays
    ## for the weight and credibility matrices; p x p matrices for the
    ## between variance-covariance matrix and total weight matrix; a
    ## vector of length p for the collective regression coefficients.
    cred <- W <- array(0, c(p, p, ncontracts))
    A <- W.s <- matrix(0, p, p)
    coll <- numeric(p)

    ## Weight matrices: we need here only the diagonal elements of
    ## X'WX, where W = diag(w_{ij}) (and not w_{ij}/w_{i.} as in the
    ## orthogonalization to keep a w_{i.} lying around). The first
    ## element is w_{i.} and the off-diagonal elements are zero by
    ## construction. Note that array W is quite different from the one
    ## in hache.origin().
    W[1, 1, ] <- weights.s
    for (i in 2:p)
        W[i, i, has.data] <-
            colSums(t(weights[has.data, ]) * x[, i]^2, na.rm = TRUE)

    ## === ESTIMATION OF THE WITHIN VARIANCE ===
    s2 <- mean(sigma2)

    ## === ESTIMATION OF THE BETWEEN VARIANCE-COVARIANCE MATRIX ===
    ##
    ## By construction, we only estimate the diagonal of the matrix.
    ## Variance components are estimated just like in the
    ## Buhlmann-Straub model (see bstraub.R for details).
    ##
    ## Should we compute the iterative estimators?
    do.iter <- method == "iterative" && diff(range(weights, na.rm = TRUE)) > .Machine$double.eps^0.5

    ## Do the computations one regression parameter at a time.
    for (i in seq_len(p))
    {
        ## Unbiased estimator
        a <- A[i, i] <- bvar.unbiased(ind[i, has.data], W[i, i, has.data],
                                      s2, eff.ncontracts)

        ## Iterative estimator
        if (do.iter)
        {
            a <- A[i, i] <-
                if (a > 0)
                    bvar.iterative(ind[i, has.data], W[i, i, has.data],
                                   s2, eff.ncontracts, start = a,
                                   tol = tol, maxit = maxit, echo = echo)
                else
                    0
        }

        ## Credibility factors and estimator of the collective
        ## regression coefficients.
        if (a > 0)
        {
            z <- cred[i, i, has.data] <- 1/(1 + s2/(W[i, i, has.data] * a))
            z. <- W.s[i, i] <- sum(z)
            coll[i] <- drop(crossprod(z, ind[i, has.data])) / z.
        }
        else
        {
            ## (credibility factors were already initialized to 0)
            w <- W[i, i, has.data]
            w. <- W.s[i, i] <- sum(w)
            coll[i] <- drop(crossprod(w, ind[i, ])) / w.
        }
    }

    ## Credibility adjusted coefficients. The coefficients of the
    ## models are replaced with these values. That way, prediction
    ## will be trivial using predict.lm().
    for (i in seq_len(ncontracts))
        fits[[i]]$coefficients <- coll + drop(cred[, , i] %*% (ind[, i] - coll))

    ## Add names to the collective coefficients vector.
    names(coll) <- rownames(ind)

    ## Results
    list(means = list(coll, ind),
         weights = list(W.s, W),
         unbiased = if (method == "unbiased") list(A, s2),
         iterative = if (method == "iterative") list(A, s2),
         cred = cred,
         nodes = list(ncontracts),
         adj.models = fits,
         transition = R)
}
