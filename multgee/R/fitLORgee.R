fitLORgee <-
function (Y, X_mat, coeffs, ncategories, id, repeated, offset, 
    link, LORterm, marpars, ipfp.ctrl, control, IM, LORem = LORem, 
    LORstr = LORstr, add) 
{
    inversematrix <- function(x) inversemat(x,IM)
    tol <- control$tolerance
    maxiter <- control$maxiter
    verbose <- control$verbose
    ncategoriesm1 <- ncategories - 1
    p <- ncol(X_mat)
    Nsubs <- max(id)
    nobs <- length(id)/ncategoriesm1
    repeatednew <- as.numeric(factor(repeated))
    noccasions <- max(repeatednew)
    noccasionpairs <- choose(noccasions, 2)
    Ti_vector <- as.numeric(lapply(split(id, id), length))/(ncategoriesm1)
    if (all(Ti_vector == 1)) 
        stop("There are no repeated responses")
    sel_mat <- lapply(seq.int(2, noccasions), function(x) which(lower.tri(diag(x * 
        ncategoriesm1))))
    if (any(Ti_vector < noccasions)) {
        times_mat <- lapply(split(repeatednew, id), unique)
    }
    code_mat <- combns(noccasions)
    if (link != "acl" & link != "bcl") {
        family <- make.link(link)
        linkinv <- family$linkinv
        mu.eta <- family$mu.eta
    }
    validmu <- function(mu) all(mu > 0) && all(mu < 1)
    if (LORstr == "fixed") {
        if (!is.matrix(LORterm)) 
            stop("'LORterm' must be inserted as matrix")
        if (ncol(LORterm) != ncategories^2) 
            stop("'LORterm' must have ", ncategories^2, " columns")
        if (nrow(LORterm) != noccasionpairs) 
            stop("'LORterm' must have ", noccasionpairs, " rows")
        LORterm <- prop.table(LORterm, 1)
    }
    if (LORstr == "independence") {
        LORterm <- matrix(1/ncategories^2, noccasionpairs, 
            ncategories^2)
    }
    formals(ipfp)$dimension <- ncategories
    formals(ipfp)$tol <- ipfp.ctrl$tol
    formals(ipfp)$maxit <- ipfp.ctrl$maxit
    index_mat <- matrix(seq.int(ncategoriesm1 * noccasions), 
        noccasions, ncategoriesm1, TRUE)
    extindex_mat <- matrix(seq.int(ncategories * noccasions), 
        noccasions, ncategories, TRUE)
    id_vector <- split(seq.int(nrow(X_mat)), id)
    extid_vector <- split(seq.int(ncategories * nobs), rep.int(seq.int(Nsubs), 
        Ti_vector * ncategories))
    deriv_mat <- if (link == "acl") 
        derivacl
    else if (link == "bcl") {
        derivbcl
    }
    else {
        derivmclm
    }
    eta <- drop(X_mat %*% coeffs) + offset
    if (link != "acl" & link != "bcl") {
        probexclude <- seq(ncategories, nobs * ncategories, ncategories)
        fitproball <- muprob(linkinv(eta), nobs, ncategoriesm1)
        if (!validmu(fitproball)) 
            stop("Please insert initial values")
        fitprob <- fitproball[-probexclude]
        dummy <- mu.eta(eta)
    }
    else {
        fitprob <- exp(matrix(eta, length(eta)/ncategoriesm1, 
            ncategoriesm1, TRUE))
        fitprob <- fitprob/(1 + .rowSums(fitprob, nobs, ncategoriesm1,FALSE))
        fitproball <- as.vector(t(cbind(fitprob, 1 - .rowSums(fitprob, 
            nobs, ncategoriesm1,FALSE))))
        if (!validmu(fitproball)) 
            stop("Please insert initial values")
        fitprob <- as.vector(t(fitprob))
        dummy <- fitprob
    }
   if (verbose) 
 cat(format("Iteration", digits = 5, trim = TRUE, 
           justify = "centre", scientific = FALSE, nsmall = 0, 
           width = 8), "\t", format("Criterion", digits = 5, 
           trim = TRUE, justify = "centre", scientific = FALSE, 
           nsmall = 5, width = 8), "\n")
    beta_mat <- matrix(coeffs)
    crit_vector <- 20
    I0_mat <- I1_mat <- matrix(0, p, p, FALSE)
    H_mat <- matrix(0, p, 1, FALSE)
        for (iter in 1:maxiter) {
            resids <- Y - fitprob
            D_mat <- deriv_mat(dummy, ncategoriesm1, X_mat)
            for (i in 1:Nsubs) {
                Ti <- Ti_vector[[i]]
                id_ind <- id_vector[[i]]
                prob <- fitprob[id_ind]
                mat1 <- diagmod(prob)
                if (Ti == 1) {
                  V_mat <- inversematrix(mat1 - tcrossprod(prob, 
                    prob))
                }
                else if (Ti == noccasions) {
                  proball <- fitproball[extid_vector[[i]]]
                  for (j in 1:(noccasions - 1)) {
                    t1 <- index_mat[j, ]
                    t2 <- extindex_mat[j, ]
                    probrow <- proball[t2]
                    for (k in (j + 1):noccasions) {
                      t3 <- extindex_mat[k, ]
                      index <- code_mat[code_mat[, 1] == j & 
                        code_mat[, 2] == k, 3]
                      ipfpfit <- ipfp(LORterm[index, ], probrow, 
                        proball[t3])
                      mat1[t1, index_mat[k, ]] <- ipfpfit[-ncategories, 
                        -ncategories]
                    }
                  }
                  mat1[sel_mat[[noccasions - 1]]] <- aperm(mat1, 
                    c(2, 1))[sel_mat[[noccasions - 1]]]
                  V_mat <- inversematrix(mat1 - tcrossprod(prob, 
                    prob))
                }
                else {
                  Tiid <- times_mat[[i]]
                  proball <- fitproball[extid_vector[[i]]]
                  for (j in 1:(Ti - 1)) {
                    t1 <- index_mat[j, ]
                    t2 <- extindex_mat[Tiid[j], ]
                    probrow <- proball[extindex_mat[j, ]]
                    for (k in (j + 1):Ti) {
                      t3 <- extindex_mat[Tiid[k], ]
                      index <- code_mat[code_mat[, 1] == Tiid[j] & 
                        code_mat[, 2] == Tiid[k], 3]
                      ipfpfit <- ipfp(LORterm[index, ], probrow, 
                        proball[extindex_mat[k, ]])
                      mat1[t1, index_mat[k, ]] <- ipfpfit[-ncategories, 
                        -ncategories]
                    }
                  }
                  mat1[sel_mat[[Ti - 1]]] <- aperm(mat1, c(2, 
                    1))[sel_mat[[Ti - 1]]]
                  V_mat <- inversematrix(mat1 - tcrossprod(prob, 
                    prob))
                }
                D_mat1 <- D_mat[id_ind, ]
                help_mat1 <- crossprod(D_mat1, V_mat)
                help_mat2 <- help_mat1 %*% resids[id_ind]
                I0_mat <- help_mat1 %*% D_mat1 + I0_mat
                H_mat <- help_mat2 + H_mat
                I1_mat <- tcrossprod(help_mat2, help_mat2) + 
                  I1_mat
            }
            naive_mat <- solve(I0_mat)
            robust_mat <- naive_mat %*% I1_mat %*% naive_mat
            if (any(eigen(robust_mat,TRUE,only.values=TRUE)$values<=0)) 
                stop("Robust covariance matrix is not positive definite")
            coeffsnew <- coeffs + naive_mat %*% H_mat
            crit <- max(abs(coeffs - coeffsnew)/pmax.int(abs(coeffs), 
                1e-06))
            crit_vector <- c(crit_vector, crit)
            beta_mat <- cbind(beta_mat, coeffsnew)
            coeffs <- coeffsnew
            eta <- drop(X_mat %*% coeffs) + offset
            if (link != "acl" & link != "bcl") {
                probexclude <- seq(ncategories, nobs * ncategories, 
                  ncategories)
                fitproball <- muprob(linkinv(eta), nobs, ncategoriesm1)
                if (!validmu(fitproball)) 
                  stop("Please insert initial values")
                fitprob <- fitproball[-probexclude]
                dummy <- mu.eta(eta)
            }
            else {
                fitprob <- exp(matrix(eta, length(eta)/ncategoriesm1, 
                  ncategoriesm1, TRUE))
                fitprob <- fitprob/(1 + .rowSums(fitprob, nobs, 
                  ncategoriesm1,FALSE))
                fitproball <- as.vector(t(cbind(fitprob, 1 - 
                  .rowSums(fitprob, nobs, ncategoriesm1,FALSE))))
                if (!validmu(fitproball)) 
                  stop("Please insert initial values")
                fitprob <- as.vector(t(fitprob))
                dummy <- fitprob
            }
            I0_mat[, ] <- I1_mat[, ] <- H_mat[, ] <- 0
            if (verbose) 
                cat(format(round(iter, 0), digits = 5, trim = TRUE, 
                  justify = "centre", scientific = FALSE, nsmall = 0, 
                  width = 8), "\t", format(round(crit, 5), digits = 5, 
                  trim = TRUE, justify = "centre", scientific = FALSE, 
                  nsmall = 5, width = 8), "\n")
            if (crit <= tol) 
                break
        }
    resids <- Y - fitprob
    D_mat <- deriv_mat(dummy, ncategoriesm1, X_mat)
    for (i in 1:Nsubs) {
        Ti <- Ti_vector[[i]]
        id_ind <- id_vector[[i]]
        prob <- fitprob[id_ind]
        mat1 <- diagmod(prob)
        if (Ti == 1) {
            V_mat <- inversematrix(mat1 - tcrossprod(prob, prob))
        }
        else if (Ti == noccasions) {
            proball <- fitproball[extid_vector[[i]]]
            for (j in 1:(noccasions - 1)) {
                t1 <- index_mat[j, ]
                probrow <- proball[extindex_mat[j, ]]
                for (k in (j + 1):noccasions) {
                  index <- code_mat[code_mat[, 1] == j & code_mat[, 
                    2] == k, 3]
                  ipfpfit <- ipfp(LORterm[index, ], probrow, 
                    proball[extindex_mat[k, ]])
                  mat1[t1, index_mat[k, ]] <- ipfpfit[-ncategories, 
                    -ncategories]
                }
            }
            mat1[sel_mat[[noccasions - 1]]] <- aperm(mat1, c(2, 
                1))[sel_mat[[noccasions - 1]]]
            V_mat <- inversematrix(mat1 - tcrossprod(prob, prob))
        }
        else {
            Tiid <- times_mat[[i]]
            proball <- fitproball[extid_vector[[i]]]
            for (j in 1:(Ti - 1)) {
                t1 <- index_mat[j, ]
                probrow <- proball[extindex_mat[j, ]]
                for (k in (j + 1):Ti) {
                  index <- code_mat[code_mat[, 1] == Tiid[j] & 
                    code_mat[, 2] == Tiid[k], 3]
                  ipfpfit <- ipfp(LORterm[index, ], probrow, 
                    proball[extindex_mat[k, ]])
                  mat1[t1, index_mat[k, ]] <- ipfpfit[-ncategories, 
                    -ncategories]
                }
            }
            mat1[sel_mat[[Ti - 1]]] <- aperm(mat1, c(2, 1))[sel_mat[[Ti - 
                1]]]
            V_mat <- inversematrix(mat1 - tcrossprod(prob, prob))
        }
        D_mat1 <- D_mat[id_ind, ]
        help_mat1 <- crossprod(D_mat1, V_mat)
        help_mat2 <- help_mat1 %*% resids[id_ind]
        I0_mat <- help_mat1 %*% D_mat1 + I0_mat
        I1_mat <- tcrossprod(help_mat2, help_mat2) + I1_mat
    }
    if (any(eigen(I0_mat,TRUE,only.values=TRUE)$values<=0)) 
        warning("'Naive' covariance matrix is not positive definite")
    naive_mat <- solve(I0_mat)
    robust_mat <- naive_mat %*% I1_mat %*% naive_mat
    if (any(eigen(robust_mat,TRUE,only.values=TRUE)$values<=0)) 
        stop("Robust covariance matrix is not positive definite")
    result <- list()
    result$beta_mat <- beta_mat
    result$naive <- naive_mat
    result$robust <- robust_mat
    result$crit <- crit_vector[-1]
    result$iter <- length(result$crit)
    ans <- diagmod(rep(0, ncategoriesm1 * noccasions))
    k <- 1
    for (i in 1:(noccasions - 1)) {
        for (j in (i + 1):noccasions) {
            ans[seq(ncategoriesm1) + ncategoriesm1 * (i - 1), 
                seq(ncategoriesm1) + ncategoriesm1 * (j - 1)] <- odds.ratio(matrix(LORterm[k, 
                ], ncategories, ncategories))
            k <- k + 1
        }
    }
    result$theta <- ans + t(ans)
    result$conv <- (result$crit[result$iter] <= tol)
    result$linear.predictor <- eta
    result$fitted.values <- fitprob
    result$residuals <- resids
    result
}

