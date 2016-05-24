#' Wrapper functions for some previous methods
#' 
#' @description These functions provide an uniform interface to three existing methods: SVA, RUV, LEAPP
#' The wrapper functions transform the data into desired forms and call the corresponding functions in the package
#' \link{sva}, \link{ruv}, \link[leapp]{leapp}
#'
#' @inheritParams cate
#' 
#' 
#' @return All functions return \code{beta.p.value} which are the p-values after adjustment. 
#' For the other returned objects, refer to \link{cate} for their meaning.
#' 
#' @examples
#' ## this is the simulation example in Wang et al. (2015).
#' n <- 100
#' p <- 1000
#' r <- 2
#' set.seed(1)
#' data <- gen.sim.data(n = n, p = p, r = r, 
#'                      alpha = rep(1 / sqrt(r), r), 
#'                      beta.strength = 3 * sqrt(1 + 1) / sqrt(n), 
#'                      Gamma.strength = c(seq(3, 1, length = r)) * sqrt(p),
#'                      Sigma = 1 / rgamma(p, 3, rate = 2),
#'                      beta.nonzero.frac = 0.05)
#' X.data <- data.frame(X1 = data$X1)
#' sva.results <- sva.wrapper(~ X1, X.data, data$Y,
#'                            r = r, sva.method = "irw")
#' ruv.results <- ruv.wrapper(~ X1, X.data, data$Y, r = r,  
#'                            nc = sample(data$beta.zero.pos, 30), ruv.method = "RUV4")
#' leapp.results <- leapp.wrapper(~ X1, X.data, data$Y, r = r)
#' cate.results <- cate(~ X1, X.data, data$Y, r = r)
#' 
#' ## p-values after adjustment
#' par(mfrow = c(2, 2))
#' hist(sva.results$beta.p.value)
#' hist(ruv.results$beta.p.value)
#' hist(leapp.results$beta.p.value)
#' hist(cate.results$beta.p.value)
#' 
#' ## type I error
#' mean(sva.results$beta.p.value[data$beta.zero.pos] < 0.05)
#' 
#' ## power
#' mean(sva.results$beta.p.value[data$beta.nonzero.pos] < 0.05)
#' 
#' ## false discovery proportion for sva
#' discoveries.sva <- which(p.adjust(sva.results$beta.p.value, "BH") < 0.2)
#' fdp.sva <- length(setdiff(discoveries.sva, data$beta.nonzero.pos)) / max(length(discoveries.sva), 1)
#' fdp.sva
#' 
#' @name wrapper
NULL

#' @rdname wrapper
#' @param sva.method parameter for \code{\link[sva]{sva}}. 
#'        whether to use an iterative reweighted algorithm (irw) or a two-step algorithm (two-step).
#' @param B parameter for \code{\link[sva]{sva}}. the number of iterations of the irwsva algorithm
#' @details The \code{beta.p.values} returned is a length \code{p} vector, each for the overall effects of all the primary variables.
#'
#' @import sva stats
#' @export
#'
sva.wrapper <- function(formula,
                        X.data = NULL,
                        Y,
                        r,
                        sva.method = c("irw", "two-step"),
                        B = 5) {
  
    method <- match.arg(sva.method, c("irw", "two-step"))
    dat <- t(Y)
    
    X <- parse.cate.formula(formula, X.data)      
    X.primary <- X$X.primary
    X.nuis <- X$X.nuis

	mod <- cbind(X.primary, X.nuis)
    mod0 <- X.nuis
    if (ncol(X.nuis) == 0)
    	mod0 <- NULL
  
    result <- sva(dat, mod, mod0, r, method = method, B= B)
    Z <- result$sv
    rownames(Z) <- rownames(X)
    modSV <- cbind(mod, Z)
    mod0SV <- cbind(mod0, Z)
    ## only can calculate a vector of p-values. It's the p-values of the effect of all the primary variables as a whole
    p.values <- f.pvalue(dat, modSV, mod0SV)
    return(list(beta.p.value = p.values, Z= Z, beta.p.post = result$pprob.b))

}

#' @rdname wrapper
#' @param ruv.method either using \code{\link[ruv]{RUV2}}, \code{\link[ruv]{RUV4}} or
#'    	         \code{\link[ruv]{RUVinv}} functions
#' @param nc parameter for \link[ruv]{ruv} functions: position of the negative controls
#' @param lambda parameter for \code{\link[ruv]{RUVinv}}
#'  		         
#' @import ruv
#' @export
#' 
ruv.wrapper <- function(formula,
                        X.data = NULL,
                        Y,
                        r,
                        nc, 
                        lambda = 1,
                        ruv.method = c("RUV2", "RUV4", "RUVinv")) {
  
    method <- match.arg(ruv.method, c("RUV2", "RUV4", "RUVinv"))

	X <- parse.cate.formula(formula, X.data)      
    X.primary <- X$X.primary
    X.nuis <- X$X.nuis
    
    X <- X.primary
    Z <- X.nuis
    if (ncol(X.nuis) == 0)
    	Z <- NULL
    
    p <- ncol(Y)
    n <- nrow(Y)
    ctl <- rep(F, p)
    ctl[nc] <- T
    if (method == "RUV2"){
        result <- RUV2(Y, X, ctl, r, Z)
    } else if (method == "RUV4") {
        result <- RUV4(Y, X, ctl, r, Z)
    } else if (method == "RUVinv") {
        ## uses p - r0 - r1 latent factors!
        if (length(nc) > n) {
            result <- RUVinv(Y, X, ctl, Z)
        } else
            result <- RUVinv(Y, X, ctl, Z, lambda = lambda)
    }
    beta <- t(result$betahat)
    Gamma <- t(result$alpha/sqrt(n))
    beta.t <- t(result$t)
    beta.p.value <- t(result$p)
    Sigma <- result$sigma2
    names(Sigma) <- colnames(Y)
    return(list(Gamma = Gamma, Sigma = Sigma, beta = beta, beta.t = beta.t,
                beta.p.value = beta.p.value, Z = result$W))
}

#' @rdname wrapper
#' @param search.tuning logical parameter for \code{\link[leapp]{leapp}}, whether using BIC to search for tuning parameter of IPOD. 
#' @param ipod.method parameter for \code{\link[leapp]{leapp}}. "hard": hard thresholding in the IPOD algorithm;
#' "soft": soft thresholding in the IPOD algorithm
#' 
#' @details Only 1 variable of interest is allowed for \code{leapp.wrapper}. The method can be slow.
#' 
#' @import corpcor
#' @export
#' 
leapp.wrapper <- function(formula,
                          X.data = NULL,
                          Y,
                          r,
                          search.tuning = F,
                          ipod.method = c("hard", "soft")) {
    method <- match.arg(ipod.method, c("hard", "soft"))

	X <- parse.cate.formula(formula, X.data)      
    X.primary <- X$X.primary
    X.nuis <- X$X.nuis

    n <- nrow(X.primary)
    data <- t(Y)
    pred.prim <- X.primary
    pred.covar <- X.nuis
 
    if (ncol(X.nuis) == 1 & (sum(X.nuis == 1) == n))
        pred.covar <- NULL
    if (search.tuning) {
        result <- leapp(data, pred.prim, pred.covar, num.fac = r, method = method)
    } else {
        result <- leapp(data, pred.prim, pred.covar, num.fac = r, method = method, length.out = 1)
    }
    beta.p.value <- result$p
    beta <- result$gamma
    names(beta) <- colnames(Y)
    Gamma <- result$uest
    Sigma <- result$sigma^2
    beta.t <- result$resOpt.scale
    return(list(Gamma = Gamma, Sigma = Sigma, beta = beta,
                beta.t = beta.t, beta.p.value = beta.p.value))

}


#' changed IPOD function
#' 
#' @importFrom leapp IPODFUN
#' @import stats
#' 
#' @keywords internal
#' 
IPOD <-
    function (X, Y, H, method = "hard", TOL = 1e-04, length.out = 50)
{
                                        #print("Using modified IPOD...")
    require
    if (is.null(X)) {
        r = 0
    }
    else {
        r = ncol(X)
    }
    N = length(Y)
    ress = NULL
    gammas = NULL
    if (is.null(X)) {
        betaInit = NULL
        lambdas = seq(round(norm(matrix(Y, N, 1), "I")/1 + 1),
            0, by = -0.1)
    }
    else {
        betaInit = rlm(Y ~ X - 1)$coefficients
        tmp = t(diag(ncol(H)) - H) %*% Y/sqrt(1 - diag(H))
        lambdas = seq(round(norm(tmp, "I")/1 + 1), 0, by = -0.1)
        if (length.out == 1)
            lambdas = round(norm(tmp, "I")/1 + 1)/2
    }
    for (sigma in lambdas) {
        sigma = sigma/sqrt(2 * log(N))
        result = IPODFUN(X, Y, H, sigma, betaInit, method = method,
            TOL = TOL)
        gammas = cbind(gammas, result$gamma)
        ress = cbind(ress, result$ress)
    }
    DF = colSums(abs(gammas) > 1e-05)
    if (!is.null(X)) {
        Q = qr.Q(qr(X), complete = TRUE)
        X.unscaled = t(Q[, (r + 1):N])
        Y.eff = X.unscaled %*% Y
        sigmaSqEsts = colSums((Y.eff %*% matrix(rep(1, ncol(gammas)),
            1, ncol(gammas)) - X.unscaled %*% gammas)^2)/(length(Y.eff) -
                                                              DF)
    }
    else {
        sigmaSqEsts = colSums((Y %*% matrix(rep(1, ncol(gammas)),
            1, ncol(gammas)) - gammas)^2)/(length(Y) - DF)
    }
    sigmaSqEsts[sigmaSqEsts < 0] = 0
    if (length.out == 1) {
    	riskEst = ((N - r) * log(sigmaSqEsts * (N - r - DF)/(N -
                                                                 r)) + (log(N - r) + 1) * (DF + 1))/(N - r)
    	optSet = which(riskEst == min(riskEst))
    	gammaOptInd = optSet[which(DF[optSet] == min(DF[optSet]))[1]]
    	gammaOpt = gammas[, gammaOptInd]
    	resOpt = ress[, gammaOptInd]
    	tau = mad(ress[gammas[, gammaOptInd] == 0, gammaOptInd])
    	resOpt.scale = resOpt/tau
    	p = 2 * pnorm(-abs(resOpt.scale))
    	return(list(p = p, resOpt.scale = resOpt.scale, gamma = gammaOpt))
    }

    if (!all(DF <= N - r)) {
        gammas = gammas[, DF <= N - r]
        sigmaSqEsts = sigmaSqEsts[DF <= N - r]
        lambdas = lambdas[DF <= N - r]
        ress = ress[, DF <= N - r]
        DF = DF[DF <= N - r]
    }
    redundantInds = which(sum(abs(gammas[, -1] - gammas[, 1:(ncol(gammas) -
                                                                 1)])) < 0.001) + 1
    if (length(redundantInds) > 0) {
        gammas = gammas[, -redundantInds]
        lambdas = lambdas[-redundantInds]
        sigmaSqEsts = sigmaSqEsts[-redundantInds]
        ress = ress[, -redundantInds]
        DF = DF[-redundantInds]
    }
    redundantInds = which(abs(sigmaSqEsts[-1] - sigmaSqEsts[1:(length(sigmaSqEsts) -
                                                                   1)]) < 0.001) + 1
    if (length(redundantInds) > 0) {
        gammas = gammas[, -redundantInds]
        lambdas = lambdas[-redundantInds]
        sigmaSqEsts = sigmaSqEsts[-redundantInds]
        ress = ress[, -redundantInds]
        DF = DF[-redundantInds]
    }
    if (!all(DF <= N/2)) {
        gammas = gammas[, DF <= N/2]
        sigmaSqEsts = sigmaSqEsts[DF <= N/2]
        lambdas = lambdas[DF <= N/2]
        ress = ress[, DF <= N/2]
        DF = DF[DF <= N/2]
    }
    riskEst = ((N - r) * log(sigmaSqEsts * (N - r - DF)/(N -
                                                             r)) + (log(N - r) + 1) * (DF + 1))/(N - r)
    optSet = which(riskEst == min(riskEst))
    gammaOptInd = optSet[which(DF[optSet] == min(DF[optSet]))[1]]
    gammaOpt = gammas[, gammaOptInd]
    resOpt = ress[, gammaOptInd]
    tau = mad(ress[gammas[, gammaOptInd] == 0, gammaOptInd])
    resOpt.scale = resOpt/tau
    p = 2 * pnorm(-abs(resOpt.scale))
    return(list(p = p, resOpt.scale = resOpt.scale, gamma = gammaOpt))
}

#' original leapp function only added the resOpt.scale(beta.t) return
#' 
#' @importFrom leapp IPODFUN ridge AlternateSVD
#' @importFrom grDevices dev.off png
#' @importFrom graphics plot
#' @keywords internal
#' 
leapp <-
    function (data, pred.prim, pred.covar = NULL, O = NULL, num.fac = "buja",
              method = "hard", sparse = TRUE, centered = FALSE, verbose = FALSE,
              perm.num = 50, TOL = 1e-04, length.out = 50)
{
    N = nrow(data)
    n = ncol(data)
    pred.prim = matrix(pred.prim, n, 1)
    if (sum(pred.prim) != 0)
        pred.prim = pred.prim - mean(pred.prim)
    if (!centered)
        data = t(apply(data, 1, scale, center = TRUE, scale = FALSE))
    if (num.fac == "buja")
        num.fac = num.sv(data, cbind(pred.prim, pred.covar),
            B = perm.num)
    if (is.null(O))
        O = qr.Q(qr(pred.prim), complete = TRUE)
    if ((O %*% pred.prim)[1] < 0)
        O = -O
    if (!is.null(pred.covar)) {
        s = ncol(pred.covar)
        pred.covar = O %*% pred.covar
        pred.covar.v1 = pred.covar[1, ]
        pred.covar.rest = pred.covar[-1, ]
    }
    else {
        s = 0
        pred.covar.rest = NULL
    }
    data = data %*% t(O)
    Y = data[, 1]
    result.svd = AlternateSVD(data[, -1], pred = pred.covar.rest,
        r = num.fac)
    sigma.all = result.svd$sigma
    X = result.svd$uest
    if (is.null(pred.covar)) {
        Y = Y/sigma.all
    }
    else {
        Y = (Y - result.svd$coef %*% t(pred.covar.v1))/sigma.all
    }
    if (!is.null(X)) {
        X = X/sigma.all
        H = X %*% solve(t(X) %*% X) %*% t(X)
    }
    else {
        H = NULL
    }
    if (!sparse) {
        result = ridge(X, Y, H, sigma = (n - s - num.fac - 1)/(n -
                                                                   s - num.fac - 3))
    }
    else {
        result = IPOD(X, Y, H, method, TOL = TOL, length.out = length.out)
    }
    RetList = list(p = result$p, vest = result.svd$vest, uest = result.svd$uest,
        gamma = as.vector(sigma.all) * result$gamma, sigma = sigma.all,
        resOpt.scale = result$resOpt.scale)
    if (verbose) {
        print(paste("number of factors:", num.fac))
        if (!is.null(X)) {
            png("outlierplot.png")
            plot(X[, 1], Y)
            dev.off()
        }
    }
    return(RetList)
}

# #' Same function in LEAPP?
# #' 
# #' @keywords internal
# #' 
# AlternateSVD <- function (x, r, pred = NULL, max.iter = 10, TOL = 1e-04)
# {
#     N = nrow(x)
#     n = ncol(x)
#     if (!is.null(pred)) {
#         R = x %*% (diag(n) - pred %*% solve(t(pred) %*% pred) %*%
#             t(pred))
#         coef = x %*% pred %*% solve(t(pred) %*% pred)
#         x = R
#     }
#     else {
#         coef = NULL
#     }
#     if (r == 0) {
#         sigma = apply(x, 1, sd)
#         return(list(sigma = sigma, u = coef, v = NULL))
#     }
#     sigma = rep(1, N)
#     iter = 0
#     sigma.old = sigma + 1
#     while (1) {
#         if (iter > max.iter | sum(abs(sigma - sigma.old))/sum(abs(sigma.old)) <
#             TOL)
#             break
#         data = as.vector(1/sigma) * x
#         result.svd = fast.svd(data, tol = 0)
#         u = result.svd$u[, 1:r] %*% diag(result.svd$d[1:r], r,
#             r)
#         v = result.svd$v[, 1:r]
#         sigma.old = sigma
#         sigma = apply(x - as.vector(sigma) * u %*% t(v), 1, sd)
#         iter = iter + 1
#     }
#     uest = as.vector(sigma) * u
#     uest = cbind(uest)
#     return(list(sigma = sigma, coef = coef, uest = uest, vest = v))
# }


# #' Same function?
# #' 
# #' @keywords internal
# #' 
# IPODFUN <- function (X, Y, H, sigma, betaInit, method = "hard", TOL = 1e-04)
# {
#     N = length(Y)
#     gamma = matrix(rep(0, N), N, 1)
#     theta = sigma * sqrt(2 * log(N))
#     if (is.null(betaInit)) {
#         r = Y
#         if (method == "hard") {
#             gamma[abs(r) > theta] = r[abs(r) > theta]
#         }
#         else if (method == "soft") {
#             gamma[r > theta] = r[r > theta] - theta[r > theta]
#             gamma[r < -theta] = r[r < -theta] + theta[r < -theta]
#         }
#     }
#     else {
#         gamma.old = gamma
#         r0 = Y - H %*% Y
#         yEff = Y
#         niter = 0
#         theta = sigma * sqrt(2 * log(N)) * sqrt(1 - diag(H))
#         while (1) {
#             niter = niter + 1
#             beta0 = ginv(t(X) %*% X) %*% t(X) %*% yEff
#             if (niter == 1) {
#                 r = Y - X %*% betaInit
#             }
#             else {
#                 r = H %*% gamma.old + r0
#             }
#             if (method == "hard") {
#                 gamma = matrix(rep(0, N), N, 1)
#                 gamma[abs(r) > theta] = r[abs(r) > theta]
#             }
#             else if (method == "soft") {
#                 gamma = matrix(rep(0, N), N, 1)
#                 gamma[r > theta] = r[r > theta] - theta[r > theta]
#                 gamma[r < -theta] = r[r < -theta] + theta[r <
#                   -theta]
#             }
#             if (max(abs(gamma - gamma.old)) < TOL) {
#                 break
#             }
#             else {
#                 gamma.old = gamma
#             }
#             yEff = Y - gamma
#         }
#     }
#     RetList = list(gamma = gamma, ress = r)
# }
