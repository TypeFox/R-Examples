### <======================================================================>
"fit.hypmv" <- function(data,
                        opt.pars = c(alpha.bar = T, mu = T, sigma = T, gamma = !symmetric),
                        symmetric = F, ...)
{
    call <- match.call(expand.dots = TRUE)

    if(!is.null(list(...)$lambda)){
        stop("Do not submit lambda! Lambda is defined as (dimension+1)/2.\n")
    }
    lambda <- (min(dim(data))+1)/2

    ghyp.object <- fit.ghypmv(data = data, lambda = lambda,
                              opt.pars = c(lambda = F, opt.pars), ...)

    ghyp.object@parametrization <- "alpha.bar"
    ghyp.object@call <- call
    return(ghyp.object)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"fit.NIGmv" <- function(data,
                        opt.pars = c(alpha.bar = T, mu = T, sigma = T, gamma = !symmetric),
                        symmetric = F, ...)
{
    call <- match.call(expand.dots = TRUE)

    if(!is.null(list(...)$lambda)){
        stop("Do not submit lambda! Lambda is defined as -0.5.\n")
    }
    ghyp.object <- fit.ghypmv(data = data, lambda = -0.5,
                              opt.pars = c(lambda = F, opt.pars),...)

    ghyp.object@parametrization <- "alpha.bar"
    ghyp.object@call <- call
    return(ghyp.object)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"fit.tmv" <- function(data, nu = 3.5,
                      opt.pars = c(lambda = T, mu = T, sigma = T, gamma = !symmetric),
                      symmetric = F, ...)
{
    call <- match.call(expand.dots = TRUE)

    if(nu < 0){
        warning("nu < 0, Variance Gamma distribution is fitted instead!\n")
    }
    if(!is.null(list(...)$alpha.bar)){
        stop("Do not submit alpha.bar! alpha.bar is defined as 0.\n")
    }
    if(!is.null(list(...)$lambda)){
        stop("Do not submit lambda! Student-t distribution works with 'nu'.")
    }
    if("nu" %in% names(opt.pars)){
        names(opt.pars)[which(names(opt.pars)=="nu")] <- "lambda"
    }
    ghyp.object <- fit.ghypmv(data = data, lambda = -nu/2, alpha.bar = 0,
                              opt.pars = c(alpha.bar = F, opt.pars), ...)

    ghyp.object@parametrization <- "alpha.bar"
    ghyp.object@call <- call
    return(ghyp.object)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"fit.VGmv" <- function(data, lambda = 1,
                       opt.pars = c(lambda = T, mu = T, sigma = T, gamma = !symmetric),
                       symmetric = F, ...)
{
    call <- match.call(expand.dots = TRUE)

    if(lambda < 0){
        warning("lambda < 0, Student-t distribution is fitted instead!\n")
    }
    if(!is.null(list(...)$alpha.bar)){
        stop("Do not submit alpha.bar! alpha.bar is defined as 0.\n")
    }
    ghyp.object <- fit.ghypmv(data = data, lambda = lambda, alpha.bar = 0,
                              opt.pars = c(alpha.bar = F, opt.pars), ...)

    ghyp.object@parametrization <- "alpha.bar"
    ghyp.object@call <- call
    return(ghyp.object)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"fit.gaussmv" <- function(data, na.rm = T, save.data = T)
{
    call <- match.call()

    data.mat <- .check.data(data, case = "mv", na.rm = na.rm, fit = TRUE)

    mu <- colMeans(data.mat, na.rm = na.rm)
    sigma <- var(data.mat, na.rm = na.rm)

    d <- length(mu)

    if(!save.data){
        data.mat <- NULL
    }
    ghyp.object <- ghyp(lambda = as.numeric(NA),
                        chi = Inf,
                        psi = Inf,
                        alpha.bar = NULL,
                        mu = mu,
                        sigma = sigma,
                        gamma = rep(0, d),
                        data = data.mat)

    ghyp.object@parametrization <- "Gaussian"
    ghyp.object@call <- call

    llh <- sum(dghyp(data.mat, ghyp.object, logvalue = TRUE))

    aic <- -2 * llh + d/2 * (d + 1) + d

    return(.fit.ghyp(ghyp.object, llh = llh, n.iter = 0,
                    converged = TRUE,
                    error.code = 0, error.message = "",
                    parameter.variance = matrix(numeric(0)),
                    fitted.params = c(mu = TRUE, sigma = TRUE), aic = aic))
}
### <---------------------------------------------------------------------->



### <======================================================================>
"fit.hypuv" <- function(data,
                        opt.pars = c(alpha.bar = T, mu = T, sigma = T, gamma = !symmetric),
                        symmetric = F, ...)
{
    call <- match.call(expand.dots = TRUE)

    if(!is.null(list(...)$lambda)){
        stop("Do not submit lambda! Lambda is defined as 1.\n")
    }
    ghyp.object <- fit.ghypuv(data = data, lambda = 1,
                              opt.pars = c(lambda = F, opt.pars),
                              symmetric = symmetric, ...)

    ghyp.object@parametrization <- "alpha.bar"
    ghyp.object@call <- call
    return(ghyp.object)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"fit.NIGuv" <- function(data,
                        opt.pars = c(alpha.bar = T, mu = T, sigma = T, gamma = !symmetric),
                        symmetric = F, ...)
{
    call <- match.call(expand.dots = TRUE)

    if(!is.null(list(...)$lambda)){
        stop("Do not submit lambda! Lambda is defined as -0.5.\n")
    }
    ghyp.object <- fit.ghypuv(data = data, lambda = -0.5,
                              opt.pars = c(lambda = F, opt.pars),
                              symmetric = symmetric, ...)

    ghyp.object@parametrization <- "alpha.bar"
    ghyp.object@call <- call
    return(ghyp.object)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"fit.tuv" <- function(data, nu = 3.5,
                      opt.pars = c(nu = T, mu = T, sigma = T, gamma = !symmetric),
                      symmetric = F, ...)
{
    call <- match.call(expand.dots = TRUE)

    if(nu < 0){
        warning("nu < 0, Variance Gamma distribution is fitted instead!\n")
    }
    if(!is.null(list(...)$alpha.bar)){
        stop("Do not submit alpha.bar! alpha.bar is defined as 0.\n")
    }
    if(!is.null(list(...)$lambda)){
        stop("Do not submit lambda! Student-t distribution works with 'nu'.")
    }
    if("nu" %in% names(opt.pars)){
        names(opt.pars)[which(names(opt.pars)=="nu")] <- "lambda"
    }
    ghyp.object <- fit.ghypuv(data = data, lambda = -nu/2, alpha.bar = 0,
                              opt.pars = c(alpha.bar = F, opt.pars),
                              symmetric = symmetric, ...)

    ghyp.object@parametrization <- "alpha.bar"
    ghyp.object@call <- call
    return(ghyp.object)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"fit.VGuv" <- function(data, lambda = 1,
                       opt.pars = c(lambda = T, mu = T, sigma = T, gamma = !symmetric),
                       symmetric = F, ...)
{
    call <- match.call(expand.dots = TRUE)

    if(lambda < 0){
        warning("lambda < 0, Student-t distribution is fitted instead!\n")
    }
    if(!is.null(list(...)$alpha.bar)){
        stop("Do not submit alpha.bar! alpha.bar is defined as 0.\n")
    }
    ghyp.object <- fit.ghypuv(data = data, lambda = lambda, alpha.bar = 0,
                              opt.pars = c(alpha.bar = F, opt.pars),
                              symmetric = symmetric, ...)

    ghyp.object@parametrization <- "alpha.bar"
    ghyp.object@call <- call
    return(ghyp.object)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"fit.gaussuv" <- function(data, na.rm = T, save.data = T)
{
    call <- match.call()
    data.vec <- .check.data(data, case = "uv", na.rm = na.rm, fit = TRUE)

    mu <- unname(mean(data.vec, na.rm = na.rm))
    sigma <- unname(as.matrix(sd(data.vec, na.rm = na.rm)))

    if(!save.data){
        data.vec <- NULL
    }

    ghyp.object <- ghyp(lambda = as.numeric(NA),
                        chi = Inf,
                        psi = Inf,
                        alpha.bar = NULL,
                        mu = mu,
                        sigma = sigma,
                        gamma = rep(0, length(mu)),
                        data = data.vec)

    ghyp.object@parametrization <- "Gaussian"
    ghyp.object@call <- call

    llh <- sum(dghyp(data.vec, ghyp.object, logvalue = TRUE))
    aic <- -2 * llh + 4
    return(.fit.ghyp(ghyp.object, llh = llh, n.iter = 0,
                    converged = TRUE,
                    error.code = 0, error.message = "",
                    parameter.variance = matrix(c(1/length(data.vec) * sigma^2, 0, 0,
                    2 * sigma^4 / (length(data.vec) - 1)), nrow = 2),
                    fitted.params = c(mu = TRUE, sigma = TRUE), aic = aic))
}
### <---------------------------------------------------------------------->



### <======================================================================>
"hyp" <- function(chi = 0.5, psi = 2, mu = 0, sigma = diag(rep(1, length(mu))),
                  gamma = rep(0, length(mu)), alpha.bar = NULL, data = NULL)
{
    call <- match.call()

    if(is.null(alpha.bar)){
        parametrization <- "chi.psi"
    }else{
        parametrization <- "alpha.bar"
    }
    ghyp.object <- ghyp(lambda = (length(mu)+1)/2, chi = chi, psi = psi, mu = mu,
                        sigma = sigma, gamma = gamma, alpha.bar = alpha.bar, data = data)
    ghyp.object@call <- call
    ghyp.object@parametrization <- parametrization
    return(ghyp.object)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"NIG" <- function(chi = 2, psi = 2, mu = 0, sigma = diag(rep(1, length(mu))),
                  gamma = rep(0, length(mu)), alpha.bar = NULL, data = NULL)
{
    call <- match.call()

    if(is.null(alpha.bar)){
        parametrization <- "chi.psi"
    }else{
        parametrization <- "alpha.bar"
    }
    ghyp.object <- ghyp(lambda = -0.5, chi = chi, psi = psi, mu = mu, sigma = sigma,
                        gamma = gamma, alpha.bar = alpha.bar, data = data)
    ghyp.object@call <- call
    ghyp.object@parametrization <- parametrization
    return(ghyp.object)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"student.t" <- function(nu = 3.5, chi = nu - 2, mu = 0, sigma = diag(rep(1, length(mu))),
                        gamma = rep(0, length(mu)), data = NULL)
{
    call <- match.call()

    if(chi == nu - 2){
        parametrization <- "alpha.bar"
    }else{
        parametrization <- "chi.psi"
    }
    ghyp.object <-  ghyp(lambda = -nu/2, chi = chi, psi = 0, mu = mu, sigma = sigma,
                         gamma = gamma, alpha.bar = NULL, data = data)
    ghyp.object@call <- call
    ghyp.object@parametrization <- parametrization
    return(ghyp.object)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"VG" <- function(lambda = 1, psi = 2*lambda, mu = 0, sigma = diag(rep(1, length(mu))),
                 gamma = rep(0, length(mu)), data = NULL)
{
    call <- match.call()

    if(psi == 2 * lambda){
        parametrization <- "alpha.bar"
    }else{
        parametrization <- "chi.psi"
    }
    ghyp.object <- ghyp(lambda = lambda, chi = 0, psi = psi, mu = mu, sigma = sigma,
                        gamma = gamma, alpha.bar = NULL, data = data)
    ghyp.object@call <- call
    ghyp.object@parametrization <- parametrization
    return(ghyp.object)
}
### <---------------------------------------------------------------------->



### <======================================================================>
### <============== alpha.delta - parametrization =========================>

"ghyp.ad" <- function(lambda = 0.5, alpha = 1.5, delta = 1, beta = rep(0, length(mu)),
                      mu = 0, Delta = diag(rep(1, length(mu))), data = NULL)
{
    call <- match.call()

    ## Test validity of parameters
    if(delta < 0){
        stop("Parameter 'delta' must be >= 0!")
    }
    if(delta == 0 && lambda <= 0){
        stop("Parameter 'lambda' must be > 0 when 'delta' = 0!")
    }
    if(length(lambda) != 1 || length(alpha) != 1 || length(delta) != 1){
        stop("Parameters 'lambda', 'alpha' and 'delta' must be of length 1!")
    }
    if(length(beta) != length(mu)){
        stop("Length of 'beta' must be equal to length of 'mu'!")
    }
    if(length(mu) != dim(as.matrix(Delta))[1] || length(mu) != dim(as.matrix(Delta))[2]){
        stop("dim(Delta) must be length(mu) x length(mu)!")
    }
    chi <- delta^2
    if(length(mu) == 1){                #Univariate case
        psi <- alpha^2 - beta^2
        if(!missing(Delta)){
            warning("Parameter 'Delta' ignored in the univariate case!")
        }
        if(alpha^2 - beta^2 < 0){
            stop("alpha^2 - beta^2 must be >= 0!")
        }
        if(isTRUE(all.equal(psi, 0, tol = 1e-8)) & (lambda >= 0)){
            stop("Student-t Distribution: lambda must be < 0 if alpha^2 - beta' Delta beta = 0!")
        }
        gamma <- beta
        ghyp.object <- ghyp(chi = chi, psi = psi, lambda = lambda,
                            mu = mu, sigma = 1, gamma = gamma, data = data)
        ghyp.object@call <- call
        ghyp.object@parametrization <- "alpha.delta"
        return(ghyp.object)
    }else{                              # Multivariate case
        if(!isTRUE(all.equal(det(Delta), 1, tol = .Machine$double.eps^0.5))){
            stop("Determinant of 'Delta' must be 1!")
        }
        psi <- as.numeric(alpha^2 - beta %*% Delta %*% beta)
        if(abs(psi) < .Machine$double.eps^0.5){
            psi <- 0
        }
        if(psi < 0){
            stop("alpha^2 - beta' * Delta * beta must be >= 0!")
        }
        if(isTRUE(all.equal(psi, 0, tol = .Machine$double.eps^0.5)) & (lambda >= 0)){
            stop("Student-t Distribution: lambda must be < 0 if alpha^2 - beta' * Delta * beta = 0!")
        }
        sigma <- Delta
        gamma <- as.numeric(Delta %*% beta)
        ghyp.object <- ghyp(chi = chi, psi = psi, lambda = lambda,
                            mu = mu, sigma = sigma, gamma = gamma, data = data)
        ghyp.object@call <- call
        ghyp.object@parametrization <- "alpha.delta"
        return(ghyp.object)
    }
}
### <---------------------------------------------------------------------->



### <======================================================================>
"hyp.ad" <- function(alpha = 1.5, delta = 1, beta = rep(0, length(mu)), mu = 0,
                     Delta = diag(rep(1, length(mu))), data = NULL)
{
    call <- match.call()
    if(length(beta) == 1 && length(mu) == 1){ # Univariate case
        if(!missing(Delta)){
            warning("Parameter 'Delta' ignored in the univariate case!")
        }
        ghyp.object <- ghyp.ad(lambda = (length(mu) + 1)/2, alpha = alpha,
                               delta = delta, mu = mu, beta = beta, data = data)
    }else{                              # Multivariate case
        ghyp.object <- ghyp.ad(lambda = (length(mu) + 1)/2, alpha = alpha, delta = delta,
                               mu = mu, Delta = Delta, beta = beta, data = data)
    }
    ghyp.object@call <- call
    return(ghyp.object)

}
### <---------------------------------------------------------------------->



### <======================================================================>
"NIG.ad" <- function(alpha = 1.5, delta = 1, beta = rep(0, length(mu)), mu = 0,
                     Delta = diag(rep(1, length(mu))), data = NULL)
{
    call <- match.call()
    if(length(beta) == 1 && length(mu) == 1){ # Univariate case
        if(!missing(Delta)){
            warning("Parameter 'Delta' ignored in the univariate case!")
        }
        ghyp.object <- ghyp.ad(lambda = -0.5, alpha = alpha, delta = delta,
                               mu = mu, beta = beta, data = data)
    }else{                              # Multivariate case
        ghyp.object <- ghyp.ad(lambda = -0.5, alpha = alpha, delta = delta,
                               mu = mu, Delta = Delta, beta = beta, data = data)
    }
    ghyp.object@call <- call
    return(ghyp.object)
}
### <---------------------------------------------------------------------->



### <======================================================================>
"student.t.ad" <- function(lambda = -2, delta = 1, beta = rep(0, length(mu)), mu = 0,
                           Delta = diag(rep(1, length(mu))), data = NULL)
{
    call <- match.call()
    if(length(beta) == 1 && length(mu) == 1){ # Univariate case
        if(!missing(Delta)){
            warning("Parameter 'Delta' ignored in the univariate case!")
        }
        ghyp.object <- ghyp.ad(lambda = lambda, alpha = sqrt(as.numeric(beta %*% as.matrix(Delta) %*% beta)),
                               delta = delta, mu = mu, beta = beta, data = data)
    }else{                              # Multivariate case
        ghyp.object <- ghyp.ad(lambda = lambda, alpha = sqrt(as.numeric(beta %*% as.matrix(Delta) %*% beta)),
                               delta = delta, beta = beta, mu = mu, Delta = Delta, data = data)
    }
    ghyp.object@call <- call
    return(ghyp.object)

}
### <---------------------------------------------------------------------->



### <======================================================================>
"VG.ad" <- function(lambda = 2, alpha = 1.5, beta = rep(0, length(mu)), mu = 0,
                    Delta = diag(rep(1, length(mu))), data = NULL)
{
    call <- match.call()
    if(length(beta) == 1 && length(mu) == 1){ # Univariate case
        if(!missing(Delta)){
            warning("Parameter 'Delta' ignored in the univariate case!")
        }
        ghyp.object <- ghyp.ad(lambda = lambda, alpha = alpha, delta = 0,
                               mu = mu, beta = beta, data = data)
    }else{                              # Multivariate case
        ghyp.object <- ghyp.ad(lambda = lambda, alpha = alpha, delta = 0,
                               mu = mu, Delta = Delta, beta = beta, data = data)
    }
    ghyp.object@call <- call
    return(ghyp.object)

}
### <---------------------------------------------------------------------->



### <======================================================================>
"gauss" <- function(mu = 0, sigma = diag(rep(1, length(mu))), data = NULL)
{
    call <- match.call()
    ghyp.object <- ghyp(chi = Inf, psi = Inf, lambda = as.numeric(NA),
                        mu = mu, sigma = sigma, gamma = rep(0, length(mu)), data = data)
    ghyp.object@call <- call
    ghyp.object@parametrization <- "Gaussian"
    return(ghyp.object)

}
### <---------------------------------------------------------------------->
