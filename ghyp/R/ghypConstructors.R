### <======================================================================>
### chi.psi and alpha.bar - parametrization
### <======================================================================>
"ghyp" <- function(lambda = 0.5, chi = 0.5, psi = 2, mu = 0, sigma = diag(rep(1, length(mu))),
                   gamma = rep(0, length(mu)), alpha.bar = NULL, data = NULL)
{
    call <- match.call()
    if(!is.null(alpha.bar)){
        tmp.abar2chipsi <- .abar2chipsi(alpha.bar,lambda)
        chi  <- tmp.abar2chipsi$chi
        psi  <- tmp.abar2chipsi$psi
        names(chi) <- "chi"
        names(psi) <- "psi"
        parametrization <- "alpha.bar"
    }else{
        alpha.bar <- sqrt(chi*psi)
        names(alpha.bar) <- "alpha.bar"
        parametrization <- "chi.psi"
    }

    .check.gig.pars(lambda, chi, psi)

    if(chi == Inf && psi == Inf){       # Gaussian case
        e.gig <- 1
        parametrization <- "Gaussian"
    }else{
        e.gig <- Egig(lambda, chi, psi, func = "x")
    }

    if(!is.numeric(gamma)){
        stop("Parameter 'gamma' must be numeric!\n")
    }
    if(!is.numeric(sigma)){
        stop("Parameter 'sigma' must be numeric!\n")
    }
    if(!is.numeric(mu)){
        stop("Parameter 'mu' must be numeric!\n")
    }

    model <- .ghyp.model(lambda, chi, psi, gamma)
    if(length(mu)==1){
        ##  univariate case

        mu <- as.vector(mu) # In case 'mu' is a 1x1 matrix
        gamma <- as.vector(gamma) # In case 'gamma' is a matrix
        
        if(length(gamma) != 1){
            stop("Parameter 'gamma' must have the same length as 'mu' !\n")
        }
        if(length(as.vector(sigma)) != 1){
            stop("Parameter 'sigma' must be scalar when 'mu' and 'gamma' are scalar!\n")
        }

        if(!is.null(data)){
            data <- .check.data(data = data, case = "uv", na.rm = FALSE, fit = FALSE)
        }else{
            data <- numeric(0)
        }

        mu <- unname(mu)
        sigma <- unname(sigma)

        if(all(gamma == 0)){
            if(is.finite(e.gig)){
                variance <- unname(e.gig * sigma^2)
            }else{
                variance <- Inf
            }
        }else{
            var.gig <- Egig(lambda, chi, psi, func = "var")
            if(is.finite(var.gig)){
                variance <- unname(var.gig * gamma^2 + e.gig * sigma^2)
            }else{
                variance <- Inf
            }
        }

    }else if(length(mu)>1){
        ## multivariate case
        if(ncol(sigma)!=length(mu)){
            stop("Dimension mismatch ( ncol(sigma)!=length(mu) )!\n")
        }
        if(length(gamma)!=length(mu)){
            stop("Dimension mismatch ( length(gamma)!=length(mu) )!\n")
        }

        mu <- as.vector(mu) # In case 'mu' is a matrix
        gamma <- as.vector(gamma) # In case 'gamma' is a matrix
        
        if(!is.null(data)){
            data <- .check.data(data = data, case = "mv", na.rm = FALSE, fit = FALSE)
        }else{
            data <- numeric(0)
        }

        if(all(gamma == 0)){
            if(is.finite(e.gig)){
                variance <- e.gig * sigma
            }else{
                variance <- matrix(Inf, ncol = length(mu), nrow = length(mu))
            }
        }else{
            var.gig <- Egig(lambda, chi, psi, func = "var")
            if(is.finite(var.gig)){
                variance <- var.gig * outer(gamma, gamma) + e.gig * sigma
            }else{
                variance <- matrix(Inf, ncol = length(mu), nrow = length(mu))
            }
        }
        dimnames(variance) <- dimnames(sigma)

    }else{
        stop("Invalid parameters!\n")
    }
    expected.value <- unname(mu + e.gig * gamma)

    if(length(mu) > 1){
        names(expected.value) <- names(mu)
    }
    return(new("ghyp", call = call,
               lambda = unname(lambda),
               chi = unname(chi),
               psi = unname(psi),
               alpha.bar = unname(alpha.bar),
               mu = mu,
               sigma = as.matrix(sigma),
               gamma= unname(gamma),
               model = model,
               dimension = length(mu),
               expected.value = expected.value,
               variance = as.matrix(variance),
               data = as.matrix(data),
               parametrization = parametrization))
}
### <---------------------------------------------------------------------->
