# lms.R
#
# created: Sep/11/2014, NU
# last mod: Aug/20/2015, NU

#--------------- main functions ---------------

# Calculate mu of multivariate normal distribution for joint vector of
# indicators (see Equations 16 and 18, 19 in Klein & Moosbrugger, 2000)
mu_lms <- function(model, z) {
    
    #check_filled(model)
    # TODO: Remove?

    matrices <- model$matrices$class1
    k    <- get_k(matrices$Omega)     # number of nonzero rows in Omega
    n    <- nrow(matrices$A)          # number of zero rows in Omega
    if (k < n & k > 0) {
        z.1 <- c(z, rep(0, n - k))   # [z_1 0]'
    } else if (k == 0) {
        z.1 <- rep(0, n - k)
    } else z.1 <- z
            # added quadratic effects, not in original equations
    A.z  <- matrices$A %*% z.1 
    mu.x <- matrices$nu.x + matrices$Lambda.x %*% (matrices$tau + A.z )
            # added tau's, not in the original equations
    mu.y <- matrices$nu.y + matrices$Lambda.y %*% (matrices$alpha +
            matrices$Gamma %*% (matrices$tau + A.z) + t(matrices$tau + A.z)
            %*% matrices$Omega %*% (matrices$tau + A.z))
    mu   <- c(mu.x, mu.y)

    mu
}

# Calculate Sigma of multivariate normal distribution for joint vector of
# indicators (see Equations 17 and 20-22 in Klein & Moosbrugger, 2000)
sigma_lms <- function(model, z) {

    #check_filled(model)

    matrices <- model$matrices$class1
    k     <- get_k(matrices$Omega)    # number of nonzero rows in Omega
    n     <- nrow(matrices$A)         # number of zero rows in Omega
    if (k < n & k > 0) {
        z.1  <- c(z, rep(0, n - k))   # [z_1 0]'
    } else if (k == 0) {
        z.1 <- rep(0, n - k)
    } else z.1 <- z
    A.z   <- matrices$A %*% z.1 
    d.mat <- get_d(n=n, k=k)
    Lx.A  <- matrices$Lambda.x %*% matrices$A
    temp  <- matrices$Gamma %*% matrices$A + t(matrices$tau + A.z) %*% matrices$Omega %*% matrices$A
            # added tau's, not in original equations
    s11   <- Lx.A %*% d.mat %*% t(Lx.A) + matrices$Theta.d 
    s12   <- Lx.A %*% d.mat %*% t(temp) %*% t(matrices$Lambda.y)
    s21   <- t(s12)
    s22   <- matrices$Lambda.y %*% temp %*% d.mat %*% t(temp) %*%
             t(matrices$Lambda.y) + matrices$Lambda.y %*% matrices$Psi %*%
             t(matrices$Lambda.y) + matrices$Theta.e
    sigma <- rbind(cbind(s11,s12), cbind(s21,s22))
    
    # check if sigma is symmetric
    if(!isSymmetric(sigma)) stop("Sigma has to be symmetric.")

    sigma
}

# Expectation step of EM-algorithm (see Klein & Moosbrugger, 2000)
estep_lms <- function(model, parameters, dat, m, ...) {

    stopifnot(count_free_parameters(model) == length(parameters))

    mod.filled <- fill_model(model=model, parameters=parameters)

    k <- get_k(mod.filled$matrices$class1$Omega)
    if (k != 0){
    quad <- quadrature(m, k)
    V <- quad$n       # matrix of node vectors m x k
    w <- quad$w       # weights
    } else {
        V <- 0
        w <- 1
        # do not need mixtures, if I do not have interactions
    }

    stopifnot(sum(w) - 1 < 1e-5)

    P <- NULL
    for(i in seq_along(w)) {
        p.ij <- w[i] * dmvnorm(dat, mean=mu_lms(model=mod.filled, z=V[i,]), 
                               sigma=sigma_lms(model=mod.filled, z=V[i,]))
        P    <- cbind(P, p.ij, deparse.level=0)
    }

    P <- P / rowSums(P)   # divide each rho_j*phi(x_i, y_i) by whole density (row)

    P
}

# log likelihood function which will be optimized in M-step (see below)
loglikelihood_lms <- function(parameters, model, dat, P, m=16, ...) {
    
    mod.filled <- fill_model(model=model, parameters=parameters)

    k <- get_k(mod.filled$matrices$class1$Omega)
    quad <- quadrature(m, k)
    V <- quad$n
    if (k == 0) V <- as.data.frame(0)

    res0 <- sapply(seq_len(nrow(V)), function(i){
              lls <- sum(dmvnorm(dat,
                            mean=mu_lms(model=mod.filled, z=V[i,]),
                            sigma=sigma_lms(model=mod.filled, z=V[i,]),
                            log=TRUE) * P[,i])
              lls
            })
    res <- sum(res0)

    -res
}

# Maximization step of EM-algorithm (see Klein & Moosbrugger, 2000)
mstep_lms <- function(parameters, model, dat, P, m, neg.hessian=FALSE,
                      optimizer=c("nlminb", "optim"), max.mstep,
                      control=list(), ...) {

    # optimizer
    optimizer <- match.arg(optimizer)

    if (optimizer == "nlminb") {

        if (is.null(control$iter.max)) {
            control$iter.max <- max.mstep
        } else warning("iter.max is set for nlminb. max.mstep will be ignored.")

        suppress_NaN_warnings(
            # See semm.R helper function
            est <- nlminb(start=parameters, objective=loglikelihood_lms, dat=dat,
                          model=model, P=P, upper=model$info$bounds$upper,
                          lower=model$info$bounds$lower, control=control,
                          ...)
        )
    } else {

        if (is.null(control$maxit)){
            control$maxit <- max.mstep
        } else warning("maxit is set for optim. max.mstep will be ignored.")

        est <- optim(par=parameters, fn=loglikelihood_lms, model=model, dat=dat,
                     P=P, upper=model$info$bounds$upper,
                     lower=model$info$bounds$lower, method="L-BFGS-B",
                     control=control, ...)
        # fit est to nlminb output
        names(est) <- gsub("value", "objective", names(est))
    }
    if (neg.hessian == TRUE){
    est$hessian <- fdHess(pars=est$par, fun=loglikelihood_lms,
                          model=model, dat=dat, P=P)$Hessian
  
    }


    est
}

#--------------- helper functions ---------------

# Count number of rows that contain at least one element which is not 0 in
# Omega
get_k <- function(Omega) {
    if (any(is.na(Omega))){
        out <- length(which(is.na(rowSums(Omega))))
    } else {
        if (any(rowSums(Omega) == 0)){
            out <- which(rowSums(Omega) == 0)[1] - 1
        } else out <- nrow(Omega)
    }
    out
}

# Get matrix D which contains of 0's and has the identity matrix in its
# lower right half (compare Equations 20-22 in Klein & Moosbrugger, 2000)
get_d <- function(n, k) {
    mat <- diag(n)
    mat[1:k, 1:k] <- 0
    mat
}

# Calculate weights and node points for mixture functions via Gauss-Hermite
# quadrature as defined in Klein & Moosbrugger (2000)
quadrature <- function(m, k) {

    one.dim         <- hermite.h.quadrature.rules(m)[[m]]
    test            <- as.matrix(expand.grid(lapply(vector("list", k), function(x) {x <- 1:m; x})))
    final.nodes     <- matrix(one.dim$x[test], ncol=k, byrow=FALSE)
    permute.weights <- matrix(one.dim$w[test], ncol=k, byrow=FALSE)
    final.weights   <- apply(permute.weights, 1, prod)

    n               <- final.nodes * sqrt(2)
    w               <- final.weights * pi^(-k/2)
  
    out <- list(n=n, w=w, k=k, m=m)
    out
}

# Convert parameters for Phi in LMS model to A 
convert_parameters_singleClass <- function(model, parameters) {

  names(parameters) <- model$info$par.names
  Phi <- matrix(0, nrow=model$info$num.xi, ncol=model$info$num.xi)
  Phi[lower.tri(Phi, diag=TRUE)] <- parameters[grep("Phi", names(parameters))]
  Phi <- fill_symmetric(Phi)
  A <- tryCatch({ t(chol(Phi)) }, 
    error=function(e) {
      warning("Starting parameters for Phi are not positive definite. Identity matrix was used instead.")
      diag(1, model$info$num.xi)}
  )
  parameters[grep("Phi", names(parameters))] <- A[lower.tri(A, diag=TRUE)]

  parameters
}
