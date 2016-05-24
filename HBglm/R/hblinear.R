###############################################################################
# Main dispatcher function for linear HB models
###############################################################################
hblinear <- function(model, family, bounds, init.vals, nsamples,
                           print.level=0, report.freq = 10)
{
    standard.bounds <- is.null(bounds$user.const)
    if (!model$has.fixed && is.null(bounds$user.const)) #NO fixed eff or bounds
      return(hblinear.bayesm(model, family, nsamples, print.level=print.level))
    else {  
      #return(hblinear.gibbs(model, family, bounds, init.vals, nsamples,
      return(mcmc(model, family, bounds, init.vals, nsamples,
          print.level = print.level, report.freq = report.freq))     
    }
}

###############################################################################
# 2-level HB estimation using rhierLinearModel() from bayesm
# Doesn't handle fixed effects or constraints
###############################################################################

hblinear.bayesm <- function(model, family, nsamp, print.level = 0)
{
    if (model$has.fixed)
        stop("Incorrect function called for fixed effects linear HB model.")
    if (print.level) 
        cat("Calling rhierLinearModel() from bayesm ... ")
    t_start <- proc.time()[3]
    Data <- list(Z = model$mat.upper, regdata = model$mat.rand.split)
    Mcmc <- list(R = nsamp, keep = 1)
    out <- rhierLinearModel(Data = Data, Mcmc = Mcmc)
    samples <- list(beta  = matrix(as.vector(out$betadraw), ncol = nsamp),
                    tau   = matrix(t(out$taudraw), ncol = nsamp),
                    theta = matrix(t(out$Deltadraw), ncol = nsamp),
                    Sigma = matrix(t(out$Vbetadraw), ncol = nsamp))
    t_samp <- proc.time()[3] - t_start
    if (print.level) 
        cat(paste0("Time taken: ", round(t_samp, 3), " seconds.\n"))
    return(samples) 
}

###############################################################################
# 2-level HB estimation handles bounds on random & fixed effects
#
# Args:
#  model, family - internal data structure of hbglm()
#  bounds        - constraints    (setup by get.box.bounds() in file model.R)
#  init.vals     - initial values (setup by get.init.vals() in file model.R)
#  nsamples      - number of sampling iterations
#  print.level   - set to greater than 0 to print iteration progress
#  report.freq   - suggested number of iterations after which to report
#                  (a value < 1 will disable reports)
#
# Returns:
#  Samples stored in a single matrix organized as follows:
#    One row per iteration. Each row has:
#      J * K slots for beta,
#      L * K slots for theta (L if no upper level covariates)
# CHECK: line above should be K if no upper level covars
#      J, K slots for tau, sigma (respectively)
#      2 slots (the last 2) for phi (only if tau is present)
#    Missing variables won't have slots - their place is taken by next in list
#  The matrix has attribute "time" giving time taken in sampling in seconds.
###############################################################################
hblinear.gibbs <- function(model, family, bounds, init.vals, nsamp,
                           print.level=0, report.freq = 10)
{
    t_start <- proc.time()[3]
    report.on <- ifelse(report.freq <= 0, FALSE, TRUE)
    report.freq <- min(report.freq, nsamp)

    # Model parameter counts
    N <- model$N
    J <- model$J
    K <- model$K
    M <- model$M
    L <- model$L

    # Allocate space to store MCMC samples
    # beta stored as array of dim: (J*K, nsamp)
    # tau stored as array of dim: (J, nsamp)
    # alpha stored as array of dim: (M, nsamp)
    # theta stored as array of dim: (L*K, nsamp)
    # Sigma stored as array  of dim: (K*K, nsamp)
    samples <- list()
    samples$beta <- matrix(rep(NA, J * K * nsamp), ncol = nsamp)
    samples$tau <- if (family$has.tau) matrix(rep(NA, J * nsamp),
                                              ncol = nsamp) else NULL
    samples$alpha <- if (model$has.fixed) matrix(rep(NA, M * nsamp),
                                                 ncol = nsamp) else NULL
    if (model$has.upper.level) {
      samples$theta <- matrix(rep(NA, L*K*nsamp), ncol = nsamp)
      samples$Sigma <- matrix(rep(NA, K*K*nsamp), ncol = nsamp)
    } else {
      samples$theta <- NULL
      samples$Sigma <- NULL
    }

    # Initialize
    beta   <- init.vals$beta
    alpha  <- init.vals$alpha
    tau    <- init.vals$tau
    theta  <- init.vals$theta
    Sigma  <- init.vals$Sigma

###TODO   
    # Precomputations to speed up
    XtX.list <- lapply(1:model$J, function(j) {
        X <- model$mat.lower.split[[j]]$X
        crossprod(X, X)
    })
    Xty.list <- lapply(1:model$J, function(j) {
        X <- model$mat.lower.split[[j]]$X
        y <- model$mat.lower.split[[j]]$y
        crossprod(X, y)
    })
    gamma = 1000
    ZA = diag(model$L) / gamma
    ZtZ <- crossprod(model$mat.upper, model$mat.upper)
    ZR  <- backsolve(chol(ZtZ + ZA), diag(model$L)) 

###TODO  
    # Store intial value as first sample
    samples$beta[ , 1] <- as.vector(beta)
    if (family$has.tau) samples$tau[ , 1] <- as.vector(tau)
    if (model$has.fixed) samples$alpha[ , 1] <- as.vector(alpha)
    if (model$has.upper.level) {
        samples$theta[ , 1] <- as.vector(theta)
        samples$Sigma[ , 1] <- Sigma
    }
    sigma <- Sigma
    nsamples <- nsamp

    # Prior parameters
    theta.bar <- matrix(rep(0, L * K), ncol = K)
    mat.A <- 0.01 * diag(L)
    nu <- 3 + K
    mat.V <- nu * diag(K)

    # Timing stuff
    t_beta  <- 0
    t_alpha <- 0
    t_tau   <- 0
    t_upper <- 0
    t0 <- proc.time()[3]
    t_prep <- t0 - t_start
 
    # Gibbs' sampling iterations
    if (print.level) cat(paste0("\nCommencing Gibbs sampling ... "))
    for (itr in 2:nsamples) {
        if (print.level && report.on && itr %% report.freq == 0) 
            cat(paste0("\n\tGibbs Iteration# ", itr, " / ", nsamples))

####TODO
        # beta sampling
        t_init <- proc.time()[3]
        beta.bar <- model$mat.upper %*% theta
        for (j in 1:model$J) {  # beta sampling
#cat(paste0("beta j = ", j, "\n"))
            R <- backsolve(chol(XtX.list[[j]] + diag(1/sigma)), 
                           diag(length(sigma)))
            beta.tilde <- R %*% crossprod(R, (Xty.list[[j]] + 
                                              beta.bar[j, ] / sigma))
            beta.samp <- beta.tilde + (R %*% rnorm(ncol(R))) * sqrt(tau[j])
            beta[j, ] <- unscale.beta(beta.samp, j, model) 
        }   # end beta sampling 
t_fin <- proc.time()[3]
t_beta <- t_beta + t_fin - t_init

t_init <- proc.time()[3]
pool.tau = F
if (F) {
        for (j in 1:model$J) {  # tau sampling
            r_j <- model$mat.lower.split[[j]]$y
                   - model$mat.lower.split[[j]]$X %*% beta[j, ]
            N_j <- length(r_j)
            shape <- N_j/2 + 0.001 #ifelse(pool.tau, phi[1], 0.001)
            scale <- sum(r_j * r_j)/2 + 0.001 # ifelse(pool.tau, phi[2], 0.001)
#browser()
            tau[j] <- rinvgamma(1, shape, scale)
        }  # end tau sampling
} else {
            for (j in 1:model$J) {  # tau sampling
              logp.tau.eval <- function(x) {
                tau[j] <- x
                lp <- lp.tauj_phi(tau[j]) +
                      lp.yj_beta.tau(model, family, beta, j, tauVec=tau) 
                return(lp)
              } 
              tau[j] <- uni.slice(tau[j], logp.tau.eval, m=10,
                                  lower = bounds$tau.lo[j], 
                                  upper = bounds$tau.hi[j])
            }  # end tau sampling
}
t_fin <- proc.time()[3]
t_tau <- t_tau + t_fin - t_init
####TODO

        # theta, Sigma sampling
        if (model$has.upper.level) {
          t_init <- proc.time()[3]
          out <- rmultireg(beta, model$mat.upper, theta.bar, mat.A, nu, mat.V)
          theta <- out$B
          Sigma <- out$Sigma
          t_upper <- t_upper + proc.time()[3] - t_init
        }

        # Store sample & prepare for next iteration
        samples$beta[ , itr] <- as.vector(beta)
        if (family$has.tau) samples$tau[ , itr] <- as.vector(tau)
        if (model$has.fixed) samples$alpha[ , itr] <- as.vector(alpha)
        if (model$has.upper.level) {
            samples$theta[ , itr] <- as.vector(theta)
            samples$Sigma[ , itr] <- Sigma
        }

    } # end Gibb's sampling 
    t_gibbs <- proc.time()[3] - t0
    t_tot <- proc.time()[3] - t_start

    if (print.level >= 1) {
        tvec <- c(t_beta, t_alpha, t_tau, t_upper)
        t_tab <- cbind(round(tvec, 2), round(100 * tvec /t_gibbs, 1))
        rownames(t_tab) <- c("beta", "alpha", "tau", "upper")
        colnames(t_tab) <- c("time(s)", "sample %")

        cat(paste0("\nTotal time(s): ", round(t_tot, 2), "\n"))
        cat(paste0("Preprocess time(s): ", round(t_prep, 2), "\n"))
        cat(paste0("\nTotal Gibbs time(s): ", round(t_gibbs, 2), "\n"))
        cat("Sampling time breakup:\n")
        print(t_tab)
    }

    return(samples)
}
###############################################################################

###############################################################################
# Log-posterior probability calculation

# log P(beta[ , k] | theta[ , k], sigmaSq[k])
lp.betak_theta.sigmaSq <- function(model, sigSqVec, thetaMat, betaMat, k)
{
  mean <- if (model$upper.reg) model$mat.upper %*% thetaMat[ , k]
          else thetaMat[ , k]     
  return(ll_normal(mean, betaMat[ , k], var = sigSqVec[k]))
}

# log P(beta[j, ] | theta, sigSq)
lp.betaj_theta.sigmaSq <- function(model, sigSqVec, thetaMat, betaMat, j)
{
    Z <- model$mat.upper
    Z.dot.theta <- if (!model$upper.reg) as.vector(thetaMat)
                   else                  as.vector(Z[j, ] %*% thetaMat)
    return(sum(sapply(c(1:model$K), 
      function(k) ll_normal(Z.dot.theta[k], betaMat[j, k], var=sigSqVec[k]))))
}

# log P(y(j) | beta[j, ], tau)  (the likelihood)
lp.yj_beta.tau <- function(model, family, betaMat, j, tauVec=NULL)
{
    grpj     <- which(model$grp.indx == model$grp.labels[j])
    X        <- model$mat.lower.split[[j]]$X 
    Xbeta    <- X %*% betaMat[j, ]
    response <- model$mat.lower.split[[j]]$y
    if (is.null(tauVec)) return(family$loglik(Xbeta, response))
    else return(family$loglik(Xbeta, response, var = tauVec[j]))
}

# log P(tau[j] | phi)
lp.tauj_phi <- function(tau.j)
{
    dinvgamma(tau.j, shape = 0.001, scale = 0.001, log = TRUE)
    #dinvgamma(tau.j, shape = phi[1], scale = phi[2], log = TRUE)
}
###############################################################################
