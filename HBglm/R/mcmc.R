###############################################################################
# Function intended to be called by hbglm() for MCMC sampling
#
# Args:
#  model, family - internal data structure of hbglm()
#  bounds        - constraints    (setup by get.box.bounds() in file model.R)
#  init.vals     - initial values (setup by get.init.vals() in file model.R)
#  nsamp         - number of sampling iterations
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
mcmc <- function(model, family, bounds, init.vals, nsamp, print.level=0,
                 report.freq = 10)
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
               
    # Store intial value as first sample
    samples$beta[ , 1] <- as.vector(beta)        
    if (family$has.tau) samples$tau[ , 1] <- as.vector(tau)        
    if (model$has.fixed) samples$alpha[ , 1] <- as.vector(alpha)
    if (model$has.upper.level) {
        samples$theta[ , 1] <- as.vector(theta)
        samples$Sigma[ , 1] <- Sigma
    }       

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
    if (print.level) cat(paste0("\nCommencing MCMC sampling ... "))
    for (itr in 2:nsamp) {
        if (print.level && report.on && itr %% report.freq == 0) 
            cat(paste0("\n\tMCMC Iteration# ", itr, " / ", nsamp))
        # Gibbs sampling sequence: beta, tau, alpha, theta, Sigma

        # Sample new beta matrix
        t_init <- proc.time()[3]
        for (j in 1:model$J) {  # loop over independent pools 
          logp.betaj <- function(betaj) {
            return(lp.yj_betaj.alpha.tauj(model, family, j, betaj, alpha,
                   tau[j]) + lp.betaj_theta.Sigma(model, family, j, betaj, 
                                                  theta, Sigma))
          }
          beta[j, ] <- MfU.Sample(beta[j, ], f = logp.betaj, 
            uni.sampler = "slice", control = MfU.Control(n = K, slice.m=10,
            slice.lower=bounds$beta.lo[j, ], slice.upper=bounds$beta.hi[j, ]))
          #beta[j, ] <- unscale.beta(beta[j, ], j , model)
        }  # end beta sampling 
        t_beta <- t_beta + proc.time()[3] - t_init

        # Sample alpha vector
        if (model$has.fixed) {
          t_init <- proc.time()[3]
          logp.alpha <- function(alpha) 
            return(lp.y_beta.alpha.tau(model, family, beta, alpha, tau))
          alpha <- MfU.Sample(alpha, f = logp.alpha, uni.sampler = "slice",
            control = MfU.Control(n = M, slice.m = 10, slice.lower = 
                bounds$alpha.lo, slice.upper=bounds$alpha.hi))
          t_alpha <- t_alpha + proc.time()[3] - t_init
        }
        
        # Sample tau vector
        if (family$has.tau) {
          t_init <- proc.time()[3]
          logp.tauj <- function(tauj) {
            return(lp.yj_betaj.alpha.tauj(model, family, j, beta[j,], alpha, 
                   tauj)) #+ lp.tauj_priors(tauj))
          }
          for (j in 1:J) {
            tau[j] <- uni.slice(tau[j], logp.tauj, m = 10, 
                                lower = bounds$tau.lo[j], 
                                upper = bounds$tau.hi[j])
          }
          #if (tau[j] < 0) browser()
          t_tau <- t_tau + proc.time()[3] - t_init
        }

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

        #cat(paste0("\nTotal time(s): ", round(t_tot, 2), "\n"))
        cat(paste0("\nTotal sampling time(s): ", round(t_gibbs, 2), "\n"))
        cat("Sampling time breakup:\n")
        print(t_tab)
    }
    return(samples)
}
###############################################################################

# Log-posterior probability calculation

# log P(beta[j, ] | theta, Sigma)
lp.betaj_theta.Sigma <- function(model, family, j, betaj, theta, Sigma)
{
    Z <- model$mat.upper
    Zj.th <- if (!model$upper.reg) as.vector(theta) else
                                   as.vector(Z[j, ] %*% theta)
    return(-0.5 * sum(Zj.th * (solve(Sigma) %*% Zj.th)))
}

# Likelihood function - evaluated for group 'j'
# log P(y[j] | betaj, alpha, tauj) (used in beta sampling) 
lp.yj_betaj.alpha.tauj <- function(model, family, j, betaj, alpha, tauj)
{
    eta <- model$mat.rand.split[[j]]$X %*% betaj
    if (model$has.fixed) {
        grpj <- which(model$grp.indx == model$grp.labels[j])
        eta <- eta + model$mat.fixed[grpj, ] %*% alpha
    }
    response <- model$mat.rand.split[[j]]$y
    if (is.null(tauj)) return(family$loglik(eta, response))
    else return(family$loglik(eta, response, var = tauj))
}

# Likelihood function - evaluated over all groups
# log P(y[j] | beta, alpha, tau)  (used when we have fixed effects)
lp.y_beta.alpha.tau <- function(model, family, beta, alpha, tau, ncores = 1)
{
    eta <- model$mat.fixed %*% alpha
    logp <- 0
    for (j in 1:model$J) {  # loop over groups
        eta.j <- eta[which(model$grp.indx == model$grp.labels[j])] +
                 model$mat.rand.split[[j]]$X %*% beta[j, ]
        response <- model$mat.rand.split[[j]]$y
        logp <- logp + if (is.null(tau)) 
                           family$loglik(eta.j, response) else
                           family$loglik(eta.j, response, var = tau[j])
    }
    return(logp)
}

# log P(tau[j] | priors)
lp.tauj_priors <- function(tau.j) 
{
    dinvgamma(tau.j, shape = 0.001, scale = 0.001, log = TRUE)
}
###############################################################################

unscale.beta <- function(beta, grpid, model)
{
    if (!model$scaled.data) return(beta)
    j <- grpid
    beta <- beta / model$mat.rand.split[[j]]$scale.facts[, 1]
    icol <- model$mat.rand.split[[j]]$intercept
    if (!is.na(icol)) {
        beta[icol] <- beta[icol] -
                      sum(beta * model$mat.rand.split[[j]]$scale.facts[, 2])
    }
    return(beta)
}
###############################################################################
