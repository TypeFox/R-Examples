##############################################################################
# Functions for post-processing
##############################################################################

##############################################################################
# Function to compute summary sample statistics
# Assumes samples are stored in a matrix. One sample per COLUMN.
# 'ignore' is the number of columns to ignore in stats calculation
get.stats <- function(samp.mat, ignore = 0)
{
    if(ignore >= ncol(samp.mat)) 
        stop("Invalid 'ignore' value given to get.stats() in summary.R")
    ignore <- ifelse(ignore <= 0, 0, ignore)
    # Discard 'ignore' columns from samp.mat before computations
    sample.mat <- samp.mat[ , (1+ignore) : ncol(samp.mat)]

    # Compute sample statistics
    stats <- cbind(apply(sample.mat, 1, mean),
                   apply(sample.mat, 1, sd),
                   round(ess(t(sample.mat))),
                   t(apply(sample.mat, 1, quantile, 
                           probs=c(0.025, 0.5, 0.975))))
    colnames(stats) <- c("mean", "sd", "ess", "2.5%", "50%", "97.5%")
    return(stats)
}

# Compute hbglm sample statistics
sample.stats <- function(samples, nburn = 0, var.names = NULL)
{
    stats <- list(beta = get.stats(samples$beta, ignore = nburn)) 
    if (!is.null(samples$tau)) stats$tau<- get.stats(samples$tau, ignore=nburn)
    if (!is.null(samples$alpha)) 
        stats$alpha <- get.stats(samples$alpha, ignore=nburn)
    if (!is.null(samples$theta))
        stats$theta <- get.stats(samples$theta, ignore=nburn)
    if (!is.null(samples$Sigma))
        stats$Sigma <- get.stats(samples$Sigma, ignore = nburn)
   
    if (!is.null(var.names))  {
        rownames(stats$beta) <- var.names$beta
        if (!is.null(stats$alpha)) rownames(stats$alpha) <- var.names$alpha 
        if (!is.null(stats$tau)) rownames(stats$tau) <- var.names$tau
        if (!is.null(stats$theta)) {
            rownames(stats$theta) <- var.names$theta
            rownames(stats$Sigma) <- var.names$Sigma
        } 
    }
    return(stats)
}

##############################################################################
# Make variable names
make.var.names <- function(model, family)
{
    beta.names <- as.vector(outer(model$grp.labels, model$rand.cov,
                            paste, sep=":"))
    alpha.names <- if(model$has.fixed) model$fixed.cov else NULL
    theta.names <- if(model$has.upper.level) {
        if(model$L == 1) model$rand.cov else
            as.vector(outer(model$upper.cov, model$rand.cov,
                            paste, sep=":"))
    } else NULL
    tau.names <- if(family$has.tau) model$grp.labels else NULL
    Sigma.names <- if(!model$has.upper.level) NULL else 
        as.vector(outer(model$rand.cov, model$rand.cov, paste, sep=":"))
    
    return(list(
        beta  = beta.names,
        theta = theta.names,
        alpha = alpha.names,
        tau   = tau.names,
        Sigma = Sigma.names
    ))
}

##############################################################################
