## Sequential normal model interface for bdpopt.

## This provides an interface to a simple group sequential model.

## Construct a sequential normal decision problem.
## n.stages is the number of stages.
## tau is the standard deviation for the prior normal distribution for theta.
## The sufficient statistic is taken to be the mean of the updated normal distribution of theta.
## sigma is the population standard deviation for each observation.

## The possible decisions are
## "c", for continue,
## "s", for stop,
## "f", for finalise, where
## "c" is available in all but the last stage,
## "s" is available in all stages, and
## "f" is only available in the last stage.
sequential.normal.dp <- function(n.stages, group.size, tau, sigma, stage.cost, final.cost, final.gain) {
    if ( !(is.numeric(n.stages) && length(n.stages) == 1 && n.stages >= 1) )
        stop("'n.stages' must be an integer greater than or equal to 1")

    if ( !(is.numeric(group.size) && length(group.size) == 1 && group.size >= 1) )
        stop("'group.size' must be an integer greater than or equal to 1")

    if ( !(is.numeric(tau) && length(tau) == 1 && tau > 0) )
        stop("'tau' must be a number greater than 0")

    if ( !(is.numeric(sigma) && length(sigma) == 1 && sigma > 0) )
        stop("'sigma' must be a number greater than 0")

    if ( !(is.numeric(stage.cost) && length(stage.cost) == 1) )
        stop("'stage.cost' must be a number")

    if ( !(is.numeric(final.cost) && length(final.cost) == 1) )
        stop("'final.cost' must be a number")

    if ( !(is.numeric(final.gain) && length(final.gain) == 1) )
        stop("'final.gain' must be a number")
    
    ## The variance of theta at a given stage
    theta.var <- function(stage) 1 / (1 / tau^2 + (stage - 1) / (sigma^2 / group.size))    
    
    ## The state in each stage is the pair (mu, var)
    ## consisting of the mean and var of the current normal distribution for the unknown
    ## population mean theta.
    post.sample <- function(stage, s, n.sims) {        
        as.list(rnorm(n.sims, mean = s, sd = sqrt(theta.var(stage))))
    }            

    ## The only possible value of d for pred.sample and update.state is "c",
    ## since these will only be called when continuing for this model.
    pred.sample <- function(stage, thetas, d) {
        Xs <- rnorm(length(thetas), mean = 0, sd = sqrt(sigma^2 / group.size))
        as.list(as.numeric(thetas) + Xs)        
    }     

    update.state <- function(stage, s, d, Xs) {
        prior.var <- theta.var(stage)
        post.var <- theta.var(stage + 1)
        obs.var <- sigma^2 / group.size

        as.list((s / prior.var + as.numeric(Xs) / obs.var) * post.var)       
    }

    term.decisions <- c(lapply(rep("s", n.stages - 1), list), list(list("s", "f")))
    term.obs.decisions <- lapply(1:n.stages, function(x) list())    
    cont.decisions <- c(lapply(rep("c", n.stages - 1), list), list(list()))

    ## Define terminal, terminal observation and continuation utility functions
    tuf <- function(d, theta) {
        if (d == "s") return(0)
        if (d == "f") return(final.gain * theta - final.cost)            
    }    
    term.utility.fun <- lapply(1:n.stages, function(x) tuf)
                               
    term.obs.utility.fun <- as.list(rep(NA, n.stages))    

    cuf <- function(d, X) -stage.cost
    cont.utility.fun <- lapply(1:n.stages, function(x) cuf)
    
    sequential.dp(n.stages,
                  post.sample,
                  pred.sample,
                  update.state,
                  term.decisions,
                  term.obs.decisions,
                  cont.decisions,                                  
                  term.utility.fun,
                  term.obs.utility.fun,
                  cont.utility.fun)
}

## Solve a sequential normal dp.
optimise.sequential.normal.eu <- function(dp, range, step.size, prior.mean = 0, n.sims = 1000, plot.results = TRUE) {

    check.dp.structure(dp)
    
    if ( !(is.numeric(range) && length(range) == 1 && range > 0) )
        stop("'range' must be a number greater than 0")

    if ( !(is.numeric(step.size) && length(step.size) == 1 && step.size > 0) )
        stop("'step.size' must be a number greater than 0")

    if ( !(is.numeric(prior.mean) && length(prior.mean)) )
        stop("'prior.mean' must be a number")
    
    solve.out <- optimise.sequential.eu(dp,
                                        prior.mean - range / 2,
                                        prior.mean + range / 2,
                                        step.size,
                                        rep(n.sims, dp$n.stages))

    if (plot.results == TRUE) {
        xs <- 1:dp$n.stages

        y.grid <- seq(prior.mean - range / 2, prior.mean + range / 2, step.size)        
        ys <- vector(mode = "numeric", length = dp$n.stages)        

        plot(x = c(min(xs), max(xs)),
             y = c(prior.mean - range / 2, prior.mean + range / 2),
             type = "n",
             xlab = "Stage",
             ylab = "Cutoff for the posterior mean")
        
        for (i in 1:dp$n.stages) {          

            j <- Position(function(y) solve.out$opt.decision(i, y) == "c" ||
                          solve.out$opt.decision(i, y) == "f", y.grid)
                
            if (identical(j, NA_integer_))
                points(xs[i], prior.mean + range / 2, pch = 4)
            else if (identical(j, 1L))
                points(xs[i], prior.mean - range / 2, pch = 4)
            else
                points(xs[i], y.grid[j], pch = 1)
        }      
    }          

    solve.out
}
