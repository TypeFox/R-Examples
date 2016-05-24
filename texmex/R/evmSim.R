evmSim <- function(o, priorParameters, prop.dist,
                    jump.const, jump.cov, iter, start,
                    thin, burn,
                    verbose, trace, theCall, ...){
    if (class(o) != "evmOpt"){
        stop("o must be of class 'evmOpt'")
    }
    # Run checks and initialize algorithm
    wh <- texmexCheckMap(o)
    jump.const <- texmexJumpConst(jump.const, o)
    seed <- initRNG()
    cov <- if (missing(jump.cov)) { o$cov } else { jump.cov }

    # Initialize matrix to hold chain
    res <- matrix(ncol=length(o$coefficients), nrow=iter)
    res[1,] <- if (is.null(start)) { o$coefficients } else { start }

    ############################# Get prior and log-likelihood
    prior <- .make.mvn.prior(priorParameters)

    evm.log.lik <- o$family$log.lik(o$data, th=o$threshold, ...)

    log.lik <- function(param) {
      evm.log.lik(param) + prior(param)
    }

    # create proposals en bloc
    proposal.fn <- switch(prop.dist,
                          gaussian=rmvnorm,
                          cauchy=.rmvcauchy,
                          function () {stop("Bad proposal distribution")})

    proposals <- proposal.fn(iter,
                             double(length(o$coefficients)),
                             cov*jump.const)

    ######################## Run the Metropolis algorithm...
    res <- texmexMetropolis(res, log.lik, proposals, verbose, trace)

    # XXX Tidy up the object below. Doesn't need any of the info in o, or acceptance
    res <- list(call=theCall, map = o,
                burn = burn, thin = thin,
                chains=res, seed=seed)

    oldClass(res) <- "evmSim"
    res <- thinAndBurn(res)
    res
}


texmexMetropolis <-
    # Metropolis algorithm. 
    # x is a matrix, initialized to hold the chain. It's first row should be the
    #    starting point of the chain.
    # proposals is a matrix of proposals
function(x, log.lik, proposals, verbose, trace){
    last.cost <- log.lik(x[1,])
    if (!is.finite(last.cost)) {
      stop("Start is infeasible.")
    }

    acc <- 0
    for(i in 2:nrow(x)){
      if( verbose){
        if(i %% trace == 0) cat(i, " steps taken\n" )
      }
      prop <- proposals[i - 1,] + x[i - 1,]
      top <- log.lik(prop)
      delta <- top - last.cost
      if (is.finite(top) && ((delta >= 0) ||
                             (runif(1) <= exp(delta)))) {
        x[i, ] <- prop
        last.cost <- top
        acc <- 1 + acc
      }
      else {
        x[i, ] <- x[i-1,]
      }
    } # Close for(i in 2:nrow

    acc <- acc / nrow(x)
    if (acc < .1) {
        warning("Acceptance rate in Metropolis algorithm is low.")
    }
    if ((trace < nrow(x)) & verbose) {
        cat("Acceptance rate:", round(acc , 3) , "\n")
    }

    attr(x, "acceptance") <- acc
    x
}

####### Support functions

# Need to check for convergence failure here. Otherwise, end up simulating
# proposals from distribution with zero variance in 1 dimension.
texmexCheckMap <- function(map){
    checkNA <- any(is.na(sqrt(diag(map$cov))))
    if (checkNA) {
        stop("MAP estimates have not converged or have converged on values for which the variance cannot be computed. Cannot proceed. Try a different prior" )
    }
    NULL
}

texmexJumpConst <- function(j, map){
    if (is.null(j)){
        j <- (2.4/sqrt(length(map$coefficients)))^2
    }
    j
}

initRNG <- function(){
    if (!exists(".Random.seed")){ runif(1)  }
    .Random.seed # Retain and add to output
}
