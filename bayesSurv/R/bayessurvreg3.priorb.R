################################################
#### AUTHOR:     Arnost Komarek             ####
####             (2005)                     ####
####                                        ####
#### FILE:       bayessurvreg3.priorb.R     ####
####                                        ####
#### FUNCTIONS:  bayessurvreg3.priorb       ####
################################################

### ======================================
### bayessurvreg3.priorb
### ======================================
##
## \item{prior.b}{a~list specifying the G-spline prior of the random intercept. It is assumed to have the same structure as
##    \code{prior} parameter of this function or the function \code{\link{bayessurvreg2}}
## }
bayessurvreg3.priorb <- function(prior.b, init, design, mcmc.par)
{  
  if (design$nrandom){

    if(length(init) == 0) ininit <- "arnost"
    else                  ininit <- names(init)

    ## Initial values of the random intercept
    ## ======================================
    tmp <- match("b", ininit, nomatch=NA)
    if(is.na(tmp)){
      init$b <- rep(0, design$ncluster)
    }
    else{
      if (length(init$b) == 0){
        init$b <- rep(0, design$ncluster)
      } 
      else{
        if (length(init$b) < design$ncluster) stop("Incorrect init$b parameter supplied.")
        init$b <- init$b[1:design$ncluster]
      }
    }  
    if (sum(is.na(init$b))) stop("Incorrect init$b parameter supplied.")

    ## Create vectors for RandomEff constructor
    ## ========================================
    bparmI <- c(1, 1, design$ncluster, design$nwithin)    
    bparmD <- as.vector(init$b)
    names(bparmI) <- c("type.prior", "nRandom", "nCluster", paste("Cl", 1:design$ncluster, sep=""))
    names(bparmD) <- paste("b.", 1:design$ncluster, sep="")        

    ## Create vectors for G-spline constructor related to the random effect vector
    ## ===========================================================================

      ## Check prior.b before passing it to give.init.Gspline
      ## ----------------------------------------------------
    if (length(prior.b) == 0) inprior <- "arnost"
    else                      inprior <- names(prior.b)
    prior.tmp <- list()
    
    tmp <- match("specification", inprior, nomatch=NA)
    if(is.na(tmp)) prior.b$specification <- 2

    tmp <- match("K", inprior, nomatch=NA)
    if(is.na(tmp)) prior.b$K <- 15

    tmp <- match("izero", inprior, nomatch=NA)
    if(is.na(tmp)) prior.b$izero <- 0

    prior.b$neighbor.system <- "uniCAR"

    tmp <- match("order", inprior, nomatch=NA)
    if(is.na(tmp)) prior.b$order <- 3

    prior.b$equal.lambda <- TRUE

    tmp <- match("prior.lambda", inprior, nomatch=NA)
    if(is.na(tmp)) stop("prior.b$prior.lambda must be given")

    prior.b$prior.gamma <- "fixed"
    prior.b$prior.intercept <- "fixed"
    if (prior.b$specification == 1){
      prior.b$prior.scale <- "fixed"
      tmp <- match("prior.sigma", inprior, nomatch=NA)
      if(is.na(tmp)) stop("prior.b$prior.sigma must be given")      
    }      
    if (prior.b$specification == 2){
      prior.b$prior.sigma <- "fixed"
      tmp <- match("prior.scale", inprior, nomatch=NA)
      if(is.na(tmp)) stop("prior.b$prior.scale must be given")
    }

    tmp <- match("c4delta", inprior, nomatch=NA)
    if(is.na(tmp)) prior.b$c4delta <- 1.5
    
    
      ## mcmc.par for the function give.init.Gspline
      ## -------------------------------------------
    if (length(mcmc.par) == 0) inmcmc <- "arnost"
    else                       inmcmc <- names(mcmc.par)
    mcmc.par.tmp <- list()
    
    tmp <- match("type.update.a.b", inmcmc, nomatch=NA)
    if(is.na(tmp)) mcmc.par$type.update.a.b <- "slice"
    mcmc.par.tmp$type.update.a <- mcmc.par$type.update.a.b

    tmp <- match("k.overrelax.a.b", inmcmc, nomatch=NA)
    if(is.na(tmp)) mcmc.par$k.overrelax.a.b <- 1
    mcmc.par.tmp$k.overrelax.a <- mcmc.par$k.overrelax.a.b
    
    tmp <- match("k.overrelax.sigma.b", inmcmc, nomatch=NA)
    if(is.na(tmp)) mcmc.par$k.overrelax.sigma.b <- 1
    mcmc.par.tmp$k.overrelax.sigma <- mcmc.par$k.overrelax.sigma.b

    tmp <- match("k.overrelax.scale.b", inmcmc, nomatch=NA)
    if(is.na(tmp)) mcmc.par$k.overrelax.scale.b <- 1
    mcmc.par.tmp$k.overrelax.scale <- mcmc.par$k.overrelax.scale.b

      ## init for the function give.init.Gspline
      ## ---------------------------------------
    if(length(init) == 0) ininit <- "arnost"
    else                  ininit <- names(init)    
    init.tmp <- list()
    
    tmp <- match("iter", ininit, nomatch=NA)
    if(is.na(tmp)) init$iter <- 0
    init.tmp$iter <- init$iter

    tmp <- match("lambda.b", ininit, nomatch=NA)
    if(is.na(tmp)) stop("Initial lambda for random intcpt G-spline must be given")
    init.tmp$lambda <- init$lambda.b

    tmp <- match("sigma.b", ininit, nomatch=NA)
    if(is.na(tmp)){
      if (prior.b$specification == 2) init$sigma.b <- 0.2
      else                            stop("Initial sigma (basis std. deviation) for random intcpt G-spline must be given")
    }
    init.tmp$sigma <- init$sigma.b

    tmp <- match("gamma.b", ininit, nomatch=NA)
    if(is.na(tmp)) init$gamma.b <- 0
    else{
      if (init$gamma.b != 0){
        init$gamma.b <- 0
        warning("init$gamma.b changed to zero")
      }  
    }  
    init.tmp$gamma <- init$gamma.b

    tmp <- match("scale.b", ininit, nomatch=NA)
    if(is.na(tmp)){  
      if (prior.b$specification == 1) init$scale.b <- 1
      else                            stop("Initial scale for random intcpt G-spline must be given")      
    }
    init.tmp$scale <- init$scale.b

    tmp <- match("intercept.b", ininit, nomatch=NA)
    if(is.na(tmp)) init$intercept.b <- 0
    else{
      if (init$intercept.b != 0){
        init$intercept.b <- 0
        warning("init$intercept.b changed to zero")
      }  
    }      
    init.tmp$intercept <- init$intercept.b

    tmp <- match("a.b", ininit, nomatch=NA)
    if(is.na(tmp)) init.tmp$a <- numeric(0)
    else           init.tmp$a <- init$a.b

    bgspl <- give.init.Gspline(prior=prior.b, init=init.tmp, mcmc.par=mcmc.par.tmp, dim=1)
    init.tmp <- attr(bgspl, "init")
    mcmc.par.tmp <- attr(bgspl, "mcmc.par")
    prior.b <- attr(bgspl, "prior")

    init$iter        <- init.tmp$iter
    init$lambda.b    <- init.tmp$lambda
    init$sigma.b     <- init.tmp$sigma
    init$gamma.b     <- init.tmp$gamma
    init$scale.b     <- init.tmp$scale
    init$intercept.b <- init.tmp$intercept
    init$a.b         <- init.tmp$a

    mcmc.par$type.update.a.b     <- mcmc.par.tmp$type.update.a
    mcmc.par$k.overrelax.a.b     <- mcmc.par.tmp$k.overrelax.a
    mcmc.par$k.overrelax.sigma.b <- mcmc.par.tmp$k.overrelax.sigma
    mcmc.par$k.overrelax.scale   <- mcmc.par.tmp$k.overrelax.scale
    
    GsplI <- bgspl$Gparmi
    GsplD <- bgspl$Gparmd
    specification <- prior.b$specification

    ## Initial allocations
    ## ====================
    tmp <- match("r.b", ininit, nomatch=NA)
    if(is.na(tmp)) init$r.b <- numeric(0)
    init$r.b <- give.init.r(init$r.b, init$b, dim=1, KK=prior.b$K,
                            gamma=init$gamma.b,      sigma=init$sigma.b, c4delta=prior.b$c4delta,
                            intcpt=init$intercept.b, scale=init$scale.b)

      ## Recalculate r into the vector of length n with entries from {1,...,total_length}
    rsm <- init$r.b + prior.b$K[1] + 1 
    
  }
  else{
    bparmI <- c(1, 0, design$n, rep(1, design$n))
    bparmD <- rep(0, design$n)
    names(bparmI) <- c("type.prior", "nRandom", "nCluster", paste("Cl", 1:design$n, sep=""))
    names(bparmD) <- paste("b.", 1:design$n, sep="")

    GsplI <- rep(0, 15)
    names(GsplI) <- c("dim", "neighbor.system", "equal.lambda", "K", "izeroR", "order",
                      "prior.for.lambda", "prior.for.gamma", "prior.for.sigma", "prior.for.intercept", "prior.for.scale",
                      "type.update.a", "k.overrelax.a", "k.overrelax.sigma", "k.overrelax.scale")
    GsplD <- rep(0, 17)
    names(GsplD) <- c("a", "lambda", "gamma", "sigma", "intercept", "scale", "c4delta", "shape.lambda", "rate.lambda",
                      "mean.gamma", "var.gamma", "shape.sigma", "rate.sigma",
                      "mean.intercept", "var.intercept", "shape.scale", "rate.scale")
    specification <- 0

    prior.b <- list()
    init$r.b <- 0
    rsm <- 0
  } 

  toreturn <- list(bparmI=bparmI, bparmD=bparmD, GsplI=GsplI, GsplD=GsplD, specification=specification, r=rsm)
  attr(toreturn, "prior.b") <- prior.b
  attr(toreturn, "init") <- init
  attr(toreturn, "mcmc.par") <- mcmc.par
  
  return(toreturn)
}
