#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       give.init.Gspline.R                 ####
####                                                 ####
#### FUNCTIONS:  give.init.Gspline                   ####
#########################################################

### ======================================
### give.init.Gspline
### ======================================
##
## Extract initial and prior information about the G-spline
##
## 16/01/2005
## ========================================================
give.init.Gspline <- function(prior, init, mcmc.par, dim)
{
  thispackage = "bayesSurv"
  #thispackage = NULL

  if(length(prior) == 0) inprior <- "arnost"
  else                   inprior <- names(prior)
  if(length(init) == 0) ininit <- "arnost"
  else                  ininit <- names(init)    
  if(length(mcmc.par) == 0) inmcmc.par <- "arnost"
  else                      inmcmc.par <- names(mcmc.par)
  
  if (dim <= 0 | dim > 2) stop("Dimension must be either 1 or 2")

  
  ## Model specification
  ## ===================
  tmp <- match("specification", inprior, nomatch=NA)
  if(is.na(tmp)) prior$specification <- 2  
  if (is.na(prior$specification)) prior$specification <- 2
  prior$specification <- prior$specification[1]
  if (prior$specification != 1 && prior$specification != 2) stop("Incorrect prior$specification given")
  
  ## Number of knots on the left and right
  ## =====================================
  tmp <- match("K", inprior, nomatch=NA)
  if(is.na(tmp)) prior$K <- rep(15, dim)
  prior$K <- prior$K[1:dim]
  if (sum(is.na(prior$K))) stop("Incorrect prior$K given")

  ## Index of the reference knot
  ## ===========================
  tmp <- match("izero", inprior, nomatch=NA)
  if(is.na(tmp)) prior$izero <- rep(0, dim)
  prior$izero <- prior$izero[1:dim]
  if (sum(is.na(prior$izero))) stop("Incorrect prior$izero given")
  if (sum(prior$izero < -prior$K) + sum(prior$izero > prior$K)) stop("Incorrect (out of range) prior$izero given")

  ## Neighbor system
  ## ===============
  tmp <- match("neighbor.system", inprior, nomatch=NA)
  if(is.na(tmp)) prior$neighbor.system <- ifelse(dim==1, "uniCAR", "eight.neighbors")
  neighbor.system <- pmatch(prior$neighbor.system, table=c("uniCAR", "eight.neighbors", "twelve.neighbors"), nomatch=0) - 1
  neighbor.system <- neighbor.system[1]
  if (neighbor.system == -1) stop("Incorrect prior$neighbor.system")
  if (dim == 1 & neighbor.system != 0) stop("Neighbor system must be 'uniCAR' for univariate data")

  ## Order of conditional autoregression
  ## ====================================
  tmp <- match("order", inprior, nomatch=NA)
  if(is.na(tmp)) prior$order <- 3
  if (neighbor.system == 1) prior$order <- 2
  else                      if (neighbor.system == 2) prior$order <- 3    
  prior$order <- prior$order[1]
  if (is.na(prior$order)) stop("Incorrect prior$order given")
  if (prior$order > 3) stop("Order of CAR higher than 3 is not implemented")

  ## Equal lambda
  ## ============
  tmp <- match("equal.lambda", inprior, nomatch=NA)
  if(is.na(tmp)) prior$equal.lambda <- TRUE
  if (dim == 1) prior$equal.lambda <- TRUE
  if (neighbor.system > 0) prior$equal.lambda <- TRUE
  prior$equal.lambda <- prior$equal.lambda[1]
  if (is.na(prior$equal.lambda)) stop("Incorrect prior$equal.lambda given")
  equal.lambda <- as.numeric(prior$equal.lambda)

  ## Prior lambda and parameters for the prior
  ## =========================================
  tmp <- match("prior.lambda", inprior, nomatch=NA)
  if(is.na(tmp)) stop("prior$prior.lambda must be given")
  if (prior$equal.lambda){
    prior$prior.lambda <- rep(prior$prior.lambda[1], dim)
  }    
  else{
    prior$prior.lambda <- prior$prior.lambda[1:dim]
  }    
  if (sum(is.na(prior$prior.lambda))) stop("Incorrect prior$prior.lambda given")
  prior.lambda <- apply(matrix(prior$prior.lambda, ncol=1), 1, pmatch, table=c("fixed", "gamma", "sduniform")) - 1  
  if (sum(prior.lambda == -1)) stop("Incorrect prior$prior.lambda given")
   
  rate.lambda <- rep(0, dim)
  shape.lambda <- rep(0, dim)
  if (sum(prior.lambda > 0)){     ## at least one prior is not 'fixed'
    if (is.null(prior$rate.lambda)) stop("Either rate for gamma prior or upper limit of uniform prior for lambda was not given")
    if (prior$equal.lambda) rate.lambda <- rep(prior$rate.lambda[1], dim)
    else                    rate.lambda <- prior$rate.lambda[1:dim]
    if (sum(is.na(rate.lambda[prior.lambda > 0]))) stop("Incorrect prior$rate.lambda given")
    if (sum(rate.lambda[prior.lambda > 0] <= 0)) stop("prior$rate.lambda must be positive")
    rate.lambda[prior.lambda == 0] <- 0
  }
  if (sum(prior.lambda == 1)){  ## at least one 'gamma' prior
      if (is.null(prior$shape.lambda)) stop("Shape for gamma prior for lambda was not given")
      if (prior$equal.lambda) shape.lambda <- rep(prior$shape.lambda[1], dim)
      else                    shape.lambda <- prior$shape.lambda[1:dim]
      if (sum(is.na(shape.lambda[prior.lambda == 1]))) stop("Incorrect prior$shape.lambda given")
      if (sum(shape.lambda[prior.lambda == 1] <= 0)) stop("prior$shape.lambda must be positive")
      shape.lambda[prior.lambda != 1] <- 0
  }      
  parms.lambda <- as.vector(rbind(shape.lambda, rate.lambda))
  prior$shape.lambda <- shape.lambda
  prior$rate.lambda <- rate.lambda  
  
  ## Prior sigma and parameters for the prior
  ## ========================================
  if (prior$specification == 2){
    prior$prior.sigma <- rep("fixed", dim)
    inprior <- names(prior)
  }
  tmp <- match("prior.sigma", inprior, nomatch=NA)
  if(is.na(tmp)) stop("prior$prior.sigma must be given")
  prior$prior.sigma <- prior$prior.sigma[1:dim]
  if (sum(is.na(prior$prior.sigma))) stop("Incorrect prior$prior.sigma given")
  prior.sigma <- apply(matrix(prior$prior.sigma, ncol=1), 1, pmatch, table=c("fixed", "gamma", "sduniform")) - 1    
  if (sum(prior.sigma == -1)) stop("Incorrect prior$prior.sigma given")
   
  rate.sigma <- rep(0, dim)
  shape.sigma <- rep(0, dim)
  if (sum(prior.sigma > 0)){        ## at least one prior is not 'fixed'
    if (is.null(prior$rate.sigma)) stop("Either rate for gamma prior or upper limit of uniform prior for sigma was not given")
    rate.sigma <- prior$rate.sigma[1:dim]
    if (sum(is.na(rate.sigma[prior.sigma > 0]))) stop("Incorrect prior$rate.sigma given")
    if (sum(rate.sigma[prior.sigma > 0] <= 0)) stop("prior$rate.sigma must be positive")
    rate.sigma[prior.sigma == 0] <- 0
  }
  if (sum(prior.sigma == 1)){    ## at least one 'gamma' prior
      if (is.null(prior$shape.sigma)) stop("Shape for gamma prior for sigma was not given")
      shape.sigma <- prior$shape.sigma[1:dim]
      if (sum(is.na(shape.sigma[prior.sigma == 1]))) stop("Incorrect prior$shape.sigma given")
      if (sum(shape.sigma[prior.sigma == 1] <= 0)) stop("prior$shape.sigma must be positive")
      shape.sigma[prior.sigma != 1] <- 0
  }      
  parms.sigma <- as.vector(rbind(shape.sigma, rate.sigma))
  prior$shape.sigma <- shape.sigma
  prior$rate.sigma <- rate.sigma  

  tmp <- match("k.overrelax.sigma", inmcmc.par, nomatch=NA)
  if (!is.na(tmp) & length(mcmc.par$k.overrelax.sigma) == 1) mcmc.par$k.overrelax.sigma <- rep(mcmc.par$k.overrelax.sigma[1], dim)
  if(is.na(tmp)) mcmc.par$k.overrelax.sigma <- rep(1, dim)
  mcmc.par$k.overrelax.sigma <- mcmc.par$k.overrelax.sigma[1:dim]
  if (sum(is.na(mcmc.par$k.overrelax.sigma))) stop("Incorrect mcmc.par$k.overrelax.sigma given")
  if (sum(mcmc.par$k.overrelax.sigma <= 0)) stop("Incorrect mcmc.par$k.overrelax.sigma given (must be all positive)")  
  
  ## Prior gamma and parameters for the prior
  ## ========================================
  if (prior$specification == 2){
    prior$prior.gamma <- rep("fixed", dim)
    inprior <- names(prior)    
  }
  tmp <- match("prior.gamma", inprior, nomatch=NA)
  if(is.na(tmp)) stop("prior$prior.gamma must be given")
  prior$prior.gamma <- prior$prior.gamma[1:dim]
  if (sum(is.na(prior$prior.gamma))) stop("Incorrect prior$prior.gamma given")
  prior.gamma <- apply(matrix(prior$prior.gamma, ncol=1), 1, pmatch, table=c("fixed", "normal")) - 1
  if (sum(prior.gamma == -1)) stop("Incorrect prior$prior.gamma given")

  mean.gamma <- rep(0, dim)
  var.gamma <- rep(0, dim)  
  if (sum(prior.gamma > 0)){     ## at least one prior is not 'fixed'
    tmp <- match("mean.gamma", inprior, nomatch=NA)
    if(is.na(tmp)) stop("prior$mean.gamma must be given")
    tmp <- match("var.gamma", inprior, nomatch=NA)
    if(is.na(tmp)) stop("prior$var.gamma must be given")
    mean.gamma <- prior$mean.gamma[1:dim]
    var.gamma <- prior$var.gamma[1:dim]
    if (sum(is.na(mean.gamma[prior.gamma > 0]))) stop("Incorrect prior$mean.gamma given")
    if (sum(is.na(var.gamma[prior.gamma > 0]))) stop("Incorrect prior$var.gamma given")    
    mean.gamma[prior.gamma == 0] <- 0
    var.gamma[prior.gamma == 0] <- 0    
  }
  parms.gamma <- as.vector(rbind(mean.gamma, var.gamma))
  prior$mean.gamma <- mean.gamma
  prior$var.gamma <- var.gamma

  ## Prior scale and parameters for the prior
  ## ========================================
  if (prior$specification == 1){
    prior$prior.scale <- rep("fixed", dim)
    inprior <- names(prior)
  }
  tmp <- match("prior.scale", inprior, nomatch=NA)
  if(is.na(tmp)) stop("prior$prior.scale must be given")
  prior$prior.scale <- prior$prior.scale[1:dim]
  if (sum(is.na(prior$prior.scale))) stop("Incorrect prior$prior.scale given")
  prior.scale <- apply(matrix(prior$prior.scale, ncol=1), 1, pmatch, table=c("fixed", "gamma", "sduniform")) - 1    
  if (sum(prior.scale == -1)) stop("Incorrect prior$prior.scale given")
   
  rate.scale <- rep(0, dim)
  shape.scale <- rep(0, dim)
  if (sum(prior.scale > 0)){        ## at least one prior is not 'fixed'
    tmp <- match("rate.scale", inprior, nomatch=NA)
    if(is.na(tmp)) stop("Either rate for gamma prior or upper limit of uniform prior for scale was not given")
    rate.scale <- prior$rate.scale[1:dim]
    if (sum(is.na(rate.scale[prior.scale > 0]))) stop("Incorrect prior$rate.scale given")
    if (sum(rate.scale[prior.scale > 0] <= 0)) stop("prior$rate.scale must be positive")
    rate.scale[prior.scale == 0] <- 0
  }
  if (sum(prior.scale == 1)){    ## at least one 'gamma' prior
      tmp <- match("shape.scale", inprior, nomatch=NA)
      if(is.na(tmp)) stop("Shape for gamma prior for scale was not given")
      shape.scale <- prior$shape.scale[1:dim]
      if (sum(is.na(shape.scale[prior.scale == 1]))) stop("Incorrect prior$shape.scale given")
      if (sum(shape.scale[prior.scale == 1] <= 0)) stop("prior$shape.scale must be positive")
      shape.scale[prior.scale != 1] <- 0
  }      
  parms.scale <- as.vector(rbind(shape.scale, rate.scale))
  prior$shape.scale <- shape.scale
  prior$rate.scale <- rate.scale  

  tmp <- match("k.overrelax.scale", inmcmc.par, nomatch=NA)
  if (!is.na(tmp) & length(mcmc.par$k.overrelax.scale) == 1) mcmc.par$k.overrelax.scale <- rep(mcmc.par$k.overrelax.scale[1], dim)  
  if(is.na(tmp)) mcmc.par$k.overrelax.scale <- rep(1, dim)
  mcmc.par$k.overrelax.scale <- mcmc.par$k.overrelax.scale[1:dim]
  if (sum(is.na(mcmc.par$k.overrelax.scale))) stop("Incorrect mcmc.par$k.overrelax.scale given")
  if (sum(mcmc.par$k.overrelax.scale <= 0)) stop("Incorrect mcmc.par$k.overrelax.scale given (must be all positive)")  

  ## Prior intercept and parameters for the prior
  ## ============================================
  if (prior$specification == 1){
    prior$prior.intercept <- rep("fixed", dim)
    inprior <- names(prior)
  }
  tmp <- match("prior.intercept", inprior, nomatch=NA)
  if(is.na(tmp)) stop("prior$prior.intercept must be given")
  prior$prior.intercept <- prior$prior.intercept[1:dim]
  if (sum(is.na(prior$prior.intercept))) stop("Incorrect prior$prior.intercept given")
  prior.intercept <- apply(matrix(prior$prior.intercept, ncol=1), 1, pmatch, table=c("fixed", "normal")) - 1
  if (sum(prior.intercept == -1)) stop("Incorrect prior$prior.intercept given")

  mean.intercept <- rep(0, dim)
  var.intercept <- rep(0, dim)  
  if (sum(prior.intercept > 0)){     ## at least one prior is not 'fixed'
    tmp <- match("mean.intercept", inprior, nomatch=NA)
    if(is.na(tmp)) stop("prior$mean.intercept must be given")
    tmp <- match("var.intercept", inprior, nomatch=NA)
    if(is.na(tmp)) stop("prior$var.intercept must be given")
    mean.intercept <- prior$mean.intercept[1:dim]
    var.intercept <- prior$var.intercept[1:dim]
    if (sum(is.na(mean.intercept[prior.intercept > 0]))) stop("Incorrect prior$mean.intercept given")
    if (sum(is.na(var.intercept[prior.intercept > 0]))) stop("Incorrect prior$var.intercept given")    
    mean.intercept[prior.intercept == 0] <- 0
    var.intercept[prior.intercept == 0] <- 0    
  }
  parms.intercept <- as.vector(rbind(mean.intercept, var.intercept))
  prior$mean.intercept <- mean.intercept
  prior$var.intercept <- var.intercept
  
  ## c4delta
  ## =======
  tmp <- match("c4delta", inprior, nomatch=NA)
  if(is.na(tmp)) prior$c4delta <- rep(1.5, dim)
  prior$c4delta <- prior$c4delta[1:dim]
  if (sum(is.na(prior$c4delta))) stop("Incorrect prior$c4delta given")
  if (sum(prior$c4delta <= 0)) stop("prior$c4delta must be positive")

  ## Type of update of a
  ## ===================
  tmp <- match("type.update.a", inmcmc.par, nomatch=NA)
  if(is.na(tmp)){
    if (dim == 1) mcmc.par$type.update.a <- "slice"
    else          mcmc.par$type.update.a <- "slice"
  }  
  tmp <- match("k.overrelax.a", inmcmc.par, nomatch=NA)
  if(is.na(tmp)) mcmc.par$k.overrelax.a <- 1
  mcmc.par$type.update.a <- mcmc.par$type.update.a[1]
  mcmc.par$k.overrelax.a <- mcmc.par$k.overrelax.a[1]
  if (is.na(mcmc.par$type.update.a)) stop("Incorrect mcmc.par$type.update.a given")
  if (is.na(mcmc.par$k.overrelax.a)) stop("Incorrect mcmc.par$k.overrelax.a given")  
  type.update.a <- pmatch(mcmc.par$type.update.a, table=c("slice", "ars.quantile", "ars.mode", "block"), nomatch=0) - 1
  if (type.update.a == -1) stop("Incorrect mcmc.par$type.update.a given")
  if (type.update.a == 0){
    if (mcmc.par$k.overrelax.a <= 0) stop("Incorrect mcmc.par$k.overrelax.a given (must be positive)")  
  }
  else{
    mcmc.par$k.overrelax.a <- 1
  }    
  
  ## Index of the first iteration
  ## =============================
  tmp <- match("iter", ininit, nomatch=NA)
  if(is.na(tmp)) init$iter <- 0
  if (is.na(init$iter)) init$iter <- 0                                 
  else                  init$iter <- init$iter[1]
  
  ## Initial values for lambda
  ## =========================
  tmp <- match("lambda", ininit, nomatch=NA)
  if(is.na(tmp)) stop("Initial lambda must be given")
  if (prior$equal.lambda) init$lambda <- rep(init$lambda[1], dim)
  else                    init$lambda <- init$lambda[1:dim]
  if (sum(is.na(init$lambda))) stop("Incorrect init$lambda given")
  if (sum(init$lambda <= 0)) stop("init$lambda must be positive")

  ## Initial values for sigma
  ## =========================
  tmp <- match("sigma", ininit, nomatch=NA)
  if(is.na(tmp)){
    if (prior$specification == 2) init$sigma <- rep(0.2, dim)
    else                          stop("Initial sigma (basis std. deviation) must be given")
  }    
  init$sigma <- init$sigma[1:dim]
  if (sum(is.na(init$sigma))) stop("Incorrect init$sigma given")
  if (sum(init$sigma <= 0)) stop("init$sigma must be positive")

  ## Initial values for gamma
  ## ========================
  tmp <- match("gamma", ininit, nomatch=NA)
  if(is.na(tmp)){
    if (prior$specification == 2) init$gamma <- rep(0, dim)
    else                          stop("Initial gamma must be given")
  }    
  init$gamma <- init$gamma[1:dim]
  if (sum(is.na(init$gamma))) stop("Incorrect init$gamma given")

  ## Initial values for scale
  ## =========================
  tmp <- match("scale", ininit, nomatch=NA)
  if(is.na(tmp)){  
    if (prior$specification == 1) init$scale <- rep(1, dim)
    else                          stop("Initial scale must be given")
  }
  if (prior$specification == 1) init$scale <- rep(1, dim)        ### Will this be temporary???
  init$scale <- init$scale[1:dim]
  if (sum(is.na(init$scale))) stop("Incorrect init$scale given")
  if (sum(init$scale <= 0)) stop("init$scale must be positive")

  ## Initial values for intercept
  ## ============================
  tmp <- match("intercept", ininit, nomatch=NA)
  if(is.na(tmp)){  
    if (prior$specification == 1) init$intercept <- rep(0, dim)
    else                          stop("Initial intercept must be given")
  }
  if (prior$specification == 1) init$intercept <- rep(0, dim)    ### Will this be temporary???
  init$intercept <- init$intercept[1:dim]
  if (sum(is.na(init$intercept))) stop("Incorrect init$intercept given")
  
  ## Initial values for a
  ## ====================
  #aconstraint <- "mean"   ### This causes problems if there are some too negative a's.
                           ### Then there must be also too positive a's and exp(a) is Inf
  aconstraint <- "reference"
  aconstraint <- pmatch(aconstraint, table=c("mean", "reference"), nomatch=0)-1
  if (aconstraint < 0) stop("Unimplemented identifiability constraint for a coefficients")
  
  ### require library(smoothSurv)!!!
    ## dim == 1, get initial 'a' by minimizing third order penalty on standardized knots
    ## dim == 2, get initial 'a' by minimizing third order penalty for each margin and then take product of marginal weights
    ##           as a joint weight (i.e. a[i,j] = a[i]+a[j] or w[i,j] = w[i]*w[j])
  tmp <- match("a", ininit, nomatch=NA)
  if(is.na(tmp) | (!is.na(tmp) & !length(init$a))){  
    acoef <- list()
    for (j in 1:dim){
      if (prior$K[j] == 0){
        acoef[[j]] <- 0
      }else{        
        delta <- prior$c4delta[j] * init$sigma[j]
        range <- delta * 2 * prior$K[j]        ## range of knots with this choice of delta
        delta2 <- (8/range) * delta          ## rescale it to knots from -4 to 4
        knots <- seq(-prior$K[j], prior$K[j], by = 1) * delta2
        sdspline <- delta2 / prior$c4delta[j]
        if (sdspline >= 0.95) sdspline <- 0.95
        minp <- minPenalty(knots = knots, sdspline = sdspline, difforder = 3, info = FALSE)
        if (minp$fail) stop("Unable to guess initial 'a' coefficients, give your own")
        acoef[[j]] <- minp$spline[, "a coef."]
        if (aconstraint == "mean"){
          acoef[[j]] <- as.vector(acoef[[j]] - mean(acoef[[j]]))
        }
        else{
          acoef[[j]] <- as.vector(acoef[[j]] - acoef[[j]][prior$izero[j] + prior$K[j] + 1])
        }
      }
    }
    if (dim == 1){ init$a <- acoef[[1]] }
    else{                 if (dim == 2){ init$a <- outer(acoef[[1]], acoef[[2]], "+")                           }
                          else         { stop("Unimplemented dimension appeared in bayesHistogram.priorInit()") }
        }
  }  
  total.length <- ifelse(dim == 1, 2*prior$K[1] + 1,
                                  (2*prior$K[1] + 1)*(2*prior$K[2] + 1))
  init$a <- init$a[1:total.length]
  if (sum(is.na(init$a))) stop("Incorrect init$a given")
  total.izero <- ifelse(dim == 1, prior$izero[1] + prior$K[1] + 1,
                                         (prior$izero[2]+prior$K[2])*(2*prior$K[1]+1) + prior$izero[1]+prior$K[1]+1)
#  if (init$a[total.izero] != 0) stop("Reference init$a should be equal to zero")
  acoef <- as.vector(init$a)
  if (dim == 2) init$a <- matrix(init$a, nrow=2*prior$K[1] + 1, ncol=2*prior$K[2] + 1)
  
  ## Put G-spline parameters to long vectors
  ## ========================================  
  Gparmi <- c(dim, neighbor.system, equal.lambda, prior$K, prior$izero, prior$order, prior.lambda, prior.gamma, prior.sigma,
              prior.intercept, prior.scale,
              type.update.a, mcmc.par$k.overrelax.a, mcmc.par$k.overrelax.sigma, mcmc.par$k.overrelax.scale,
              aconstraint)
  names(Gparmi) <- c("dim", "neighbor.system", "equal.lambda", paste("K", 1:dim, sep=""),
                      paste("izeroR", 1:dim, sep=""), "order",
                      paste("prior.for.lambda", 1:dim, sep=""),
                      paste("prior.for.gamma", 1:dim, sep=""),
                      paste("prior.for.sigma", 1:dim, sep=""),
                      paste("prior.for.intercept", 1:dim, sep=""),
                      paste("prior.for.scale", 1:dim, sep=""),
                      "type.update.a", "k.overrelax.a",
                      paste("k.overrelax.sigma", 1:dim, sep=""),
                      paste("k.overrelax.scale", 1:dim, sep=""),
                      "aconstraint")
  
  Gparmd <- c(acoef, init$lambda, init$gamma, init$sigma, init$intercept, init$scale,
              prior$c4delta, parms.lambda, parms.gamma, parms.sigma, parms.intercept, parms.scale)
  names(Gparmd) <- c(paste("a", 1:total.length, sep=""), paste("lambda", 1:dim, sep=""),
                     paste("gamma", 1:dim, sep=""), paste("sigma", 1:dim, sep=""),
                     paste("intercept", 1:dim, sep=""), paste("scale", 1:dim, sep=""),                     
                     paste("c4delta", 1:dim, sep=""),
                     paste(c("shape.lambda", "rate.lambda"), rep(1:dim, rep(2, dim)), sep=""),
                     paste(c("mean.gamma", "var.gamma"), rep(1:dim, rep(2, dim)), sep=""),
                     paste(c("shape.sigma", "rate.sigma"), rep(1:dim, rep(2, dim)), sep=""),
                     paste(c("mean.intercept", "var.intercept"), rep(1:dim, rep(2, dim)), sep=""),
                     paste(c("shape.scale", "rate.scale"), rep(1:dim, rep(2, dim)), sep=""))                     

  toreturn <- list(Gparmi = Gparmi, Gparmd = Gparmd, specification=prior$specification)
  attr(toreturn, "init") <- init
  attr(toreturn, "prior") <- prior
  attr(toreturn, "mcmc.par") <- mcmc.par
  
  return(toreturn)
}  

