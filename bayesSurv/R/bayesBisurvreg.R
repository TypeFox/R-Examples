#################################################
#### AUTHOR:     Arnost Komarek              ####
####             (2005)                      ####
####                                         ####
#### FILE:       bayesBisurvreg.R            ####
####                                         ####
#### FUNCTIONS:  bayesBisurvreg              ####
#################################################

### ======================================
### bayesBisurvreg
### ======================================
## 21/12/2004: start working on it
##
bayesBisurvreg <- function
(  formula,
   formula2,
   data = parent.frame(),
   na.action = na.fail,
   onlyX = FALSE,
   nsimul = list(niter = 10, nthin = 1, nburn = 0, nwrite = 10),
   prior,
   prior.beta,
   init = list(iter = 0),
   mcmc.par = list(type.update.a = "slice", k.overrelax.a = 1, k.overrelax.sigma = 1, k.overrelax.scale = 1),
   prior2,
   prior.beta2,
   init2,
   mcmc.par2 = list(type.update.a = "slice", k.overrelax.a = 1, k.overrelax.sigma = 1, k.overrelax.scale = 1),                           
   store = list(a = FALSE, a2 = FALSE, y = FALSE, y2 = FALSE, r = FALSE, r2 = FALSE),
   dir = getwd())
{
  thispackage = "bayesSurv"
  #thispackage = NULL
  store <- bayesBisurvreg.checkStore(store)
  nsimul <- bayessurvreg.checknsimul(nsimul)
  
  transform = function(t){log(t)}
  dtransform = function(t){1/t}
  transform2 = function(t){t}
  dtransform2 = function(t){1}
  
  ## Give a function call to be recorded in a resulting object.
  call <- match.call(expand.dots = TRUE)

  ## Extract all the design information from the function call
  ## Do not transform the response at this moment
  doubly <- ifelse(missing(formula2), FALSE, TRUE)
  m <- match.call(expand.dots = FALSE)   
  des <- bayessurvreg.design(m=m, formula=formula, data=data, transform=transform2, dtransform=transform2)  
  if (doubly) des2 <- bayessurvreg.design(m=m, formula=formula2, data=data, transform=transform2, dtransform=transform2)
  else        des2 <- list(nX = 0, n = des$n, nrandom = 0)

  if (onlyX){
    if (doubly) return(list(X=des$X, X2=des2$X))
    else        return(des$X)
  } 

  ## Find out dimension (univariate/bivariate) and perform some checks
  nP <- des$ncluster
  if (sum(table(des$cluster) == 1) == nP){       dim <- 1 }
  else{ if (sum(table(des$cluster) == 2) == nP){ dim <- 2 }
        else                                   { stop("Incorrect cluster indicator or missing data indicated") }
      }  

  if (doubly){
    nP2 <- des2$ncluster
    if (sum(table(des2$cluster) == 1) == nP2)      dim2 <- 1
    else if (sum(table(des2$cluster) == 2) == nP2) dim2 <- 2
         else stop("Incorrect cluster indicator in formula2 or missing data indicated")
    if (dim != dim2) stop("Inconsistent formula and formula2 indicated")
    if (nP != nP2) stop("Inconsistent formula and formula2 indicated")
  }

  ## Priors and inits for beta parameters
  if (missing(init))       init <- list()
  if (missing(init2))       init2 <- list()
  
  if (missing(prior.beta)) prior.beta <- list()
  betadi <- bayesBisurvreg.priorBeta(prior.beta, init, des)
  init$beta  <- attr(betadi, "init")
  prior.beta <- attr(betadi, "prior.beta")
  
  if (missing(prior.beta2)) prior.beta2 <- list()
  betadi2 <- bayesBisurvreg.priorBeta(prior.beta2, init2, des2)
  init2$beta  <- attr(betadi2, "init")  
  prior.beta2 <- attr(betadi2, "prior.beta")  

  ## Priors and inits for G-spline, (censored) observations and allocations
  if (!doubly){
    prior2 <- list()
    mcmc.par2 <- list()
  }
  prinit <- bayesBisurvreg.priorInit(dim, prior, init, des, mcmc.par, prior2, init2, des2, mcmc.par2, doubly)
  init      <- attr(prinit, "init")
  prior     <- attr(prinit, "prior")  
  mcmc.par  <- attr(prinit, "mcmc.par")
  init2     <- attr(prinit, "init2")
  prior2    <- attr(prinit, "prior2")  
  mcmc.par2 <- attr(prinit, "mcmc.par2")
  
  ## Compute quantities to determine the space needed to be allocated
  ##   and numbers of iterations in different phases
  if (nsimul$nburn >= nsimul$niter) nsimul$nburn <- nsimul$niter - 1
  if (nsimul$nburn < 0) nsimul$nburn <- 0
 
  if (nsimul$nburn == 0) nruns <- 1
  else                   nruns <- 2

  nrun <- numeric(2)
  nrun[2] <- nsimul$niter - nsimul$nburn
  nrun[1] <- nsimul$nburn

  nwrite.run <- nrun
  nwrite.run[nsimul$nwrite <= nrun] <- nsimul$nwrite   
  max.nwrite <- max(nwrite.run)

  ## Write headers to files with stored values
  bayesBisurvreg.writeHeaders(dir, dim, nP, doubly, prinit, store, des, des2)
  
  ## Combine similar parameters into one vector
  dims <- c(nP, as.numeric(doubly))
  storeV <- c(store$a, store$y, store$r, store$a2, store$y2, store$r2)
  nsimul.run1 <- c(nrun[1], nsimul$nthin, nwrite.run[1])
  nsimul.run2 <- c(nrun[2], nsimul$nthin, nwrite.run[2])
  
  cat("Simulation started on                       ", date(), "\n", sep = "")
  fit <- .C("bayesBisurvreg", as.character(dir),
                              dims = as.integer(dims),
                              X1 = as.double(if(des$nX) t(des$X) else 0),
                              X2 = as.double(if(des2$nX) t(des2$X) else 0),
                              y1.left = as.double(prinit$y.left),
                              y1.right = as.double(prinit$y.right),
                              status1 = as.integer(prinit$status),
                              t2.left = as.double(prinit$t2.left),
                              t2.right = as.double(prinit$t2.right),
                              status2 = as.integer(prinit$status2),
                              Ys1 = as.double(prinit$y),
                              Ys2 = as.double(prinit$y2),
                              r1 = as.integer(prinit$r),
                              r2 = as.integer(prinit$r2),
                              specif = as.integer(prinit$specification),
                              GsplI1 = as.integer(prinit$Gparmi),
                              GsplD1 = as.double(prinit$Gparmd),
                              GsplI2 = as.integer(prinit$Gparmi2),
                              GsplD2 = as.double(prinit$Gparmd2),
                              priorBetaI1 = as.integer(betadi$parmI),
                              priorBetaD1 = as.double(betadi$parmD),
                              priorBetaI2 = as.integer(betadi2$parmI),
                              priorBetaD2 = as.double(betadi2$parmD),
                              iter = as.integer(prinit$iter),
                              nsimul = as.integer(nsimul.run1),
                              store = as.integer(storeV),
                              mainSimul = as.integer(0),
                              err = integer(1),
           PACKAGE = thispackage)

  if (fit$err != 0) stop ("Something went wrong during the simulation.")
  cat("Burn-up finished on                         ", date(), "   (iteration ", fit$iter, ")", "\n", sep = "")

  ## Rewrite sampled values by new files
  bayesBisurvreg.writeHeaders(dir, dim, nP, doubly, prinit, store, des, des2)
  
  ## Main simulation
  fit <- .C("bayesBisurvreg", as.character(dir),
                              dims = as.integer(dims),
                              X1 = as.double(fit$X1),
                              X2 = as.double(fit$X2),
                              y1.left = as.double(fit$y1.left),
                              y1.right = as.double(fit$y1.right),
                              status1 = as.integer(fit$status1),
                              t2.left = as.double(fit$t2.left),
                              t2.right = as.double(fit$t2.right),
                              status2 = as.integer(fit$status2),
                              Ys1 = as.double(fit$Ys1),
                              Ys2 = as.double(fit$Ys2),
                              r1 = as.integer(fit$r1),
                              r2 = as.integer(fit$r2),
                              specif = as.integer(fit$specif),
                              GsplI1 = as.integer(fit$GsplI1),
                              GsplD1 = as.double(fit$GsplD1),
                              GsplI2 = as.integer(fit$GsplI2),
                              GsplD2 = as.double(fit$GsplD2),
                              priorBetaI1 = as.integer(fit$priorBetaI1),
                              priorBetaD1 = as.double(fit$priorBetaD1),
                              priorBetaI2 = as.integer(fit$priorBetaI2),
                              priorBetaD2 = as.double(fit$priorBetaD2),
                              iter = as.integer(fit$iter),
                              nsimul = as.integer(nsimul.run2),
                              store = as.integer(fit$store),
                              mainSimul = as.integer(1),            
                              err = integer(1),
          PACKAGE = thispackage)  
  if (fit$err != 0) stop ("Something went wrong during the simulation.")  
  cat("Simulation finished on                      ", date(), "   (iteration ", fit$iter, ")", "\n", sep = "")     

  toreturn <- fit$iter
  attr(toreturn, "call") <- call

  attr(toreturn, "init") <- init
  attr(toreturn, "prior") <- prior
  attr(toreturn, "prior.beta") <- prior.beta
  attr(toreturn, "mcmc.par") <- mcmc.par
  
  if (doubly){
    attr(toreturn, "init2") <- init2
    attr(toreturn, "prior2") <- prior2
    attr(toreturn, "prior.beta2") <- prior.beta2
    attr(toreturn, "mcmc.par2") <- mcmc.par2
  }

  class(toreturn) <- "bayesBisurvreg"
  return(toreturn)    
}
