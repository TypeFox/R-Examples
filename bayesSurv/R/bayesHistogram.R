#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       bayesHistogram.R                    ####
####                                                 ####
#### FUNCTIONS:  bayesHistogram                      ####
#########################################################

### ======================================
### bayesHistogram
### ======================================
## 26/10/2004: start working on it
##
bayesHistogram <- function(
      y1,
      y2,
      nsimul = list(niter = 10, nthin = 1, nburn = 0, nwrite = 10),
      prior,
      init = list(iter = 0),
      mcmc.par = list(type.update.a = "slice", k.overrelax.a = 1, k.overrelax.sigma = 1, k.overrelax.scale = 1),
      store = list(a = FALSE, y = FALSE, r = FALSE),
      dir = getwd())
{
  thispackage = "bayesSurv"
  #thispackage = NULL
  store <- bayesHistogram.checkStore(store)
  nsimul <- bayessurvreg.checknsimul(nsimul)
  
  ## Give a function call to be recorded in a resulting object.
  call <- match.call(expand.dots = TRUE)

  ## Extract all the design information from the function call
  des <- bayesHistogram.design(y1, y2)

  ## Manipulate with prior/initial information
  prior.init <- bayesHistogram.priorInit(prior, init, mcmc.par, des)
  
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
  bayesHistogram.writeHeaders(dir, des, prior.init, store)

  ## Combine similar parameters into one vector
  dims <- c(length(des$status)/des$dim)
  storeV <- c(store$a, store$y, store$r)
  nsimul.run1 <- c(nrun[1], nsimul$nthin, nwrite.run[1])
  nsimul.run2 <- c(nrun[2], nsimul$nthin, nwrite.run[2])

#print(prior.init$Gparmi)
#print(prior.init$Gparmd)

  ## Burn-in
  keep.init <- prior.init
  cat("Simulation started on                       ", date(), "\n", sep = "")
  fit <- .C("bayesHistogram", as.character(dir),
                              dims = as.integer(dims),
                              y.left = as.double(des$y.left),
                              y.right = as.double(des$y.right),
                              status = as.integer(des$status),
                              r = as.integer(prior.init$r),
                              Ys = as.double(prior.init$y),
                              iter = as.integer(prior.init$iter),
                              specif = as.integer(prior.init$specification),
                              GsplineI = as.integer(prior.init$Gparmi),
                              GsplineD = as.double(prior.init$Gparmd),
                              nsimul = as.integer(nsimul.run1),
                              store = as.integer(storeV),
                              mainSimul = as.integer(0),
                              err = integer(1),
            PACKAGE = thispackage)

  if (fit$err != 0) stop ("Something went wrong during the simulation.")
  cat("Burn-up finished on                         ", date(), "   (iteration ", fit$iter, ")", "\n", sep = "")

  ## Give new initials
  prior.init$r <- fit$r
  prior.init$y <- fit$Ys
  prior.init$iter <- fit$iter

  nGparmi <- names(prior.init$Gparmi)
  nGparmd <- names(prior.init$Gparmd)
  prior.init$Gparmi <- fit$GsplineI
  prior.init$Gparmd <- fit$GsplineD
  names(prior.init$Gparmi) <- nGparmi
  names(prior.init$Gparmd) <- nGparmd  
  
  ## Rewrite sampled values by new files
  bayesHistogram.writeHeaders(dir, des, prior.init, store)
  
  ## Main simulation
  fit <- .C("bayesHistogram", as.character(dir),
                              dims = as.integer(dims),
                              y.left = as.double(des$y.left),
                              y.right = as.double(des$y.right),
                              status = as.integer(des$status),
                              r = as.integer(prior.init$r),
                              Ys = as.double(prior.init$y),
                              iter = as.integer(prior.init$iter),
                              specif = as.integer(prior.init$specification),
                              GsplineI = as.integer(prior.init$Gparmi),
                              GsplineD = as.double(prior.init$Gparmd),
                              nsimul = as.integer(nsimul.run2),
                              store = as.integer(storeV),
                              mainSimul = as.integer(1),
                              err = integer(1),
            PACKAGE = thispackage)
  if (fit$err != 0) stop ("Something went wrong during the simulation.")  
  cat("Simulation finished on                      ", date(), "   (iteration ", fit$iter, ")", "\n", sep = "")     

  toreturn <- fit$iter
  attr(toreturn, "call") <- call
  attr(toreturn, "prior") <- attr(keep.init, "prior")
  attr(toreturn, "init") <- attr(keep.init, "init")
  attr(toreturn, "mcmc.par") <- attr(keep.init, "mcmc.par")  
  class(toreturn) <- "bayesHistogram"
  
  return(toreturn)
}  
