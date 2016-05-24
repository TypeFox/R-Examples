## -----------------------------------------------------------------------------
## Pseudorandom Search Optimisation Routine
## -----------------------------------------------------------------------------

pseudoOptim <- function (f, p, ..., lower, upper, control = list() ) {

  ## check input
  npar  <- length(p)
  if (npar == 1)
      stop("number of parameters to estimate should be > 1 in pseudoOptim")

  if (! all(is.finite(lower))) stop("lower cannot be Inf or -Inf")
  if (! all(is.finite(upper))) stop("upper cannot be Inf or -Inf")
  if (length(lower) != npar & length(lower)!= 1)
    stop("length of 'lower' should be either 1 or equal to number of parameters")
  if (length(upper) != npar & length(upper)!= 1)
    stop("length of 'upper' should be either 1 or equal to number of parameters")

  ## Initialisation
  con <- list(npop     = max(5*length(p), 50), # nr elements in population
              numiter  = 10000,                # number of iterations
              centroid = 3,                    # number of points in centroid
              varleft  = 1e-8,                 # relative variation upon stopping
              verbose  = FALSE)
  nmsC <- names(con)

  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]) > 0)
    stop("unknown names in control: ", paste(noNms, collapse = ", "))

  npop     <- con$npop
  numiter  <- con$numiter
  centroid <- con$centroid
  varleft  <- con$varleft

  cost     <- function (par) f(par, ...)
  tiny     <- 1e-8
  varleft  <- max(tiny,varleft)
  rsstrace <- NULL

  populationpar <- matrix(nrow = npop, ncol = npar, byrow = TRUE,
    data = lower + runif(npar*npop) * rep((upper - lower), npop))
  colnames(populationpar) <- names(p)
  populationpar[1,] <- p

  populationcost <- apply(populationpar, FUN = cost, MARGIN = 1)
  iworst         <- which.max(populationcost)
  worstcost      <- populationcost[iworst]

  ## Hybridisation phase
  iter <- 0
  lastbest  <- -Inf
  while (iter < numiter && (max(populationcost) - min(populationcost))
                                  > (min(populationcost)*varleft)) {
    iter <- iter + 1

    selectpar <- sample(1:npop, size = centroid)  # for cross-fertilisation
    mirrorpar <- sample(1:npop, size = 1)         # for mirroring

    newpar    <- colMeans(populationpar[selectpar,])    # centroid
    newpar    <- 2*newpar - populationpar[mirrorpar,]   # mirroring

    newpar    <- pmin(pmax(newpar, lower), upper)

    newcost   <- cost(newpar)
    if(!is.nan(newcost) & !is.na(newcost) & !is.null(newcost)) {
      if (newcost < worstcost) {
        populationcost[iworst]  <- newcost
        populationpar [iworst,] <- newpar
        iworst     <- which.max(populationcost) # new worst member
        worstcost  <- populationcost[iworst]
        bestcost   <- min(populationcost)
        if (bestcost != lastbest)
          rsstrace  <- rbind(rsstrace, c(iter, min(populationcost)))
        lastbest <- bestcost
      }
    }
  } # end while loop

  ibest    <- which.min(populationcost)
  bestpar  <- populationpar[ibest,]
  bestcost <- populationcost[ibest]

  res <- list(par = bestpar, cost = bestcost, iterations = iter)
  if (con$verbose) {
    res$poppar   <- populationpar
    res$popcost  <- populationcost
    res$rsstrace <- rsstrace
  }
  return (res)
}

