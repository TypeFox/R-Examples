DEoptim.control <- function(VTR = -Inf, strategy = 2, bs = FALSE, NP = 50,
                            itermax = 200, CR = 0.5, F = 0.8, trace = TRUE,
                            initialpop = NULL, storepopfrom = itermax + 1,
                            storepopfreq = 1, checkWinner = FALSE,
                            avWinner = TRUE, p = 0.2) {
  if (itermax <= 0) {
    warning("'itermax' <= 0; set to default value 200\n", immediate. = TRUE)
    itermax <- 200
  }
  if (NP < 4) {
    warning("'NP' < 4; set to default value 50\n", immediate. = TRUE)
    NP <- 50
  }
  if (F < 0 | F > 2) {
    warning("'F' not in [0,2]; set to default value 0.8\n", immediate. = TRUE)
    F <- 0.8
  }
  if (CR < 0 | CR > 1) {
    warning("'CR' not in [0,1]; set to default value 0.5\n", immediate. = TRUE)
    CR <- 0.5
  }
  if (strategy < 1 | strategy > 6) {
    warning("'strategy' not in {1,...,6}; set to default value 2\n",
            immediate. = TRUE)
    strategy <- 2
  }

  bs <- (bs > 0)

  if ( trace < 0 ) {
    warning("'trace' cannot be negative; set to 'TRUE'")
    trace <- TRUE
  }

  storepopfreq <- floor(storepopfreq)
  if (storepopfreq > itermax)
    storepopfreq <- 1

  if (p <= 0 || p > 1) {
    warning("'p' not in (0,1]; set to default value 0.2\n", immediate. = TRUE)
    p <- 0.2
  }

  list(VTR = VTR, strategy = strategy, NP = NP, itermax = itermax, CR
       = CR, F = F, bs = bs, trace = trace, initialpop = initialpop,
       storepopfrom = storepopfrom, storepopfreq = storepopfreq,
       checkWinner = checkWinner, avWinner = avWinner, p = p)
}

DEoptim <- function(fn, lower, upper, control = DEoptim.control(), ...) {
  ##fn1  <- function(par) fn(par, ...)
  if (length(lower) != length(upper))
    stop("'lower' and 'upper' are not of same length")
  if (!is.vector(lower))
    lower <- as.vector(lower)
  if (!is.vector(upper))
    upper <- as.vector(upper)
  if (any(lower > upper))
    stop("'lower' > 'upper'")
  if (any(lower == "Inf"))
    warning("you set a component of 'lower' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(lower == "-Inf"))
    warning("you set a component of 'lower' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(upper == "Inf"))
    warning("you set a component of 'upper' to 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (any(upper == "-Inf"))
    warning("you set a component of 'upper' to '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  if (!is.null(names(lower)))
    nam <- names(lower)
  else if (!is.null(names(upper)) & is.null(names(lower)))
    nam <- names(upper)
  else
    nam <- paste("par", 1:length(lower), sep = "")

  env <- new.env()

  ctrl <- do.call(DEoptim.control, as.list(control))
  ctrl$npar <- length(lower)
  if (ctrl$NP < 4) {
    warning("'NP' < 4; set to default value 50\n", immediate. = TRUE)
    ctrl$NP <- 50
  }
  if (ctrl$NP < 10*length(lower))
    warning("For many problems it is best to set 'NP' (in 'control') to be at least ten times the length of the parameter vector. \n", immediate. = TRUE)
  if (!is.null(ctrl$initialpop)) {
    ctrl$specinitialpop <- TRUE
    if(!identical(as.numeric(dim(ctrl$initialpop)), c(ctrl$NP, ctrl$npar)))
      stop("Initial population is not a matrix with dim. NP x length(upper).")
  }
  else {
    ctrl$specinitialpop <- FALSE
    ctrl$initialpop <- matrix(0,1,1)    # dummy matrix
  }
  ##
  ctrl$trace <- as.numeric(ctrl$trace)
  ctrl$specinitialpop <- as.numeric(ctrl$specinitialpop)

  outC <- DEoptim_impl(lower, upper, fn, ctrl, env)
  ##
  if (length(outC$storepop) > 0) {
    nstorepop <- floor((outC$iter - ctrl$storepopfrom) / ctrl$storepopfreq)
    ## storepop <- list()
    ## cnt <- 1
    ## for(i in 1:nstorepop) {
    ##   idx <- cnt:((cnt - 1) + (ctrl$NP * ctrl$npar))
    ##   storepop[[i]] <- matrix(outC$storepop[idx], nrow = ctrl$NP, ncol = ctrl$npar,
    ##                      byrow = TRUE)
    ##   cnt <- cnt + (ctrl$NP * ctrl$npar)
    ##   dimnames(storepop[[i]]) <- list(1:ctrl$NP, nam)
    ## }
    storepop <- outC$storepop[1:nstorepop]
    for (i in 1:length(storepop)) dimnames(storepop[[i]]) <- list(1:ctrl$NP, nam)
  }
  else {
    storepop = NULL
  }

  ## optim
  bestmem <- as.numeric(outC$bestmem)
  names(bestmem) <- nam
  bestval <- as.numeric(outC$bestval)
  nfeval <- as.numeric(outC$nfeval)
  iter <- as.numeric(outC$iter)

  ## member
  names(lower) <- names(upper) <- nam
  #bestmemit <- matrix(outC$bestmemit, nrow = iter, ncol = ctrl$npar, byrow = TRUE)
  bestmemit <- outC$bestmemit

  dimnames(bestmemit) <- list(1:iter, nam)
  bestvalit <- as.numeric(outC$bestvalit[1:iter])
  #pop <- matrix(outC$pop, nrow = ctrl$NP, ncol = ctrl$npar, byrow = TRUE)
  pop <- outC$pop
  storepop <- as.list(storepop)

  outR <- list(optim = list(
              bestmem = bestmem,
              bestval = bestval,
              nfeval = nfeval,
              iter = iter),
            member = list(
              lower = lower,
              upper = upper,
              bestmemit = bestmemit,
              bestvalit = bestvalit,
              pop = pop,
              storepop = storepop))

  attr(outR, "class") <- "DEoptim"
  return(outR)
}

