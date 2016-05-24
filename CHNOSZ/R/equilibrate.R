# CHNOSZ/equilibrate.R
# functions to calculation logarithm of activity
# of species in (metastable) equilibrium

equilibrate <- function(aout, balance=NULL, loga.balance=NULL, 
  ispecies=1:length(aout$values), normalize=FALSE, as.residue=FALSE,
  method=c("boltzmann", "reaction")) {
  ### set up calculation of equilibrium activities of species from the affinities 
  ### of their formation reactions from basis species at known activities
  ### split from diagram() 20120925 jmd
  ## number of possible species
  nspecies <- length(aout$values)
  ## get the balancing coefficients
  n.balance <- balance(aout, balance)
  ## take selected species in 'ispecies'
  if(length(ispecies)==0) stop("the length of ispecies is zero")
  # take out species that have NA affinities
  ina <- sapply(aout$values, function(x) all(is.na(x)))
  ispecies <- ispecies[!ina[ispecies]]
  if(length(ispecies)==0) stop("all species have NA affinities")
  if(!identical(ispecies, 1:nspecies)) {
    msgout(paste("equilibrate: using", length(ispecies), "of", nspecies, "species\n"))
    aout$species <- aout$species[ispecies, ]
    aout$values <- aout$values[ispecies]
    n.balance <- n.balance[ispecies]
  }
  ## number of species that are left
  nspecies <- length(aout$values)
  ## say what the balancing coefficients are
  if(length(n.balance) < 100) msgout(paste("equilibrate: n.balance is", c2s(n.balance), "\n"))
  ## logarithm of total activity of the balance
  if(is.null(loga.balance)) {
    # sum up the activities, then take absolute value
    # in case n.balance is negative
    sumact <- abs(sum(10^aout$species$logact * n.balance))
    loga.balance <- log10(sumact)
  }
  msgout(paste0("equilibrate: loga.balance is ", loga.balance, "\n"))
  ## normalize -- normalize the molar formula by the balance coefficients
  m.balance <- n.balance
  isprotein <- grepl("_", as.character(aout$species$name))
  if(normalize | as.residue) {
    if(any(n.balance < 0)) stop("one or more negative balancing coefficients prohibit using normalized molar formulas")
    n.balance <- rep(1, nspecies)
    if(as.residue) msgout(paste("equilibrate: using 'as.residue' for molar formulas\n"))
    else msgout(paste("equilibrate: using 'normalize' for molar formulas\n"))
  } else m.balance <- rep(1, nspecies)
  ## Astar: the affinities/2.303RT of formation reactions with
  ## formed species in their standard-state activities
  Astar <- lapply(1:nspecies, function(i) { 
    # 'starve' the affinity of the activity of the species,
    # and normalize the value by the nor-molar ratio
    (aout$values[[i]] + aout$species$logact[i]) / m.balance[i] 
  })
  ## chose a method and compute the equilibrium activities of species
  if(missing(method)) {
    if(all(n.balance==1)) method <- method[1]
    else method <- method[2]
  }
  msgout(paste("equilibrate: using", method[1], "method\n"))
  if(method[1]=="boltzmann") loga.equil <- equil.boltzmann(Astar, n.balance, loga.balance)
  else if(method[1]=="reaction") loga.equil <- equil.reaction(Astar, n.balance, loga.balance)
  ## if we normalized the formulas, get back to activities to species
  if(normalize & !as.residue) {
    loga.equil <- lapply(1:nspecies, function(i) {
      loga.equil[[i]] - log10(m.balance[i])
    })
  }
  ## put together the output
  out <- c(aout, list(balance=balance, m.balance=m.balance, n.balance=n.balance,
    loga.balance=loga.balance, Astar=Astar, loga.equil=loga.equil))
  # done!
  return(out)
}

equil.boltzmann <- function(Astar, n.balance, loga.balance) {
  # 20090217 new "abundance" function
  # return equilibrium logarithms of activity of species
  # works using Boltzmann distribution
  # A/At = e^(Astar/n.balance) / sum(e^(Astar/n.balance))
  # where A is activity of the ith residue and
  # At is total activity of residues
  # advantages over equil.reaction
  # 2) loops over species only - much faster
  # 3) no root finding - those games might fail at times
  # disadvantage:
  # 1) only works for per-residue reactions
  # 2) can create NaN logacts if the Astars are huge/small
  if(any(n.balance!=1)) stop("won't run equil.boltzmann for balance <> 1")
  # initialize output object
  A <- Astar
  # remember the dimensions of elements of Astar (could be vector,matrix)
  Astardim <- dim(Astar[[1]])
  Anames <- names(Astar)
  # first loop: make vectors
  A <- palply("", 1:length(A), function(i) as.vector(A[[i]]))
  # second loop: get the exponentiated Astars (numerators)
  # need to convert /2.303RT to /RT
  #A[[i]] <- exp(log(10)*Astar[[i]]/n.balance[i])/n.balance[i]
  A <- palply("", 1:length(A), function(i) exp(log(10)*Astar[[i]]/n.balance[i]))
  # third loop: accumulate the denominator
  # initialize variable to hold the sum
  At <- A[[1]]
  At[] <- 0
  for(i in 1:length(A)) At <- At + A[[i]]*n.balance[i]
  # fourth loop: calculate log abundances and replace the dimensions
  A <- palply("", 1:length(A), function(i) loga.balance + log10(A[[i]]/At))
  # fifth loop: replace dimensions
  for(i in 1:length(A)) dim(A[[i]]) <- Astardim
  # add names and we're done!
  names(A) <- Anames
  return(A)
}

equil.reaction <- function(Astar, n.balance, loga.balance) {
  # to turn the affinities/RT (A) of formation reactions into 
  # logactivities of species (logact(things)) at metastable equilibrium
  # 20090217 extracted from diagram and renamed to abundance.old
  # 20080217 idea / 20120128 cleaned-up strategy
  # for any reaction stuff = thing,
  #   A = logK - logQ 
  #     = logK - logact(thing) + logact(stuff)
  # given Astar = A + logact(thing),
  # given Abar = A / n.balance,
  #   logact(thing) = Astar - Abar * n.balance  [2]
  # where n.balance is the number of the balanced quanitity
  # (conserved component) in each species
  # equilibrium values of logact(thing) satifsy:
  # 1) Abar is equal for all species
  # 2) log10( sum of (10^logact(thing) * n.balance) ) = loga.balance  [1]
  #
  # because of the logarithms, we can't solve the equations directly
  # instead, use uniroot() to compute Abar satisfying [1]

  # we can't run on one species
  if(length(Astar)==1) stop("at least two species needed for reaction-based equilibration")
  # remember the dimensions and names
  Adim <- dim(Astar[[1]])
  Anames <- names(Astar)
  # make a matrix out of the list of Astar
  Astar <- list2array(lapply(Astar, c))
  # that produces the same result (other than colnames) and is much faster than
  #Astar <- sapply(Astar, c)  
  # also, latter has NULL nrow for length(Astar[[x]])==1
  # some function definitions
  # to calculate log of activity of balanced quantity from logact(thing) of all species [1]
  logafun <- function(logact) log10(sum(10^logact * n.balance))
  # to calculate logact(thing) from Abar for the ith condition [2]
  logactfun <- function(Abar, i) Astar[i,] - Abar * n.balance
  # to calculate difference between logafun and loga.balance for the ith condition
  logadiff <- function(Abar, i) loga.balance - logafun(logactfun(Abar, i))
  # to calculate a range of Abar that gives negative and positive values of logadiff for the ith condition
  Abarrange <- function(i) {
    # starting guess of Abar (min/max) from range of Astar / n.balance
    Abar.range <- range(Astar[i, ] / n.balance)
    # diff(Abar.range) can't be 0 (dlogadiff.dAbar becomes NaN)
    if(diff(Abar.range)==0) {
      Abar.range[1] <- Abar.range[1] - 0.1
      Abar.range[2] <- Abar.range[2] + 0.1
    }
    # the range of logadiff
    logadiff.min <- logadiff(Abar.range[1], i)
    logadiff.max <- logadiff(Abar.range[2], i)
    # we're out of luck if they're both infinite
    if(is.infinite(logadiff.min) & is.infinite(logadiff.max))
      stop("FIXME: there are no initial guesses for Abar that give
        finite values of the differences in logarithm of activity
        of the conserved component")
    # if one of them is infinite we might have a chance
    if(is.infinite(logadiff.min)) {
      # decrease the Abar range by increasing the minimum
      Abar.range[1] <- Abar.range[1] + 0.99 * diff(Abar.range)
      logadiff.min <- logadiff(Abar.range[1], i)
      if(is.infinite(logadiff.min)) stop("FIXME: the second initial guess for Abar.min failed")
    }
    if(is.infinite(logadiff.max)) {
      # decrease the Abar range by decreasing the maximum
      Abar.range[2] <- Abar.range[2] - 0.99 * diff(Abar.range)
      logadiff.max <- logadiff(Abar.range[2], i)
      if(is.infinite(logadiff.max)) stop("FIXME: the second initial guess for Abar.max failed")
    }
    iter <- 0
    while(logadiff.min > 0 | logadiff.max < 0) {
      # the change of logadiff with Abar
      # it's a weighted mean of the n.balance
      dlogadiff.dAbar <- (logadiff.max - logadiff.min) / diff(Abar.range)
      # change Abar to center logadiff (min/max) on zero
      logadiff.mean <- mean(c(logadiff.min, logadiff.max))
      Abar.range <- Abar.range - logadiff.mean / dlogadiff.dAbar
      # one iteration is enough for the examples in the package
      # but there might be a case where the range of logadiff doesn't cross zero
      # (e.g. for the carboxylic acid example in revisit.Rd)
      logadiff.min <- logadiff(Abar.range[1], i)
      logadiff.max <- logadiff(Abar.range[2], i)
      iter <- 1
      if(iter > 5) {
        stop("FIXME: we seem to be stuck! This function (Abarrange() in
          equil.reaction()) can't find a range of Abar such that the differences
          in logarithm of activity of the conserved component cross zero")
      }
    }
    return(Abar.range)
  }
  # to calculate an equilibrium Abar for the ith condition
  Abarfun <- function(i) {
    # get limits of Abar where logadiff brackets zero
    Abar.range <- Abarrange(i)
    # now for the real thing: uniroot!
    Abar <- uniroot(logadiff, interval=Abar.range, i=i)$root
    return(Abar)
  }
  # calculate the logact(thing) for each condition
  logact <- palply("", 1:nrow(Astar), function(i) {
    # get the equilibrium Abar for each condition
    Abar <- Abarfun(i)
    return(logactfun(Abar, i))
  }) 
  # restore the dimensions and names
  if(length(Adim)==1) logact <- list2array(logact)
  else logact <- sapply(logact, c)
  logact <- lapply(1:nrow(logact), function(i) {
    thisla <- list(logact[i,])[[1]]
    dim(thisla) <- Adim
    return(thisla)
  })
  names(logact) <- Anames
  # all done!
  return(logact)
}

balance <- function(aout, balance=NULL) {
  ## generate n.balance from user-given or automatically identified basis species
  ## extracted from equilibrate() 20120929
  # 'balance' can be:
  #   NULL                    - autoselect using which.balance
  #   name of basis species   - balanced on this basis species
  #   "length"                   - balanced on sequence length of proteins 
  #                             (default if balance is missing and all species are proteins)
  #   1                       - balanced on one mole of species
  #   numeric vector          - user-defined n.balance
  #   "volume"                - standard-state volume listed in thermo$obigt
  # the index of the basis species that might be balanced
  ibalance <- numeric()
  # deal with proteins
  isprotein <- grepl("_", as.character(aout$species$name))
  if(is.null(balance) & all(isprotein)) balance <- "length"
  # try to automatically find a balance
  if(is.null(balance)) {
    ibalance <- which.balance(aout$species)
    # no shared basis species - balance on one mole of species
    if(length(ibalance) == 0) balance <- 1
  } 
  if(is.numeric(balance[1])) {
    # a numeric vector
    n.balance <- rep(balance, length.out=length(aout$values))
    msgout("balance: from numeric argument value\n")
  } else {
    # "length" for balancing on protein length
    if(identical(balance, "length")) {
      if(!all(isprotein)) stop("length was the requested balance, but some species are not proteins")
      n.balance <- protein.length(aout$species$name)
      msgout("balance: from protein length\n")
    } else if(identical(balance, "volume")) {
      n.balance <- info(aout$species$ispecies, check.it=FALSE)$V
      msgout("balance: from volume")
    } else {
      # is the balance the name of a basis species?
      if(length(ibalance)==0) {
        ibalance <- match(balance, rownames(aout$basis))
        if(is.na(ibalance)) stop("basis species (", balance, ") not available to balance transformations")
      }
      # the name of the basis species (need this if we got ibalance which which.balance, above)
      balance <- colnames(aout$species)[ibalance[1]]
      msgout(paste("balance: from moles of", balance, "in formation reactions\n"))
      # the balance vector
      n.balance <- aout$species[, ibalance[1]]
      # we check if that all formation reactions contain this basis species
      if(any(n.balance==0)) stop("some species have no ", balance, " in the formation reaction")
    }
  }
  return(n.balance)
}
