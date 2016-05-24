probabilities.function <- function(hypothesis, cons, doR=FALSE) {
  # Creates a probability function.
  #
  # The probability function computes all probabilities associated with a
  # specific set of genotypes. 
  #
  # In practice, this function makes is easy to substitute C vs R
  # implementations, as well as specialize the C implementations. 

  if((!hypothesis$doDropin) && !doR) {
    probabilities.no.dropin.C <- function(res, vDoseDropout, csp, unc, ...) {
        if(length(res) != ncol(vDoseDropout))
          stop("output vector and vDoseDropout have incompatible sizes.")
        if(length(csp) != nrow(vDoseDropout))
          stop("csp and vDoseDropout have incompatible sizes.")
        if(length(unc) != nrow(vDoseDropout))
          stop("unc and vDoseDropout have incompatible sizes.")
        if(any(dim(vDoseDropout) != dim(cons$zeroAll)))
          stop("expected vDoseDropout and zeroAll to match.")
        .Call(.cpp.probabilitiesNoDropin, res, vDoseDropout, !csp & !unc, csp,
              cons$zeroAll, PACKAGE="likeLTD") 
    }
    return(probabilities.no.dropin.C) 
  } else if(hypothesis$doDropin && !doR) {
    probabilities.with.dropin.C <- function(res, vDoseDropout, csp, unc, rate) {
       if(length(res) != ncol(vDoseDropout))
         stop("output vector and vDoseDropout have incompatible sizes.")
       if(length(csp) != nrow(vDoseDropout))
         stop("csp and vDoseDropout have incompatible sizes.")
       if(length(unc) != nrow(vDoseDropout))
         stop("unc and vDoseDropout have incompatible sizes.")
       if(any(dim(vDoseDropout) != dim(cons$zeroAll)))
         stop("expected vDoseDropout and zeroAll to match.")
       if(any(dim(cons$freqMat) != dim(cons$zeroAll)))
         stop("expected vDoseDropout and zeroAll to match.")
       .Call(.cpp.probabilitiesWithDropin, res, vDoseDropout, !csp & !unc, csp,
             cons$zeroAll, cons$freqMat, rate, PACKAGE="likeLTD")
     }
    return(probabilities.with.dropin.C)
  }
  # Otherwise, return an R function. 
  function(res, vDoseDropout, csp, unc, rate) {
     res <- res                                             *
            selective.col.prod(!csp & !unc, vDoseDropout)   *
            selective.col.prod(csp, 1 - vDoseDropout)  
    if(hypothesis$doDropin) 
      res <- res                                                           *
             selective.col.prod(csp & cons$zeroAll, rate * cons$freqMat)   *
             selective.col.prod(!csp & !unc & cons$zeroAll,
                                1 - rate * cons$freqMat)
    res
  }
}


empty.alleles = function(genotypes, dropoutPresence, nUnknowns) {
  # Mask over genotype matrix indicating empty/null elements.
  #
  # Parameters:
  #   genotypes: Matrix of indices into dosage matrix.
  #   dropoutPresence: Profiles with dropout.
  # Returns: A matrix indicating which alleles will be NA in dosage.
  if(!is.matrix(dropoutPresence)) stop("Expected matrix in input.")
  if(ncol(genotypes) == 0) return(matrix(nrow=nrow(dropoutPresence), ncol=0))
  if(nrow(dropoutPresence) != 0) knownZero = rowSums(dropoutPresence) == 0
  else knownZero = array(FALSE, ncol(genotypes))
  if(nUnknowns == 0) 
    return(matrix(knownZero, nrow=length(knownZero), ncol=ncol(genotypes)))
  # Performs following call in C:
  # iota = 1:length(knownZero)
  # apply(genotypes, 2, function(n) (!iota %in% n) & knownZero) 
  # A new object is returned. 
  .Call(.cpp.emptyAlleles, genotypes, knownZero, PACKAGE="likeLTD")
}

relatedness.factors <- function(input, genotypes, alleleDb, queriedProfile,
                                relatedness=c(0, 0)) {
  # Creates relatedness factor 
  #
  # The relatedness factor takes into account the ancestry of the queried
  # subject and modifies allele frequencies accordingly.

  if(any(abs(relatedness) > 1e-8)) {
    indices = c(which(rownames(alleleDb) %in% queriedProfile[[1]][[1]]), 
                which(rownames(alleleDb) %in% queriedProfile[[1]][[2]]) )
    frequencies = alleleDb[indices, 1]
    
    .Call(.cpp.relatednessFactors, input, relatedness, genotypes, indices,
          frequencies, PACKAGE="likeLTD")
  }
  input 
  # defactored = function(n) {
  #   result = 1.0 - sum(relatedness)
  #   hasFirst = n %in% indices[[1]]
  #   hasSecond = n %in% indices[[2]]
  #   if(any(hasFirst)) {
  #     factor = relatedness[[1]] * 0.5 / alleleDb[indices[[1]], 1]
  #     if(!all(hasFirst)) factor = factor * 0.5
  #     result = result + factor
  #   }
  #   if(any(hasSecond)) {
  #     factor = relatedness[[1]] * 0.5 / alleleDb[indices[[2]], 1]
  #     if(!all(hasSecond)) factor = factor * 0.5
  #     result = result + factor
  #   }
  #   if(any(hasFirst) & any(hasSecond)) {
  #     factor = relatedness[[2]] / alleleDb[indices[[1]], 1] /
  #              alleleDb[indices[[2]], 1] 
  #     if(indices[[1]] != indices[[2]]) factor = factor * 0.5
  #     result = result + factor
  #   }
  #   return(result)
  # }
  # apply(genotypes[1:2, ], 2, defactored) 
}

genotype.factors <- function(genotypes, alleleDb, nUnknowns, doDropin,
                             queriedProfile, relatedness=c(0, 0)) {
  # Fraction cum heterozygote factors.
  #

  if(nUnknowns == 0 & !doDropin) {
    het = 1
    fractions <- matrix(alleleDb[genotypes, 1], ncol=ncol(genotypes))
    result <- apply(fractions, 2, prod) * het
  }
  else {
    # The following is performed directly in C. 
    # n = nrow(genotypes) / 2
    # odd = 2*1:n - 1; even = 2*1:n
    # het <- 1 + matrix(genotypes[odd, ] < genotypes[even, ], nrow=length(odd))
    # het <- apply(het, 2, prod)
    # fractions <- matrix(alleleDb[genotypes, 1], ncol=ncol(genotypes))
    # result <- apply(fractions, 2, prod) * het
    result <- .Call(.cpp.fractionsAndHet, genotypes, alleleDb[, 1],
                    PACKAGE="likeLTD")
  }

  relatedness.factors(result, genotypes, alleleDb, queriedProfile, relatedness)
}

likelihood.constructs.per.locus = function(hypothesis) {
  # Creates the locus-specific data needed by each likehood function.
  #
  # Parameters:
  #   hypothesis: A hypothesis, for instance one returned by
  #               prosecution.hypothesis(...) or defence.hypothesis(...)
  alleles = rownames(hypothesis$alleleDb)
  if(is.null(alleles)) stop("Could not figure out alleles names.")
  alleles.vector = function(n) alleles %in% unlist(n)
  cspPresence     = apply(hypothesis$cspProfile, 1, alleles.vector)
  dropoutPresence = apply(hypothesis$dropoutProfs, 1, alleles.vector)
  uncPresence     = apply(hypothesis$uncProf, 1, alleles.vector)
  if(!is.matrix(cspPresence))
    cspPresence = matrix(ncol=0, nrow=length(alleles))
  if(!is.matrix(dropoutPresence))
    dropoutPresence = matrix(ncol=0, nrow=length(alleles))
  if(!is.matrix(uncPresence))
    uncPresence = matrix(ncol=0, nrow=length(alleles))

  missingReps = apply(hypothesis$cspProfile, 1, is.na)

  genotypes <- compatible.genotypes(cspPresence, dropoutPresence, alleles,
                                    hypothesis$nUnknowns, hypothesis$doDropin,
                                    missingReps)
  zeroAll = empty.alleles(genotypes, dropoutPresence, hypothesis$nUnknowns) 
  # Only take into account relatedness between Q and X unedr Hd 
  if(hypothesis$hypothesis=="defence")
	{
  	factors = genotype.factors(genotypes, hypothesis$alleleDb,
                             hypothesis$nUnknowns, hypothesis$doDropin,
                             hypothesis$queriedProfile,
                             hypothesis$relatedness) 
	} else {
  	factors = genotype.factors(genotypes, hypothesis$alleleDb,
                             hypothesis$nUnknowns, hypothesis$doDropin,
                             hypothesis$queriedProfile,
                             c(0,0))
	}

  list(cspPresence=cspPresence, dropoutPresence=dropoutPresence,
       uncPresence=uncPresence, missingReps=missingReps,
       genotypes=genotypes, zeroAll=zeroAll, factors=factors,
       freqMat=hypothesis$alleleDb[, 1])
}

create.likelihood.per.locus <- function(hypothesis, addAttr=FALSE, likeMatrix = FALSE) {
  # Creates a likelyhood function for a given hypothesis and locus
  #
  # A hypothesis is given by the number of unknown contributors, whether to model
  # dropin, so on and so forth.

  cons = likelihood.constructs.per.locus(hypothesis)
  doR = !is.null(hypothesis$doR) && hypothesis$doR == TRUE
  probabilities = probabilities.function(hypothesis, cons, doR=doR)

  result.function <- function(locusAdjustment, power, dropout,
                              degradation=NULL, rcont=NULL, dropin=NULL, ...) {
    # Likelyhood function for a given hypothesis and locus
    #
    # This function is specific to the hypothesis for which it was created.
    # It hides everything except the nuisance parameters over which to
    # optimize.
    #
    # Parameters:
    #   locusAdjustment: a scalar floating point value.
    #   power: a scalar floating point value.
    #   dropout: the dropout rate for each replicate.
    #   degradation: relative degradation from each profiled individual in this
    #                hypothesis
    #   rcont: relative contribution from each profiled individual and unknown
    #          contributor in this hypothesis If it is one less than that, then
    #          an additional 1 is inserted at a position given by
    #          hypothesis$refIndiv. In general, this means either the queried
    #          individual (if subject to dropout and prosecution hypothesis) or
    #          the first individual subject to droput is the reference individual.
    #   ...: Any other parameter, e.g. for penalty functions. 
    #        These parameters are ignored here.
    # Returns: A scalar value giving the likelihood for this locus and
    #          hypothesis
    if(is.null(degradation)) degradation = c()
    if(is.null(rcont)) rcont = c()
    if(length(rcont) + 1 == hypothesis$nUnknowns
                            + ncol(cons$dropoutPresence)) {
      if(hypothesis$refIndiv == 1) rcont = c(1, rcont)
      else if(hypothesis$refIndiv > length(rcont)) rcont = c(rcont, 1)
      else rcont = c(rcont[1:hypothesis$refIndiv-1], 1,
                     rcont[hypothesis$refIndiv:length(rcont)])
    }
    if(length(rcont) != hypothesis$nUnknowns + ncol(cons$dropoutPresence))
      stop(sprintf("rcont should be %d long.",
                   hypothesis$nUnknowns + ncol(cons$dropoutPresence)))
    if(length(degradation) != hypothesis$nUnknowns +
                              ncol(cons$dropoutPresence))
      stop(sprintf("degradation should be %d long.",
                   hypothesis$nUnknowns + ncol(cons$dropoutPresence)))
    if(length(dropout) != length(cons$missingReps))
      stop(sprintf("dropout should be %d long.", ncol(cons$dropoutPresence)))
    if(any(rcont < 0)) stop("found negative relative contribution.")
    if(any(degradation < 0)) stop("found negative degradation parameter.")
    if(hypothesis$doDropin && is.null(dropin)) 
      stop("Model requires missing argument 'dropin'")
    else if(is.null(dropin)) dropin = 0
    if(hypothesis$doDropin && dropin < 0) 
      stop("Dropin rate should be between 0 and 1 (included).")
    if(power >= 0)
      stop("Power parameter cannot be positive or null.")
    if(any(dropout < 0) || any(dropout > 1)) 
      stop("Dropout rates must be between 0 and 1 (included).")
    if(length(locusAdjustment) != 1)
      stop("locusAdjustment should be a scalar")
    if(locusAdjustment < 0) stop("locusAdjustment must be positive.")

    allEPG <- all.epg.per.locus(rcont, degradation, cons$dropoutPresence,
                                hypothesis$alleleDb[, 2], cons$genotypes,
                                hypothesis$nUnknowns > 0)
    # Goes to C to do exponentiation. This is quite expensive an operation.
    # Doing it in C is only interesting for larger matrices and when compiled
    # with OpenMP. 
    # allEPG = (allEPG * locusAdjustment)^power
    # Nothing is returned by this function. Works directly on allEPG. 
    .Call(.cpp.powerAdjustment, allEPG, cons$zeroAll, locusAdjustment,
          power, PACKAGE="likeLTD")
          
    # If any allEPG are infinite change to zero
    # so unlikely that probability might as well be zero
    allEPG[which(is.infinite(allEPG))] = 0

    # res: (Partial) Likelihood per allele.
    res = array(1, length(cons$factors))
    # Loop over replicates.
    for(i in 1:length(cons$missingReps)) {
      if(!cons$missingReps[i]) {
        csp = cons$cspPresence[, i]
        unc = cons$uncPresence[, i]

        # Goes to C to do fraction. Avoids intermediate step.
        # For smaller matrices, the C code seems a bit slower. For larger
        # matrices, it can be twice as fast. 
        # vDoseDropout = allEPG * dropout[i]
        # vDoseDropout = vDoseDropout / (vDoseDropout + 1 - dropout[i])
        vDoseDropout <- .Call(.cpp.doseFraction, allEPG, cons$zeroAll, dropout[[i]],
                              PACKAGE="likeLTD")
        # When there is no reference individual the dropin rate is just the dropin probability
	      if(nrow(hypothesis$dropoutProfs)+hypothesis$nUnknowns==0) 
		      {		
		      dirate = dropin
		      } else {
		      dirate = dropin * (1 - dropout[i])
		      }
        res <- probabilities(res, vDoseDropout, csp, unc,
                             dirate)
      } # End of if(any(csp)) 
    } # End of loop over replicates.

    # Figure out likelihood for good and return.
    if(likeMatrix==FALSE)
	{
    	return(sum(res * cons$factors))
	} else {
	return(res * cons$factors)
	}
  }

  if(addAttr) {
    attr(result.function, "hypothesis") <- hypothesis
    attr(result.function, "constructs") <- cons
  }
  result.function
}


# Penalties to apply to the likelihood.
# Documentation is in man directory.
penalties <- function(locusAdjustment, power, dropout, degradation=NULL,
                      rcont=NULL, dropin=NULL, locusAdjPenalty=50,
                      dropinPenalty=2, degradationPenalty=50, bemn=-4.35,
                      besd=0.38, ...) {
  result = 1
  # Normalizes by number of loci so product of penalties same as in old code.
  # Some penalties are per locus and other for all locus. Hence this bit of run
  # around.
  normalization = 1.0 / max(length(locusAdjustment), 1)

  if(!missing(dropin) & !is.null(dropin))
    result = result * exp(-dropin * dropinPenalty * normalization)
  if(!missing(degradation) & !is.null(degradation))
    result = result * exp(-sum(degradation) * degradationPenalty 
                           * normalization )
  if(!missing(power) & !is.null(power))
    result = result * dnorm(power, bemn, besd) ^ normalization
  if(!missing(locusAdjustment) & !is.null(locusAdjustment))
    result = result * dgamma(as.vector(unlist(locusAdjustment)),
                             locusAdjPenalty, locusAdjPenalty)

  return(result)
}

# Creates a likelihood function from the input hypothesis
# Documentation is in man directory.
create.likelihood.vectors <- function(hypothesis, addAttr=FALSE, likeMatrix=FALSE, ...) {

  hypothesis = add.args.to.hypothesis(hypothesis, ...)
  sanity.check(hypothesis) # makes sure hypothesis has right type.
  locusCentric = transform.to.locus.centric(hypothesis)
  functions <- mapply(create.likelihood.per.locus, locusCentric,
                      MoreArgs=list(addAttr=addAttr, likeMatrix=likeMatrix))

  if(is.null(hypothesis$locusAdjPenalty)) hypothesis$locusAdjPenalty = 50 
  if(is.null(hypothesis$dropinPenalty)) hypothesis$dropinPenalty = 2 
  if(is.null(hypothesis$degradationPenalty)) hypothesis$degradationPenalty = 50 
  if(is.null(hypothesis$bemn)) hypothesis$bemn = -4.35
  if(is.null(hypothesis$besd)) hypothesis$besd = 0.38 

  likelihood.vectors <- function(locusAdjustment, power, dropout,
                                 degradation=NULL, rcont=NULL, dropin=NULL,
                                 locusAdjPenalty=hypothesis$locusAdjPenalty, 
				 dropinPenalty=hypothesis$dropinPenalty,
                                 degradationPenalty=hypothesis$degradationPenalty, 
				 bemn=hypothesis$bemn,
                                 besd=hypothesis$besd, ...) {
    # Call each and every function in the array.
    arguments = list(power=power, dropout=dropout,
                     degradation=degradation, rcont=rcont,
                     dropinPenalty=dropinPenalty,
                     degradationPenalty=degradationPenalty, bemn=bemn,
                     besd=besd, dropin=dropin)
    callme <- function(objective, adj) {
      args = append(arguments, list(locusAdjustment=adj))
      do.call(objective, args)
    }
    if(length(locusAdjustment) == 1)
      locusAdjustment = rep(locusAdjustment, length(functions))
    if(setequal(names(locusAdjustment), colnames(hypothesis$cspProfile)))
      locusAdjustment <- locusAdjustment[colnames(hypothesis$cspProfile)]
    objectives = mapply(callme, functions, locusAdjustment)
    arguments = append(arguments, list(locusAdjustment=locusAdjustment))
    arguments = append(arguments, list(...))
    pens <- do.call(penalties, arguments)
    list(objectives=objectives, penalties=pens)
  }
  if(addAttr) {
    attr(likelihood.vectors, "hypothesis") <- hypothesis
    attr(likelihood.vectors, "functions") <- functions
  }
  return(likelihood.vectors)
}


# Creates a likelihood function from the input hypothesis
# Documentation is in man directory.
create.likelihood <- function(hypothesis, addAttr=FALSE, ...) {

  vecfunc <- create.likelihood.vectors(hypothesis, addAttr, ...)

  likelihood.scalar <- function(...) { 
    result <- vecfunc(...) 
    prod(result$objectives * result$penalties)
  }

  attributes(likelihood.scalar) <- attributes(vecfunc)
  return(likelihood.scalar)
}

# Creates a likelihood function from the input hypothesis
# Documentation is in man directory.
create.likelihood.log <- function(hypothesis, addAttr=FALSE, ...) {

  vecfunc <- create.likelihood.vectors(hypothesis, addAttr, ...)
  likelihood.log <- function(...) { 
    result <- vecfunc(...) 
    sum(log(result$objectives)) + sum(log(result$penalties))
  }
  attributes(likelihood.log) <- attributes(vecfunc)
  return(likelihood.log)
}
