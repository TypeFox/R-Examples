# CHNOSZ/swap.basis.R
# functions related to swapping basis species
# extracted from basis() 20120114 jmd

# return the current basis elements
basis.elements <- function(basis = get("thermo")$basis) {
  if(is.null(basis)) stop("basis species are not defined")
  return(as.matrix(basis[, 1:nrow(basis), drop=FALSE]))
}

# calculate chemical potentials of elements from logarithms of activity of basis species
element.mu <- function(basis = get("thermo")$basis, T = 25) {
  # matrix part of the basis definition
  basis.mat <- basis.elements(basis)
  # the standard Gibbs energies of the basis species
  if(T==25) G <- get("thermo")$obigt$G[basis$ispecies]
  else G <- unlist(subcrt(basis$ispecies, T=T, property="G")$out)
  # chemical potentials of the basis species
  species.mu <- G - convert(basis$logact, "G", T=convert(T, "K"))
  # chemical potentials of the elements
  element.mu <- solve(basis.mat, species.mu)
  # give them useful names
  names(element.mu) <- colnames(basis.mat)
  return(element.mu)
}

# calculate logarithms of activity of basis species from chemical potentials of elements
basis.logact <- function(emu, basis = get("thermo")$basis, T = 25) {
  # matrix part of the basis definition
  basis.mat <- basis.elements(basis)
  # elements in emu can't be less than the number in the basis
  if(length(emu) < ncol(basis.mat)) stop("number of elements in 'emu' is less than those in basis")
  # sort names of emu in order of those in basis.mat
  ielem <- match(names(emu), colnames(basis.mat))
  # check that elements of basis.mat and emu are identical
  if(any(is.na(ielem))) stop(paste("element(s)", paste(names(emu)[is.na(ielem)], collapse=" "), "not found in basis"))
  # the standard Gibbs energies of the basis species
  if(T==25) G <- get("thermo")$obigt$G[basis$ispecies]
  else G <- unlist(subcrt(basis$ispecies, T=T, property="G")$out)
  # the chemical potentials of the basis species in equilibrium
  # with the chemical potentials of the elements
  basis.mu <- colSums((t(basis.mat)*emu)) - G
  # convert chemical potentials to logarithms of activity
  basis.logact <- -convert(basis.mu, "logK", T=convert(T, "K"))
  # give them useful names
  names(basis.logact) <- rownames(basis.mat)
  return(basis.logact)
}

ibasis <- function(species) {
  # get the index of a basis species from a species index, name or formula
  basis <- basis()
  if(is.numeric(species)) ib <- match(species, basis$ispecies)
  else {
    # character: first look for formula of basis species
    ib <- match(species, rownames(basis))
    # if that doesn't work, look for name of basis species
    if(is.na(ib)) ib <- match(species, get("thermo")$obigt$name[basis$ispecies])
  }
  return(ib)
}

# swap in one basis species for another
swap.basis <- function(species, species2, T = 25) {
  # before we do anything, remember the old basis definition
  oldbasis <- get("thermo")$basis
  # and the species definition
  ts <- species()
  # the delete the species
  species(delete=TRUE)
  if(is.null(oldbasis)) 
    stop("swapping basis species requires an existing basis definition")
  # both arguments must have length 1
  if(missing(species) | missing(species2))
    stop("two species must be identified")
  if(length(species) > 1 | length(species2) > 2)
    stop("can only swap one species for one species")
  # arguments are good, now find the basis species to swap out
  ib <- ibasis(species)
  if(is.na(ib)) stop(paste("basis species '",species,"' is not defined",sep=""))
  # find species2 in the thermodynamic database
  if(is.numeric(species2)) ispecies2 <- species2
  else ispecies2 <- suppressMessages(info(species2))
  if(is.na(ispecies2) | is.list(ispecies2))
    stop(paste("a species matching '",species2,"' is not available in thermo$obigt",sep=""))
  # try to load the new basis species
  ispecies <- oldbasis$ispecies
  ispecies[ib] <- ispecies2
  newbasis <- put.basis(ispecies)
  # if put.basis didn't stop with an error, we're good to go!
  # what were the original chemical potentials of the elements?
  emu <- element.mu(oldbasis, T=T)
  # the corresponding logarithms of activities of the new basis species
  bl <- basis.logact(emu, newbasis, T=T)
  # update the basis with these logacts
  mb <- mod.basis(ispecies, logact=bl)
  # restore species if they were defined
  if(!is.null(ts)) {
    suppressMessages(species(ts$ispecies))
    suppressMessages(species(1:nrow(get("thermo")$species), ts$logact))
  }
  # all done, return the new basis definition
  return(mb)
}

