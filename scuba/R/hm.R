#
# 	hm.R
#
#	$Revision: 1.24 $	$Date: 2013/12/09 06:00:43 $
#
################################################################
#
#  Haldane-type models
#
#

hm <- local({

  # main function
  hm <- function(HalfT, M0=NULL, dM=NULL, ...,
                 N2 = list(HalfT=HalfT, M0=M0, dM=dM),
                 He = NULL,
                 title="user-defined model",
                 cnames=NULL,
                 mixrule=NULL) {
    if(length(list(...)) > 0)
      warning("Some unrecognised arguments were ignored")
    
    N2 <- validate(N2, cnames)
    pars <- list(N2=N2)
    nc <- nrow(N2)
    if(!is.null(He)) {
      He <- validate(He, cnames)
      if(nrow(He) != nc)
        stop("Data for He and N2 have different numbers of compartments")
      pars <- append(pars, list(He = He))
    }
    
    out <- list(title=title, pars=pars)

    if(is.null(mixrule)) 
      mixrule <- if(is.null(He)) "N2" else "interpolate"
    out$mixrule <- match.arg(mixrule, c("N2", "interpolate"))

    class(out) <- c("hm", class(out))
    return(out)
  }

  # helper function to validate entries
  validate <- function(X, cnames, validnames=c("HalfT", "M0", "dM")) {
    Xname <- deparse(substitute(X))
    isnull <- unlist(lapply(X, is.null))
    X <- X[!isnull]
    ok <- unlist(lapply(X, function(z) { is.numeric(z) && all(z > 0) }))
    if(!all(ok))
      stop(paste("All", Xname, "data should be positive numbers"))
    namesX <- names(X)
    if(length(namesX) == 0)
      stop(paste("The components of", sQuote(Xname), "must be labelled",
                 commasep(sQuote(validnames))))
    else if(!all(nzchar(namesX)))
      stop(paste("some components of", sQuote(Xname), "are not labelled"))
    else if(!all(ok <- namesX %in% validnames)) {
      bad <- !ok
      stop(paste("Unrecognised",
                 ngettext(sum(bad), "label", "labels"),
                 commasep(dQuote(namesX[bad])),
                 "in list", sQuote(Xname), "- recognised names are",
                 commasep(dQuote(validnames))))
    }
    if(with(X, is.null(M0) && !is.null(dM)))
      stop("If dM is provided, then M0 must be provided")
    X <- as.data.frame(X)
    if(!is.null(cnames)) {
      if(length(cnames) != nrow(X))
        stop(paste("length of cnames = ", length(cnames),
                   "!=", nrow(X), "= number of compartments"))
      rownames(X) <- cnames
    }
    return(X)
  }

  hm
})
            
print.hm <- function(x, ...) {
  stopifnot(inherits(x, "hm"))
  cat("Haldane type decompression model\n")
  y <- summary(x)
  cat(paste("Name:", y$title, "\n"))
  cat(paste(y$nc, ngettext(y$nc, "compartment\n", "compartments\n")))
  species <- y$species
  cat(paste("inert", ngettext(length(species), "gas:", "gases:"),
            paste(species, collapse=", "), "\n"))
  print(as.data.frame(x))
  return(invisible(NULL))
}

summary.hm <- function(object, ...) {
  stopifnot(inherits(object, "hm"))
  pars <- object$pars
  # names of inert gas species 
  species <- names(pars)
  # compartments
  cnames <- rownames(pars[[1]])
  nc     <- nrow(pars[[1]])
  # modelling capability for each species: 1 = halftimes, 2 = ndl, 3 = deco
  capabil <- lapply(pars, ncol)
  # rule for deriving M-values for mixed gases
  mixrule <- object$mixrule
  out <- list(title   = object$title,
              species = species,
              cnames  = cnames,
              nc      = nc,
              capabil = capabil,
              mixrule = mixrule,
              df = as.data.frame(object))
  class(out) <- c("summary.hm", class(out))
  return(out)
}

print.summary.hm <- function(x, ...) {
  cat("Haldane type decompression model\n")
  cat(paste("Name:", x$title, "\n"))
  cat(paste(x$nc, ngettext(x$nc, "compartment\n", "compartments\n")))
  cat(paste("inert gas species:",
            paste(x$species, collapse=", "), "\n"))
  for(i in seq(length(x$capabil))) {
    cat(paste("Data for ", x$species[i], ":\t ", sep="" ))
    switch(x$capabil[[i]],
         cat("Halftimes only\n"),
         cat("Halftimes and Surfacing M-values\n"),
         cat("Halftimes and M-values for any depth\n"))
  }
  cat(paste("Rule for generating M-values of mixed gases:",
            sQuote(x$mixrule), "\n"))
  print(x$df, ...)
  return(invisible(NULL))
}

as.data.frame.hm <- function(x, ...) {
  as.data.frame(x$pars, ...)
}

#################
#

# extracting properties/data

capable <- function(model, g="N2", what="HalfT") {
  stopifnot(inherits(model, "hm"))
  y <- summary(model)
  if(is.gas(g)) {
    ok <- TRUE
    if(g$fN2 > 0) ok <- capable(model, "N2", what)
    if(g$fHe > 0) ok <- ok && capable(model, "He", what)
    return(ok)
  }
  if(is.dive(g)) 
    return(capable(model, allspecies(g, inert=FALSE), what=what))
  if(!is.character(g))
    stop("Unrecognised format for argument g")
  # character string or character vector
  if(length(g) > 1) {
    g <- as.list(g)
    ok <- all(unlist(lapply(g, function(x, ...) { capable(g=x, ...) },
                            model=model, what=what)))
    return(ok)
  }
  # single character string
  if(g == "O2")
    return(TRUE)
  if(!(g %in% y$species))
    return(FALSE)
  k <- y$capabil[[g]]
  switch(what,
         HalfT=return(k >= 1),
         M0   =return(k >= 2),
         dM   =return(k >= 3))
  return(NA)
}

param <- function(model, species="N2", what="HalfT") {
  stopifnot(inherits(model, "hm"))
  spec <- model$pars[[species]]
  result <- spec[[what]]
  if(is.null(result)) return(NULL)
  names(result) <- row.names(spec)
  return(result)
}

########

# compute surfacing M-value for a mixture of inert gases

M0mix <- function(model, fN2, fHe, mixrule=NULL) {
  if(is.character(model)) model <- pickmodel(model)
  if(!inherits(model, "hm")) 
    stop(paste("model must be an object of class \'hm\'",
               "or a string giving the name of a Haldane model"))
  if(is.null(mixrule))
    mixrule <- summary(model)$mixrule
  if(is.gas(g <- fN2)) {
    fN2 <- g$fN2
    fHe <- g$fHe
  }
  # ensure compatible lengths
  df <- data.frame(fN2 = fN2, fHe=fHe)
  fN2 <- df$fN2
  fHe <- df$fHe
  ntimes <- length(fN2)
  M0.N2 <- param(model, "N2", "M0")
  one <- rep(1, ntimes)
  if(all(fHe == 0))
    return(outer(one, M0.N2, "*"))
  # Helium is present
  if(!capable(model, "He", "M0"))
    stop("Model does not provide surfacing M-values for Helium")
  switch(mixrule,
         N2 = { 
           M0 <- outer(one, M0.N2, "*")
         },
         interpolate={
           dM.N2 <- param(model, "N2", "dM")
           M0.He <- param(model, "He", "M0")
           dM.He <- param(model, "He", "dM")
           # Buehlmann equivalent parameters
           aN2 <- M0.N2 - dM.N2
           bN2 <- 1/dM.N2
           aHe <- M0.He - dM.He
           bHe <- 1/dM.He
           # mixture fraction (time-dependent)
           fIG <- fN2+fHe
           denom <- ifelse(fIG > 0, fIG, 1)
           z <- fN2/denom
           # apply to Buehlmann parameters
           a <- outer(z, aN2, "*") + outer(1-z, aHe, "*")
           b <- outer(z, bN2, "*") + outer(1-z, bHe, "*")
           # convert back to M0
           M0 <- a + 1/b
         },
         stop("Unrecognised mixture rule")
         )
  return(M0)
}

# compute M-value at depth for a mixture of inert gases

Mmix <- function(model, depth, fN2, fHe, mixrule=NULL) {
  if(is.character(model)) model <- pickmodel(model)
  if(!inherits(model, "hm")) 
    stop(paste("model must be an object of class \'hm\'",
               "or a string giving the name of a Haldane model"))
  if(!capable(model, "N2", "dM"))
    stop("Model does not provide M-value gradients (dM) for Nitrogen")
  if(is.null(mixrule))
    mixrule <- summary(model)$mixrule
  if(is.gas(g <- fN2)) {
    fN2 <- g$fN2
    fHe <- g$fHe
  }
  # ensure all lengths are the same
  df <- data.frame(depth=depth, fN2=fN2, fHe=fHe)
  depth <- df$depth
  fN2   <- df$fN2
  fHe   <- df$fHe
  ntimes <- length(depth)
  M0.N2 <- param(model, "N2", "M0")
  dM.N2 <- param(model, "N2", "dM")
  # pressure components
  one <- rep(1, ntimes)
  extraP <- depth/10
  if(all(fHe == 0))
    return(outer(one, M0.N2, "*") + outer(extraP, dM.N2, "*"))
  # Helium is present
  if(!capable(model, "He", "M0"))
    stop("Model does not provide surfacing M-values for Helium")
  if(!capable(model, "He", "dM"))
    stop("Model does not provide M-value gradients (dM) for Helium")
  switch(mixrule,
         N2 = { 
           M <- outer(one, M0.N2, "*") + outer(extraP, dM.N2, "*")
         },
         interpolate={
           M0.He <- param(model, "He", "M0")
           dM.He <- param(model, "He", "dM")
           # Buehlmann equivalent parameters
           aN2 <- M0.N2 - dM.N2
           bN2 <- 1/dM.N2
           aHe <- M0.He - dM.He
           bHe <- 1/dM.He
           # mixture fraction (time-dependent)
           fIG <- fN2+fHe
           denom <- ifelse(fIG > 0, fIG, 1)
           z <- fN2/denom
           # apply to Buehlmann parameters
           a <- outer(z, aN2, "*") + outer(1-z, aHe, "*")
           b <- outer(z, bN2, "*") + outer(1-z, bHe, "*")
           # compute M in Buehlmann form
           Pmat <- matrix(1 + extraP, length(extraP), ncol(b))
           M <- a + Pmat/b
         },
         stop("Unrecognised mixture rule")
         )
  return(M)
}

deco.ceiling <- function(state, model, what=c("depth", "pressure")) {
  what <- match.arg(what)
  if(is.character(model)) model <- pickmodel(model)
  if(!inherits(model, "hm")) 
    stop(paste("model must be an object of class \'hm\'",
               "or a string giving the name of a Haldane model"))
  sm <- summary(model)
  ncomp <- sm[["nc"]]
  compnames <- sm[["cnames"]]
  if(is.vector(state) && length(state) == ncomp) {
    # assume state contains N2 tensions at a single time point
    if(!capable(model, "N2", "dM"))
      stop("Model does not allow calculation of decompression ceiling")
    species <- "N2"
    profile <- array(state, dim=c(1, ncomp, 1))
    dimnames(profile) <- list(NULL, compnames, "N2")
  } else if((is.matrix(state) || is.data.frame(state))
            && nrow(state) == ncomp) {
    state <- as.matrix(state)
    if(ncol(state) > length(sm[["species"]]))
      stop("Format of state is not understood")
    # match columns to gases
    cn <- colnames(state)
    if(length(cn) == ncol(state) && all(nzchar(cn))) {
      m <- match(cn, sm[["species"]])
      if(any(is.na(m)))
        stop("Unrecognised names of gases")
      state <- state[,m]
      species <- colnames(state)
    } else {
      species <- sm[["species"]][1:ncol(state)]
      colnames(state) <- species
    }
    profile <- array(state, dim=c(1,dim(state)),
                     dimnames=append(list(NULL), dimnames(state)))
  } else if(is.array(state) && dim(state)[2] == ncomp) {
    nplanes <- dim(state)[3]
    if(nplanes > length(sm$species))
      stop("Format of state is not understood")
    # match planes to gases
    pn <- dimnames(state)[[3]]
    if(length(pn) == nplanes && all(nzchar(pn))) {
      m <- match(pn, sm$species)
      if(any(is.na(m)))
        stop("Unrecognised names of gases")
      state <- state[,,m,drop=FALSE]
      species <- dimnames(state)[[3]][m]
    } else {
      species <- sm$species[1:nplanes]
      dimnames(state)[[3]] <- species
    }
    profile <- state
  } else stop("Format of state is not understood")

  result <- array(, dim=dim(profile), dimnames=dimnames(profile))
  ntimes <- dim(profile)[1]
  
  for(g in species) {
    if(!capable(model, g, "dM"))
      stop("Model does not allow calculation of decompression ceiling")
    tensions <- profile[ , , g]
    M0 <- param(model, g, "M0")
    dM <- param(model, g, "dM")
    M0 <- matrix(M0, nrow=ntimes, ncol=ncomp, byrow=TRUE)
    dM <- matrix(dM, nrow=ntimes, ncol=ncomp, byrow=TRUE)
    Pambtol <- pmax(0, (tensions - M0)/dM + 1)
    result[,,g] <- switch(what,
                          pressure = Pambtol,
                          depth = pmax(0, (Pambtol - 1) * 10))
  }
  return(result[,,,drop=TRUE])
}

  
###################################################################
#
#  Haldane type models that are installed
#

pickmodel <- function(model) {
  if(missing(model)) {
    cat(paste("Available options are:",
              paste(sQuote(names(.scuba.models)), collapse=", "), "\n"))
    return(invisible(NULL))
  }
  
  stopifnot(is.character(model))

  k <- pmatch(model, names(.scuba.models))
  if(is.na(k))
    stop(paste("Unrecognised model", sQuote(model)))
  return(.scuba.models[[k]])
}

################################################################
#
#  Standard models
#
#

.buehlmannL16A <-
  list(tN2=c(4, 5, 8, 12.5, 18.5, 27, 38.3, 54.3, 77,
             109, 146, 187, 239, 305, 390, 498, 635),
       aN2=c(
         1.2599,
         1.1696,
         1.0000,
         0.8618,
         0.7562,	
         0.6667,	
         0.5933,
         0.5282,
         0.4701,	
         0.4187,
         0.3798,
         0.3497,
         0.3223,
         0.2971,
         0.2737,
         0.2523,
         0.2327),
       bN2=c(
         0.5050,
         0.5578,
         0.6514,
         0.7222,
         0.7825,
         0.8126,
         0.8434,
         0.8693,
         0.8910,
         0.9092,
         0.9222,
         0.9319,
         0.9403,
         0.9477,
         0.9544,
         0.9602,
         0.9653),
       tHe=c(
         1.51,
         1.88,
         3.02,
         4.72,
         6.99,
         10.21,
         14.48,
         20.53,
         29.11,
         41.20,
         55.19,
         70.69,
         90.34,
         115.29,
         147.42,
         188.24,
         240.03),
       aHe=c(
         1.7424,
         1.6189,
         1.3830,
         1.1919,
         1.0458,
         0.9220,
         0.8205,
         0.7305,
         0.6502,
         0.5950,
         0.5545,
         0.5333,
         0.5189,
         0.5181,
         0.5176,
         0.5172,
         0.5119),
       bHe=c(
         0.4245,
         0.4770,
         0.5747,
         0.6527,
         0.7223,
         0.7582,
         0.7957,
         0.8279,
         0.8553,
         0.8757,
         0.8903,
         0.8997,
         0.9073,
         0.9122,
         0.9171,
         0.9217,
         0.9267))

.scuba.models <-
  list(Haldane = hm(
         c(5, 10, 20, 40, 75),
         rep(2 * 0.79, 5),
         rep(2 * 0.79, 5),
         title="Haldane"),
       USN = hm(
         c(5, 10, 20, 40, 80, 120),
         c(104, 88, 72, 56, 52, 51)/32.646,
         c(2.27, 2.01, 1.67, 1.34, 1.26, 1.19),
         title="USN"),
       DSAT = hm(
         c(5, 10, 20, 30, 40, 60, 80, 120),
         c(99.08, 82.68, 66.89, 59.74, 55.73, 51.44, 49.21, 46.93)/32.646,
         title="DSAT"),
       Workman65 = hm(
         c(5, 10, 20, 40, 80, 120, 160, 200, 240),
         c(31.5, 26.7, 21.8, 17.0, 16.4, 15.8, 15.5, 15.5, 15.2)/10,
         c(1.8, 1.6, 1.5, 1.4, 1.3, 1.2, 1.2, 1.1, 1.1),
         title="Workman65"),
       "ZH-L16A"= {
           hm(with(.buehlmannL16A, tN2),
              with(.buehlmannL16A, aN2 + 1/bN2),
              with(.buehlmannL16A, 1/bN2),
              title="Buehlmann ZH-L16A",
              cnames=c("1", "1b", paste(2:16)),
              He=list(HalfT=with(.buehlmannL16A, tHe),
                      M0=with(.buehlmannL16A, aHe + 1/bHe),
                      dM=with(.buehlmannL16A, 1/bHe)),
              mixrule="interpolate")
         })
