"overlap" <-
  function(modes, dv, nmodes=20) {
    if(missing(modes))
      stop("overlap: 'modes' must be prodivded")
    if(missing(dv))
      stop("overlap: 'dv' must be prodivded")

    if ("pca" %in% class(modes)) {
      ev <- modes$U
      mass <- NULL
      first.mode <- 1
    }
    else if("nma" %in% class(modes)) {
      ev <- modes$modes
      mass <- modes$mass
      first.mode <- modes$triv.modes+1
      nmodes <- modes$triv.modes + nmodes
    }
    else {
      if(class(modes)!="matrix" && class(modes)!="pca.loadings")
        stop("overlap: 'modes' must be an object of type 'pca', 'nma', or 'matrix'")
      ev <- modes
      mass <- NULL
      first.mode <- 1
    }

    if (nrow(ev)!=length(dv))
      stop("overlap: unequal vector lengths")
    
    if ( ncol(ev) < nmodes ) {
      nmodes <- dim(ev)[2L]
      warning("nmodes larger than dimensions of 'modes'")
    }

    inds <- seq(first.mode, nmodes)
    ev <- ev[,inds]

    ## Normalize vectors - mass-weighted if normal modes are
    ev <- normalize.vector(ev, mass)
    dvn <- normalize.vector(dv, mass)
    
    overlap.values <- inner.prod(ev, dvn, mass)**2
    cum <- cumsum(overlap.values)
    
    out <- list(overlap=overlap.values, overlap.cum=cum)
    
    return(out)
  }
