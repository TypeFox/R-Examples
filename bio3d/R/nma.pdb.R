"nma.pdb" <- function(pdb, inds=NULL, ff='calpha', pfc.fun=NULL, mass=TRUE,
                      temp=300.0, keep=NULL, hessian=NULL, outmodes=NULL, ... ) {
  
  ## Log the call
  cl <- match.call()

  if(!is.pdb(pdb))
    stop("please provide a 'pdb' object as obtained from 'read.pdb()'")

  if(!is.null(outmodes) & !is.select(outmodes))
    stop("provide 'outmodes' as obtained from function atom.select()")
  
  ## Prepare PDB
  ## Take only first frame of multi model PDB files
  if(nrow(pdb$xyz)>1) {
    warning("multimodel PDB file detected - using only first frame")
    pdb$xyz=pdb$xyz[1,, drop=FALSE]
  }

  ## Trim to only CA atoms
  if(is.null(inds)) {
    ca.inds <- atom.select(pdb, "calpha", verbose=FALSE)
    pdb.in <- trim.pdb(pdb, ca.inds)
  }

  ## or to user selection
  else {
    pdb.in <- trim.pdb(pdb, inds)
    if(!all(pdb.in$atom$elety=="CA"))
      stop("non-CA atoms detected")
  }

  ## Indices for effective hessian
  if(is.select(outmodes)) {
    ## re-select since outmodes indices are based on input PDB
    inc.inds <- .match.sel(pdb, pdb.in, outmodes)
    pdb.out <- trim.pdb(pdb.in, inc.inds)
  }
  else {
    pdb.out <- pdb.in
    inc.inds <- atom.select(pdb.in, "all", verbose=FALSE)
  }

  ## fetch number of atoms and sequence
  natoms.in  <- ncol(pdb.in$xyz)/3
  natoms.out <- ncol(pdb.out$xyz)/3
  sequ <- pdb.in$atom$resid
  
  if (natoms.in<3)
    stop("nma: insufficient number of atoms")

  ## check structure connectivity
  conn <- inspect.connectivity(pdb.in$xyz)
  if(!conn) {
    warning("Possible multi-chain structure or missing in-structure residue(s) present\n", 
            "  Fluctuations at neighboring positions may be affected.")
  }
  
  ## Process input arguments
  init <- .nma.init(ff=ff, pfc.fun=pfc.fun, sequ=sequ, ...)
  
  ## Use aa2mass to fetch residue mass
  if (mass) {
    masses.in <- do.call('aa2mass', c(list(pdb=sequ, inds=NULL), init$am.args))
    masses.out <- masses.in[ inc.inds$atom ]
  }
  
  ## No mass-weighting
  else {
    masses.out <- NULL;
  }

  ## NMA hessian
  hessian <- .nma.hess(pdb.in$xyz, init=init,
                       hessian=hessian, inc.inds=inc.inds)

  ## mass weight hessian
  if(!is.null(masses.out))
    hessian <- .nma.mwhessian(hessian, masses=masses.out)

  ## diagaonalize - get eigenvectors
  ei <- .nma.diag(hessian)

  ## make a NMA object
  m <- .nma.finalize(ei, xyz=pdb.out$xyz, temp=temp, masses=masses.out,
                        natoms=natoms.out, keep=keep, call=cl)
  return(m)
}


".nma.init" <- function(ff=NULL, pfc.fun=NULL, ...) {

    ## Arguments to functions build.hessian and aa2mass
    bh.names <- names(formals( build.hessian ))
    am.names <- names(formals( aa2mass ))
    
    dots <- list(...)
    bh.args <- dots[names(dots) %in% bh.names]
    am.args <- dots[names(dots) %in% am.names]

    ## Define force field
    if (is.null(pfc.fun)) {
      ff <- load.enmff(ff)
    }
    else {
      ## Use customized force field
      if(!is.function(pfc.fun))
        stop("'pfc.fun' must be a function")
      bh.args <- bh.args[ !('pfc.fun' %in% names(bh.args)) ]
      ff <- pfc.fun
    }

    ## Check for optional arguments to pfc.fun
    ff.names <- names(formals( ff ))
    ff.args  <- dots[names(dots) %in% ff.names]
    
    ## Redirect them to build.hessian
    bh.args  <- c(bh.args, ff.args)
    
    ## Arguments without destination
    all.names <- unique(c(bh.names, am.names, ff.names))
    if(!all(names(dots) %in% all.names)) {
      oops <- names(dots)[!(names(dots) %in% all.names)]
      stop(paste("argument mismatch:", oops))
    }

    if(length(bh.args)==0)
      bh.args=NULL
    if(length(am.args)==0)
      am.args=NULL
    #if(length(ff.args)==0)
    #  ff.args=NULL
        
    out <- list(pfcfun=ff, bh.args=bh.args, am.args=am.args)
    return(out)
  }



## extract effective hessian
".nma.trim.hessian" <- function(hessian, inc.inds=NULL) {
  if(!is.matrix(hessian))
    stop("hessian must be a matrix")
  if(is.null(inc.inds))
    stop("indices must be provided")
  
  kaa     <- hessian[inc.inds, inc.inds]
  kqq.inv <- solve(hessian[-inc.inds, -inc.inds])
  kaq     <- hessian[inc.inds, -inc.inds]
  kqa     <- t(kaq)
  k <- kaa - ((kaq %*% kqq.inv) %*% kqa)
  return(k)
}

## mass-weight hessian
".nma.mwhessian" <- function(hessian, masses=NULL) {
  if(!is.matrix(hessian))
    stop("hessian must be a matrix")
  if(is.null(masses))
    stop("masses must be provided")

  #cat(" Mass weighting Hessian...")
  #ptm <- proc.time()
   
  dims <- dim(hessian)
  natoms <- dims[1] / 3

  if(length(masses)!=natoms)
    stop("dimension mismatch")
  
  masses <- sqrt(masses)
  inds <- rep(1:natoms, each=3)
  col.inds <- seq(1, ncol(hessian), by=3)

  for ( i in 1:natoms ) {
    m <- col.inds[i]
    hessian[,m:(m+2)] <- hessian[,m:(m+2)] * (1/masses[i])
    hessian[,m:(m+2)] <- hessian[,m:(m+2)] * (1/masses[inds])
  }

  #t <- proc.time() - ptm
  #cat("\tDone in", t[[3]], "seconds.\n")

  return(hessian)
}

## wrapper for generating the hessian matrix
".nma.hess" <- function(xyz, init=NULL,
                        hessian=NULL, inc.inds=NULL) {
  
  natoms <- ncol(as.xyz(xyz))/3
  if(nrow(xyz)>1)
    xyz=xyz[1,,drop=FALSE]

  ## Build the Hessian Matrix
  if(is.null(hessian)) {
    cat(" Building Hessian...")
    ptm <- proc.time()
    H <- do.call('build.hessian',
                 c(list(xyz=xyz, pfc.fun=init$pfcfun), init$bh.args))
    t <- proc.time() - ptm
    cat("\t\tDone in", t[[3]], "seconds.\n")
  }
  else {
    H <- hessian
  }
  
  ## Effective Hessian
  if(!is.null(inc.inds)) {
    if(ncol(xyz)>length(inc.inds$xyz)) {
      cat(" Extracting effective Hessian..")
      ptm <- proc.time()
      H <- .nma.trim.hessian(H, inc.inds=inc.inds$xyz)
      t <- proc.time() - ptm
      cat("\tDone in", t[[3]], "seconds.\n")
    }
  }
  return(H)
}

## diagonalize hessian
".nma.diag" <- function(H) {
  
  ## Diagonalize matrix
  cat(" Diagonalizing Hessian...")
  ptm <- proc.time()
  ei <- eigen(H, symmetric=TRUE)
  t <- proc.time() - ptm
  cat("\tDone in", t[[3]], "seconds.\n")
  return(ei)
}

## build a NMA object
".nma.finalize" <- function(ei, xyz, temp, masses, natoms, keep, call) {
  if(length(masses)>0)
    mass <- TRUE
  else
    mass <- FALSE

  xyz=as.xyz(xyz)
  dims <- dim(ei$vectors)
  dimchecks <- c(ncol(xyz)/3==natoms,
                 ifelse(mass, length(masses)==natoms, TRUE),
                 dims[1]/3==natoms,
                 dims[2]/3==natoms)
  
  if(!all(dimchecks))
    stop(paste("dimension mismatch when generating nma object\n",
               paste(dimchecks, collapse=", ")))
  
  ## Raw eigenvalues
  ei$values <- round(ei$values, 6)

  ## Trivial modes first - sort on abs(ei$values)
  sort.inds  <- order(abs(ei$values))
  ei$values  <- ei$values[sort.inds]
  ei$vectors <- ei$vectors[, sort.inds]
  
  ## hard code 6 trivial modes
  triv.modes <- seq(1, 6)

  ## keep only a subset of modes - including trivial modes
   if(!is.null(keep)) {
    if(keep>ncol(ei$vectors))
      keep <- ncol(ei$vectors)
    keep.inds <- seq(1, keep)
    ei$vectors <- ei$vectors[,keep.inds]
    ei$values <- ei$values[keep.inds]
  }

    ## Frequencies are given by
    if (mass)  {
      pi <- 3.14159265359
      freq <- sqrt(abs(ei$values)) / (2 * pi)
      force.constants <- NULL
    } else {
      freq <- NULL
      force.constants <- ei$values
    }

    ## Raw unmodified eigenvectors:
    ## ei$vectors

    ## V holds the eigenvectors converted to unweighted Cartesian coords:
    V <- ei$vectors

    ## Change to non-mass-weighted eigenvectors
    if(mass) {
      wts.sqrt <- sqrt(masses)
      tri.inds <- rep(1:natoms, each=3)
      V <- apply(V, 2, '*', 1 / wts.sqrt[tri.inds])
    }

    ## Temperature scaling
    kb <- 0.00831447086363271
    if ( !is.null(temp) ) {
      if (!is.null(freq)) {
        amplitudes <- sqrt(2* temp * kb) / (2* pi * freq[ -triv.modes ])
        amplitudes <- c(rep(1,length(triv.modes)), amplitudes)
      }
      else if(!is.null(force.constants)) {
        amplitudes <- sqrt((2* temp * kb) /
                           force.constants[ -triv.modes ])
        amplitudes <- c(rep(1,length(triv.modes)), amplitudes)
      }
    } else {
      amplitudes <- rep(1, times=3*natoms)
    }

    ## Temperature scaling of eigenvectors
    for ( i in (length(triv.modes)+1):ncol(V) ) {
      V[,i] <- (V[,i] * amplitudes[i])
    }

    ## Check if first modes are zero-modes
    if(any(ei$values<0)) {
      warning("Negative eigenvalue(s) detected! \
              This can be an indication of an unphysical input structure.")
    }

    ## Output to class "nma"
    nma <- list(modes=V,
                frequencies=NULL,
                force.constants=NULL,
                fluctuations=NULL,
                U=ei$vectors, L=ei$values,

                xyz=xyz,
                mass=masses,
                temp=temp,
                triv.modes=length(triv.modes),
                natoms=natoms,
                call=call)

    if(mass) {
      class(nma) <- c("VibrationalModes", "nma")
      nma$frequencies <- freq
    }
    else {
      class(nma) <- c("EnergeticModes", "nma")
      nma$force.constants <- force.constants
    }

    ## Calculate mode fluctuations
    nma$fluctuations <- fluct.nma(nma, mode.inds=NULL)

    ## Notes:
    ## U are the raw unmodified eigenvectors
    ## These mode vectors are in mass-weighted coordinates and not
    ## scaled by the thermal amplitudes, so they are orthonormal.

    ## V holds the eigenvectors converted to unweighted Cartesian
    ## coordinates.  Unless you set temp=NULL, the modes are
    ## also scaled by the thermal fluctuation amplitudes.

    return(nma)
  }

".match.sel" <- function(a, b, inds) {
  ## a= original pdb
  ## b= trimmed pdb
  ## inds= indices of pdb 'a' to keep
  ## find corresponding atoms in b
  
  names.a <- paste(a$atom[inds$atom, "chain"],
                   a$atom[inds$atom, "resno"],
                   a$atom[inds$atom, "elety"],
                   a$atom[inds$atom, "eleno"], sep="-")
  
  names.b <- paste(b$atom[, "chain"],
                   b$atom[, "resno"], 
                   b$atom[, "elety"],
                   b$atom[, "eleno"], sep="-")
  
  inds <- which(names.b %in% names.a)
  out <- list(atom=inds, xyz=atom2xyz(inds))
  class(out) <- "select"
  return(out)
}
