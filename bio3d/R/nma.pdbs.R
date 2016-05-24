## Two options - both calculating modes on the FULL structure:
## 1 -  use k <- kaa - ((kaq %*% kqq.inv) %*% kqa) to derive hessian for core atoms
## 2 - return the full objects

"nma.pdbs" <- function(pdbs, fit=TRUE, full=FALSE, subspace=NULL,
                       rm.gaps=TRUE, varweight=FALSE, 
                       outpath = NULL, ncore=1, ...) {
 
  
  if(!inherits(pdbs, "pdbs"))
    stop("input 'pdbs' should be a list object as obtained from 'read.fasta.pdb'")

  ## Log the call
  cl <- match.call()
  
  if(!is.null(outpath))
    dir.create(outpath, FALSE)

  ## Parallelized by parallel package
  ncore <- setup.ncore(ncore, bigmem = TRUE)
  prev.warn <- getOption("warn")
  
  if(ncore>1) {
    mylapply <- mclapply
    options(warn=1)
  }
  else
    mylapply <- lapply

  ## Passing arguments to functions aa2mass and nma
  am.names <- names(formals( aa2mass ))
  nm.names <- names(formals( nma.pdb ))
  
  dots <- list(...)
  am.args <- dots[names(dots) %in% am.names]
  nm.args <- dots[names(dots) %in% nm.names]

  ## Limiting input 
  if("mass" %in% names(nm.args))
    mass <- nm.args$mass
  else
    mass <- TRUE
  if("ff" %in% names(nm.args))
    ff <- nm.args$ff
  else
    ff <- 'calpha'
  if("temp" %in% names(nm.args))
    temp <- nm.args$temp
  else
    temp <- 300
  if("keep" %in% names(nm.args))
    nm.keep <- nm.args$keep
  else
    nm.keep <- NULL
  
  if(!all((names(nm.args) %in% c("mass", "ff", "temp", "keep")))) {
    war <- paste(names(nm.args)[! names(nm.args) %in% c("mass", "ff", "temp", "keep") ], collapse=", ")
    warning(paste("ignoring arguments:", war))
  }
  
  ## Force field
  pfc.fun <- load.enmff(ff)
  
  ## Check for optional arguments to pfc.fun
  ff.names <- names(formals( pfc.fun ))
  ff.args  <- dots[names(dots) %in% ff.names]
  
  ## Set indicies
  gaps.res <- gap.inspect(pdbs$resno)
  gaps.pos <- gap.inspect(pdbs$xyz)

  ## check for missing masses before we start calculating
  if(any(pdbs$ali=="X") & mass==TRUE) {
    resnames <- c(bio3d::aa.table$aa3, names(am.args$mass.custom))
    
    ops.inds <- which(pdbs$ali=="X", arr.ind=TRUE)
    unknowns <- c()
    
    for(i in 1:nrow(ops.inds)) {
      j <- ops.inds[i, ]
      
      resid <- pdbs$resid[j[1], j[2]]
      if(any(!(resid %in% resnames)))
        unknowns <- c(unknowns, resid[!resid%in%resnames])
    }
    
    if(length(unknowns)>0) {
      options(warn=prev.warn)
      unknowns <- paste(unique(unknowns), collapse=", ")
      stop(paste0("Unknown mass for amino acid(s): ", unknowns,
                  "\n  Provide mass with argument 'mass.custom=list(", unknowns[1], "=100.00)',", 
                  "\n  or ommit mass weighting with argument 'mass=FALSE'."))
    }
  }

  ## Check connectivity
  con <- inspect.connectivity(pdbs, cut=4.05)
  if(!all(con)) {
    warning(paste(paste(basename(pdbs$id[which(!con)]), collapse=", "),
                  "might have missing residue(s) in structure:\n",
                  "  Fluctuations at neighboring positions may be affected."))
  }
            
  ## Use for later indexing
  pdbs$inds <- matrix(NA, ncol=ncol(pdbs$resno), nrow=nrow(pdbs$resno))
 
  ## Number of modes to store in U.subspace
  if(is.null(subspace)) {
    keep <- length(gaps.pos$f.inds)-6
  }
  else {
    keep <- subspace
    if (length(gaps.pos$f.inds) < (keep+6))
      keep <- length(gaps.pos$f.inds)-6
  }
  
  ## Coordiantes - fit or not
  if(fit) {
    xyz <- fit.xyz(fixed = pdbs$xyz[1, ], mobile = pdbs,
                   fixed.inds = gaps.pos$f.inds, mobile.inds = gaps.pos$f.inds,
                   ncore = ncore)
  }
  else
    xyz <- pdbs$xyz
  
  ## Fluctuations for each structure
  if(rm.gaps)
    flucts <- matrix(NA, nrow=nrow(gaps.res$bin), ncol=length(gaps.res$f.inds))
  else
    flucts <- matrix(NA, nrow=nrow(gaps.res$bin), ncol=ncol(gaps.res$bin))

  ## store residue numbers (same as pdbs$inds[,gaps$f.inds])
  ##resnos <- flucts

  ## List object to store each modes object
  if(full)
    all.modes <- list()
  else
    all.modes <- NULL

  ## 3D array- containing the modes vectors for each structure
  if(rm.gaps)
    modes.array <- array(NA, dim=c(length(gaps.pos$f.inds), keep, nrow(gaps.res$bin)))
  else
    modes.array <- array(NA, dim=c(ncol(pdbs$xyz), keep, nrow(gaps.res$bin)))

  ## store eigenvalues of the first modes
  L.mat <- matrix(NA, ncol=keep, nrow=nrow(gaps.res$bin))

  if(is.null(outpath))
    fname <- tempfile(fileext = "pdb")

  
  ##### Start calculation of variance weighting  #####
  wts <- NULL
  if(is.logical(varweight)) {
    if(varweight)
      wts <- var.xyz(xyz, weights=TRUE)
  }
  else {
    dims.vw <- dim(varweight)
    if(all(dims==ncol(xyz)/3)) {
      wts <- varweight
      varweight <- TRUE
    }
    else
      stop("incompatible length of varweight vector")
  }

  ### Memory usage ###
  dims <- dim(modes.array)
  mem.usage <- sum(c(as.numeric(object.size(modes.array)),
                     as.numeric(object.size(L.mat)),
                     as.numeric(object.size(flucts)),
                     as.numeric(object.size(matrix(NA, ncol=dims[3], nrow=dims[3]))) ))*2
                   
  if(full) {
    if(is.null(nm.keep))
      tmpncol <- dims[2]
    else
      tmpncol <- nm.keep
    
    size.mat <- object.size(matrix(0.00000001, ncol=tmpncol, nrow=dims[1]))
    size.vec <- object.size(vector(length=dims[1], 'numeric'))

    tot.size <- ((size.mat * 2) + (size.vec * 4)) * length(pdbs$id)
    mem.usage <- mem.usage+tot.size
  }
  mem.usage=round(mem.usage/1048600,1)
  

  #### Print overview of scheduled calcualtion ####
  cat("\nDetails of Scheduled Calculation:\n")
  cat(paste("  ...", length(pdbs$id), "input structures", "\n"))
  if(keep>0)
    cat(paste("  ...", "storing", keep, "eigenvectors for each structure", "\n"))
  if(keep>0)
    cat(paste("  ...", "dimension of x$U.subspace: (",
              paste(dims[1], dims[2], dims[3], sep="x"), ")\n"))
  
  if(fit)
    cat(paste("  ...", "coordinate superposition prior to NM calculation", "\n"))

  if(varweight)
    cat(paste("  ...", "weighting force constants based on structural variance", "\n"))

  if(full)
    cat(paste("  ... individual complete 'nma' objects will be stored", "\n"))
  
  if(rm.gaps)
    cat(paste("  ... aligned eigenvectors (gap containing positions removed) ", "\n"))

  if(mem.usage>0)
    cat(paste("  ...", "estimated memory usage of final 'eNMA' object:", mem.usage, "Mb \n"))
  
  cat("\n")


  ##### Start modes calculation #####
  pb <- txtProgressBar(min=0, max=length(pdbs$id), style=3)
  
  ## shared memory to follow progress bar
  if(ncore>1)
    iipb <- bigmemory::big.matrix(1, length(pdbs$id), init=NA)

  
  ## call .calcAlnModes for each structure in 'pdbs'
  alnModes <- mylapply(1:length(pdbs$id), .calcAlnModes,
                       pdbs, xyz, gaps.res,
                       mass, am.args, nm.keep, temp, keep, wts,
                       rm.gaps, full, 
                       pfc.fun, ff, ff.args, outpath, pb, ncore, env=environment())
  close(pb)
  
  ##### Collect data #####
  for(i in 1:length(alnModes)) {
    tmp.modes <- alnModes[[i]]
    
    modes.array[,,i] = tmp.modes$U
    L.mat[i, ]       = tmp.modes$L
    flucts[i, ]      = tmp.modes$flucts

    if(full)
      all.modes[[i]] = tmp.modes$modes
  }

  remove(alnModes)
  invisible(gc())
  
  ##### RMSIP ######
  rmsip.map <- NULL
  if(rm.gaps) {
    rmsip.map <- .calcRMSIP(modes.array, ncore=ncore)
    rownames(rmsip.map) <- basename(rownames(pdbs$xyz))
    colnames(rmsip.map) <- basename(rownames(pdbs$xyz))

    if(!fit)
      warning("rmsip calculated on non-fitted structures:
               ensure that your input coordinates are pre-fitted.")
  }

  if(ncore>1) {
#    rm(iipb, pos = ".GlobalEnv")  ## remove global iipb variable
    options(warn=prev.warn)       ## restore warning option
  }

  rownames(flucts) <- basename(rownames(pdbs$xyz))
  out <- list(fluctuations=flucts, rmsip=rmsip.map,
              U.subspace=modes.array, L=L.mat, full.nma=all.modes,
              xyz=xyz, call=cl)
  class(out) = "enma"
  return(out)
}

.calcRMSIP <- function(x, ncore=1) {
  if(ncore>1)
    mylapply <- mclapply
  else
    mylapply <- lapply
  
  n <- dim(x)[3]
  mat <- matrix(NA, n, n)
  inds <- rbind(pairwise(n),
                matrix(rep(1:n,each=2), ncol=2, byrow=T))

  mylist <- mylapply(1:nrow(inds), function(row) {
    return(list(
      rmsip=rmsip(x[,,inds[row,1]], x[,,inds[row,2]])$rmsip,
      i=inds[row,1], j=inds[row,2])
           )
  })
  
  for ( i in 1:length(mylist)) {
    tmp.rmsip <- mylist[[i]]
    mat[tmp.rmsip$i, tmp.rmsip$j] <- tmp.rmsip$rmsip
  }

  mat[ inds[,c(2,1)] ] = mat[ inds ]
  return(round(mat, 4))
}


## Calculate 'aligned' normal modes of structure i in pdbs
.calcAlnModes <- function(i, pdbs, xyz, gaps.res,
                          mass, am.args, nm.keep, temp, keep, wts,
                          rm.gaps, full, 
                          pfc.fun, ff, ff.args, outpath, pb, ncore, env=NULL) {

  ## Set indices for this structure only
  f.inds <- NULL
  f.inds$res <- which(gaps.res$bin[i,]==0)
  f.inds$pos <- atom2xyz(f.inds$res)

  ## similar to $resno but sequential indices
  pdbs$inds[i, f.inds$res] <- seq(1, length(f.inds$res))
  
  ## Indices to extract from Hessian
  inds.inc <- pdbs$inds[i, gaps.res$f.inds]
  inds.exc <- pdbs$inds[i, gaps.res$t.inds][ !is.na(pdbs$inds[i, gaps.res$t.inds]) ]

  inds.inc.xyz <- atom2xyz(inds.inc)
  inds.exc.xyz <- atom2xyz(inds.exc)
  
  ## Generate content of PDB object
  tmp.xyz <- xyz[i, f.inds$pos]
  resno   <- pdbs$resno[i,f.inds$res]
  chain   <- pdbs$chain[i,f.inds$res]

  ## Fix for missing chain IDs
  chain[is.na(chain)] <- ""
  
  ## 3-letter AA code is provided in the pdbs object
  ## avoid using aa123() here (translates TPO to THR)
  resid <- pdbs$resid[i,f.inds$res]
  sequ  <- resid
  
  ## Build a dummy PDB to use with function nma.pdb()
  pdb.in <- as.pdb.default(xyz=tmp.xyz, elety="CA", 
                           resno=resno, chain=chain, resid=resid, verbose=FALSE)

  if(!is.null(outpath)) {
    fname <- file.path(outpath, basename(pdbs$id[i]))
    write.pdb(pdb.in, file=fname)
  }
  
  if(mass) {
    masses <- try(
      do.call('aa2mass', c(list(pdb=sequ, inds=NULL), am.args)),
      silent=TRUE
      )
    
    if(inherits(masses, "try-error")) {
      hmm <- attr(masses,"condition")
      cat("\n\n")
      stop(paste(hmm$message, "in file", basename(pdbs$id[i])))
    }
    masses.in <- masses
    masses.out <- masses
  }
  else {
    masses.in <- NULL
    masses.out <- NULL
  }

  natoms.in <- nrow(pdb.in$atom)
  natoms.out <- natoms.in
  
  if(rm.gaps) {
    ## Second PDB - containing only the aligned atoms
    sele <- list(atom=inds.inc, xyz=inds.inc.xyz)
    class(sele) <- "select"
    
    pdb.out <- trim.pdb(pdb.in, sele)
    natoms.out <- nrow(pdb.out$atom)

    if(mass)
      masses.out <- masses.in[ inds.inc ]

    inc.inds <- list(xyz=inds.inc.xyz)
  }
  else {
    pdb.out <- pdb.in
    inc.inds <- NULL
  }
  
  ## Build effective hessian
  bh.args <- c(list(sequ=sequ, fc.weights=wts[f.inds$res, f.inds$res]), ff.args)
  init <- list(pfcfun=pfc.fun, bh.args=bh.args)
  ##print(init$bh.args$fc.weights)

  invisible(capture.output( hessian <-
                           .nma.hess(pdb.in$xyz, init=init,
                                     hessian=NULL, inc.inds=inc.inds) ))

  ## Mass-weight hessian
  if(!is.null(masses.out))
    invisible(capture.output( hessian <- .nma.mwhessian(hessian, masses=masses.out)))
  
  ## Diagonalize
  invisible(capture.output( ei <- .nma.diag(hessian) ))
  
  ## Build an NMA object
  invisible(capture.output( modes <-
                           .nma.finalize(ei, xyz=pdb.out$xyz, temp=temp,
                                         masses=masses.out,
                                         natoms=natoms.out, keep=nm.keep, call=NULL) ))

  if(rm.gaps)
    modes.mat <- matrix(NA, ncol=keep, nrow=nrow(modes$U))
  else
    modes.mat <- matrix(NA, ncol=keep, nrow=ncol(pdbs$xyz))
  
  j <- 1
  for(k in 7:(keep+6)) {
    if(rm.gaps)
      modes.mat[, j] <- modes$U[,k]
    else
      modes.mat[f.inds$pos, j] <- modes$U[,k]
    j <- j+1
  }
  
  if(rm.gaps) {
    flucts <- modes$fluctuations
  }
  else {
    flucts <- rep(NA, length=ncol(pdbs$resno))
    flucts[f.inds$res] <- modes$fluctuations
  }

  ## Progress bar
  if(ncore>1) {
    iipb <- get("iipb", envir = env)
    iipb[1,i] <- 1
    j <- length(which(!is.na(bigmemory::as.matrix(iipb))))
  }
  else
    j <- i
  setTxtProgressBar(pb, j)

  L <- modes$L[seq(7, keep+6)]
  if(!full) {
    remove(modes)
    modes <- NULL
  }
  
  out <- list(modes=modes,
              U=modes.mat,
              L=L,
              flucts=flucts)
              ##f.inds=f.inds)
  
  invisible(gc())
  return(out)
}
  
  



