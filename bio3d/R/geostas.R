#.getverbose <- function(...) {
#  dots <- list(...)
#  if(length(dots[names(dots) %in% "verbose"])==0)
#    verbose <- FALSE
#  else
#    verbose <- (dots[names(dots) %in% "verbose"])$verbose
#  return(verbose)
#}

geostas <- function(...)
  UseMethod("geostas")

geostas.nma <- function(nma, m.inds=7:11, verbose=TRUE, ...) {
  if(verbose) cat("  .. generating trajectory from", length(m.inds), "modes\n")
  
  trj <- NULL
  for(i in m.inds) {
    trj <- rbind(trj, mktrj.nma(nma, mode=i))
  }
  
  gs <- geostas.xyz(trj, verbose=verbose, ...)
  return(gs)
}

geostas.enma <- function(enma, pdbs=NULL, m.inds=1:5, verbose=TRUE, ...) {
  if(!inherits(enma, "enma"))
    stop("provide an 'enma' object as obtained by function 'nma.pdbs()'")
  if(!inherits(pdbs, "pdbs"))
    stop("provide an 'pdbs' object as obtained by function 'pdbaln' or 'read.fasta.pdb'")

  if(verbose) cat("  .. generating trajectory from", length(m.inds), "modes\n")

  trj <- mktrj.enma(enma, pdbs, m.inds=m.inds, rock=FALSE, ...)
  gs <- geostas.xyz(trj, verbose=verbose, ...)
  return(gs)
}

geostas.pdb <- function(pdb, inds=NULL, verbose=TRUE, ...) {
  if(!is.pdb(pdb))
    stop("provide a 'pdb' object as obtained from function 'read.pdb'")
  if(!nrow(pdb$xyz)>2)
    stop("provide a multi-model (>2) 'pdb' with more than ")

  if(is.null(inds)) {
    inds <- atom.select(pdb, "calpha")
    if(verbose) cat("  ..", length(inds$atom), "'calpha' atoms selected\n")
  }
  
  xyz <- pdb$xyz[,inds$xyz]
  gs <- geostas.xyz(xyz, verbose=verbose, ...)

  ## map back so that indices matches input 'pdb'
  if(verbose) {
    cat("  .. converting indices to match input 'pdb' object \n")
    cat("     (additional attribute 'atomgrps' generated) \n")
  }

  gs$fit.inds <- inds$xyz[gs$fit.inds]
  resid <- paste(pdb$atom$resid, pdb$atom$resno, pdb$atom$chain, sep="-")
  
  grps <- rep(NA, nrow(pdb$atom))
  for(i in 1:length(gs$inds)) {
    gs$inds[[i]] <- as.select(inds$atom[gs$inds[[i]]$atom])
    gs$inds[[i]]$call <- NA

    tmp.inds <- which(resid %in% resid[ gs$inds[[i]]$atom ])
    grps[ tmp.inds ] <- i
  }

  gs$atomgrps <- grps
  return(gs) 
}

geostas.pdbs <- function(pdbs, verbose=TRUE, ...) {
  if(!inherits(pdbs, "pdbs"))
    stop("provide an 'pdbs' object as obtained by function 'pdbaln' or 'read.fasta.pdb'")

  ## identify non-gap regions
  gaps.res <- gap.inspect(pdbs$ali)
  gaps.pos <- gap.inspect(pdbs$xyz)

  if(verbose) cat("  ..", length(gaps.pos$f.inds), "non-gap positions selected\n")
  
  xyz <- pdbs$xyz[, gaps.pos$f.inds]
  gs <- geostas.xyz(xyz, verbose=verbose, ...)

  ## map back so that indices matches input 'pdbs'
  gs$fit.inds <- gaps.pos$f.inds[gs$fit.inds]
  grps <- rep(NA, ncol(pdbs$ali))
  for(i in 1:length(gs$inds)) {
    gs$inds[[i]] <- as.select(gaps.res$f.inds[ gs$inds[[i]]$atom ])
    gs$inds[[i]]$call <- NA

    grps[ gs$inds[[i]]$atom ] <- i
  }
  
  return(gs) 
}

geostas.default <- function(...)
  geostas.xyz(...)

geostas.xyz <- function(xyz, amsm=NULL, k=3, pairwise=TRUE, clustalg="kmeans",
                        fit=TRUE, ncore=NULL, verbose=TRUE, ...) {
  cl <- match.call()
  
  xyz <- as.xyz(xyz)
  if(!nrow(xyz)>2)
    stop("provide a trajectory (e.g xyz object) with multiple (>2) frames")

  if(verbose) cat("  .. 'xyz' coordinate data with", nrow(xyz), "frames \n")
  
  if(!clustalg %in% c("hclust", "kmeans"))
    stop("'clustalg' should be 'kmeans' or 'hclust'")

  if(k<2)
    stop("provide 'k>1'")
  
  if(fit & is.null(amsm)) {
    if(verbose) cat("  .. 'fit=TRUE': running function 'core.find'\n")
    
    invisible(capture.output( core <- core.find(xyz) ))
    fit.inds <- core$xyz

    xyz <- fit.xyz(xyz[1,], xyz, fixed.inds=fit.inds, mobile.inds=fit.inds, ncore=ncore)
    if(is.null(fit.inds))
      warning("core indices not found. fitting to all atoms")
    
    if(verbose) cat("  .. coordinates are superimposed to core region\n")
  }
  else {
    if(verbose) cat("  .. coordinates are not superimposed prior to geostas calculation\n")
    fit.inds <- NULL
  }
  
  if(is.null(amsm)) {
    if(verbose) cat("  .. calculating atomic movement similarity matrix ('amsm.xyz()') \n")
    amsm <- amsm.xyz(xyz, ncore=ncore)
    dims <- dim(amsm)
    if(verbose) cat("  .. dimensions of AMSM are ", dims[1], "x", dims[2], "\n", sep="")
  }
  else {
    if(!all(dim(amsm)==ncol(xyz)/3))
      stop("dimension mismatch ('xyz' and 'amsm')")
  }

  if(pairwise) {
    cm <- 1-amsm
  }
  else {
    cm.tmp  <- normalize.vector(amsm)
    cm <- 1 - apply(cm.tmp, 2, function(x,y) x %*% y, cm.tmp)
  }
  
  ## hierarchical clustering
  if(clustalg=="hclust") {
    if(verbose) cat("  .. clustering AMSM using 'hclust' \n")
    dis    <- as.dist(cm)
    hc     <- hclust(dis, ...)
    grps   <- cutree(hc, k=k)
  }
  
  ## k-means clustering
  if(clustalg=="kmeans") {
    if(verbose) cat("  .. clustering AMSM using 'kmeans' \n")
    grps <- kmeans(cm, centers=k, ...)$cluster
  }

  ## return indices for the identified domains
  inds <- list()
  for(i in 1:length(unique(grps))) {
    inds[[i]] <- as.select(grps==i)
  }
  
  out <- list(call=cl, amsm=amsm, fit.inds=fit.inds, grps=grps, inds=inds)
  class(out) <- "geostas"
  return(out)
}

