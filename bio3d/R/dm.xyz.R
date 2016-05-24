dm.xyz <- function(xyz, grpby=NULL, scut=NULL, mask.lower=TRUE, ncore=1, ...) {
  ## Parallelized by parallel package
  ncore <- setup.ncore(ncore, bigmem = FALSE)

  if(ncore > 1) {
    mcparallel <- get("mcparallel", envir = getNamespace("parallel"))
    mccollect <- get("mccollect", envir = getNamespace("parallel"))
  }

  ## function for multicore calculation of dmats
  calcdm <- function(r.inds, core.id, xyz) {
    j <- 1
    out <- vector("list", length=length(r.inds))
    for(i in r.inds) {
      dmi <- .dm.xyz1(xyz[i,], grpby=grpby, scut=scut, mask.lower=mask.lower)
      out[[j]] <- dmi
      j <- j+1
    }
    return(out)
  }
  
  xyz <- as.xyz(xyz)
  if(nrow(xyz)>1) {
    ## dimensions of final array
    d3 <- nrow(xyz)
    
    if(!is.null(grpby))
      d1 <- length(unique(grpby))
    else
      d1 <- ncol(xyz)/3
    
    dms <- array(data=0.00, dim=c(d1, d1, d3))

    ## multicore setup
    if(ncore>1) {
      jobs <- list()
    }

    ## run calcdm() for each core
    core.ids <- sort(rep(1:ncore, length.out=d3))
    for( i in 1:ncore ) {
      r.inds <- which(core.ids==i)

      if(ncore>1) {
        q <- mcparallel(calcdm(r.inds, i, xyz))
        jobs[[i]] <- q
      }
      else {
        dm.list <- calcdm(r.inds, i, xyz)
      }
    }

    ## Collect all jobs
    if(ncore>1) 
      res <- mccollect(jobs, wait=TRUE)
    else
      res <- list(dm.list)

    ## list to array
    i <- 1
    for ( job in res ) {
      for(mat in job) {
        dms[,,i] <- mat
        i <- i+1
      }
    }  
  }
  else {
    dms <- .dm.xyz1(xyz, grpby=grpby, scut=scut, mask.lower=mask.lower)
  }
  return(dms)
}

.dm.xyz1 <-
function(xyz, grpby=NULL, scut=NULL, mask.lower=TRUE) {
  ##-- New distance matrix function with 'grpby' option
  ##  dm(pdb$xyz, grpby=pdb$atom[,"resno"], scut=3)
  
  xyz=as.xyz(xyz)
  if(dim(xyz)[1L]>1) {
    warning("multiple frames detected - using only the first row in xyz matrix")
    xyz = xyz[1,, drop=FALSE]
  }
  xyz=as.vector(xyz)
  
  ##- Full Distance matrix (could use 'dm' or 'dist.xyz')
  d <- as.matrix(dist(matrix(xyz, ncol = 3, byrow = TRUE)))
  
  ##- Mask lower.tri  
  if( mask.lower )
    d[lower.tri(d)] = NA

  ##- Mask concetive atoms
  if( is.null(grpby) ) {
    if (!is.null(scut)) {
      d[diag.ind(d, n = scut)] = NA
      if(!mask.lower) 
         d[lower.tri(d)] = t(d)[lower.tri(d)]
    }

    return(d)
    
  } else {

    ##- Group by concetive numbers in 'grpby'
    if( length(xyz) != (length(grpby)*3) )
      stop("dimension miss-match in 'xyz' and 'grpby', check lengths")

    ##- Bounds of 'grpby' numbers
    inds <- bounds(grpby, dup.inds=TRUE)
    nres <- nrow(inds)
    
    ##- Per-residue matrix
    m <- matrix(, ncol=nres, nrow=nres)
    ij <- pairwise(nres)
    
    ##  Ignore concetive groups (residues)
    if (!is.null(scut))
      ij <- ij[ij[,2]-ij[,1] > (scut-1),]
    
    ##- Min per residue
    for(k in 1 : nrow(ij) ) {
      m[ij[k,1],ij[k,2]] <-
        min( d[ (inds[ij[k,1],"start"]:inds[ij[k,1],"end"]),
                  (inds[ij[k,2],"start"]:inds[ij[k,2],"end"])],
            na.rm=TRUE )
    }
    if( !mask.lower )
      m[lower.tri(m)] = t(m)[lower.tri(m)]
    
    return(m)
  
  }
}

