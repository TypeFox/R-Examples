"deformation.nma" <- function(nma, mode.inds=NULL, pfc.fun=NULL, ncore=NULL) {
  if(!inherits(nma, "nma"))
    stop("provide input of class 'nma' as obtained from function 'nma'")
  
  if(is.null(mode.inds)) {
    nmodes <- 20
    if(length(nma$L) < (nmodes+nma$triv.modes))
      nmodes <- length(nma$L) - nma$triv.modes
    mode.inds <- seq(nma$triv.modes+1, nma$triv.modes+nmodes)
  }
  else {
    if(length(nma$L) < (length(mode.inds)+nma$triv.modes)) {
      nmodes <- length(nma$L) - nma$triv.modes
      mode.inds <- seq(nma$triv.modes+1, nma$triv.modes+nmodes)
      warning("'mode.inds' was modified to include all modes")
    }
  }
    
  ## Check for multiple cores
  ncore <- setup.ncore(ncore)

  ## calcualte deformation energies for a specific mode
  def.mode <- function(mode.id, nma, xyz, ff, natoms) {
    mode.vec <- nma$modes[,mode.id]

    ## for normalization
    norm <- sqrt(sum(mode.vec**2) / natoms)

    ## better work with a matrix 
    mode.vec <- matrix(mode.vec, ncol=3, byrow=T)
  
    def.e <- rep(0, length=natoms)
    for ( i in 1:(natoms)) {
      
      ## distance vectors between a and b
      rs <- t(apply(xyz, 1, "-", xyz[i,]))
      
      ## differences mode vectors
      vs <- t(apply(mode.vec, 1, "-", mode.vec[i,]))
      
      ##rl <- apply(rs, 1, function(x) sqrt(sum(x^2)))
      rsq <- rowSums(rs^2)
      rl <- sqrt(rsq)
      
      ##l <- inner.prod(t(rs), t(vs)) / norm
      l <- colSums(t(rs)*t(vs)) / norm
            
      k <- ff(rl)
      l2 <- k*l*l / rsq
      
      l2[i] <- 0
      
      def.e <- def.e + (0.5*l2)
    }
    return(def.e)
  }

  if(is.null(pfc.fun))
    ff     <- load.enmff("calpha")
  else {
    if (!is.function(pfc.fun))
      stop("'pfc.fun' must be a function")
    ff <- pfc.fun
  }

  natoms <- length(nma$xyz)/3      
  xyz    <- matrix(nma$xyz, ncol=3, byrow=T)
  
  if(ncore>1)
    hm <- mclapply(mode.inds, def.mode, nma, xyz, ff, natoms, mc.cores=ncore)
  else
    hm <- lapply(mode.inds, def.mode, nma, xyz, ff, natoms)

  ## Collect results and make a matrix to store results
  hm <- t(do.call(rbind, hm))
  sums <- colSums(hm)

  out <- list(ei=hm, sums=sums)
  return(out)
}
