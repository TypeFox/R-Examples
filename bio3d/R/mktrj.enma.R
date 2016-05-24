"mktrj.enma" <- function(enma=NULL,   # enma data structure
                         pdbs=NULL,   # pdbs object 
                         s.inds=NULL, # structure ids
                         m.inds=NULL, # modes ids
                         mag=10,      # magnification factor
                         step=1.25,   # step size
                         file=NULL,   # output pdb file
                         rock=TRUE, ncore=NULL,
                         ... ) {      # args for write.pdb

  ## make a trjactory of atomic displacments along a given mode
  if(!inherits(enma, "enma"))
    stop("mktrj.enma: must supply 'enma' object, i.e. from 'nma.pdbs'")

  ## Parallelized by parallel package
  ncore <- setup.ncore(ncore, bigmem = FALSE)
  
  if(ncore>1)
    mylapply <- mclapply
  else
    mylapply <- lapply
  
  if(is.null(s.inds))
    s.inds <- 1:nrow(enma$fluctuations)

  if(is.null(m.inds))
    m.inds <- 1:5
  
  if(is.null(file) & length(s.inds)==1 & length(m.inds)==1)
    file <- paste("mode_", m.inds+6, "-s", s.inds, ".pdb", sep="")
  
  if(is.null(enma$call$rm.gaps))
    rm.gaps <- TRUE
  else if(enma$call$rm.gaps=="T" || enma$call$rm.gaps=="TRUE")
    rm.gaps <- TRUE
  else
    rm.gaps <- FALSE

  if(!rm.gaps & length(s.inds)>1 & length(m.inds)>1)
    stop(paste("enma object must be calculated with argument rm.gaps=TRUE", "\n",
               "for trajectory generation of multiple structures and modes"))

  if(any(enma$L[s.inds, m.inds]<=0))
    warning("Mode with eigenvalue <=0 detected. Check 'mode' index.")

  nstep <- c(seq(step, to=mag, by=step))
  zcoor <- cbind(1) %*% nstep
  scor  <- function(x,u,m) { return(x*u+m) }
  
  myMktrj <- function(i) {
    coor <- NULL
    ind <- s.inds[i]
    
    for(j in 1:length(m.inds)) {
      mode <- m.inds[j]
      
      u.inds <- which(!is.na(enma$U.subspace[,mode,ind]))
      if(rm.gaps)
        xyz.inds <- gap.inspect(enma$xyz)$f.inds
      else
        xyz.inds <- u.inds
      
      plus  <- sapply(c(zcoor), scor, u=enma$U.subspace[u.inds,mode,ind], m=enma$xyz[ind,xyz.inds])
      minus <- sapply(c(-zcoor), scor, u=enma$U.subspace[u.inds,mode,ind], m=enma$xyz[ind,xyz.inds])
      
      if(rock) {
        tmp  <- cbind(pdbs$xyz[ind,xyz.inds],
                      plus,  plus[,rev(1:ncol(plus))],
                      enma$xyz[ind,xyz.inds],
                      minus, minus[,rev(1:ncol(minus))])
      }
      else {
        tmp  <- cbind(plus[,rev(1:ncol(plus))],
                      enma$xyz[ind,xyz.inds],
                      minus)
      }
      
      coor <- rbind(coor, t(tmp))
    }
    return(coor)
  }

  ## do the calc
  coor <- mylapply(1:length(s.inds), myMktrj)
  coor <- do.call(rbind, coor)
  class(coor) <- "xyz"
  
  if(!is.null(file)) {
    if(rm.gaps)
      xyz.inds <- gap.inspect(enma$xyz)$f.inds
    else
      xyz.inds <- which(!is.na(enma$U.subspace[,m.inds[1],s.inds[1]]))

    if(is.null(pdbs)) 
      write.pdb(xyz=coor, file=file, ...)
    else {
      write.pdb(xyz=coor, file=file, 
                chain=pdbs$chain[s.inds[1], xyz2atom(xyz.inds)],
                resno=pdbs$resno[s.inds[1], xyz2atom(xyz.inds)],
                resid=pdbs$resid[s.inds[1], xyz2atom(xyz.inds)],
                b=enma$fluctuations[s.inds[1], !is.gap(enma$fluctuations[s.inds[1],])],
                ...)
    }
    invisible(coor)
  }
  else {
    return(coor)
  }
}
