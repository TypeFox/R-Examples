"mktrj.pca" <- function(pca=NULL,   # pca data structure
                        pc=1,       # which pc to move along
                        mag=1,      # magnification factor
                        step=0.125, # step size
                        file=NULL,  # output pdb file
                        ... ) {     # args for write.pdb

  ## make a trjactory of atomic displacments along a given pc
  if(class(pca)!="pca") {
    stop("input should be a list object of class 'pca' (from 'pca.xyz')")
  }

  if(is.null(file))
    file <- paste("pc_", pc, ".pdb", sep="")
  
  nstep <- c(seq(step, to=mag, by=step))
  zcoor  <- cbind(sqrt(pca$L[pc])) %*% nstep

  ##- Bug fix: Fri Jun 15 14:49:24 EDT 2012
  ##  plus  <- apply(zcoor, 2, pca.z2xyz, pca)
  ##  minus <- apply( (-(zcoor)), 2, pca.z2xyz, pca)

  scor  <- function(x,u,m) { return(x*u+m) }
  plus  <- sapply(c(zcoor), scor, u=pca$U[,pc], m=pca$mean)
  minus <- sapply(c(-zcoor), scor, u=pca$U[,pc], m=pca$mean)

  coor  <- t(cbind(pca$mean,
                   plus, plus[,rev(1:ncol(plus))],
                   pca$mean,
                   minus, minus[,rev(1:ncol(minus))]))

  write.pdb(xyz=coor, file=file, ...)
  invisible(coor)
}
