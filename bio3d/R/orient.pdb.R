"orient.pdb" <-
function (pdb, atom.subset = NULL, verbose = TRUE ) {

  ## x <- c(rep(10,3), rep(0,3), rep(-10,3))
  ## write.pdb(xyz=x, file="t1.pdb")
  ## write.pdb(xyz=orient.pdb(x), file="t2.pdb")

  if (missing(pdb)) {
    stop("pdb.orient: must supply 'pdb' object, e.g. from 'read.pdb'")
  }
  
  if(is.list(pdb)) { xyz <- pdb$xyz
  } else {
    if (!is.vector(pdb)) {
      stop("pdb.orient: input 'pdb' should NOT be a matrix")
    }
    xyz <- pdb
  }
  xyz <- matrix( xyz, ncol=3, byrow=TRUE )   
  
  if (is.null(atom.subset)) atom.subset <- c(1:nrow(xyz))
  if (length(atom.subset) > nrow(xyz)) {
    stop("pdb.orient: there are more 'atom.subset' inds than there atoms")
  }
  

  ## Center on mean xyz positions
  xyz.bar <- apply(xyz[atom.subset, ], 2, mean)
  xyz <- sweep(xyz, 2, xyz.bar)
  
  ## Determine principal axis
  S <- var(xyz[atom.subset, ])
  prj <- eigen(S, symmetric = TRUE)
  
  ## Mke rotation explicitly rh system
  ## z <- xyz %*% (prj$vectors)
  A <- prj$vectors
  b <- A[,1]; c <- A[,2]
  A[1,3] <- (b[2] * c[3]) - (b[3] * c[2])
  A[2,3] <- (b[3] * c[1]) - (b[1] * c[3])
  A[3,3] <- (b[1] * c[2]) - (b[2] * c[1])
  
  ## Rotate
  z <- xyz %*% (A)
  
  if (verbose) {
    cat("Dimensions:", "\n")
    cat(" x  min=",  round(min(z[, 1]), 3),
        "   max=",   round(max(z[, 1]), 3),
        "   range=", round(max(z[, 1]) - min(z[, 1]), 3), "\n")
    
    cat(" y  min=", round(min(z[, 2]), 3),
        "   max=", round(max(z[, 2]), 3),
        "   range=", round(max(z[, 2]) - min(z[, 2]), 3), "\n")
    
    cat(" z  min=", round(min(z[, 3]), 3),
        "   max=", round(max(z[, 3]), 3),
        "   range=", round(max(z[, 3]) - min(z[, 3]), 3), "\n")
    
  }
  z <- round(as.vector(t(z)),3)
  z <- as.xyz(z)
  invisible(z)
}
