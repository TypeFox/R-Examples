rmsip <- function(...)
  UseMethod("rmsip")

rmsip.enma <- function(enma, ncore=NULL, subset=10, ...) {
  if(!inherits(enma, "enma"))
    stop("provide a 'enma' object as obtain from function 'nma.pdbs()'")
  if(any(is.na(enma$fluctuations)))
    stop("provide 'enma' object calculated with argument 'rm.gaps=TRUE'")
  
  ncore <- setup.ncore(ncore, bigmem = FALSE)
  
  if(ncore>1)
    mylapply <- mclapply
  else
    mylapply <- lapply
  
  dims <- dim(enma$fluctuations)
  m <- dims[1]

  mat <- matrix(NA, m, m)
  ##inds <- pairwise(m)
  inds <- rbind(pairwise(m),
                matrix(rep(1:m,each=2), ncol=2, byrow=T))
  
  mylist <- mylapply(1:nrow(inds), function(row) {
    i <- inds[row,1]; j <- inds[row,2];
    r <- rmsip.default(enma$U.subspace[,,i],
                      enma$U.subspace[,,j],
                      subset=subset)
    out <- list(val=r$rmsip, i=i, j=j)
    cat(".")
    return(out)
  })

  for ( i in 1:length(mylist)) {
    tmp <- mylist[[i]]
    mat[tmp$i, tmp$j] <- tmp$val
  }
  
  mat[ inds[,c(2,1)] ] = mat[ inds ]
  ##diag(mat) <- rep(1, n)
  colnames(mat) <- basename(rownames(enma$fluctuations))
  rownames(mat) <- basename(rownames(enma$fluctuations))

  cat("\n")
  return(round(mat, 6))
}

rmsip.default <- function(modes.a, modes.b, subset = 10,
                         row.name="a", col.name="b", ...) {
    
    if(missing(modes.a))
      stop("rmsip: 'modes.a' must be prodivded")
    if(missing(modes.b))
      stop("rmsip: 'modes.b' must be prodivded")
    
    m1 <- .fetchmodes(modes.a, subset=subset)
    m2 <- .fetchmodes(modes.b, subset=subset)

    dims.a <- dim(m1$U)
    dims.b <- dim(m2$U)
    subset <- dims.a[2]
    
    if( dims.a[1] != dims.b[1] )
      stop("dimension mismatch")
    if( dims.a[2] != dims.b[2] )
      stop("dimension mismatch")

    mass.a <- NULL; mass.b <- NULL;
    x <- normalize.vector(m1$U, mass.b)
    y <- normalize.vector(m2$U, mass.b)

    if(is.null(mass.a))
      o <- ( t(x) %*% y ) **2
    else
      o <- t(apply(x, 2, inner.prod, y, mass.a) **2)

    if (!is.null(row.name)) {
      rownames(o) <- paste(row.name, c(1:subset), sep="")
    }
    
    if (!is.null(col.name)) {
      colnames(o) <- paste(col.name, c(1:subset), sep="")
    }
    
    rmsip <- sqrt(sum(o)/subset)
    out <- list(overlap=round(o,3), rmsip=rmsip)

    class(out) <- "rmsip"
    return( out )
  }


.fetchmodes <- function(x, subset=NULL) {
  if (inherits(x, "pca")) {
    U <- x$U; L <- x$L;
    first.mode <- 1
  }
  else if (inherits(x, "nma")) {
    U <- x$U; L <- x$L;
    mass <- x$mass
    first.mode <- x$triv.modes+1
  }
  else {
    if( !(inherits(x, "matrix") | inherits(x, "pca.loadings")) )
      stop("provide an object of type 'pca', 'nma', or 'matrix'")
    U <- x; L <- NULL;
    first.mode <- 1
  }

  dims <- dim(U)
  if(is.null(subset)) {
    n <- dims[2]
  }
  else {
    n <- subset + first.mode - 1
    if(n>dims[2])
      n <- dims[2]
  }

  U <- U[, first.mode:n, drop=FALSE]
  L <- L[first.mode:n]
  
  out <- list(U=U, L=L, mass=NULL)
  return(out)
}
