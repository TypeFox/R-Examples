## class DistanceCalculator --> calculate()
"amsm.xyz" <- function(xyz, ncore=NULL) {
  if(!is.matrix(xyz))
    stop("'xyz' must be a trajectory matrix")
  
  natoms <- ncol(xyz) / 3
  atom.pairs <- natoms * (natoms -1) / 2

  ## Distance Calculator
  ## returns a list of 4x4 matrices
  M <- .amsm.distCalc(xyz, ncore=ncore)

  ## we use 'list' in this version of the code
  if(class(M)=="list")
    Mlist <- TRUE
  else
    Mlist <- FALSE

  ## assign the atom mov. sim. matrix
  matrixCorr <- matrix(0, nrow=natoms, ncol=natoms)
  maxEigenVal <- 0
  
  ## solve eigenvalue problems
  ij <- combn(natoms,2)
  for ( i in 1:atom.pairs ) {
    atom.inds <- ij[,i]

    if(Mlist) {
      ev      <- eigen(M[[i]])
      maxDist <- M[[i]][1,1]
    }
    else {
      ev      <- eigen(M[,,i])
      maxDist <- M[1,1,i]
    }
    
    maxEigenVal <- max(ev$values)
    hei <- 1 - (sqrt(maxDist / maxEigenVal))
    matrixCorr[atom.inds[1], atom.inds[2]] <- hei
    matrixCorr[atom.inds[2], atom.inds[1]] <- hei

    ##matrixCorr[atom.inds[1], atom.inds[2]] <- maxDist
    ##matrixCorr[atom.inds[2], atom.inds[1]] <- maxDist
    
    ##if(max(ev$values)>maxEigenVal) 
    ##  maxEigenVal <- max(ev$values)
  }

  ##matrixCorr <- 1 - sqrt(matrixCorr / maxEigenVal)
  diag(matrixCorr) <- 1
  return(matrixCorr)
}


## class DistanceCalculator --> DistanceCalculator()
".amsm.distCalc" <- function(xyz, ncore=NULL) {

  if(!is.matrix(xyz))
    stop("'xyz' must be a trajectory matrix")

  ## Parallelized by package 'parallel'
  ncore <- setup.ncore(ncore, bigmem = FALSE)
    
  ## used for vectProdSum
  vectPS <- function(xyz.ab) {
    a <- xyz.ab[1:3]; b <- xyz.ab[4:6];
    m <- (a[2] * b[3]) - (a[3] * b[2])
    n <- (a[3] * b[1]) - (a[1] * b[3])
    o <- (a[1] * b[2]) - (a[2] * b[1])
    return(c(m,n,o))
  }

  ## used for matrixPSum
  matrPS <- function(xyz.ab) {
    a <- xyz.ab[1:3]; b <- xyz.ab[4:6];
    m <- a[1] * b + b[1] * a
    n <- a[2] * b + b[2] * a
    o <- a[3] * b + b[3] * a
    return(c(m,n,o))
  }

  for.atompair2 <- function(i, xyz, ij, M) {
    ij <- ij[,i]
    
    inds <- rep(ij*3,each=3) - c(2,1,0)
    xyz.ab <- xyz[, inds]

    vectProdSum <- rowSums(apply(xyz.ab, 1, vectPS))
    matrixPSum  <- rowSums(apply(xyz.ab, 1, matrPS))
    suma <- sum((xyz.ab[,1:3]-xyz.ab[,4:6])^2)

    ##M <- matrix(0, ncol=4, nrow=4)
    M[1,1] = suma
    M[1,2] = vectProdSum[1];
    M[1,3] = vectProdSum[2];
    M[1,4] = vectProdSum[3];
    
    M[2,2] = matrixPSum[1];
    M[2,3] = matrixPSum[2];
    M[2,4] = matrixPSum[3];
    
    M[3,3] = matrixPSum[5];
    M[3,4] = matrixPSum[6];
    M[4,4] = matrixPSum[9];
    
    ## symmerty
    M[2,1] = M[1,2]
    M[3,1] = M[1,3]
    M[4,1] = M[1,4]
    
    M[3,2] = M[2,3]
    M[4,2] = M[2,4]
    M[4,3] = M[3,4]
    
    return(M)
  }

  ##atom.pairs <- natoms * (natoms -1) / 2
  natoms <- ncol(xyz) / 3
  M      <- matrix(0, ncol=4, nrow=4)
  ij     <- combn(natoms,2)
  
  if(ncore==1) {
    all.Ms <- lapply(1:ncol(ij), for.atompair2, xyz, ij, M)
  }
  else {
    all.Ms <- mclapply(1:ncol(ij), for.atompair2, xyz, ij, M,
                       mc.cores=ncore)
  }
  return(all.Ms)
}
