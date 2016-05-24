## Set the view of the rgl scene

viewAxis <- function(V1, V2){
  if(is.list(V1)|is.list(V2))
    stop("'V1' and 'V2' can not be lists")
  if(length(V1)!=3|length(V2)!=3)
    stop("'V1' and 'V2' must be of length 3")
  if(vectNorm(V1)==0|vectNorm(V2)==0)
    stop("'V1' and 'V2' can not be null vectors")
  V1 <- matrix(c(V1,0))
  V2 <- matrix(c(V2,0))
  A <- acos(V1[1]/sqrt(sum(V1[c(1,2)]^2)))
  A <- ifelse(V1[2] > 0, -A, A)
  R1 <- rotationMatrix(A ,0 , 0, 1)
  V1 <- R1%*%V1
  V2 <- R1%*%V2
  A <- acos(V1[3]/sqrt(sum(V1[c(1,3)]^2)))
  A <- ifelse(V1[1] > 0, -A, A)
  R2 <- rotationMatrix(A ,0 , 1, 0)
  Vz <- R2%*%V1
  V2 <- R2%*%V2
  A <- acos(V2[1]/sqrt(sum(V2[c(1,2)]^2)))
  A <- ifelse(V2[2] > 0, -A, A)
  R3 <- rotationMatrix(A ,0 , 0, 1)
  par3d(userMatrix=R3%*%R2%*%R1)
}

## Helper functions

viewXY <- function()
  viewAxis(c(1,0,0),c(0,1,0))

viewYZ <- function()
  viewAxis(c(0,1,0),c(0,0,1))

viewZX <- function()
  viewAxis(c(0,0,1),c(1,0,0))

viewAB <- function(cryst1){
  if(missing(cryst1))
    stop("Please specify 'cryst1' to defined the lattice vectors used to set the view")
  cell <- cell.coords.cryst1(cryst1)
  viewAxis(cell[,"a"],cell[,"b"])
}

viewBC <- function(cryst1){
  if(missing(cryst1))
    stop("Please specify 'cryst1' to defined the lattice vectors used to set the view")
  cell <- cell.coords.cryst1(cryst1)
  viewAxis(cell[,"b"],cell[,"c"])
}

viewCA <- function(cryst1){
  if(missing(cryst1))
    stop("Please specify 'cryst1' to defined the lattice vectors used to set the view")
  cell <- cell.coords.cryst1(cryst1)
  viewAxis(cell[,"c"],cell[,"a"])
}

viewInertia <- function(x, m = NULL){
  M <- diag(4)
  M[1:3,1:3] <- eigen(inertia(x, m))$vectors
  par3d(userMatrix=M)
}
