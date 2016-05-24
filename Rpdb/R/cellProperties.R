#  Compute the Cartesian coordinates of lattice vectors.

cell.coords <- function(...)
  UseMethod("cell.coords")

cell.coords.default <- function(abc, abg = c(90,90,90), digits = 3, ...)
{
  if(missing(abc)) stop("Please provide at list a 'abc' vector containing the length of the lattice vectors")
  if(length(abc) != 3) stop("'abc' must be a vector of length 3")
  if(length(abg) != 3) stop("'abg' must be a vector of length 3")
  
  abg <- abg*pi/180
  
  M <- matrix(ncol=3,nrow=3)
  M[ ,1] <- c(abc[1],0,0)
  M[ ,2] <- c(abc[2]*cos(abg[3]),abc[2]*sin(abg[3]),0)
  M[1,3] <-   abc[3]*(cos(abg[2]))
  M[2,3] <-   abc[3]*(cos(abg[1])-cos(abg[2])*cos(abg[3]))/sin(abg[3])
  M[3,3] <-   abc[3]*sqrt(1+2*cos(abg[1])*cos(abg[2])*cos(abg[3])-(cos(abg[1]))^2-(cos(abg[2]))^2-(cos(abg[3]))^2)/sin(abg[3])
  M <- round(M, digits=3)
  dimnames(M) <- list(c("x","y","z"), c("a","b","c"))
  
  return(M)
}

cell.coords.cryst1 <- function(x, digits = 3, ...)
{
  if(!is.cryst1(x)) stop("'x' must be an object of class 'cryst1'")
  
  M <- cell.coords.default(x$abc, x$abg, digits)
  return(M)
}

cell.coords.pdb <- function(x, digits = 3, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'atoms'")
  if(is.null(x$cryst1)) stop("'x' must contained a 'cryst1' object")
  
  M <- cell.coords.cryst1(x$cryst1, digits = 3)
  return(M)
}

#  Compute the volume of a unit cell.

cell.volume <- function(...)
  UseMethod("cell.volume")

cell.volume.cryst1 <- function(x, ...)
{
  if(!is.cryst1(x)) stop("'x' must be an object of class 'cryst1'")
  
  V <-  prod(x$abc)*sqrt(1 - sum(cos(x$abg*pi/180)^2) + 2*prod(cos(x$abg*pi/180)))
  #   attr(V, "unit") <- "AngtromCube"
  return(V)
}

cell.volume.pdb <- function(x, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'atoms'")
  if(is.null(x$cryst1)) stop("'x' must contained a 'cryst1' object")
  
  V <- cell.volume.cryst1(x$cryst1)
  return(V)
}

#  Compute the density of a unit cell.

cell.density <- function(...)
  UseMethod("cell.density")

cell.density.default <- function(masses, volume, ...) {
  Na <- universalConstants["Na","Value"]
  d <- sum(masses)/(volume*1E-24*Na)
  #   attr(d, "unit") <- "g.cm-3"
  return(d)
}

cell.density.pdb <- function(x, ...) {
  M <- masses(x)
  V <- cell.volume(x)
  d <- cell.density(M, V)
  return(d)
}
