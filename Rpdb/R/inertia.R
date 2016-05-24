## Compute the inertia tensor of a molecular system

inertia <- function(...)
  UseMethod("inertia")

inertia.coords <- function(x, m = NULL, ...){
  if(!is.coords(x))
    stop("'x' must be an object of class coords")
  if(is.null(m))
    stop("Please specify the masses")
  if(any(is.na(m))|any(is.na(x)))
    stop("NA values not permetted")
  Ixx<-sum(m*(x$x2^2+x$x3^2))
  Iyy<-sum(m*(x$x1^2+x$x3^2))
  Izz<-sum(m*(x$x1^2+x$x2^2))
  Ixy<-sum(m*(x$x1*x$x2))
  Ixz<-sum(m*(x$x1*x$x3))
  Iyz<-sum(m*(x$x2*x$x3))
  I<-matrix(c(Ixx,-Ixy,-Ixz,-Ixy,Iyy,-Iyz,-Ixz,-Iyz,Izz),ncol=3)
  return(I)
}

inertia.atoms <- function(x, m = NULL, ...){
  if(is.null(m))
    m <- masses(toSymbols(x$elename))
  inertia(coords(x), m)
}

inertia.pdb <- function(x, m = NULL, ...)
  inertia(x$atoms, m)
