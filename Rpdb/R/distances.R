distances <- function(...)
  UseMethod("distances")

distances.default <- function(dx1 = numeric(0), dx2 = numeric(0), dx3 = numeric(0), basis = "xyz", ...){
  if(!is.numeric(dx1) | !is.numeric(dx2) | !is.numeric(dx3))
    stop("'dx1', 'dx2' and 'dx3' must be numeric")
  if(is.null(dim(dx1))){
    if(!is.null(dim(dx2)) & !is.null(dim(dx3)))
      stop("'dx1', 'dx2' and 'dx3' must have the same dimensions")
    if(length(dx1) != length(dx2) | length(dx1 != length(dx3)))
      stop("'dx1', 'dx2' and 'dx3' must have the same dimensions")
  }
  else{
    if(any(dim(dx1) != dim(dx2)) | any(dim(dx1) != dim(dx3)))
      stop("'dx1', 'dx2' and 'dx3' must have the same dimensions")    
  }
  if(!basis %in% c("xyz","abc"))
    stop("'basis' must be equal to 'xyz' or 'abc'")
  
  to.return <- list(dx1 = dx1, dx2 = dx2, dx3 = dx3)
  class(to.return) <- c("distances", "list")
  attr(to.return, "basis") <- basis
  return(to.return)
}

distances.coords <- function(x, sel1, sel2, ...){
  if(missing(sel1) | missing(sel2))
    stop("Please specify 'sel1' and 'sel2'")
  if(is.numeric(sel1)){
    if(any(round(sel1) != sel1))
      stop("'sel1' must be a logical or integer vector")
    if(any(sel1 > natom(x)))
      stop("'sel1' contains indices out of range")
  }
  if(is.numeric(sel2)){
    if(any(round(sel2) != sel2))
      stop("'sel2' must be a logical or integer vector")
    if(any(sel2 > natom(x)))
      stop("'sel2' contains indices out of range")
  }
  if(is.logical(sel1) & length(sel1) != natom(x))
    stop("'sel1' length must be equal to natom(x)")
  if(is.logical(sel2) & length(sel2) != natom(x))
    stop("'sel2' length must be equal to natom(x)")  

  xyz1 <- x[sel1,]
  xyz2 <- x[sel2,]
  
  dx1 <- t(outer(xyz2$x1, xyz1$x1, "-"))
  dx2 <- t(outer(xyz2$x2, xyz1$x2, "-"))
  dx3 <- t(outer(xyz2$x3, xyz1$x3, "-"))
  
  dimnames(dx1) <- list(sel1 = NULL, sel2 = NULL)
  dimnames(dx2) <- list(sel1 = NULL, sel2 = NULL)
  dimnames(dx3) <- list(sel1 = NULL, sel2 = NULL)
  
  to.return <- distances.default(dx1, dx2, dx3, basis = basis(x))
  return(to.return)
}

distances.atoms <- function(x, sel1, sel2, ...)
  distances.coords(coords(x), sel1 = sel1, sel2 = sel2, ...)

distances.pdb <- function(x, sel1, sel2, ...)
  distances.atoms(x$atoms, sel1 = sel1, sel2 = sel2, ...)

is.distances <- function(x)
  any(class(x) == "distances")

norm <- function(...)
  UseMethod("norm")

norm.distances <- function(x, type = "xyz", ...){
  if(basis(x) == "abc")
    stop("Please provide Cartesian coordinates. See 'abc2xyz'.")
  
  to.return <- switch(type,
    x   = x$dx1,
    y   = x$dx2,
    z   = x$dx3,
    xy  = sqrt(x$dx1^2+x$dx2^2),
    yz  = sqrt(x$dx2^2+x$dx3^2),
    zx  = sqrt(x$dx3^2+x$dx1^2),
    xyz = sqrt(x$dx1^2+x$dx2^2+x$dx3^2))

  return(to.return)
}

