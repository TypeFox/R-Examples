## Mirror operation on atomic coordinates

mirror <- function(...)
  UseMethod("mirror")

mirror.coords <- function(x, p1, p2 = NULL, p3 = NULL, mask = TRUE, cryst1 = NULL, ...){
  if(missing(p1))
    stop("Please specify at least 'p1'")
  if(is.null(p2) & is.null(p3)){
    if(ncol(p1)!=3 | nrow(p1)!=3)
      stop("When 'p2' and 'p3' are not specifyed, 'p1' must be a 3x3 matrix or data.frame")
    p3 <- p1[3,]
    p2 <- p1[2,]
    p1 <- p1[1,]
  } else {
    if(length(p3) != 3 | length(p2) != 3 | length(p1) != 3)
      stop("'p1', 'p2' and 'p3' must be vectors of length 3")
  }
  if(all(p1==p2)|all(p1==p3)|all(p2==p3))
    stop("'p1', 'p2' and 'p3' must be different to define the mirror")
  if(length(mask) != natom(x)){
    if(length(mask) != 1)
      warning("'mask' has been recycled")
    mask <- rep(mask, length = natom(x))
  }
  basis.ori <- basis(x)
  if(basis.ori != "xyz"){
    if(is.null(cryst1))
      stop("Please specify a 'cryst1' obj to convert your fractional into Cartesian coordinates")
    x <- abc2xyz(x, cryst1 = cryst1)
  }
  v12 <- p2 - p1  
  v23 <- p3 - p2

  vn <- vectProd(v12,v23)
  vn <- vn/vectNorm(vn)
  x <- Txyz(x,  p1[1],  p1[2],  p1[3], mask=mask)
  rotM <- diag(3) - as.matrix(
    rbind(
      c(2*vn[1]*vn[1],2*vn[1]*vn[2],2*vn[1]*vn[3]),
      c(2*vn[2]*vn[1],2*vn[2]*vn[2],2*vn[2]*vn[3]),
      c(2*vn[3]*vn[1],2*vn[3]*vn[2],2*vn[3]*vn[3])
      )
    )
  x[mask,] <- coords(as.matrix(x[mask,])%*%rotM, basis = "xyz")
  x <- Txyz(x, -p1[1], -p1[2], -p1[3], mask=mask)
  return(x)
}

mirror.pdb <- function(x, p1, p2 = NULL, p3 = NULL, mask = TRUE, cryst1 = x$cryst1, ...){
  coords(x) <- mirror(coords(x), p1=p1, p2=p2, p3=p3, mask=mask, cryst1=cryst1, ...)
  return(x)
}
