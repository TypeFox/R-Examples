## Rotation of atomic coordinates

R <- function(...)
  UseMethod("R")

R.coords <- function(obj, angle = 0, x = 0, y = 0, z = 1, mask = TRUE, cryst1 = NULL, ...){
  if(!is.coords(obj)) stop("'object' must be an obj of class 'coords'")

  if(length(mask) != natom(obj)){
    if(length(mask) != 1)
      warning("'mask' has been recycled")
    mask <- rep(mask, length = natom(obj))
  }
  
  basis.ori <- basis(obj)
  if(basis.ori != "xyz"){
    if(is.null(cryst1))
      stop("Please specify a 'cryst1' obj to convert your fractional into Cartesian coordinates")
    obj <- abc2xyz(obj, cryst1 = cryst1)
  }
  M <- rotationMatrix(angle=angle*pi/180, x = x, y = y, z = z)[1:3,1:3]
  obj[mask,] <- coords(as.matrix(obj[mask,])%*%M, basis = "xyz")
  if(basis.ori != "xyz")
    obj <- xyz2abc(obj, cryst1 = cryst1)
  return(obj)
}

R.pdb <- function(obj, angle = 0, x = 0, y = 0, z = 1, mask = TRUE, cryst1 = obj$cryst1, ...){
  if(!is.pdb(obj)) stop("'object' must be an obj of class 'pdb'")
  
  coords(obj) <- R(coords(obj), angle = angle, x = x, y = y, z = z, mask = mask, cryst1 = cryst1, ...)
  
  return(obj)
}
