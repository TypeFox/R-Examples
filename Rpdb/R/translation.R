## Translation of atomic coordinates

## Translation of Cartesian coordinates
Txyz <- function(...)
  UseMethod("Txyz")

Txyz.coords <- function(obj, x = 0, y = 0, z = 0, mask = TRUE, thickness = NULL, cryst1 = NULL, ...){
  if(!is.coords(obj)) stop("'object' must be an obj of class 'coords'")
  
  if(length(mask) != natom(obj)){
    if(length(mask) != 1)
      warning("'mask' has been recycled")
    mask <- rep(mask, length = natom(obj))
  }
  
  v <- coords(x,y,z, basis = "xyz")
  T <- coords(0,0,0, basis = "xyz")
  if(basis(obj) != "xyz"){
    if(is.null(cryst1))
      stop("Please specify a 'cryst1' obj to convert your fractional into Cartesian coordinates")
    v <- xyz2abc(v, cryst1 = cryst1)
    T <- xyz2abc(T, cryst1 = cryst1)
  }

  vn <- coords(0,0,0, basis = "xyz")
  if(sqrt(sum(v^2)) != 0) vn <- v/sqrt(sum(v^2))
  
  if(!is.null(thickness)) {
    if(length(thickness) != 1) stop("'thickness must be a single element numeric vector'")
    T <- as.matrix(obj[mask,])%*%t(vn)
    T <- diff(range(T))*vn*thickness
  }
  
  obj$x1[mask] <- obj$x1[mask] + v$x1 + T$x1
  obj$x2[mask] <- obj$x2[mask] + v$x2 + T$x2
  obj$x3[mask] <- obj$x3[mask] + v$x3 + T$x3

  return(obj)
}

Txyz.pdb <- function(obj, x = 0, y = 0, z = 0, mask = TRUE, thickness = NULL, cryst1 = obj$cryst1, ...){
  if(!is.pdb(obj)) stop("'object' must be an obj of class 'pdb'")
  
  coords(obj) <- Txyz(coords(obj), x = x, y = y, z = z, mask = mask, thickness = thickness, cryst1 = cryst1, ...)
  
  return(obj)
}

## Translation of fractional coordinates
Tabc <- function(...)
  UseMethod("Tabc")

Tabc.coords <- function(obj, a = 0, b = 0, c = 0, mask = TRUE, thickness = NULL, cryst1 = NULL, ...){  
  if(!is.coords(obj)) stop("'object' must be an obj of class 'coords'")

  if(length(mask) != natom(obj)){
    if(length(mask) != 1)
      warning("'mask' has been recycled")
    mask <- rep(mask, length = natom(obj))
  }
  
  v <- coords(a,b,c, basis = "abc")
  T <- coords(0,0,0, basis = "abc")
  if(basis(obj) != "abc"){
    if(is.null(cryst1))
      stop("Please specify a 'cryst1' obj to convert your Cartesian into fractional coordinates")
    v <- abc2xyz(v, cryst1 = cryst1)
    T <- abc2xyz(T, cryst1 = cryst1)
  }

  vn <- coords(0,0,0, basis = "abc")
  if(sqrt(sum(v^2)) != 0) vn <- v/sqrt(sum(v^2))

  if(!is.null(thickness)) {
    if(length(thickness) != 1) stop("'thickness must be a single element numeric vector'")
    T <- as.matrix(obj[mask,])%*%t(vn)
    T <- diff(range(T))*vn*thickness
  }
  
  obj$x1[mask] <- obj$x1[mask] + v$x1 + T$x1
  obj$x2[mask] <- obj$x2[mask] + v$x2 + T$x2
  obj$x3[mask] <- obj$x3[mask] + v$x3 + T$x3
  
  return(obj)
}

Tabc.pdb <- function(obj, a = 0, b = 0, c = 0, mask = TRUE, thickness = NULL, cryst1 = obj$cryst1, ...){
  if(!is.pdb(obj)) stop("'object' must be an obj of class 'pdb'")
  
  coords(obj) <- Tabc(coords(obj), a = a, b = b, c = c, mask = mask, thickness = thickness, cryst1 = cryst1, ...)
  
  return(obj)
}
