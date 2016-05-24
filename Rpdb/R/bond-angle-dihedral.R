bond <- function(...)
  UseMethod("bond")

bond.coords <- function(x, sel1, sel2, ...){
  if(!is.coords(x))
    stop("'x' must be an object of class 'coords'")
  if(basis(x) != "xyz")
    stop("'x' coordinates must be Cartesian coordinates")
  if(missing(sel1) | missing(sel2))
    stop("Please specify 'sel1' and 'sel2'")
  if(length(sel1) != length(sel2))
    stop("'sel1' and 'sel3' must have the same length")
  
  if(is.logical(sel1))
    sel1 <- which(sel1)
  if(is.logical(sel2))
    sel2 <- which(sel2)
  
  ill.bonds <- any(sel1==sel2)
  if(ill.bonds)
    warning("Ill defined bonds")
  
  at1 <- x[sel1,]
  at2 <- x[sel2,]
  if(natom(at1) != natom(at2))
    stop("'sel1' and 'sel2' must select the same number of atoms")
  
  B <- cbind(
    at2$x1 - at1$x1,
    at2$x2 - at1$x2,
    at2$x3 - at1$x3)
  B <- apply(B, 1, vectNorm)
  
  return(B)
}

bond.pdb <- function(x, sel1, sel2, ...){
  if(!is.pdb(x))
    stop("'x' must be an object of class 'pdb'")
  if(missing(sel1) | missing(sel2))
    stop("Please specify 'sel1' and 'sel2'")
  if(length(sel1) != length(sel2))
    stop("'sel1' and 'sel2' must have the same length")
  if(basis(x)=="abc"){
    if(is.null(x$cryst1))
      stop("'x' contains fractional coordinates but not 'cryst1' component for convertion to Cartesien coordinates")
    x <- abc2xyz(x)    
  }
  A <- bond.coords(x$atoms, sel1, sel2, ...)
  return(A)
}

angle <- function(...)
  UseMethod("angle")

angle.coords <- function(x, sel1, sel2, sel3, ...){
  if(!is.coords(x))
    stop("'x' must be an object of class 'coords'")
  if(basis(x) != "xyz")
    stop("'x' coordinates must be Cartesian coordinates")
  if(missing(sel1) | missing(sel2) | missing(sel3))
    stop("Please specify 'sel1','sel2' and 'sel3'")
  if(length(sel1) != length(sel2) | length(sel1) != length(sel3))
    stop("'sel1', 'sel2' and 'sel3' must have the same length")
  
  if(is.logical(sel1))
    sel1 <- which(sel1)
  if(is.logical(sel2))
    sel2 <- which(sel2)
  if(is.logical(sel3))
    sel3 <- which(sel3)
  
  ill.angles <- any(sel1==sel2) | any(sel2==sel3)
  if(ill.angles)
    warning("Ill defined angles")
  
  at1 <- x[sel1,]
  at2 <- x[sel2,]
  at3 <- x[sel3,]
  if(natom(at1) != natom(at2) | natom(at1) != natom(at3))
    stop("'sel1', 'sel2' and 'sel3' must select the same number of atoms")
  
  U <- cbind(
    at1$x1 - at2$x1,
    at1$x2 - at2$x2,
    at1$x3 - at2$x3)
  Un <- apply(U, 1, vectNorm)
  U <- U/Un
  
  V <- cbind(
    at3$x1 - at2$x1,
    at3$x2 - at2$x2,
    at3$x3 - at2$x3)
  Vn <- apply(V, 1, vectNorm)
  V <- V/Vn
  
  A <- (180/pi)*acos(rowSums(U*V))
  return(A)
}

angle.pdb <- function(x, sel1, sel2, sel3, ...){
  if(!is.pdb(x))
    stop("'x' must be an object of class 'pdb'")
  if(missing(sel1) | missing(sel2) | missing(sel3))
    stop("Please specify 'sel1','sel2' and 'sel3'")
  if(length(sel1) != length(sel2) | length(sel1) != length(sel3))
    stop("'sel1', 'sel2' and 'sel3' must have the same length")
  if(basis(x)=="abc"){
    if(is.null(x$cryst1))
      stop("'x' contains fractional coordinates but not 'cryst1' component for convertion to Cartesien coordinates")
    x <- abc2xyz(x)    
  }
  A <- angle.coords(x$atoms, sel1, sel2, sel3, ...)
  return(A)
}

dihedral <- function(...)
  UseMethod("dihedral")

dihedral.coords <- function(x, sel1, sel2, sel3, sel4, ...){
  if(!is.coords(x))
    stop("'x' must be an object of class 'coords'")
  if(basis(x) != "xyz")
    stop("'x' coordinates must be Cartesian coordinates")
  if(missing(sel1) | missing(sel2) | missing(sel3) | missing(sel4))
    stop("Please specify 'sel1','sel2', 'sel3' and 'sel4'")
  if(length(sel1) != length(sel2) |
       length(sel1) != length(sel3) |
       length(sel1) != length(sel4))
    stop("'sel1', 'sel2', 'sel3' and 'sel4' must have the same length")
  
  if(is.logical(sel1))
    sel1 <- which(sel1)
  if(is.logical(sel2))
    sel2 <- which(sel2)
  if(is.logical(sel3))
    sel3 <- which(sel3)
  if(is.logical(sel4))
    sel4 <- which(sel4)
  
  at1 <- x[sel1,]
  at2 <- x[sel2,]
  at3 <- x[sel3,]
  at4 <- x[sel4,]
  if(natom(at1) != natom(at2) |
       natom(at1) != natom(at3) |
       natom(at1) != natom(at4))
    stop("'sel1', 'sel2', 'sel3' and 'sel4' must select the same number of atoms")
  
  ill.dihedrals <- any(sel1==sel2) | any(sel2==sel3) |
    any(sel2==sel3) | any(sel3==sel4)
  if(ill.dihedrals)
    warning("Ill defined dihedrals")
  
  U <- cbind(
    at2$x1 - at1$x1,
    at2$x2 - at1$x2,
    at2$x3 - at1$x3)
  U <- U/apply(U, 1, vectNorm)
  
  V <- cbind(
    at3$x1 - at2$x1,
    at3$x2 - at2$x2,
    at3$x3 - at2$x3)
  V <- V/apply(V, 1, vectNorm)
  
  W <- cbind(
    at4$x1 - at3$x1,
    at4$x2 - at3$x2,
    at4$x3 - at3$x3)
  W <- W/apply(W, 1, vectNorm)
  
  N1 <- apply(cbind(U,V), 1,
              function(x){
                x <- vectProd(x[1:3], x[4:6])
                x <- x/vectNorm(x)
              })
  N2 <- apply(cbind(V,W), 1,
              function(x){
                x <- vectProd(x[1:3], x[4:6])
                x <- x/vectNorm(x)
              })
  
  D <- (180/pi)*acos(colSums(N1*N2))
  
  return(D)
}

dihedral.pdb <- function(x, sel1, sel2, sel3, sel4, ...){
  if(!is.pdb(x))
    stop("'x' must be an object of class 'pdb'")
  if(missing(sel1) | missing(sel2) | missing(sel3) | missing(sel4))
    stop("Please specify 'sel1','sel2', 'sel3' and 'sel4'")
  if(length(sel1) != length(sel2) |
       length(sel1) != length(sel3) |
       length(sel1) != length(sel4))
    stop("'sel1', 'sel2', 'sel3' and 'sel4' must have the same length")
  if(basis(x)=="abc"){
    if(is.null(x$cryst1))
      stop("'x' contains fractional coordinates but not 'cryst1' component for convertion to Cartesien coordinates")
    x <- abc2xyz(x)    
  }
  A <- dihedral.coords(x$atoms, sel1, sel2, sel3, sel4, ...)
  return(A)
}

