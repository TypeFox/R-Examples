"com.xyz" <- function(xyz, mass=NULL, ...) {
  xyz <- as.xyz(xyz)
  natoms <- ncol(xyz)/3
  
  if(is.null(mass))
    mass <- rep(1, times=natoms)

  if (natoms != length(mass))
    stop("com.xyz: length of input vector 'mass' uequal to number of atoms (ncol(xyz)/3)")

  com1 <- function(x) {
    xyz <- matrix(x, ncol=3, byrow=T)
    com <- colSums(xyz * mass) / sum(mass)
    return(com)
  }

  com <- t(apply(xyz, 1, com1))
  colnames(com) <- c("x", "y", "z")
  return(com)
}
