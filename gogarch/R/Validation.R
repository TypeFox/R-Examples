##
## Validity function for objects of class "Goinit"
##
validGoinitObject <- function(object){
  m <- nrow(object@V)
  if(all.equal(diff(dim(object@V)), 0, check.attributes = FALSE)){
    TRUE
  } else {
    print("\nObject 'V' is not a square matrix.\n")
  }
  if(all.equal(diff(dim(object@P)), 0, check.attributes = FALSE)){
    TRUE
  } else {
    print("\nObject 'P' is not a square matrix.\n")
  }
  if(all.equal(diff(dim(object@Dsqr)), 0, check.attributes = FALSE)){
    TRUE
  } else {
    print("\nObject 'Dsqr' is not a square matrix.\n")
  }
  if(all.equal(det(object@Dsqr), prod(diag(object@Dsqr)), check.attributes = FALSE)){
    TRUE
  } else {
    print("\nObject 'Dsqr' is not a diagonal matrix.\n")
  }  
  if(all.equal(object@V, object@P %*% object@Dsqr^2 %*% t(object@P), check.attributes = FALSE)){ 
    TRUE
  } else {
    print("\nCovariance matrix cannot be replicated from singular values.\n")
  }
}
## Setting as validity function
setValidity("Goinit", validGoinitObject)
##
## Validity function for objects of class "Orthom"
##
validOrthomObject <- function(object){
  m <- nrow(object@M)
  Id <- diag(m)
  if(diff(dim(object@M)) == 0){
    TRUE
  } else {
    print("\nObject is not a square matrix.\n")
  }
  if(isTRUE(all.equal(1, abs(det(object@M))))){
    TRUE
  } else {
    print("\nAbsolute value of Determinant of object is not equal to 1.\n")
  }
  if(isTRUE(all.equal(Id, crossprod(object@M), check.attributes = FALSE))){
    TRUE
  } else {
    print("\nThe cross product of the object is not the Identity matrix.\n")
  }
}
## Setting as validity function
setValidity("Orthom", validOrthomObject)
