# =========== Convergence checks for weight functions ===========

#'Convergence checks
#'
#'The convergence check functions compare two weight matrices and returns the value of the
#'convergence criterion describing the difference between the two weight matrices.
#'
#'@param Wnew Current weight matrix
#'@param Wold Weight matrix from the previous iteration
#'
#'@return Scalar value of the convergence criterion
#'
#'@name convergenceCheck
NULL


#'@describeIn convergenceCheck  maximum of relative differences between
#'weights from two iterations
#'
#'@export
convCheck.relative <- function(Wnew, Wold){
  max(abs((Wold[Wnew != 0]-Wnew[Wnew != 0])/Wnew[Wnew != 0]))
}

#'@describeIn convergenceCheck maximum of squared absolute 
#'differences between weights from two iterations.
#'
#'@export

convCheck.square <- function(Wnew, Wold){
  max((Wold-Wnew)^2)
}

#'@describeIn convergenceCheck maximum of absolute differences
#'between weights from two iterations
#'
#'@export


convCheck.absolute <- function(Wnew, Wold){
  max(abs(Wold-Wnew))
}
