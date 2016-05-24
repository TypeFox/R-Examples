### functions for QR or sparseQR class
qrrank<-function(QR, tol=1e-8){
  if(class(QR)=="sparseQR")
    return(sum(abs(diag(QR@R))>tol))
  else
    return(QR$rank)
}
qrRM<-function(QR){
  if(class(QR)=="sparseQR")
    return(qrR(QR))
  else
    return(qr.R(QR))
}

# invTwobyTwoSymMatrix: inverts two-by-two symmetric matrices
# M: symmetric matrix
# OUTPUT: a two-by-two symmetric matrix
invTwobyTwoSymMatrix <- function(M) {
  if(missing(M)) stop("invTwobyTwoSymMatrix: Supply M!") 
  if(nrow(M) != 2 || ncol(M) != 2) stop("invTwobyTwoSymMatrix: The matrix is not 2 by 2")
  if(M[1,2] != M[2,1]) stop("invTwobyTwoSymMatrix: The matrix is not symmetric!")
  detM = M[1,1] * M[2,2] - M[1,2]^2
  Minverse = M
  Minverse[2,2] = M[1,1]
  Minverse[1,1] = M[2,2]
  Minverse[1,2] = Minverse[2,1] = -1 * M[1,2]
  return(1/detM * Minverse)
}

# quadSolver: solves the quadratic equation ax^2 + bx + c = 0
#             Can also be used to solve eigenvalues of two-by-two matrices
#             E.x. if [d,e; f,g] is the two-by-two matrix
#                  eigenvalues are defined as the solution to det( [d,e; f,g] - lambda * I) = 0
#                  Or det([d - lambda, e; f, g - lambda]) = 0
#                  Or (d-lambda) * (g -lambda) - e*f = 0
#                  Or lambda^2 - (d +g) *lambda + d*g - e*f = 0
quadSolver = function(a,b,c) {
  if(missing(a) || missing(b) || missing(c)) stop("Supply the coefficients!")
  if(length(a) != 1 || length(b) != 1 || length(c) != 1) stop("Supply only scalar values to a,b, and c")
  if(a == 0) stop("Something is wrong! Leading coefficient is zero")
  
  discriminant = b^2 - 4 * a * c
  if(discriminant < 0) return(c(NA,NA))
  if(discriminant == 0) return( rep(-b/(2*a),2))
  if(discriminant > 0) {
    if(a > 0) {
	  lowerRoot = -b / (2*a) - sqrt(discriminant) / (2*a)
	  upperRoot = -b / (2*a) + sqrt(discriminant) / (2*a)
	  return(c(lowerRoot,upperRoot))
	} else {
	  lowerRoot = -b / (2*a) + sqrt(discriminant) / (2*a)
	  upperRoot = -b / (2*a) - sqrt(discriminant) /(2*a)
	  return(c(lowerRoot,upperRoot))
	}
  }
}

