#' Rotate X to match Y
#'
#' Function to rotate X to match Y via singular value decomposition
#'
#' @param X matrix to be rotated
#' @param Y objective matrix 
#' @return rotated object \code{Xrot}, and the rotation matrix \code{R}
#' @export

rotXtoY<-function(X, Y){
	SVD <- svd(t(Y) %*% X) 
	R <- SVD$v %*% t(SVD$u)
	return(list(Xrot = X%*%R, R = R))
	}