# This file contains all the functions for relative importance.


# Internal Functions ----------------------------------------------------------- 


#' Relative weights
#' 
#' Function to implement Johnson's (2000) relative weight computation.
#'
#' @param rxx A matrix of predictor intercorrelations.
#' @param rxy A vector of predictor, criterion correlations.
#' @return DO THIS JEFF
#' @author Jeff Jones and Allen Goebl
#' @references Johnson, J. (2000). A heuristic method for estimating the
#'             relative weight of predictor variables in multiple regression.
#'             \emph{Multivariate Behavioral Research, 35}, 1--19.
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.relWt
.relWt <- function(rxx, rxy) {
    #Define Terms
    VLV <- eigen(rxx)
    V <- VLV$vectors
    L <- diag(VLV$values)
    #Equations
    lambda_star <- V %*% sqrt(L) %*% t(V)
    beta_star <- solve(lambda_star, rxy) #NOTE: Test qr.solve
    eps <- data.frame(EPS = (lambda_star^2) %*% (beta_star^2), stringsAsFactors=FALSE)
    #Return Values
    rownames(eps) <- rownames(rxx)
    return(list(eps = eps, beta_star=beta_star, lambda_star=lambda_star))
}


# External functions -----------------------------------------------------------


#' Relative weights
#' 
#' Function to implement Johnson's (2000) relative weight computation.
#'
#' @param r_mat A correlation matrix.
#' @param y_col A vector of columns representing criterion variables.
#' @param x_col A vector of columns representing predictor variables.
#' @return  A list containing the objects eps, beta_star, and lambda_star. The object eps contains 
#'          the vector of relative weights of the predictors whose sum is equivalent to the model \eqn{R^2}
#'          (see Johnson, 2000, ps 8 - 9). The object beta_star contains the regression weights from
#'          regressing the criterion on Z, the 'best fitting orthogonal approximation' of the predictor 
#'          variables (see Johnson, 2000, p. 5). The object lambda_star contains the regression coefficients 
#'          from regressing Z on the predictor variables (see Jonhson, 2000, p. 8).
#' @author Jeff Jones and Allen Goebl
#' @references Johnson, J. (2000). A heuristic method for estimating the
#'             relative weight of predictor variables in multiple regression.
#'             \emph{Multivariate Behavioral Research, 35}, 1--19.
#' @examples
#' Rs <- matrix(c(1.0, 0.2,  0.3, 0.4, -0.4,
#'                0.2, 1.0,  0.5, 0.1,  0.1,
#'                0.3, 0.5,  1.0, 0.2, -0.3,
#'                0.4, 0.1,  0.2, 1.0,  0.4,
#'               -0.4, 0.1, -0.3, 0.4,  1.0), 5, 5)
#' ys <- 5
#' xs <- 1:4
#' 
#' relWt(Rs, ys, xs)
#' @export
relWt <- function(r_mat, y_col, x_col) {
    #Check Input
    .checkIndex(r_mat=r_mat, y_col=y_col, x_col=x_col)
    #Call Functions
    term <- .indexMat(r_mat=r_mat, y_col=y_col, x_col=x_col)
    out <- .relWt(rxx=term$rxx, rxy=term$rxy)
    #Format Output
    return(out)
}
