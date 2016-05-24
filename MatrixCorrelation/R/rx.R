#' @aliases r2 r3 r4
#' @title Correlational Measures for Matrices
#'
#' @description Matrix similarity as described by Ramsey et al. (1984).
#'
#' @param X1 first \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param X2 second \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param ncomp1 (GCD) number of subspace components from the first \code{matrix} (default: full subspace).
#' @param ncomp2 (GCD) number of subspace components from the second \code{matrix} (default: full subspace).
#'
#' @details Details can be found in Ramsey's paper:
#' \itemize{
#'  \item{r1:}{ inner product correlation}
#'  \item{r2:}{ orientation-independent inner product correlation}
#'  \item{r3:}{ spectra-independent inner product correlations (including orientation)}
#'  \item{r4:}{ Spectra-Independent inner product Correlations}
#'  \item{GCD:}{ Yanai's GCD Measure. To reproduce the original GCD, use all components.}
#' }
#'
#' @return A single value measuring the similarity of two matrices.
#'
#' @author Kristian Hovde Liland
#'
#' @references Ramsay, JO; Berg, JT; Styan, GPH; 1984. "Matrix Correlation". Psychometrica 49(3): 403-423.
#'
#' @seealso \code{\link{SMI}}, \code{\link{RV}} (RV2/RVadj).
#'
#' @examples
#' X1  <- matrix(rnorm(100*300),100,300)
#' usv <- svd(X1)
#' X2  <- usv$u[,-3] %*% diag(usv$d[-3]) %*% t(usv$v[,-3])
#'
#' r1(X1,X2)
#' r2(X1,X2)
#' r3(X1,X2)
#' r4(X1,X2)
#' GCD(X1,X2)
#' GCD(X1,X2, 5,5)
#'
#' @export
r1 <- function(X1, X2){
  Trace(crossprod(X1,X2))/sqrt(Trace(crossprod(X1))*Trace(crossprod(X2)))
}

#' @rdname r1
#' @export
r2 <- function(X1, X2){
  udv_x1 <- svd(X1)
  udv_x2 <- svd(X2)
  r1(udv_x1$u %*% diag(udv_x1$d),udv_x2$u %*% diag(udv_x2$d))
}

#' @rdname r1
#' @export
r3 <- function(X1, X2){
  udv_x1 <- svd(X1)
  udv_x2 <- svd(X2)
  r1(tcrossprod(udv_x1$u, udv_x1$v),tcrossprod(udv_x2$u, udv_x2$v))
}

#' @rdname r1
#' @export
r4 <- function(X1, X2){
  udv_x1 <- svd(X1)
  udv_x2 <- svd(X2)
  r1(udv_x1$u, udv_x2$u)
}

#' @rdname r1
#' @export
GCD <- function(X1, X2, ncomp1 = Rank(X1)-1, ncomp2 = Rank(X2)-1){
  udv_x1 <- svd(X1)
  udv_x2 <- svd(X2)
  r1(tcrossprod(udv_x1$u[,1:ncomp1]), tcrossprod(udv_x2$u[,1:ncomp2]))
}
