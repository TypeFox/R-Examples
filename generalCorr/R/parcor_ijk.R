#' Generalized partial correlation coefficient between Xi and Xj removing the
#' effect of all other columns using the matrix R* of generalized correlation coefficients.
#'
#' This function uses a generalized correlation matrix R* as input to compute
#' generalized partial correlations between \eqn{X_i} and one of the remaining
#' variables after removing the effect of all other variables in the matrix.
#'
#' @param x {has a p by p matrix R* of generalized corr coefficients.}
#' @param i {A column number identifying the first variable.}
#' @param j {A column number identifying the second variable.}
#' @return 
#' \item{ouij}{partial correlation Xi on Xj (=cause) after removing all other X's}
#' \item{ouji}{partial correlation Xj on Xi (=cause) after removing all other X's}
#' \item{myk}{list of column numbers whose effect has been removed}
#' @note This function calls \code{\link{minor}}, and \code{\link{cofactor}} and is called 
#'   by \code{parcor_ridge}.
#' @examples 
#' 
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
#' colnames(x)=c("V1", "v2", "V3")
#' gm1=gmcmtx0(x)
#' parcor_ijk(gm1, 2,3)
#' }#' 
#' @export

parcor_ijk=function(x,i,j){
 n=NROW(x)
 p=NCOL(x)
 if (n<i) stop("n<i, parcor undefined")
if (p<j) stop("p<j, parcor undefined")
if (i<=0 | j<=0) stop("i OR j <=0, parcor undefined")
 myn=1:n
 myp=1:p
 myk=myp[c(-i,-j)]
 numij=det(cofactor(x,i,j))
 numji=det(cofactor(x,j,i))
deni=abs(det(cofactor(x,i,i)))
denj=abs(det(cofactor(x,j,j)))
ouij=(numij)/sqrt(deni*denj)
ouji=(numji)/sqrt(deni*denj)
list(ouij=ouij, ouji=ouji, myk=myk)}