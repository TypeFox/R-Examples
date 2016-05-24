#' Find large correlations
#' 
#' Find any correlations for which their absolute value exceeds a
#' specified amount (rho). Returns a dataframe with row and column names and correlation
#' from lower triangular matrix.
#' 
#' @param x correlation matrix
#' @param rho abolute value for lower bound of correlation
#' @return a dataframe with 3 columns var1=row name,
#' var2= column name or number, Value of matrix element. Only contains rows
#' in which matrix element satisfies logical expression.
#' @export
#' @author Jeff Laake
find_large_cor=function(x,rho=0.9)
{
    if(!is.matrix(x))stop("x must be a matrix\n")
	if(nrow(x)!=ncol(x))stop("x must be square\n")
	if(any(t(x)[lower.tri(t(x))]!=x[lower.tri(x)]))stop("x must be symmetric\n")
	if(any(abs(x)>1.00001))stop("x does not appear to be a correlation matrix\n")
	if(is.null(rownames(x)))rownames(x)=1:nrow(x)
	if(is.null(colnames(x)))colnames(x)=1:ncol(x)
	data.frame(var1=rownames(x)[row(x)[abs(x)>rho&lower.tri(x)]],
			var2=colnames(x)[col(x)[abs(x)>rho&lower.tri(x)]],
			cor=x[abs(x)>rho&lower.tri(x)])
}
