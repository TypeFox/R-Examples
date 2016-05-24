#' Function to do compute the minor of a matrix defined by row r and column c.
#' 
#' @param x {the input matrix}
#' @param r {the row number}
#' @param c {the column number}
#' @return the appropriate `minor' matrix from the original matrix.
#' 
#' @note this is needed by the cofactor function
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @examples
#' 
#' \dontrun{
#'  x=matrix(1:20,ncol=4)
#' minor(x,1,2)}
#' 
#' @export

minor <-
function(x,r,c){
 #x is n by p matrix we want its minor
 #after eliminating rth row and cth column
 n=nrow(x)
 p=ncol(x)
 myn=1:n
 myp=1:p
if (n<r) stop("n<r, minor undefined")
if (p<c) stop("p<c, minor undefined")
if (c<=0) stop("c<=0, minor undefined")
 newr=myn[-r]
 newc=myp[-c]
 out=x[newr,newc]
return(out)}
