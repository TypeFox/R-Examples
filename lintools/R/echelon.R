#' Reduced row echelon form
#' 
#' Transform the equalities in a system of linear (in)equations or Reduced Row Echelon form (RRE)
#' 
#' 
#' @section Details:
#' 
#' The parameters \code{A}, \code{b} and \code{neq} describe a system of the form \code{Ax<=b}, where
#' the first \code{neq} rows are equalities. The equalities are transformed to RRE form.
#' 
#' 
#' A system of equations is in \href{https://en.wikipedia.org/wiki/Row_echelon_form}{reduced row echelon} form when
#' \itemize{
#' \item{All zero rows are below the nonzero rows}
#' \item{For every row, the leading coefficient (first nonzero from the left) is always right of the leading coefficient of the row above it.}
#' \item{The leading coefficient equals 1, and is the only nonzero coefficient in its column.}
#' }
#' 
#' 
#' @param A \code{[numeric]} matrix
#' @param b \code{[numeric]} vector
#' @param neq \code{[numeric]} The first \code{neq} rows of \code{A}, \code{b} are treated as equations.
#' @param nleq [\code{numeric}] The \code{nleq} rows after \code{neq} are treated as
#'   inequations of the form \code{a.x<=b}. All remaining rows are treated as strict inequations
#'   of the form \code{a.x<b}.
#' @param eps \code{[numeric]} Values of magnitude less than \code{eps} are considered zero (for the purpose of handling
#' machine rounding errors).
#' 
#' @return 
#' A list with the following components:
#' \itemize{
#'   \item{\code{A}: the \code{A} matrix with equalities transformed to RRE form.}
#'   \item{\code{b}: the constant vector corresponding to \code{A}}
#'   \item{\code{neq}: the number of equalities in the resulting system.}
#'   \item{\code{nleq}}: the number of inequalities of the form \code{a.x <= b}. This will only
#'   be passed to the output.
#' }
#' 
#' @examples 
#' echelon(
#'  A = matrix(c(
#'     1,3,1,
#'     2,7,3,
#'     1,5,3,
#'     1,2,0), byrow=TRUE, nrow=4)
#'  , b = c(4,-9,1,8)
#'  , neq=4
#' )
#' 
#' 
#' @export
echelon <- function(A, b, neq=nrow(A), nleq=0, eps=1e-8){
    check_sys(A=A,b=b,neq=neq,eps=eps)
  
    Ab <- cbind(A,b)
    
    ineq <- Ab[neq + seq_len(nrow(A)-neq),,drop=FALSE]
    Ab <- Ab[seq_len(neq),,drop=FALSE]
    
    k <- min(ncol(Ab),nrow(Ab))
    I <- seq_len(nrow(Ab))
    for ( i in 1:k ){
        I1 <- which(I >= i)
        ip <- I1[which.max(abs(Ab[I1,i]))]
        p <- Ab[ip,]
        if ( abs(p[i]) < eps ) next
        if ( ip > i ) Ab[c(ip,i),] <- Ab[c(i,ip),]
        Ab[-i,] <- Ab[-i,] - outer(Ab[-i,i],p/p[i])
    }
    
    d <- diag(Ab)
    id <- abs(d) > eps
    Ab[id,] <- Ab[id,]/d[id]
    I0 <- rowSums(abs(Ab) < eps) == ncol(Ab)
    Ab <- rbind(Ab[!I0,,drop=FALSE],Ab[I0,,drop=FALSE])
    i <- rowSums(abs(Ab) > eps) > 0 
    Ab <- rbind(Ab[i,,drop=FALSE],ineq)
    list(
      A = as.matrix(Ab[,1:ncol(A),drop=FALSE])
      , b = Ab[,ncol(A)+1,drop=TRUE]
      , neq = sum(i)
      , nleq=nleq
    )
}


