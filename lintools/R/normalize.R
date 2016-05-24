
#' Bring a system of (in)equalities in a standard form
#' 
#' @section Details:
#' For this package, a set of equations is in normal form when 
#' \itemize{
#' \item The first \code{neq} rows represent linear equalities.
#' \item The next \code{nleq} rows represent inequalities of the form \code{a.x <= b}
#' \item All other rows are strict inequalities of the form \code{a.x < b}
#' }
#' 
#' If \code{unit>0}, the strict inequalities \code{a.x < b} are replaced with 
#' inequations of the form \code{a.x <= b-unit}, where \code{unit} represents
#' the precision of measurement.
#' 
#' 
#' @param A \code{[numeric]} Matrix
#' @param b \code{[numeric]} vector
#' @param operators \code{[character]} operators in \code{{<,<=,==,>=,>}}.
#' @param unit \code{[numeric]} (nonnegative) Your unit of measurement. This is used to
#' replace strict inequations of the form \code{a.x < b} with \code{a.x <= b-unit}.
#' Typically, \code{unit} is related to the units in which your data 
#' is measured.  If unit is 0, inequations are not replaced.
#' 
#' @return A \code{list} with the folowing components
#' 
#' \itemize{
#'   \item{\code{A}: the \code{A} corresponding to the normalized sytem.}
#'   \item{\code{b}: the constant vector corresponding to the normalized system}
#'   \item{\code{neq}: the number of equations}
#'   \item{\code{nleq}}: the number of non-strict inequations (<=)
#'   \item{\code{order}: the index vector used to permute the original rows of \code{A}.}
#' }
#' 
#' @examples 
#' 
#' A <- matrix(1:12,nrow=4)
#' b <- 1:4
#' ops <- c("<=","==","==","<")
#' normalize(A,b,ops)
#' normalize(A,b,ops,unit=0.1)
#' 
#' @export
normalize <- function(A, b, operators, unit=0 ){
  stopifnot(unit >= 0)
  geq <- operators  == ">="
  gt <- operators == ">"
  A[geq|gt,] <- -A[geq|gt,,drop=FALSE]
  b[geq|gt] <- -b[geq|gt]
  if (unit > 0){
    b[gt] <- b[gt] - unit
    operators[geq|gt] <- "<="
    lt <- operators == "<"
    b[lt] <- b[lt] - unit
    operators[lt] <- "<="
  }
   
  ii <- order(operators,decreasing=TRUE)
    
  list(A = A[ii,,drop=FALSE], b=b[ii], neq=sum(operators=="=="), nleq = sum(operators=="<="),order=ii)
}

