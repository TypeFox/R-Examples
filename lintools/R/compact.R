#' Remove spurious variables and restrictions
#' 
#' A system of linear (in)equations can be compactified by removing
#' zero-rows and zero-columns (=variables). Such rows and columns may
#' arise after substitution (see \code{\link{subst_value}}) or eliminaton
#' of a variable (see \code{\link{eliminate}}).
#' 
#' @section Details:
#' It is assumend that the system of equations is in normalized form (see \code{link{normalize}}).
#'  
#'    
#' 
#' @param A [\code{numeric}] matrix
#' @param b [\code{numeric}] vector
#' @param x [\code{numeric}] vector
#' @param neq [\code{numeric}] The first \code{neq} rows in \code{A} and
#'   \code{b} are treated as linear equalities. 
#' @param nleq [\code{numeric}] The \code{nleq} rows after \code{neq} are treated as
#'   inequations of the form \code{a.x<=b}. All remaining rows are treated as strict inequations
#'   of the form \code{a.x<b}.
#' @param eps [\code{numeric}] Anything with absolute value \code{< eps} is considered zero.
#' @param remove_columns [\code{logical}] Toggle remove spurious columns from \code{A} and variables from \code{x}
#' @param remove_rows [\code{logical}] Toggle remove spurious rows
#' @param implied_equations [\code{logical}] replace cases of \code{a.x<=b} and \code{a.x>=b} with
#'        \code{a.x==b}.
#'
#' 
#' 
#' @return A \code{list} with the following elements.
#' \itemize{
#'   \item{\code{A}: The compactified version of input \code{A}}
#'   \item{\code{b}: The compactified version of input \code{b}}
#'   \item{\code{x}: The compactified version of input \code{x}}
#'   \item{\code{neq}: number of equations in new system}
#'   \item{\code{nleq}: number of inequations of the form \code{a.x<=b} in the new system}
#'   \item{\code{cols_removed}: [\code{logical}] indicates what elements of \code{x} (columns of \code{A}) have been removed}
#' }
#'  
#'   
#'  
#' @export
compact <- function(A, b, x=NULL, neq=nrow(A), nleq=0, eps=1e-8
    , remove_columns=TRUE, remove_rows=TRUE, implied_equations=TRUE){
  check_sys(A=A,b=b,neq=neq,eps=eps) 
 
  if ( nrow(A)==0 | ncol(A)== 0){
    return(list(
      A = matrix()[0,0,drop=FALSE]
      , b=numeric(0)
      , neq=0
      , nleq=0
      , cols_removed=!logical(ncol(A))
    ))
  }
  
  Ai <- abs(A) > eps
  
  ops <- rep("<",nrow(A))
  ops[seq_len(neq)] <- "=="
  ops[neq + seq_len(nleq)] <- "<="
  
  cols_removed <- logical(ncol(A))
  if (remove_columns){
    cols_removed <- colSums(Ai) == 0
    A <- A[,!cols_removed,drop=FALSE]
    if( !is.null(x) ) x <- x[!cols_removed]
  }
  
  if ( remove_rows ){
    I <- rowSums(Ai) == 0
    rows_removed <- (ops != "<" & I & abs(b) < eps ) | (ops == "<" & I & b < -eps)
    A <- A[!rows_removed,,drop=FALSE]
    b <- b[!rows_removed]
    neq <- neq - sum(ops == "==" & rows_removed)
    nleq <- nleq - sum(ops=="<=" & rows_removed)
  }
  
  ineqs <- neq + seq_len(nleq)
  if (implied_equations && length(ineqs) > 0 ){
    Ab <- cbind(A[ineqs,,drop=FALSE])
    
    # 
    bi <- abs(b[ineqs])
    bi[bi < eps] <- 1  # avoid dividing by zero
    Ab <- Ab/bi        # normalize coefficients
    # compare all rows (avoid redundancies)
    I <- rep(seq_len(nleq),times=nrow(A))
    J <- rep(seq_len(nleq),each=nrow(A))
    ii <- I < J
    I <- I[ii]
    J <- J[ii]
    # compute implied equations
    ieqn <- abs( rowSums(Ab[I,,drop=FALSE] + Ab[J,,drop=FALSE]) ) < eps 
    if ( any(ieqn) ){
      keep <-  I[ieqn]
      throw <- J[ieqn]
      ##
      A <- rbind(
        A[seq_len(neq),,drop=FALSE]            # original equalities
        , A[ineqs[keep],,drop=FALSE]           # combined equalities
        , A[ineqs[-c(throw,keep)],,drop=FALSE] # remaining inequalities
        , A[-seq_len(neq+nleq),,drop=FALSE]    # strict inequalities
        )
      b <- c(
        b[seq_len(neq)]            # original equalities
        , b[ineqs[keep]]           # combined equalities
        , b[ineqs[-c(throw,keep)]] # remaining inequalities
        , b[-seq_len(neq+nleq)]    # strict inequalities
      )
      neq <- neq + length(ieqn)
      nleq <- length(ineqs) - length(throw) - length(keep)
    }
  } 
  
  list(A=A, b=b, x=x, neq=neq, nleq=nleq, cols_removed=cols_removed)
}

