#' Check wether a matrix is totally unimodular.
#'
#'
#' A matrix for which the determinant of every square submatrix equals \eqn{-1},
#' \eqn{0} or \eqn{1} is called
#' \href{https://en.wikipedia.org/wiki/Unimodular_matrix}{totally unimodular}. 
#' This function tests wether a matrix with coefficients in \eqn{\{-1,0,1\}} is
#' totally unimodular. It tries to reduce the matrix using the reduction method
#' described in Scholtus (2008). Next, a test based on Heller and Tompkins
#' (1956) or Raghavachari is performed.
#'
#' @title Test for total unimodularity of a matrix.
#'
#' @param A An object of class \code{\link{matrix}}.
#' @return logical
#' 
#' @example ../examples/unimodular.R
#' @references 
#' Heller I and Tompkins CB (1956). An extension of a theorem of Danttzig's In kuhn HW and Tucker AW (eds.),
#' pp. 247-254. Princeton University Press.
#'
#' Raghavachari M (1976). A constructive method to recognize the total
#' unimodularity of a matrix. _Zeitschrift fur operations research_,
#' *20*, pp. 59-61.
#'
#' Scholtus S (2008). Algorithms for correcting some obvious
#' inconsistencies and rounding errors in business survey data. Technical
#' Report 08015, Netherlands.
#'
#' @export
is_totally_unimodular <- function(A) {
    
    # A matrix with elements not in {-1,0,1} cannot be totally unimodular.
    if ( !all(A %in% c(-1,0,1)) ){
        value <- FALSE
    }

    # After reduction, A has no rows or columns containing less than 2 elements.
    A <- reduceMatrix(A)
    if ( length(A) == 0 ){ # if reduced to nothingness, we are ready.
        value <- TRUE
    # HT-criterium, by rows or columns
    } else if (max(colSums(abs(A))) == 2){ 
        value <- hellerTompkins(A)
    } else if (max(rowSums(abs(A))) == 2){
        value <- hellerTompkins(t(A))
    # raghavachari criterium recurses over columns. Minimize effort by
    # transposition when possible.    
    } else {                        
        if ( nrow(A) >= ncol(A) ){ 
            value <- raghavachari(A)
        } else {
            value <- raghavachari(t(A))
        }
    }
    return(value)
}


#' Apply reduction method from Scholtus (2008)
#'
#' Apply the reduction method in the appendix of Scholtus (2008) to a matrix.
#' Let \eqn{A} with coefficients in \eqn{\{-1,0,1\}}. If, after a possible 
#' permutation of columns it can be written 
#' in the form \eqn{A=[B,C]} where each column in \eqn{B} has at most 1 nonzero
#' element, then \eqn{A} is totally unimodular if and only if \eqn{C} is totally
#' unimodular. By transposition, a similar theorem holds for the rows of A. This
#' function iteratively removes rows and columns with only 1 nonzero element
#' from \eqn{A} and returns the reduced result.
#'
#' @param A An object of class matrix in \eqn{\{-1,0,1\}^{m\times n}}. 
#' @return The reduction of A.
#' @seealso \code{\link{is_totally_unimodular}}
#' 
#' @references
#'
#' Scholtus S (2008). Algorithms for correcting some obvious
#' inconsistencies and rounding errors in business survey data. Technical
#' Report 08015, Netherlands.
#'
#' @keywords internal
reduceMatrix <- function(A){
    d1 <- c(0,0)
    d <- dim(A)
    while ( !all(d1==d) ){
        A <- A[, colSums(abs(A)) >= 2,drop=FALSE]
        A <- A[rowSums(abs(A)) >= 2, ,drop=FALSE]
        d1 <- d
        d <- dim(A)
    }
    return(A)
}

#' Determine if a matrix is totally unimodular using Heller and Tompkins criterium.
#'
#' This function is \code{deducorrect} internal
#' 
#'
#' @param A An object of class matrix in \eqn{\{-1,0,1\}^{m\times n}}. 
#'  Each column  must have exactly 2 nonzero elements. (This is tested by 
#'      \code{\link{is_totally_unimodular}}).
#'
#' @return \code{TRUE} if matrix is unimodular, otherwise \code{FALSE}
#' @seealso \code{\link{is_totally_unimodular}}
#' @keywords internal
hellerTompkins <- function(A){
    # If the matrix has columns with two elements, and those elements differ in
    # sign for all those columns, the matrix is unimodular. 
    if ( !any(abs(colSums(A))==2) ){ 
        return(TRUE)
    }
    # Loop over ways to split a matrix in 2 by rows.
    # Return TRUE when HT criterium is met, FALSE otherwise.
    for ( m in 1:(nrow(A)%/%2) ){
        I <- combn(1:nrow(A), m)
        for ( i in 1:ncol(I)){
            M1 <- A[I[ , i], ,drop=FALSE]
            M2 <- A[-I[ , i], ,drop=FALSE]
            if ( !any(abs(colSums(M1)) == 2) && !any(abs(colSums(M2)) == 2) ){
                return(TRUE)
            }
        }
    }
    return(FALSE)
}

#' Test if a list of matrices are all unimodular
#'
#' Helper function for \code{\link{raghavachari}}
#' @param L A list of objects of class matrix.
#' @return logical vector of length \code{length(L)}
#' @seealso \code{\link{is_totally_unimodular}}
#' @keywords internal
allTotallyUnimodular <- function(L){
    for ( i in 1:length(L) ){
        if ( !is_totally_unimodular(L[[i]]) ){
            return(FALSE)
        }
    }
    return(TRUE)
}

#' Determine if a matrix is unimodular using recursive Raghavachari criterium
#' 
#' This function is \code{lintools} internal
#'
#' @param A An object of class Matrix in \eqn{\{-1,0,1\}^{m\times n}}. 
#' @return \code{TRUE} or \code{FALSE}
#' @seealso \code{\link{is_totally_unimodular}}
#' @keywords internal
raghavachari <- function(A){
    J <- colSums(abs(A))>=2
    j <- which.max(colSums(abs(A[ ,J,drop=FALSE])))
    j <- which(J)[j]
    a_j <- A[ , j, drop=FALSE]
    i_j <- which(a_j != 0)
    L <- lapply(seq_along(i_j), function(i, A){
            irow <- A[i_j[i], , drop=FALSE] 
            A[i_j[-i],] <- A[i_j[-i], , drop=FALSE] - (A[i_j[-i], j, drop=FALSE]*irow[j]) %*% irow
            A[, -j, drop=FALSE]
        }, A)
    return(allTotallyUnimodular(L))
}


