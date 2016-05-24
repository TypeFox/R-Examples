#' Rao-Stirling diversity index based on proportions of cited disciplines.
#'
#' This function calculates the Rao-Stirling diversity index of a single publication, based on its proportions of citations to different disciplines.
#' Unlike the function \code{\link{RaoProportions}}, the arguments of this function are not validated.
#'
#' @param p A vector of proportions of citations to different disciplines of a single publication.
#' @param sim A positive semi-definite matrix that encodes the similarity between disciplines, as explained in Porter and Rafols (2009).
#' The dimensions of this matrix are \emph{n} x \emph{n}, being \emph{n} the total number of disciplines.
#' The number of rows and the number columns of this matrix need to be equal to the length of \code{counts}.
#' The self-similarities (i.e. the diagonal elements) have to be 1.
#' @return The Rao-Stirling diversity index of a publication.
#' @references
#' Porter, A. and Rafols, I. (2009) Is science becoming more interdisciplinary? Measuring and mapping six research fields over time. Scientometrics, Vol. 81, No. 3 (719-745). DOI:10.1007/s11192-008-2197-2
#' @keywords internal
RaoProportionsNoValidation <- function(p, sim) {
  return(drop( 1 - (p %*% sim %*% p) ))
}

#' Rao-Stirling diversity index based on counts of cited disciplines.
#'
#' This function calculates the Rao-Stirling diversity index of a single publication, based on its count of citations to different disciplines.
#' Unlike the function \code{\link{RaoCounts}}, the arguments of this function are not validated.
#'
#' @param c A vector of counts of citations to different disciplines of a single publication.
#' @param sim A positive semi-definite matrix that encodes the similarity between disciplines, as explained in Porter and Rafols (2009).
#' The dimensions of this matrix are \emph{n} x \emph{n}, being \emph{n} the total number of disciplines.
#' The number of rows and the number columns of this matrix need to be equal to the length of \code{counts}.
#' The self-similarities (i.e. the diagonal elements) have to be 1.
#' @return The Rao-Stirling diversity index of a publication.
#' @references
#' Porter, A. and Rafols, I. (2009) Is science becoming more interdisciplinary? Measuring and mapping six research fields over time. Scientometrics, Vol. 81, No. 3 (719-745). DOI:10.1007/s11192-008-2197-2
#' @keywords internal
RaoCountsNoValidation <- function(c, sim) {
  return(RaoProportionsNoValidation(c / sum(c), sim))
}

#' Rao-Stirling diversity index based on proportions of cited disciplines.
#'
#' This function calculates the Rao-Stirling diversity index of a single publication, based on its proportions of citations to different disciplines.
#'
#' @param proportions A vector of proportions of citations to different disciplines of a single publication.
#' Since the elements of this vector are proportions, they need to be within the interval [0,1].
#' @param similarity A positive semi-definite matrix that encodes the similarity between disciplines, as explained in Porter and Rafols (2009).
#' The dimensions of this matrix are \emph{n} x \emph{n}, being \emph{n} the total number of disciplines.
#' The number of rows and the number columns of this matrix need to be equal to the length of \code{counts}.
#' The self-similarities (i.e. the diagonal elements) have to be 1.
#' @return The Rao-Stirling diversity index of a publication.
#' @references
#' Porter, A. and Rafols, I. (2009) Is science becoming more interdisciplinary? Measuring and mapping six research fields over time. Scientometrics, Vol. 81, No. 3 (719-745). DOI:10.1007/s11192-008-2197-2
#' @keywords internal
RaoProportions <- function(proportions, similarity) {

  n <- length(proportions)
  # Error handling.
  if (n <= 1 || n != nrow(similarity) || n != ncol(similarity)) {
    stop("Parameters 'proportions' and 'similarity' have incompatible sizes.")
  }
  if (any(is.nan(proportions)) || any(proportions < 0) || any(proportions > 1)) {
    stop("Elements of 'proportions' are out of range.")
  }
  if (all(proportions == 0)) {
    stop("Elements of 'proportions' cannot be all 0.")
  }

  return(RaoProportionsNoValidation(proportions, similarity))
}

#' Rao-Stirling diversity index based on the counts of cited disciplines.
#'
#' This function calculates the Rao-Stirling diversity index of a single publication, based on the count of citations to different disciplines.
#' Note that this function returns the correct result if a vector of proportions is supplied as input. Since a normalization by (approximately) 1.0 is performed, the result is identical to the result of the function \code{\link{RaoProportions}} up to numerical precision.
#'
#' @param counts A vector of counts of citations to different disciplines of a single publication.
#' The length of the vector should be equal to the total number of disciplines.
#' @param similarity A positive semi-definite matrix that encodes the similarity between disciplines, as explained in Porter and Rafols (2009).
#' The dimensions of this matrix are \emph{n} x \emph{n}, being \emph{n} the total number of disciplines.
#' The number of rows and the number columns of this matrix need to be equal to the length of \code{counts}.
#' The self-similarities (i.e. the diagonal elements) have to be 1.
#' @return The Rao-Stirling diversity index of a publication.
#' @references
#' Porter, A. and Rafols, I. (2009) Is science becoming more interdisciplinary? Measuring and mapping six research fields over time. Scientometrics, Vol. 81, No. 3 (719-745). DOI:10.1007/s11192-008-2197-2
#' @keywords internal
RaoCounts <- function(counts, similarity) {

  n <- length(counts)
  # Error handling.
  if (n <= 1 || n != nrow(similarity) || n != ncol(similarity)) {
    stop("Parameters 'counts' and 'similarity' have incompatible sizes.")
  }
  if (any(is.nan(counts)) || any(counts < 0)) {
    stop("Elements of 'counts' are out of range.")
  }
  if (all(counts == 0)) {
    stop("Elements of 'counts' cannot be all 0.")
  }
  return(RaoCountsNoValidation(counts, similarity))
}
