#' Local retention of a connectivity matrix
#' 
#' Local retention is defined as the diagonal elements of the connectivity matrix.
#' 
#' @param conn.mat A square connectivity matrix.
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @example tests/test.retentionStats.R
#' @export
localRetention <- function(conn.mat)
  diag(conn.mat)

#' Relative local retention of a connectivity matrix
#' 
#' Relative local retention is defined as the diagonal elements of the
#' connectivity matrix divided by the sum of the corresponding column of the
#' connectivity matrix.
#' 
#' @param conn.mat A square connectivity matrix.
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @example tests/test.retentionStats.R
#' @export
relativeLocalRetention <- function(conn.mat) {
  ss = apply(conn.mat,2,sum)
  ss[ ss == 0 ] = 1
  localRetention(conn.mat) / ss  
}

#' Self recruitment of a connectivity matrix
#' 
#' If egg production is uniform over sites, then self recruitment is defined as 
#' the diagonal elements of the connectivity matrix divided by the sum of the 
#' corresponding row of the connectivity matrix.  If not, then the elements of 
#' the dispersal matrix must be weighted by the number of eggs produced.
#' 
#' @param conn.mat A square connectivity matrix.
#' @param eggs A vector of egg production values for each site.  Defaults to
#'   \code{NULL}, equivalent to assuming all sites have equal egg production.
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @example tests/test.retentionStats.R
#' @export
selfRecruitment <- function(conn.mat,eggs=NULL) {
  if (!is.null(eggs))
    conn.mat = conn.mat %*% diag(eggs)
  
  ss = apply(conn.mat,1,sum)
  ss[ ss == 0 ] = 1
  localRetention(conn.mat) / ss  
}

