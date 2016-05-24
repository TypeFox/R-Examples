#' In-migration Coefficient of Variation
#'
#' As "the coefficient of variation is defined as the standard deviation to mean ratio of a distribution", the In-migration Coefficient of Variation is computed by dividing the standard deviation (with the nominator being \eqn{n} instead of \eqn{n-1}) of the in-migration flows by the mean.
#' @param m migration matrix
#' @return A numeric vector of standardized values where a higher (\eqn{\neq 0}) shows more spatial focus.
#' @references \itemize{
#'   \item Andrei Rogers and Stuart Sweeney (1998) Measuring the Spatial Focus of Migration Patterns. \emph{The Professional Geographer} \bold{50}, 232--242
#' }
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.cv.in(migration.hyp)    # 0.2000000 0.5000000 0.3333333
#' migration.cv.in(migration.hyp2)   # 0.2000000 0.0000000 0.4285714
#' }
#' @export
#' @seealso \code{\link{migration.cv.out}} \code{\link{migration.acv.in}} \code{\link{migration.acv.out}} \code{\link{migration.acv}}
migration.cv.in <- function(m) {
    diag(m) <- NA
    n       <- ncol(m) - 1
    apply(m, 2, function(m.row) (sd(m.row, na.rm = TRUE) * sqrt((n - 1) / n)) / mean(m.row, na.rm = TRUE))
}


#' Out-migration Coefficient of Variation
#'
#' As "the coefficient of variation is defined as the standard deviation to mean ratio of a distribution", the Out-migration Coefficient of Variation is computed by dividing the standard deviation (with the nominator being \eqn{n} instead of \eqn{n-1}) of the out-migration flows by the mean.
#' @param m migration matrix
#' @return A numeric vector of standardized values where a higher (\eqn{\neq 0}) shows more spatial focus.
#' @references \itemize{
#'   \item Andrei Rogers and Stuart Sweeney (1998) Measuring the Spatial Focus of Migration Patterns. \emph{The Professional Geographer} \bold{50}, 232--242
#' }
#' @examples \dontrun{
#' data(migration.hyp)
#' migration.cv.out(migration.hyp)    # 0 0 0
#' migration.cv.out(migration.hyp2)   # 0.00 0.25 0.00
#' }
#' @export
#' @seealso \code{\link{migration.cv.in}} \code{\link{migration.acv.in}} \code{\link{migration.acv.out}} \code{\link{migration.acv}}
migration.cv.out <- function(m) {
    diag(m) <- NA
    n       <- ncol(m) - 1
    apply(m, 1, function(m.row) (sd(m.row, na.rm = TRUE) * sqrt((n - 1) / n)) / mean(m.row, na.rm = TRUE))
}


#' Aggregated In-migration Coefficient of Variation
#'
#' The Aggregated In-migration Coefficient of Variation is the weighted average of the In-migration Coefficient of Variation (\code{\link{migration.cv.in}}).
#' @param m migration matrix
#' @return A number where a higher (\eqn{\neq 0}) shows more spatial focus.
#' @references \itemize{
#'   \item Andrei Rogers and Stuart Sweeney (1998) Measuring the Spatial Focus of Migration Patterns. \emph{The Professional Geographer} \bold{50}, 232--242
#' }
#' @examples
#' data(migration.hyp)
#' migration.acv.in(migration.hyp)    # 0.3333333
#' migration.acv.in(migration.hyp2)   # 0.25
#' @export
#' @seealso \code{\link{migration.cv.in}} \code{\link{migration.cv.out}} \code{\link{migration.acv.out}} \code{\link{migration.acv}}
migration.acv.in <- function(m) {
    diag(m) <- NA
    n       <- ncol(m) - 1
    sum(colSums(m, na.rm = TRUE) / sum(m, na.rm = TRUE) * migration.cv.in(m))
}


#' Aggregated Out-migration Coefficient of Variation
#'
#' The Aggregated Out-migration Coefficient of Variation is the weighted average of the Out-migration Coefficient of Variation (\code{\link{migration.cv.out}}).
#' @param m migration matrix
#' @return A number where a higher (\eqn{\neq 0}) shows more spatial focus.
#' @references \itemize{
#'   \item Andrei Rogers and Stuart Sweeney (1998) Measuring the Spatial Focus of Migration Patterns. \emph{The Professional Geographer} \bold{50}, 232--242
#' }
#' @examples
#' data(migration.hyp)
#' migration.acv.out(migration.hyp)    # 0
#' migration.acv.out(migration.hyp2)   # 0.125
#' @export
#' @seealso \code{\link{migration.cv.in}} \code{\link{migration.cv.out}} \code{\link{migration.acv.in}} \code{\link{migration.acv}}
migration.acv.out <- function(m) {
    diag(m) <- NA
    n       <- ncol(m) - 1
    sum(rowSums(m, na.rm = TRUE) / sum(m, na.rm = TRUE) * migration.cv.out(m))
}


#' Aggregated System-wide Coefficient of Variation
#'
#' The Aggregated System-wide Coefficient of Variation is simply the sum of the Aggregated In-migration (\code{\link{migration.acv.in}}) and the Aggregated Out-migration Coefficient of Variation (\code{\link{migration.acv.out}}).
#' @param m migration matrix
#' @return A number where a higher (\eqn{\neq 0}) shows more spatial focus.
#' @references \itemize{
#'   \item Andrei Rogers and Stuart Sweeney (1998) Measuring the Spatial Focus of Migration Patterns. \emph{The Professional Geographer} \bold{50}, 232--242
#' }
#' @examples
#' data(migration.hyp)
#' migration.acv(migration.hyp)    # 0.3333333
#' migration.acv(migration.hyp2)   # 0.375
#' @export
#' @seealso \code{\link{migration.cv.in}} \code{\link{migration.cv.out}} \code{\link{migration.acv.in}} \code{\link{migration.acv.out}}
migration.acv <- function(m)
    migration.acv.in(m) + migration.acv.out(m)
