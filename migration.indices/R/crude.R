#' Crude Migration Rate
#'
#'
#' @param m migration matrix
#' @param PAR population at risk (estimated average population size)
#' @param k scaling constant (set to \code{100} by default to result in percentage)
#' @return percentage (when \code{k=100})
#' @references \itemize{
#'   \item Philip Rees, Martin Bell, Oliver Duke-Williams and Marcus Blake (2000) Problems and Solutions in the Measurement of Migration Intensities: Australia and Britain Compared. \emph{Population Studies} \bold{54}, 207--222
#' }
#' @examples
#' data(migration.world)
#' migration.cmr(migration.world, 6e+9)
#' @export
migration.cmr <- function(m, PAR, k = 100) {

    check.migration.matrix(m)
    if (missing(PAR))
        stop('Estimated population size was not provided!')

    k*sum(m)/PAR

}
