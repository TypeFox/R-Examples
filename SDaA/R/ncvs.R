#' Victimization Incidents in the July-December 1989 NCVS
#' 
#' Selected variables for victimization incidents in the 
#' July-December 1989 NCVS. Note that some variables were 
#' recoded from the original data file.
#' @name ncvs
#' @docType data
#' @format Data frame with the following seven variables: 
#' \describe{
#'   \item{wt}{incident weight}
#'   \item{sex}{factor with levels \code{male} and \code{female}}
#'   \item{violent}{violent crime? factor with levels \code{no} and \code{yes}}
#'   \item{injury}{did the victim have injuries? factor with levels \code{no}
#'     and \code{yes}}
#'   \item{medcare}{factor with levels \code{yes} if the victim received medical 
#'     care and \code{no} otherwise}
#'   \item{reppol}{was the incident reported to the police? factor with levels
#'     \code{yes} and \code{no}}
#'   \item{numoff}{number of offenders involved in crime; factor with levels
#'     \code{one}, \code{more} (more than one) and \code{dontknow}}
#' }
#' @source Incident-level concatenated file, NCS8864I, in NCJ-130915, U.S.
#'   Department of Justice 1991.
#' @references Lohr (1999). Sampling: Design and Analysis, Duxbury, p. TODO and
#'   443.
#' @export
NULL
