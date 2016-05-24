#' Egg Size from Coots 
#'
#' Selected information on egg size from coots, from a study
#' by Arnold (1991). Data courtesy of Todd Arnold.
#' @name coots
#' @docType data
#' @format Data frame with the following 11 variables: 
#' \describe{
#'   \item{clutch}{clutch number from which eggs were subsampled}
#'   \item{csize}{number of eggs in clutch (Mi)}
#'   \item{length}{length of egg (mm)}
#'   \item{breadth}{maximum breadth of egg (mm)}
#'   \item{volume}{calculated as 0.00507 x length x breadth^2}
#'   \item{tmt}{received supplemental feeding? factor with levels
#'     \code{no} and \code{yes}}
#' }   
#' @note Not all observations are used for this data set, so
#' results may not agree with those in Arnold (1991)
#' @source Arnold, T.W. (1991). Intraclutch variation in egg size
#' of American Coots, \emph{The Condor}, 93: 19--27
#' @references Lohr (1999). Sampling: Design and Analysis, Duxbury, p. TODO and
#'   440. 
#' @export
NULL
