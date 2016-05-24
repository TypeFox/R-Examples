#' @title Single census mark-recapture data with lengths for Brown Trout from Valley Creek, MN.
#' 
#' @description Single censuse mark-recapture data for Brown Trout (\emph{Salmo trutta}) from Valley Creek, MN captured in April, 1988.  Length of trout was recorded so that abundance estimated can be made by length categories.
#' 
#' @name BrownTroutVC1
#' 
#' @docType data
#' 
#' @format A data frame with 1014 observations on the following 3 variables.
#' \describe{
#'   \item{len}{A numeric vector of total length measurements (cm)}
#'   \item{sample}{A factor variable representing the sample in which the fish was captured.  The marking run is labelled with \code{first} and the recapture run is labelled with \code{second}}
#'   \item{recap}{A factor variable representing whether the fish was a \dQuote{recap}ture  in the second sample (\code{YES}) or not (\code{NO})}
#' }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Population Size
#'    \item Abundance
#'    \item Mark-Recapture
#'    \item Capture-Recapture
#'    \item Petersen
#'  }
#' 
#' @concept Abundance 'Population Size' 'Mark-Recapture' 'Capture-Recapture' 'Petersen'
#' 
#' @source Obtained directly from Tom Kwak, North Carolina Cooperate Unit at North Carolina State University and part of the data published in Kwak, T.J. and T.F. Waters.  1997.  Trout production dynamics and water quality in Minnesota streams.  Transactions of the American Fisheries Society, 126:35-48.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(BrownTroutVC1)
#' str(BrownTroutVC1)
#' head(BrownTroutVC1)
#' hist(BrownTroutVC1$len,main="")
#' 
NULL
