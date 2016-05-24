#' Stress-Anxiety data
#'
#' A second example is a data from a study that investigates the relationship 
#' between stress and anxiety..
#'
#' @format A data frame with 166 rows and 2 variables:
#' \describe{
#'   \item{Anxiety}{Scores on Anxiety subscale}
#'   \item{Stress}{Scores on Stress subscale}
#' }
#' @source \url{https://dl.dropboxusercontent.com/u/1857674/betareg/betareg.html}
"AnxStrData"


#' Juror data
#'
#' Juror Judgment Study.
#'
#' @format A data frame with 104 rows and 3 variables:
#' \describe{
#'   \item{crc99}{The ratings of confidence levels with rescaling into the (0, 1) interval to avoide 1 and 0 values.}
#'   \item{vert}{ was the dummy variable for coding the conditions of verdict types, whereas }
#'   \item{confl}{ was the dummy variable for coding the conflict conditions}
#' }
#' @source \url{https://dl.dropboxusercontent.com/u/1857674/betareg/betareg.html}
"JurorData"


#' IPCC data-set
#'
#' The IPCC data-set comprises the lower, best, and upper estimates 
#' for the phrases "likely" and "unlikely" in six IPCC report sentences.
#'
#' @format A data frame with 4014 rows and 8 variables:
#' \describe{
#'   \item{subj}{Subject ID number}
#'   \item{treat}{Experimental conditions}
#'   \item{valence}{Valence of the sentences}
#'   \item{prob}{raw probability estimates} 
#'   \item{probm}{Linear transformed prob into (0, 1) interval} 
#'   \item{mid}{Distinguish lower, best and upper estiamtes }
#'   \item{high}{Distinguish lower, best and upper estiamtes } 
#'   \item{Question}{IPCC question number} 
#' }
"IPCC"