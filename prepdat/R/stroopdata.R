#' Reaction-times and accuracy for color naming in a Stroop task (e.g., Stroop,
#' 1935).
#'
#' A dataset containing reaction-times, accuracy, and other attributes of 5400
#' experimental trials.
#' 
#' @usage data(stroopdata)
#'
#' @format A data frame with 5401 rows and 10 columns:
#' \describe{
#'   \item{subject}{Case identifier, in numerals}
#'   \item{block}{Percent of congruent target_type trials in a block.
#'    1 means 80 percent congruent, 2 means 20 percent congruent}
#'   \item{age}{Age of subject, in integers}
#'   \item{gender}{Gender of subject, in integers. 1 means male, 2 means
#'   female}
#'   \item{order}{Order of blocks, in integers.
#'    1 means subject did 80 percent congruent block first and 20 percent
#'    congruent block second.
#'    2 means subject did 20 percent congruent block first and 80 percent
#'    congruent block second.}
#'   \item{font_size}{Font size of the stimulus, in integers}
#'   \item{trial_num}{Trial number, in integers}
#'   \item{target_type}{Type of stimulus for a given trial. 1 means
#'   congruent stimulus, 2 means incongruent stimulus}
#'   \item{rt}{Reaction time, in milliseconds}
#'   \item{ac}{Accuracy, 1 means correct, 0 means incorrect}
#' }
#' 
#' @references Stroop, J. R. (1935). Studies of interference in serial verbal
#' reactions. \emph{Journal of experimental psychology, 18}(6), 643.
#' 
#' @examples data(stroopdata)
#' head(stroopdata)
"stroopdata"