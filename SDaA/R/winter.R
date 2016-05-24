#' ASU Winter Closure Survey
#' 
#' Selected variables from the Arizona State University Winter Closure
#' Survey, taken in January 1995. This survey was taken to investigate
#' the attitudes and opinions of university employees toward the closing
#' of the university between December 25 and January 1.
#'
#' @name winter
#' @docType data
#' @format data frame with the following 6 variables:
#' \describe{
#'   \item{class}{stratum number; factor with levels \code{faculty}, 
#'     \code{classstaff} (classified staff), \code{admstaff} (administrative 
#'     staff) and \code{acprof} (academic professional)}
#'   \item{yearasu}{factor with levels \code{1} (1-2 years), \code{2} (3-4 years),
#'     \code{3} (5-9 years), \code{4} (10-14 years) and \code{5} (15 or more years)}
#'   \item{vacation}{In the past, have you \emph{usually} taken vacation days in 
#'     the entire period between December 25 and January 1? factor with levels
#'     \code{no} and \code{yes}}
#'   \item{work}{Did you work on campus during Winter Break Closure? factor with
#'     levels \code{no} and \code{yes}}
#'   \item{havediff}{Did the Winter Break Closure cause you any difficulty/concerns?
#'     factor with levels \code{no} and \code{yes}}
#'   \item{negaeffe}{Did the Winter Break Closure \emph{negatively} affect your work
#'     productivity? factor with levels \code{no} and \code{yes}}
#'   \item{ownsupp}{I was unable to obtain staff support in my department/office. 
#'     factor with levels \code{yes} and \code{no}}
#'   \item{othersup}{I was unable to obtain staff support in other departments/offices.
#'     factor with levels \code{yes} and \code{no}}
#'   \item{utility}{I was unable to access computers, copy machine, etc. in my 
#'     department/office. factor with levels \code{yes} and \code{no}}
#'   \item{environ}{I was unable to endure environmental conditions - e.g., not
#'     properly climatized. factor with levels \code{yes} and \code{no}}
#'   \item{uniserve}{I was unable to access university services necessary to my 
#'     work; factor with levels \code{yes} and \code{no}}
#'   \item{workelse}{I was unable to work on my assignments because I work in 
#'     another department/office; factor with levels \code{yes} and \code{no}}
#'   \item{offclose}{I was unable to work on my assignments because my office 
#'     was closed; factor with levels \code{yes} and \code{no}}
#'   \item{treatsta}{compared to other departments/offices, I feel staff in my 
#'      department/office were treated fairly; factor with levels \code{strongagr}
#'      (strongly agree), \code{agree}, \code{undecided}, \code{disagree}, 
#'      \code{strdisagr} (strongly disagree)}
#'   \item{treatme}{compared to other people working in my department/office, I
#'      feel I was treated fairly; factor with levels \code{strongagr}
#'      (strongly agree), \code{agree}, \code{undecided}, \code{disagree}, 
#'      \code{strdisagr} (strongly disagree)}
#'   \item{process}{How satisfied are you with the process used to inform staff 
#'      about Winter Closure? factor with levels \code{verysat} (very satisfied),
#'      \code{satisfied}, \code{undecided}, \code{dissatisfied} and \code{verydissat}
#'      (very dissatisfied)}
#'   \item{satbreak}{How satisfied are you with the fact that ASU had a Winter Break
#'      Closure this year? factor with levels \code{verysat} (very satisfied),
#'      \code{satisfied}, \code{undecided}, \code{dissatisfied} and \code{verydissat}
#'      (very dissatisfied)}
#'   \item{breakaga}{Would you want to have Winter Break Closure again? factor with
#'      levels \code{no} and \code{yes}}
#' }
#' @source courtesy of the ASU Office of University Evaluation.
#' @references Lohr (1999). Sampling: Design and Analysis, Duxbury, p. TODO and
#' 447--448.
#' @export
NULL
