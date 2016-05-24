#' Data from the 1994 Survey of ASA Membership on Certification 
#'
#' @name certify
#' @docType data
#' @format Data frame with the following 11 variables: 
#' \describe{
#'   \item{certify}{should the ASA develop some form of certification? factor
#'     with levels \code{yes}, \code{possibly}, \code{noopinion}, 
#'     \code{unlikely} and \code{no}}
#'   \item{approve}{would you approve of a certification program similar to 
#'     that described in the July 1993 issue of \emph{Amstat News}? factor
#'     with levels \code{yes}, \code{possibly}, \code{noopinion}, 
#'     \code{unlikely} and \code{no}}
#'   \item{speccert}{Should there be specific certification programs for
#'     statistics subdisciplines? factor with levels \code{yes}, 
#'     \code{possibly}, \code{noopinion}, \code{unlikely} and \code{no}}
#'   \item{wouldyou}{If the ASA developed a certification program, would you
#'     attempt to become certified? factor with levels \code{yes}, 
#'     \code{possibly}, \code{noopinion}, \code{unlikely} and \code{no}}
#'   \item{recert}{If the ASA offered certification, should recertification
#'     be required every several years? factor with levels \code{yes}, 
#'     \code{possibly}, \code{noopinion}, \code{unlikely} and \code{no}}
#'   \item{subdisc}{Major subdiscipline; factor with levels \code{BA} (Bayesian),
#'     \code{BE} (business and economic), \code{BI} (biometrics), \code{BP}
#'     (biopharmaceutical), \code{CM} (computing), \code{EN} (environment), 
#'     \code{EP} (epidemiology), \code{GV} (government), \code{MR} (marketing),
#'     \code{PE} (physical and engineering), \code{QP} (quality and productivity),
#'     \code{SE} (statistical education), \code{SG} (statistical graphics), 
#'     \code{SP} (sports), \code{SR} (survey research), \code{SS} (social statistics),
#'     \code{TH} (teaching  statistics in health sciences), \code{O} (other)}
#'   \item{college}{Highest collegiate degree; factor with levels \code{B} (BS or BA),
#'     \code{M} (MS), \code{N} (none), \code{P} (PhD) and \code{O} (other)}
#'   \item{employ}{Employment status; factor with levels \code{E} (employed), 
#'     \code{I} (in school), \code{R} (retired), \code{S} (self-employed), 
#'     \code{U} (unemployed) and \code{O} (other)}
#'   \item{workenv}{Primary work environment; factor with levels \code{A} (academia),
#'     \code{G} (government), \code{I} (industry), \code{O} (other)
#'   \item{workact}{Primary work activity; factor with levels \code{C} (consultant), 
#'     \code{E} (educator), \code{P} (practitioner), \code{R} (researcher), 
#'     \code{S} (student) and \code{O} (other)}
#'   \item{yearsmem}{For how many years have you been a member of ASA?}
#' }   
#' @note The full dataset is on Statlib
#' @references Lohr (1999). Sampling: Design and Analysis, Duxbury, p. TODO and
#'   439. 
#' \url{http://lib.stat.cmu.edu/asacert/certsurvey}
#' @export
NULL
