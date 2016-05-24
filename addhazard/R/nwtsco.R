#' Dataset from the National Wilms Tumor Study (NWTS) 
#' @format A data frame with 3915 rows and 12 variables:
#' \describe{
#' \item{trel}{Time to relapse orlast date seen (yr), continuous}
#' \item{tsur}{Time to death or last date seen (yr), continuous}
#' \item{relaps}{Indicator of relapse,
#'                0 = Alive no prior relapse when last seen,
#'                1 = Relapsed after trel years}
#' \item{dead}{Indicator of death, 
#'               0 = Alive when last seen,
#'               1 = Died after tsur years}
#' \item{study}{NWTS study, 
#'              3 = NWTS-3,
#'              4 = NWTS-4}
#'\item{stage}{Stage of disease, 
#'             1=I,
#'             2=II,
#'             3=III,
#'             4=IV}
#' \item{histol}{Central Path histology, 
#'              0 = Favorable (FH),
#'              1 = Unfavorable (UH)}
#' \item{instit}{Institutional histology,
#'               0 = Favorable (FH),
#'              1 = Unfavorable (UH)}
#' \item{age}{Age at diagnosis (yr), continuous}
#' \item{yr}{Year of diagnosis, calendar year}
#' \item{specwgt}{Weight of tumor bearing specimen, in grams (continuous)}
#' \item{tumdiam}{Diameter of tumor, in centimeters (continuous)}
#' }
#' @source \url{http://faculty.washington.edu/norm/datasets/nwts-share.txt}
#'  Originally used by M. Kulich and D.Y. Lin:
#'  Improving the efficiency of relative-risk estimation in case-cohort studies. J Amer Statis 
#'  Assoc 99:832-844, 2004.
#' @docType data
#' @keywords datasets
#' @name nwtsco
"nwtsco"

