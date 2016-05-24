#' Cohort study on breast cancer patients from the Netherlands.
#'
#' This dataset is borrowed from
#' "Flexible parametric survival analysis using Stata: beyond the Cox model" (Roystone and Lambert, 2011).
#'  It contains follow-up data on 2982 woman with breast cancer who have gone through breast surgery.
#'  The women are followed from the time of surgery until death, relapse or censoring.
#'
#' @docType data
#' @name rott2
#' @usage data(rott2)
#' @format The dataset \code{rott2} contains the following variables:
#' \describe{
#'   \item{pid}{patient ID number.}
#'   \item{year}{year of breast surgery (i.e. year of enrollment into the study), between the years 1978-1993.}
#'   \item{rf}{relapse free interval measured in months.}
#'   \item{rfi}{relapse indicator.}
#'   \item{m}{metastasis free.}
#'   \item{mfi}{metastasis status.}
#'   \item{os}{overall survival}
#'   \item{osi}{overall survival indicator}
#'   \item{age}{age at surgery measured in years.}
#'   \item{meno}{menopausal status with levels "\code{pre}" and "\code{post}".}
#'   \item{size}{tumor size in three classes: \code{<=20mm, >20-50mmm} and \code{>50mm}.}
#'   \item{grade}{differentiation grade with levels 2 or 3.}
#'   \item{pr}{progesterone receptors, fmol/l.}
#'   \item{er}{oestrogen receptors, fmol/l.}
#'   \item{nodes}{the number of positive lymph nodes.}
#'   \item{hormon}{hormonal therapy with levels "\code{no}" and "\code{yes}".}
#'   \item{chemo}{categorical variable indicating whether the patient recieved chemotheraphy or not, with levels "\code{no}" and "\code{yes}".}
#'   \item{recent}{a numeric indicator of whether the tumor was discovered recently with levels "\code{1978-87}" and "\code{1988-93}".}
#'   \item{no.chemo}{a numerical indicator of whether the patient did not recieved chemotherapy. Recoded version of "\code{chemo}" where "\code{yes}" is recoded as 0 and "\code{no}" is recoded as 1.}
#' }
#'
#' \strong{The following changes have been made to the original data in Roystone and Lambert (2011):}
#'
#' - The variable "\code{chemo}" is recoded into the numeric indicator variable "\code{no.chemo}":
#'
#' \code{rott22$no.chemo <- as.numeric(rott2$chemo == "no")}
#'
#' The follwing variables have been removed from the original dataset: \code{enodes, pr_1, enodes_1, _st, _d, _t, _t0}
#' since they are recodings of some existing variables which are not used in this analysis.
#' @references Royston, Patrick & Lambert, Paul. C (2011). \emph{Flexible parametric survival analysis using Stata: beyond the Cox model}. College Station, Texas, U.S, Stata press.
#' @references \url{http://www.stata-press.com/data/fpsaus.html}
NULL
