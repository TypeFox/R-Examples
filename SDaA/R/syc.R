#' Survey of Youth in Custody, 1987
#' 
#' The 1987 Survey of Youth in Custody sampled juveniles and young adults in 
#' long-term, state-operated juvenile institutions. Residents of facilities
#' at the end of 1987 were interviewed about family background, previous
#' criminal history, and drug and alcohol use. Selected variables from the
#' survey are contained in the syc data frame.
#' @name syc
#' @docType data
#' @format \describe{
#'   \item{stratum}{stratum number}
#'   \item{psu}{psu (facility) number}
#'   \item{psusize}{number of eligible residents in psu}
#'   \item{initwt}{initial weight}
#'   \item{finalwt}{final weight}
#'   \item{randgrp}{random group number}
#'   \item{age}{age of resident}
#'   \item{race}{race of resident: factor with levels \code{1} (white),
#'      \code{2} (black), \code{3} (Asian/Pacific Islander),
#'      \code{4} (American Indian, Aleut, Eskimo}, \code{5} (other)}
#'   \item{ethnicty}{ethnicity; factor with levels \code{hispanic} and
#'      \code{notHispanic}}
#'   \item{educ}{highest grade before sent to correctional institution; factor
#'      with levels \code{0} (never attended), \code{1}-\code{12} (highest grade
#'      attended), \code{13} (GED), \code{14} (other)}
#'   \item{sex}{factor with levels \code{male} and \code{female}}
#'   \item{livewith}{factor with levels \code{1} (mother only), \code{2} (father only), 
#'     \code{3} (both mother and father), \code{4} (grandparents), \code{5} (other relatives),
#'     \code{6} (friends), \code{7} (foster home), \code{8} (agency or institution), 
#'     \code{9} (someone else)}
#'   \item{famtime}{Has anyone in your family, such as your mother, father, brother, sister,
#'     ever served time in jail or prison? factor with levels \code{yes} and \code{no}}
#'   \item{crimtype}{most serious crime in current offense; one of \code{violent} (e.g. murder,
#'     rape, robbery, assault), \code{property} (e.g. burglary, larceny, arson, fraud, motor
#'     vehicle theft), \code{drug} (drug possession or trafficking), \code{publicorder} 
#'     (weapons violation, perjury, failure to appear in court), \code{juvenile} (juvenile-status
#'     offense, e.g. truancy, running away, incorrigible behavior)}
#'   \item{everviol}{Ever put on probation or sent to correctional institution for violent
#'     offense? factor with levels \code{no} and \code{yes}}
#'   \item{numarr}{number of times arrested (integer)}
#'   \item{probtn}{number of times on probation}
#'   \item{corrinst}{number of times previously committed to correctional institution}
#'   \item{evertime}{Prior to being sent here, did you ever serve time in a correctional
#'     institution? factor with levels \code{yes} and \code{no}}
#'   \item{prviol}{previously arrested for violent offense; factor with levels 
#'      \code{no} and \code{yes}}
#'   \item{prprop}{previously arrested for property offense; factor with levels 
#'      \code{no} and \code{yes}}
#'   \item{prdrug}{previously arrested for drug offense; factor with levels 
#'      \code{no} and \code{yes}}
#'   \item{prpub}{previously arrested for public-order offense; factor with levels 
#'      \code{no} and \code{yes}}
#'   \item{prjuv}{previously arrested for juvenile-status offense; factor with levels 
#'      \code{no} and \code{yes}}
#'   \item{agefirst}{age first arrested (integer)}      
#'   \item{usewepn}{Did you use a weapon... for this incident? factor with levels
#'      \code{yes} and \code{no}}
#'   \item{alcuse}{Did you drink alcohol at all during the year before being sent 
#'      here this time? factor with levels \code{yes}, \code{noduringyear}, \code{noatall}}
#'   \item{everdrug}{Ever used illegal drugs? factor with levels \code{no}, \code{yes}} 
#' }
#' @source Inter-University Consortium on Political and Social Research,
#'   NCJ-130915, U.S. Department of Justice 1989.
#' @references Lohr (1999). Sampling: Design and Analysis, Duxbury, p. 235--239 and
#'   445.
#' @export
NULL
