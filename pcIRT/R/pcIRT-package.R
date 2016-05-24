#' Data set B5PO-R personality questionnaire
#'
#' This object contains data from the B5PO-R (Kubinger and Siptroth, in prep.), the revised version of the B5PO (Holocher-Ertl, Kubinger and Menghin, 2003). The items are from the extraversion scale.
#'#'
#' @name extraversion
#' @docType data
#' @format workspace
#' @source Study
#' @keywords datasets
#' @references Holocher-Ertl, S., Kubinger, K. D. and Menghin, S. (2003). Big Five Plus One Persoenlichkeitsinventar (B5PO). Moedling: Wiener Testsystem/Schuhfried.
NULL





#' Data set META reasoning test.
#'
#' This object contains data from the reasoning test 'META' by Gatternig and Kubinger (1994). The test includes 11 encoding tasks.
#'#'
#' @name reason.test
#' @docType data
#' @format A matrix with 22 variables and 380 observations. Variables 'I1' to 'I11' contain the responses to the eleven items, 'BT1' to 'BT11' the response times for each item in seconds.
#' @source Study
#' @keywords datasets
#' @references Gatternig, J. and Kubinger, K. D. (1994). Erkennen von Metaregeln. Frankfurt: Swets.
#'
NULL





#' IRT Models for Polytomous and Continuous Item Responses
#'
#' The multidimensional polytomous Rasch model (Rasch,
#' 1961) can be estimated with pcIRT. It provides functions to set linear restrictions on the item category
#' parameters of this models. With this functions it is possible to test
#' whether item categories can be collapsed or set as linear dependent. Thus it
#' is also possible to test whether the multidimensional model can be reduced
#' to a unidimensional model that is whether item categories represent a
#' unidimensional continuum. For this case the scoring parameter of the
#' categories is estimated.
#'
#' This package estimates the Continuous Rating Scale model by Mueller (1987).
#' It is an extension of the Rating Scale Model by Andrich (1978) on continuous
#' responses (e.g. taken by a visual analog scale).
#'
#' \tabular{ll}{ Package: \tab pcIRT\cr Type: \tab Package\cr Version: \tab
#' 0.1\cr Date: \tab 2013-11-13\cr License: \tab GPL-3\cr }
#'
#' @name pcIRT-package
#' @aliases pcIRT-package pcIRT
#' @docType package
#' @author Christine Hohensinn Maintainer: Christine Hohensinn
#' <christine.hohensinn@@univie.ac.at>
#' @seealso \code{\link{MPRM}} \code{\link{CRSM}}
#' @references Andersen, E. B. (1995). Polytomous Rasch models and their
#' estimation. In G. H. Fischer and I. Molenaar (Eds.). Rasch Models -
#' Foundations, Recent Developements, and Applications. Springer.
#'
#' Fischer, G. H. (1974). Einfuehrung in die Theorie psychologischer Tests
#' [Introduction to test theory]. Bern: Huber.
#'
#' Mueller, H. (1987). A Rasch model for continuous ratings. Psychometrika, 52,
#' 165-181.
#'
#' Rasch, G. (1961). On general laws and the meaning of measurement in
#' psychology, Proceedings Fourth Berekely Symposium on Mathematical
#' Statistiscs and Probability 5, 321-333.
#' @keywords package IRT Item Response Theory psychometrics multidimensional
#' polytomous Rasch model Continuous Rating Scale Model
#' 
#' @import stats graphics grDevices utils methods
#' 
#' @examples
#'
#' #simulate data set according to the multidimensional polytomous Rasch model (MPRM)
#' simdat <- simMPRM(rbind(matrix(c(-1.5,0.5,0.5,1,0.8,-0.3, 0.2,-1.2), ncol=4),0), 500)
#'
#' #estimate MPRM item parameters
#' res_mprm <- MPRM(simdat$datmat)
#'
#' summary(res_mprm)
#'
#'
NULL



