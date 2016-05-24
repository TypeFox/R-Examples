#' @title Performing Multivariate Dose-Response Meta-Analysis 
#' 
#' @name dosresmeta-package
#' 
#' @description The package \code{dosresmeta} consists of a collection of functions to estimate a dose-response relation from either a single or multiple summarized dose-response data.
#' The method was first formalized by Greenland and Longnecker (1992); the authors described how to approximate the covariances of reported log relative risks and how use them
#' to efficiently estimate an exposure-disease relation. The study specific estimates are combined through multivariate random-effect meta-analytical model, to obtaind
#' a pooled dose-response association.
#' 
#' 
#' @details
#' \tabular{ll}{
#' Package: \tab dosresmeta\cr
#' Type: \tab Package\cr
#' Version: \tab 1.3.2\cr
#' Date: \tab 2015-08-11\cr
#' License: \tab GPL-2\cr
#' }
#' @import mvmeta aod
#' @importFrom stats BIC delete.response lm model.frame model.matrix model.response optim pchisq pnorm qnorm symnum
#' @docType package
#' @author Alessio Crippa, \email{alessio.crippa@@ki.se}
#' @references 
#' Greenland, S., Longnecker, M. P. (1992). Methods for trend estimation from summarized dose-response data, with applications to meta-analysis. American journal of epidemiology, 135(11), 1301-1309.
#' 
#' Orsini, N., Bellocco, R.,  Greenland, S. (2006). Generalized least squares for trend estimation of summarized dose-response data. Stata Journal, 6(1), 40.
#' 
#' Hamling, J., Lee, P., Weitkunat, R., Ambuhl, M. (2008). Facilitating meta-analyses by deriving relative effect and precision estimates for alternative comparisons from a set of estimates presented
#' by exposure level or disease category. Statistics in medicine, 27(7), 954-970.
#'
#' Orsini, N., Li, R., Wolk, A., Khudyakov, P., Spiegelman, D. (2012). Meta-analysis for linear and nonlinear dose-response relations: examples, an evaluation of approximations, and software. 
#' American journal of epidemiology, 175(1), 66-73.
#' 
#' Gasparrini, A., Armstrong, B.,  Kenward, M. G. (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. Statistics in Medicine, 31(29), 3821-3839.
#' 
#' Liu, Q., Cook, N. R., Bergstrom, A., Hsieh, C. C. (2009). A two-stage hierarchical regression model for meta-analysis of epidemiologic nonlinear dose-response data.
#' Computational Statistics & Data Analysis, 53(12), 4157-4167.
#' 
#' 
#' @seealso \code{\link{dosresmeta}}, \code{\link{mvmeta}}
#' 
NULL



#' Case-control data on alcohol and breast cancer risk
#'
#' @name cc_ex
#' @description The dataset reports the summarized dose-response results from a case-control study
#' on alcohol and breast cancer, first presented by Rohan and McMichael.
#' @docType data
#' @format A data frame with 4 observations on the following 10 variables:
#' \tabular{ll}{
#' \code{gday} \tab label for exposure levels.\cr
#' \code{dose} \tab assigned dose level.\cr
#' \code{case} \tab number of cases for each exposure level.\cr
#' \code{control} \tab number of controls for each exposure level.\cr
#' \code{n} \tab total number of subjects for each exposure level.\cr
#' \code{crudeor} \tab unadjusted odds ratios for each exposure level.\cr
#' \code{adjrr} \tab adjusted odds ratios for each exposure level.\cr
#' \code{lb} \tab lower bound for the confidence limit of the adjusted odds ratios.\cr
#' \code{ub} \tab upper bound for the confidence limit of the adjusted odds ratios.\cr
#' \code{logrr} \tab natural logarithm of adjusted odds ratios.\cr
#' } 
#' 
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' @import Matrix aod mvmeta
#' @export coef.dosresmeta vcov.dosresmeta predict.dosresmeta print.dosresmeta summary.dosresmeta print.summary.dosresmeta 
#' 
#' @references 
#' Rohan, T. E., McMichael, A. J. (1988). Alcohol consumption and risk op breast cancer. International journal of cancer, 41(5), 695-699.
#' 
#' Greenland, S.,  Longnecker, M. P. (1992). Methods for trend estimation from summarized dose-response data, with applications to meta-analysis. American journal of epidemiology, 135(11), 1301-1309.
#' 
#' @keywords data
NULL


#' Six published studies on the relation between alcohol intake and vascular disease risk.
#' 
#' @name alcohol_cvd
#' @description The dataset reports the summarized dose-response results from six observational
#' studies on the relation between alcohol intake and vascular disease risk (Qin Liu 2009). Four are case-control studies, 
#' two prospective (cumulative-incidence data). 
#' @docType data
#' @format A data frame with 25 observations on the following 8 variables:
#' \tabular{ll}{
#' \code{id} \tab id of the studies included in the analysis.\cr
#' \code{author} \tab names of the first author of the study.\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose level.\cr
#' \code{case} \tab number of cases for each exposure level.\cr
#' \code{n} \tab total number of subjects for each exposure level.\cr
#' \code{logrr} \tab natural logarithm of adjusted "relative risks".\cr
#' \code{se} \tab standard errornatural for the logarithm of adjusted "relative risks".\cr
#' }
#' 
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references 
#' Liu, Q., Cook, N. R., Bergstrom, A., Hsieh, C. C. (2009). A two-stage hierarchical regression model for meta-analysis of epidemiologic nonlinear dose-response data. Computational Statistics & Data Analysis, 53(12), 4157-4167.
#'   
#' @keywords data
NULL


#' Eight published studies on the relation between alcohol intake and colon-rectal cancer.
#' 
#' @name alcohol_crc
#' @description The dataset reports the summarized dose-response results from eight prospective
#' #' studies on the relation between alcohol intake and colon-rectal risk (Orsini 2012).
#' @docType data
#' @format A data frame with 48 observations on the following 7 variables:
#' \tabular{ll}{
#' \code{id} \tab label for author's names (id variable).\cr
#' \code{type} \tab code for study design.\cr
#' \code{dose} \tab assigned dose level.\cr
#' \code{cases} \tab number of cases for each exposure level.\cr
#' \code{peryears} \tab amount of person-time for each exposure level.\cr
#' \code{logrr} \tab natural logarithm of adjusted "relative risks".\cr
#' \code{se} \tab standard errornatural for the logarithm of adjusted "relative risks".\cr
#' }
#' 
#' @author Alessio Crippa, <\email{alessio.crippa@@ki.se}>
#' 
#' @references 
#' Orsini, N., Li, R., Wolk, A., Khudyakov, P., Spiegelman, D. (2012). Meta-analysis for linear and nonlinear dose-response relations: examples, an evaluation of approximations, and software. 
#' American journal of epidemiology, 175(1), 66-73.
#'   
#' @keywords data
NULL