#' @name rddtools
#' @docType package
#' @import KernSmooth
#' @import np
#' @import ggplot2
#' @title Regression Discontinuity Design 


utils::globalVariables(c("x", "y", "position", "cutpoint", "LATE", "CI_low", "CI_high", "sd", "quantile", "ks.test", "t.test", "coef", "density"))
utils::globalVariables(c("abline", "as.formula", "coef density", "df.residual", "fitted", "glm", "hist", "ksmooth",
"lines", "lm", "model.frame", "model.matrix", "na.pass", "par", "pnorm", "points", "poly",
"predict", "printCoefmat", "qnorm", "qt", "rbeta", "residuals", "rnorm", "segments", "title", "var", "vcov"))

#' @name indh
#' @docType data
#' @title INDH data set
#' @description Data from the Initiative Nationale du Development Humaine, collected as the part of the SNSF project "Development Aid and Social Dynamics"
#' @format A data frame with two variables with 720 observations each
#' @source Development Aid and social Dyanmics website: \url{http://qua.st/Development-Aid-Social-Dynamics}
#' @references Arcand, Rieger, and Nguyen (2015) 'Development Aid and Social Dyanmics Data Set'
#' @examples 
#' # load the data
#' data(indh)
#' 
#' # construct rdd_data frame
#' rdd_dat_indh <- rdd_data(y=choice_pg, x=poverty, data=indh, cutpoint=30)
#' 
#' # inspect data frame
#' summary(rdd_dat_indh)
#' 
#' # perform non-parametric regression
#' ( reg_np_indh <- rdd_reg_np(rdd_dat_indh) )
#' plot(reg_np_indh)
NULL
#' @name house
#' @docType data
#' @title Dataset used in Lee (2008)
#' @description Randomized experiments from non-random selection in U.S. House elections
#' @description Dataset described used in Imbens and Kalyamaran (2012), and probably the same dataset used in Lee (2008),
#' @format A data frame with 6558 observations and two variables:
#' \describe{
#' \item{x}{Vote at election t-1}
#' \item{y}{Vote at election t}
#' }
#' @source Guido Imbens webpage: \url{http://scholar.harvard.edu/imbens/scholar_software/regression-discontinuity}
#' @references Imbens, Guido and Karthik Kalyanaraman. (2012) 'Optimal Bandwidth Choice for the regression discontinuity estimator,' 
#' Review of Economic Studies (2012) 79, 933-959
#' @references   Lee, D. (2008) Randomized experiments from non-random selection in U.S. House elections, 
#' \emph{Journal of Econometrics}, 142, 675-697
#' @examples 
#' data(house)
#' rdd_house <- rdd_data(x=x, y=y, data=house, cutpoint=0)
#' summary(rdd_house)
#' plot(rdd_house)
NULL
#' @name STAR_MHE
#' @docType data
#' @title Transformation of the STAR dataset as used in Angrist and Pischke (2008)
#' @description Transformation of the STAR dataset as used in Table 8.2.1 of Angrist and Pischke (2008) 
#' @usage STAR_MHE
#' @seealso \code{\link[AER]{STAR}} for the original dataset.
#' @format A data frame containing 5743 observations and 6 variables. The first variable is from the original dataset, 
#' all other are created by Angrist and Pischke STAT code.
#' \describe{
#' \item{schidkn}{School ID in kindergarden (original variable, schoolidk in \code{\link[AER]{STAR}})}
#' \item{pscore}{The propensity score  (computed by A & P)}
#' \item{classid}{The id of the class (computed by A & P)}
#' \item{cs}{Class size (computed by A & P)}
#' \item{female, nwhite}{Various covariates (computed by A & P)}
#' }
#' @details ). This is a transformation of the dataset from the project STAR (Student/Teacher Achievement Ratio. 
#' The full dataset is described and available in package AER, \code{\link[AER]{STAR}}. 
#' The transformed data was obtained using the STATA script krueger.do, obtained from Joshua Angrist website 
#' (\url{http://economics.mit.edu/faculty/angrist/data1/mhe/krueger}), on the webstar.dta.
#' @references Krueger, A. (1999) 'Experimental Estimates Of Education Production Functions,' 
#' \emph{The Quarterly Journal of Economics}, Vol. 114(2), pages 497-532, May.
#' @references Angrist, A. ad  Pischke J-S (2008) \emph{Mostly Harmless Econometrics: An Empiricist's Companion}, 
#' Princeton University press
#' @source Data obtained using the script krueger.do on data webstar.rda, found on J. Angrist website 
#' \url{http://economics.mit.edu/faculty/angrist/data1/mhe/krueger}, retrieved on 26 November 2012.
#' @examples 
#' data(STAR_MHE)
#' 
#' # Compute the group means:
#' STAR_MHE_means <- aggregate(STAR_MHE[, c('classid', 'pscore', 'cs')],
#'                             by=list(STAR_MHE$classid), mean)
#' 
#' # Regression of means, with weighted average:
#' reg_krug_gls <- lm(pscore~cs, data=STAR_MHE_means, weights=cs)
#' coef(summary(reg_krug_gls))[2,2]
NULL