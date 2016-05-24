##
## This is the PACKAGE documentation
##

#' Constrained inference for linear mixed models.
#' 
#' @docType package
#' @name CLME-package
#' @rdname CLME-package
#' 
#' @description
#' Constrained inference on linear fixed and mixed models using residual bootstrap. 
#' Covariates and random effects are permitted but not required.
#'
#' Appropriate credit should be given when publishing results obtained using \pkg{CLME}, or when 
#' developing other programs/packages based off of this one. Use \code{citation(package="CLME")} 
#' for Bibtex information.
#' 
#' The work was produced in part with funding from the Intramural Research Program of the NIH, 
#' National Institute of Environmental Health Sciences (Z01 ES101744).
#' 
#' @author Casey M. Jelsema <\email{casey.jelsema@@gmail.com}>
#' @author Shyamal D. Peddada
#' 
#' @note
#' Unless otherwise noted, Casey M. Jelsema wrote the functions in this package.
#' 
#' @details
#' \tabular{ll}{
#' Package: \tab CLME\cr
#' Type:    \tab Package\cr
#' Version: \tab 2.0.4\cr
#' Date:    \tab 2015-05-18\cr
#' License: \tab GLP-3  \cr
#' }
#' 
#' This package was introduced in Jelsema and Peddada (2015, under review). The method which 
#' is implemented is the constrained linear mixed effects model described in 
#' Farnan, Ivanova, and Peddada (2014). See that paper for more details regarding the method. 
#' Here we give a brief overview of the assumed model:
#' 
#' \deqn{ Y = X_{1}\theta_{1} + X_{2}\theta_{2} + U\xi + \epsilon }{Y = X1*theta1 + X2*theta2 + U*xi + e}
#' 
#' where
#' 
#' \itemize{
#' \item \eqn{X_1}{X1} is a \eqn{N \times p_1}{N x p1} design matrix.
#' \item \eqn{\theta_1}{theta1} are the coefficients (often treatment effects).
#' \item \eqn{X_2}{X2} is a \eqn{N \times p_2}{N x p2} matrix of fixed covariates.
#' \item \eqn{\theta_1}{theta2} are the coefficients for the covariates.
#' \item \eqn{U}{U} is a \eqn{N \times c}{N x c} matrix of random effects.
#' \item \eqn{\xi}{xi} is a zero-mean random vector with covariance \eqn{T}{T} (see below).
#' \item \eqn{\epsilon}{e} is a zero-mean random vector with covariance \eqn{\Sigma}{Sigma} (see below).
#' }
#' 
#' Neither covariates (\eqn{X_2}{X2}) nor random effects (\eqn{U}{U}) are required by the model or \pkg{CLME}. The covariance matrix of \eqn{\xi}{xi} is given by:
#'   
#' \deqn{ T = diag\left( \tau^{2}_{1}I_{c_{1}}, \tau^{2}_{2}I_{c_{2}}  , \dots , \tau^{2}_{q}I_{c_{q}} \right) }{ T = diag( tau1^2 I_c1, tau2^2 I_c2 , ... , tauq^2 I_cq) }
#' 
#' The first \eqn{c_{1}}{c1} random effects will share a common variance, \eqn{\tau^{2}_{1}}{tau1^2}, the next \eqn{c_{2}}{c2} random effects will share a common variance, and so on. Note that \eqn{c = \sum_{i=1}^{q} c_i}{c = SUM(ci), i=1,...q}. Homogeneity of variances in the random effects can be induced by letting \eqn{q=1}{q=1} (hence \eqn{c_{1}=c=ncol(U)}{c1=c=ncol(U)}).
#' 
#' Similarly, the covariance matrix of \eqn{\epsilon}{e} is given by:
#'   
#' \deqn{ \Sigma = diag\left( \sigma^{2}_{1}I_{n_{1}}, \sigma^{2}_{2}I_{n_{2}}  , \dots , \sigma^{2}_{q}I_{n_{k}} \right) }{ Sigma = diag( sigma1^2 I_n1, sigma2^2 I_n2 , ... , sigmak^2 I_nk)}
#' 
#' Again, the first \eqn{n_{1}}{n1} observations will share a common variance, \eqn{\sigma^{2}_{1}}{sigma1^2}, the next \eqn{n_{2}}{n2} will share a common variance, and so on. Note that \eqn{N = \sum_{i=1}^{k} n_i}{N = SUM(n_i), i=1,...k}. Homogeneity of variances in the residuals can be induced by letting \eqn{k=1}{k=1}.
#' 
#' The order constraints are defined by the matrix \eqn{A}{A}. This is an \eqn{r \times p}{r x p} matrix where \eqn{r}{r} is the number of constraints, and \eqn{p=p_{1}+p_{2}}{p = p1 + p2} is the dimension of \eqn{ \theta = ( \theta_{1}' , \theta_{2}')'}{ theta = ( theta1' , theta2')'}. Formally the hypothesis being tested is:
#'   
#' \deqn{ H_{a}: A\theta > 0 }{Ha: A*theta > 0 }
#' 
#' For several default orders (simple, umbrella, simple tree) the \eqn{A}{A} matrix can be automatically generated. Alternatively, the user may define a custom \eqn{A}{A} matrix to test other patterns among the elements of \eqn{\theta}{theta}. See \code{\link{create.constraints}} and \code{\link{clme}} for more details.
#' 
#' For computational reasons, the implementation is not identical to the model expressed. Particularly, the fixed-effects matrix (or matrices) and the random effects matrix are assumed to be columns in a data frame, not passed as matrices. The \eqn{A}{A} matrix is not \eqn{r\ times p}{r x p}, but \eqn{r\ times 2}{r x 2}, where each row gives the indices of the constrained coefficients. See \code{\link{create.constraints}} for further explanation.
#' 
#' The primary function of \pkg{CLME} is \code{\link{clme}}. The other functions in this package may be run separately, but in general are designed for use by \code{\link{clme}}.
#' 
#' The creation of this package \pkg{CLME}, this manual, and the vingette were all supported by the Intermural Research Program of the United States' National Institutes of Health (Z01 ES101744).
#' 
#' @references
#' Farnan, L., Ivanova, A., and Peddada, S. D. (2014).
#' Linear Mixed Efects Models under Inequality Constraints with Applications.
#' \emph{PLOS ONE}, 9(1). e84778. doi: 10.1371/journal.pone.0084778
#' \url{http://www.plosone.org/article/info:doi/10.1371/journal.pone.0084778}
#' 
#' Jelsema, C.M. and Peddada, S.D. (2015).
#' CLME: An R Package for Linear Mixed Effects Models under Inequality Constraints.
#' \emph{(under review)}.
#' 
#' 
#' @import methods
#' @import stats
#'
#' 
#' 
NULL



#' @title 
#' Fibroid Growth Study
#' 
#' @description
#' This data set contains a subset of the data from the Fibroid Growth Study.
#'
#' \tabular{rll}{
#' [,1] \tab  ID  \tab ID for subject. \cr
#' [,2] \tab fid  \tab ID for fibroid (each women could have multiple fibroids). \cr
#' [,3] \tab lfgr \tab log fibroid growth rate. See details. \cr
#' [,4] \tab age  \tab age category Younger, Middle, Older. \cr
#' [,5] \tab loc  \tab location of fibroid, corpus, fundus, or lower segment. \cr
#' [,6] \tab bmi  \tab body mass index of subject. \cr
#' [,7] \tab preg \tab parity, whether the subject had delivered a child. \cr
#' [,8] \tab race \tab race of subject (Black or White only). \cr
#' [,9] \tab vol  \tab initial volume of fibroid. \cr
#' }
#' 
#' @details
#' The response variable \code{lfgr} was calculated as the change in log fibroid volume, 
#' divided by the length of time between measurements. The growth rates were averaged to produce
#'  a single value for each fibroid, which was scaled to represent a 6-month percent change in volume.
#' 
#' @references
#' Peddada, Laughlin, Miner, Guyon, Haneke, Vahdat, Semelka, Kowalik, Armao, Davis, and Baird(2008).
#'  Growth of Uterine Leiomyomata Among Premenopausal Black and White Women.
#'   Proceedings of the National Academy of Sciences of the United States of America, 105(50),
#'    19887-19892. URL \url{http://www.pnas.org/content/105/50/19887.full.pdf}.
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name fibroid
#' @usage data(fibroid)
#' @format A frame containing 240 observations on 9 variables.
NULL



#' @title 
#' Experiment on mice
#' 
#' @description
#' This data set contains the data from an experiment on 24 Sprague-Dawley rats from Cora et al (2012).
#'
#' \tabular{rll}{
#' [,1]  \tab id   \tab ID for rat (factor). \cr
#' [,2]  \tab time \tab time period (in order, 0 , 6, 24, 48, 72, 96 hours). \cr
#' [,3]  \tab temp \tab storage temperature reference (\code{''Ref''}) vs. room temperature (\code{''RT''}). \cr
#' [,4]  \tab sex  \tab sex, male (\code{''Male''}) vs. female (\code{''Female''}). Coded as \code{''Female''=1}. \cr
#' [,5]  \tab wbc  \tab white blood cell count (\eqn{10^3 / \mu L}{10^3 / mu L}). \cr
#' [,6]  \tab rbc  \tab red blood cell count )\eqn{10^6 / \mu L}{10^6 / mu L}). \cr
#' [,7]  \tab hgb  \tab hemoglobin concentration (g/dl). \cr
#' [,8]  \tab hct  \tab hematocrit (\%). \cr
#' [,9]  \tab spun \tab (HCT \%). \cr
#' [,10] \tab mcv  \tab MCV, a measurement of erythrocyte volume (fl). \cr
#' [,11] \tab mch  \tab mean corpuscular hemoglobin (pg). \cr     % ????
#' [,12] \tab mchc \tab mean corpuscular hemoglobin concentration (g/dl). \cr
#' [,13] \tab plts \tab platelet count (\eqn{10^3 / \mu L}{10^3 / mu L}). \cr
#' }
#' 
#' @details
#' The response variable \code{lfgr} was calculated as the change in log fibroid volume, 
#' divided by the length of time between measurements. The growth rates were averaged to produce
#'  a single value for each fibroid, which was scaled to represent a 6-month percent change in volume.
#' 
#' @references
#' Cora M, King D, Betz L, Wilson R, and Travlos G (2012). 
#' Artifactual changes in Sprauge-Dawley rat hematologic parameters after storage of samples at 3 C and 21 C.
#'  Journal of the American Association for Laboratory Animal Science, 51(5), 616-621. 
#'  URL \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3447451/}.
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name rat.blood
#' @usage data(rat.blood)
#' @format A frame containing 241 observations on 13 variables.
NULL














