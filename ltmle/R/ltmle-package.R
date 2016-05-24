

#' Targeted Maximum Likelihood Estimation for Longitudinal Data
#' 
#' Targeted Maximum Likelihood Estimation (TMLE) of treatment/censoring
#' specific mean outcome or marginal structural model for point-treatment and
#' longitudinal data. Also provides Inverse Probability of Treatment/Censoring
#' Weighted estimate (IPTW) and maximum likelihood based G-computation estimate
#' (G-comp). Can be used to calculate additive treatment effect, risk ratio,
#' and odds ratio.
#' 
#' 
#' @name ltmle-package
#' @docType package
#' @author Joshua Schwab, Samuel Lendle, Maya Petersen, and Mark van der Laan,
#' with contributions from Susan Gruber
#' 
#' Maintainer: Joshua Schwab \email{joshuaschwab@@yahoo.com}
#' @seealso \code{\link{ltmle}}
#' @references 
#' Bang, Heejung, and James M. Robins. "Doubly robust estimation in missing data
#' and causal inference models." Biometrics 61.4 (2005): 962-973.
#' 
#' Lendle, Samuel, Schwab, Joshua, Petersen, Maya and and van der
#' Laan, Mark J "ltmle: An R Package Implementing Targeted Minimum Loss-based
#' Estimation for Longitudinal Data", Forthcoming
#' 
#' Petersen, Maya, Schwab, Joshua and van der Laan, Mark J, "Targeted Maximum
#' Likelihood Estimation of Marginal Structural Working Models for Dynamic
#' Treatments Time-Dependent Outcomes", Journal of Causal Inference, 2014
#' \url{http://www.ncbi.nlm.nih.gov/pubmed/25909047}
#' 
#' Robins JM, Sued M, Lei-Gomez Q, Rotnitsky A. (2007). Comment: Performance of 
#' double-robust estimators when Inverse Probability weights are highly 
#' variable. Statistical Science 22(4):544-559.
#' 
#' van der Laan, Mark J. and Gruber, Susan, "Targeted Minimum Loss Based
#' Estimation of an Intervention Specific Mean Outcome" (August 2011). U.C.
#' Berkeley Division of Biostatistics Working Paper Series. Working Paper 290.
#' \url{http://biostats.bepress.com/ucbbiostat/paper290}
#' 
#' van der Laan, Mark J. and Rose, Sherri, "Targeted Learning: Causal Inference
#' for Observational and Experimental Data" New York: Springer, 2011.
#' @keywords package
#' @examples
#' 
#' ## For examples see examples(ltmle)
#' 
NULL




#' Sample data, regimes, and summary measures
#' 
#' Sample data for use with ltmleMSM. Data: n=1000: male age CD4_1 A1 Y1 CD4_2
#' A2 Y2 CD4_3 A3 Y3 A1..A3 are treatment nodes, Y1..Y3 are death, CD4_1..CD4_3
#' are time varying covariates. We are interested in static regimes where a
#' patient switches at some time. In summary.measures, switch.time is first
#' time where At is 1 (4 if never switch), time is the horizon.
#' 
#' regimes: 200 x 3 x 4 [n x numACnodes x numRegimes] summary.measures: 4 x 2 x
#' 3 [numRegimes x numSummaryMeasures x numFinalYnodes]
#' 
#' @name sampleDataForLtmleMSM
#' @docType data
#' @format List with three components: data, regimes, summary.measures
#' @source simulated data
#' @keywords datasets
#' @examples
#' 
#' data(sampleDataForLtmleMSM)
#' 
NULL




