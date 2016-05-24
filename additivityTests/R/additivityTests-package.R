#' Testing for interaction in the two way ANOVA with single sub-class numbers.
#' 
#' @name additivityTests-package
#' @aliases additivityTests
#' @docType package
#' @title Additivity tests in the two way ANOVA with single sub-class numbers.
#' @author Petr Simecek <simecek@@gmail.com>
#' @description In many applications of statistical methods, it is assumed that the response variable is 
#' a sum of several factors and a random noise. In a real world this may not be an appropriate model. 
#' For example, some patients may react differently to the same drug treatment or the effect of fertilizer 
#' may be influenced by the type of a soil. There might exist an interaction between factors.
#' 
#' If there is more than one observation per cell then standard ANOVA techniques may be applied. Unfortunately, 
#' in many cases it is infeasible to get more than one observation taken under the same conditions. 
#' For instance, it is not logical to ask the same student the same question twice.
#' 
#' Six tests of additivity hypothesis (under various alternatives) are included in this package: 
#' Tukey test, modified Tukey test, Johnson-Graybill test, LBI test, Mandel test and Tussel test.
#' 
#' @keywords package
NULL
#' @name Boik
#' @title Multi-headed Machine Data
#' @description Performance of a multiple-headed machine used to fill bottles. Weights for six heads on five occasions were recorded.
#' @docType data
#' @usage data(Boik)
#' @source Robert J. Boik: A comparison of three invariant tests of additivity in two-way classifications with no replications, Computational Statistics \& Data Analysis, 1993.
#' @keywords datasets
NULL