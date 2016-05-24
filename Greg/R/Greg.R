#' Regression Helper Functions
#' 
#' This R-package provides functions that primarily aimed at helping 
#' you work with regression models. While much of the data presented by the
#' standard regression output is useful and important - there is often a
#' need for further simplification prior to publication. The methods implemented
#' in this package are inspired by some of the top journals such as NEJM, BMJ,
#' and other medical journals as this is my research field. 
#' 
#' @section Output functions:
#' 
#' The package has function that automatically prints the crude unadjusted estimates
#' of a function next to the adjusted estimates, a common practice for medical
#' publications. 
#' 
#' The forestplot wrappers allows for easily displaying regression
#' estimates, often convenient for models with a large number of variables. 
#' There is also functionality that can help you comparing different models,
#' e.g. subsets of patients or compare different regression types.
#' 
#' @section Time-splitter:
#' 
#' When working with Cox regressions the proportional hazards can sometimes be violated.
#' As the \code{tt()} approach doesn't lend itself that well to big datasets I often
#' rely on time-splitting the dataset and then using the start time as an interaction
#' term. See the function \code{\link{timeSplitter}} and the associated \code{vignette("timeSplitter")}.
#' 
#' @section Other regression functions:
#' 
#' In addition to these funciton the package has some extentions to linear regression
#' where it extends the functionality by allowing for robust covariance matrices. 
#' by integrating the sandwich-package for rms::ols().
#' 
#' @section Important notice:
#' 
#' This package has an extensive test-set for ensuring that everything behaves as expected.
#' Despite this I strongly urge you to check that the values make sense. I commonly use
#' the regression methods available in the \code{rms}-package and in the \code{stats} package.
#' In addition I use the \code{\link[survival]{coxph}} in many of my analyses and should
#' also be safe. Please send me a notice if you are using the package with some other 
#' regression models, especially if you have some tests verifying the functionality.
#' 
#' @author Max Gordon
#' @name Greg-package
#' @docType package
NULL