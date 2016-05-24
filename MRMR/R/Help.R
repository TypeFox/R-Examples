#' Multivariate Regression Models for Reserving
#' @description
#'
#' MRMR allows an actuary to create sets of loss data and forecast liabilities. It uses a set of 3 S4 objects to store data, models and predictions.
#' 
#' @details
#' \strong{Triangle}
#' 
#' A Triangle is a collection of aggregate loss data. All triangles must have a defined set of OriginPeriods, a defined set of DevelopmentIntervals 
#' and data along those axes. A triangle may carry additional descriptive information such as line of business, geographic region and so on.
#' 
#' \strong{TriangleModel}
#' 
#' A TriangleModel is a statistical model fit to triangle data. The formula may be defined by the user and will generally 
#' be a linear or generalized linear model. A triangle may have more than one model. It usually will.
#' 
#' \strong{TriangleProjection}
#' 
#' A TriangleProjection is a prediction based on a TriangleModel. A TriangleModel may have more than one projection.
#' 
#'
#' @docType package
#' @name MRMR
#' @aliases MRMR MRMR-package
#' @importFrom(plyr, daply)
#' @importFrom(plyr, dlply)
NULL