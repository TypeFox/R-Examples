#' Unified interface to your data sources.
#'
#' This package provides several S4 classes that make it easier to
#' collect and work with datasets. The package is inspired
#' by the \href{https://bitbucket.org/ScraperWiki/scraperwiki}{scraperwiki project}, 
#' which provides a webbased service for data collection. Also inspiring 
#' are \href{http://reference.wolfram.com/mathematica/ref/CountryData.html}{Mathematica's xxxData functions}, 
#' which provide in-built parametrizable datasets.
#' 
#' You can specify web resources with the \code{urldata} and the \code{xsparql} functions. For working with locally saved data,
#' see the \code{internalData} and the \code{csvdata} function. The objects instantiated with these functions can than be passed
#' to the generic \code{query} along with some parameters to get to the data.
#'
#' You can combine several resources with the \code{datamart} function. 
#'
#' Besides parameterized queries ("read" operations), the package also aims to support "write" operations.
#' For this purpose, some functions (currently \code{mdreport}, \code{swvreport}) for defining targets
#' as well as some functions (currently \code{blogger} and \code{dirloc}) for defining locations 
#' are provided. The generic \code{put} then builds the target and puts it at the defined location.
#'
#' Some examples aim to proof the concept, for instance \code{dbpedia}, \code{sourceforge}, \code{expenditures},
#' and \code{city_coords}.
#' 
#' The package is highly experimental, and likely to change heavily without backward compatiblity.
#'
#' @references Karsten Weinert, \href{http://factbased.blogspot.com/search/label/datamart}{factbased blogspot.}
#' @docType package
#' @name datamart-package
NULL
