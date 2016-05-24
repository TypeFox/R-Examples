#' This packages provides functions to estimate and visualize multilvel propensity
#' score analysis.
#' 
#' This package extends the principles put forth by the \code{PSAgraphics} 
#' (Helmreich, Pruzek, & Xiong, 2010) for multilevel, or clustered, data.
#'
#' Propensity score analyses are typically done in two phases. In phase I, a
#' statistical model prediciting treatment using the available individual covaraites 
#' is estimated. This package currently currently provides functions to perform
#' propensity score estimates using logistic regression (see \code{\link{mlpsa.logistic}})
#' and conditional inference trees (see \code{\link{mlpsa.ctree}}). The latter method
#' provides explicit stratifications as defined by each leaf node. The former however,
#' results in a numerical value ranging from zero to one (i.e. the fitted values).
#' A common approach is to then create stratificaitons using quintiles. However,
#' other approaches such as Loess regression are also provided.
#'
#' Phase II of typical pronsity score analyses concerns with the comparison of an
#' outcome between the treatment and comparison groups. The \code{\link{mlpsa}}
#' method will perform this analysis in a multilevel, or clustered, fashion. That
#' is, the results of the \code{\link{mlpsa}} procedure produce summary results
#' at level one (i.e. each strata within each cluster), level two (i.e. overall results
#' for each cluster), and overall (i.e. overall results across all clusters).
#'
#' This package also provides a number of visualizaions that provide a critical
#' part in presenting, understanding, and interpreting the results. See
#' \code{\link{plot.mlpsa}} for details.
#' 
#' @name multilevelPSA-package
#' @aliases multilevelPSA
#' @docType package
#' @title Multilevel Propensity Score Analysis
#' @author Jason Bryer \email{jason@@bryer.org}
#' @references \url{http://cran.r-project.org/web/packages/PSAgraphics/PSAgraphics.pdf}
#' 		\url{http://www.jstatsoft.org/v29/i06/}
#' @keywords propensity score analysis psa multilevel graphics
#' @seealso \code{PSAgraphics}
#' @import plyr
#' @import PSAgraphics
#' @import ggplot2
#' @import party
#' @importFrom grid grob gTree editGrob vpPath viewport vpTree grid.layout 
#'                  getGrob grobWidth grobHeight unit.c pushViewport grid.draw
#'                  upViewport grid.newpage
#' @importFrom reshape melt cast melt.data.frame
#' @importFrom psych describeBy
#' @importFrom MASS stepAIC
#' @importFrom stats binomial density fitted glm median model.matrix na.omit qt quantile sd var
#' @importFrom utils capture.output object.size setTxtProgressBar txtProgressBar
NA

#' North American (i.e. Canada, Mexico, and United States) student results of the 2009
#' Programm of International Student Assessment.
#'
#' Student results from the 2009 Programme of International Student Assessment (PISA)
#' as provided by the Organization for Economic Co-operation and Development (OECD).
#' See \url{http://www.pisa.oecd.org/} for more information including the code book.
#'
#' Note that missing values have been imputed using the 
#' \href{mice}{http://cran.r-project.org/web/packages/mice/index.html} package.
#' Details on the specific procedure are in the \code{pisa.impute} function
#' in the \href{http://github.com/jbryer/pisa}{\code{pisa} package}.
#' 
#' @references Organisation for Economic Co-operation and Development (2009).
#'             Programme for International Student Assessment (PISA). 
#'             \url{http://www.pisa.oecd.org/}
#' @name pisana
#' @docType data
#' @format a data frame with 66,548 ovservations of 65 variables.
#' @source Organization for Economic Co-operation and Development
#' @keywords datasets
NULL

#' Mapping of variables in \code{\link{pisana}} with full descriptions.
#' 
#' This data frame provides three variables, \code{Variable} corresponding to the
#' column names in \code{\link{pisana}}, \code{ShortDesc} providing a short
#' description of the variable as a valid R object name, and \code{Desc} providing
#' a longer description of the variable.
#' 
#' @name pisa.colnames
#' @docType data
#' @format a data frame with 50 rows of 3 variables.
#' @keywords datasets
NULL

#' Character vector representing the list of covariates used for estimating
#' propensity scores.
#' 
#' @name pisa.psa.cols
#' @docType data
#' @format a character vector with covariate names for estimating propensity scores.
#' @keywords datasets
NULL

#' Data frame mapping PISA countries to their three letter abbreviation.
#' 
#' This data frame has two columns, \code{CNT3} for the three letter abbreviation
#' of each country and \code{Country} that provides the full country name in English.
#' 
#' @name pisa.countries
#' @docType data
#' @format data frame with 65 rows of 2 variables.
#' @keywords datasets
NULL
