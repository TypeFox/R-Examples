#' County Wikipedia lookup table
#' 
#' Table of links to wikipedia articles on U.S. counties.  Adapted from table here: \url{http://en.wikipedia.org/wiki/List_of_United_States_counties_and_county_equivalents}.
#' 
#' @section Variables:
#' 
#' \itemize{
#'  \item \code{fips}: FIPS county code
#'  \item \code{county}: county name
#'  \item \code{state}: state abbreviation
#'  \item \code{pop2013}: 2013 population (from wikipedia table)
#'  \item \code{href}: web link to wikipedia page about the county
#' }
#' @docType data
#' @name wikiCounty
#' @usage wikiCounty
#' @format A data frame with 3143 rows and 5 columns.
#' @examples
#' head(wikiCounty)
NULL

