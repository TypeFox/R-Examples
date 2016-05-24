#' postGIStools: Tools for interacting with PostgreSQL / PostGIS Databases
#'
#' postGIStools facilitates the import/export of data tables between R
#' and PostgreSQL, in particular those with associated geometries (PostGIS
#' extension) and hstore (key-value pairs) type columns.
#'
#' @section Key Functions:
#' \code{\link{get_postgis_query}} works like \code{\link[DBI]{dbGetQuery}},
#' with the additional benefit of parsing hstore types (as a list-column in the
#' resulting R data frame) and geometry types (producing a spatial data frame
#' in R).
#'
#' \code{\link{\%->\%}} reproduces the behavior of the PostgreSQL
#' \code{hstore -> key} operator.
#'
#' \code{\link{postgis_insert}} and \code{\link{postgis_update}} respectively
#' insert new rows or update existing rows in a PostgreSQL table based on the
#' contents of a R data frame. For spatial data and hstore columns, they
#' performs the same conversions as \code{get_postgis_query}, in reverse.
#'
#' @docType package
#' @name postGIStools-package
#' @aliases postGIStools
#'
#' @import sp
#' @importFrom methods is
NULL
