#' @title Read European Fluxes CSV Files
#'
#' @description The European Eddy Fluxes Database Cluster distributes fluxes of different Green House Gases measured mainly using the eddy covariance technique acquired in sites involved in EU projects but also single sites in Europe, Africa and others continents that decided to share their measurements in the database (cit. http://gaia.agraria.unitus.it ). The package provides two functions to load and row-wise bind CSV files distributed by the database. Currently only L3 and L4 (L=Level), half-hourly and daily (aggregation) files are supported.
#'
#' @name efreadr-package
#' @docType package
#' @author Marco Bascietto \email{marco.bascietto@@crea.gov.it}
#' @keywords package
#' @references Source code is hosted at GitHub (\url{https://github.com/mbask/efreadr})
NULL

#' @author Marco Bascietto \email{marco.bascietto@@crea.gov.it}
`: dataframe_with_filename_and_siteid` <- ensures_that(
  is.data.frame(.),
  sum(c("efreadr_year", "efreadr_file_name", "efreadr_site_id") %in% colnames(.)) == 3,
  err_desc = "Something wrong with the returned dataframe, are any fluxes files present in the directories?")

globalVariables(c(".", "starts_with", "funs"))