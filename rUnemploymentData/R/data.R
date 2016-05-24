#' A data.frame containing US State Unemployment Rates  
#'
#' Contains annualized data from 2000-2013. Data comes from the US Bureau
#' of Labor Statistics (BLS) Local Area Unemployment Statistics (LAUS) webpage: http://www.bls.gov/lau/.
#' @name df_state_unemployment
#' @usage data(df_state_unemployment)
#' @docType data
#' @references Created via build_state_df() on January 4, 2015.
#' @keywords data
#' @author Ari Lamstein
#' @examples
#' data(df_state_unemployment)
#' 
#' head(df_state_unemployment)
#' boxplot(df_state_unemployment[, -1],
#'         main="USA State Unemployment Data",
#'         xlab="Year", 
#'         ylab="Percent Unemployment")
#'         
#' \dontrun{
#' state_unemployment_choropleth(year=2013)
#' }
NULL

#' A data.frame Containing US County Unemployment Rates  
#'
#' Contains annualized data from 1990-2013. Data comes from the US Bureau
#' of Labor Statistics (BLS) Local Area Unemployment Statistics (LAUS) webpage: http://www.bls.gov/lau/.
#' The "region" column contains the numeric version of the County FIPS Code.
#' @name df_county_unemployment
#' @usage data(df_county_unemployment)
#' @docType data
#' @references Created via build_county_df() on January 4, 2015.
#' @keywords data
#' @author Ari Lamstein
#' @examples
#' 
#' data(df_county_unemployment)
#' 
#' head(df_county_unemployment)
#' boxplot(df_county_unemployment[, c(-1, -2, -3)],
#'         main="USA County Unemployment Data",
#'         xlab="Year", 
#'         ylab="Percent Unemployment")
#'         
#' \dontrun{
#' county_unemployment_choropleth(year=2013)
#' }
NULL