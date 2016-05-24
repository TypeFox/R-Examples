#' U.S. Census Region and Demographic Data
#'
#' The R package \code{noncensus} provides a collection of various regional
#' information determined by the U.S. Census Bureau along with demographic
#' data. We have also included scripts to download, process, and load the
#' original data from their sources.
#'
#' @docType package
#' @name noncensus
#' @aliases noncensus package-noncensus
NULL

#' Demographic Data and Census Regions of U.S. States and Territories
#' 
#' A dataset containing demographic information and the census regions of each
#' U.S. state as defined by the U.S. Census Bureau. Also included are the
#' U.S. territories, such as Puerto Rico and Guam.
#'
#' The variables included are:
#'
#' \itemize{
#'   \item state. State abbreviation
#'   \item name. State name
#'   \item region. Region as defined by the U.S. Census Bureau
#'   \item division. Subregion as defined by the U.S. Census Bureau
#'   \item capital. Capital city.
#'   \item area. Land area in square miles
#'   \item population. Population from 2010 Census
#' }
#' 
#' The U.S. is divided into four regions:
#'
#' \enumerate{
#'   \item Midwest
#'   \item Northeast
#'   \item South
#'   \item West
#' }
#'
#' Within each region, states are further partitioned into divisions. For more
#' details about census regions, see:
#' \url{http://en.wikipedia.org/wiki/List_of_regions_of_the_United_States#Census_Bureau-designated_regions_and_divisions}
#'
#' Much of the state data was extracted from
#' \url{http://www.census.gov/popest/data/state/totals/2013/index.html}
#' 
#' @docType data
#' @keywords datasets
#' @name states
#' @usage data(states)
#' @format A data frame with 56 rows and 7 variables
NULL

#' Data for U.S. Counties and County-Equivalent Entities
#'
#' Data containing state and county FIPS codes for U.S. counties and
#' county-equivalent entities (CEE) along with county-level demographic
#' data. The CEE includes non-state locations, such as Puerto Rico (PR) and Guam
#' (GU).
#'
#' \itemize{
#'   \item county_name. County Name and Legal/Statistical Area Description
#'   \item state. State Postal Code
#'   \item state_fips. State FIPS Code
#'   \item county_fips. County FIPS Code
#'   \item fips_class. FIPS Class Code
#'   \item CSA. Combined Statistical Area
#'   \item CBSA. Core-based Statistical Area
#'   \item population. County population from 2010 Census
#' }
#'
#' The U.S. Census Bureau groups counties into CSAs and CBSAs primarily based on
#' county population. We provide listings of both in
#' \code{\link[noncensus]{combined_areas}} and
#' \code{\link[noncensus]{corebased_areas}}.
#'
#' For a detailed description, Wikipedia has excellent discussions of both
#' areas: \url{http://en.wikipedia.org/wiki/Combined_Statistical_Area} and
#' \url{http://en.wikipedia.org/wiki/Core_Based_Statistical_Area}.  Also, the
#' following map from Wikipedia is excellent to visualize the areas:
#' \url{http://upload.wikimedia.org/wikipedia/commons/7/7b/Combined_statistical_areas_of_the_United_States_and_Puerto_Rico.gif}
#'
#' NOTE: Not all counties are members of a CSA or CBSA.
#'
#' The following details about FIPS Class Codes have been blatantly taken from
#' the Census Bureau's website:
#'
#' \itemize{
#'   \item H1. Identifies an active county or statistically equivalent entity that does not qualify under subclass C7 or H6.
#'   \item H4. Identifies a legally defined inactive or nonfunctioning county or statistically equivalent entity that does not qualify under subclass H6.
#'   \item H5. Identifies census areas in Alaska, a statistical county equivalent entity.
#'   \item H6. Identifies a county or statistically equivalent entity that is areally coextensive or governmentally consolidated with an incorporated place, part of an incorporated place, or a consolidated city.
#'   \item C7: Identifies an incorporated place that is an independent city; that is, it also serves as a county equivalent because it is not part of any county, and a minor civil division (MCD) equivalent because it is not part of any MCD.
#' }
#'
#' For more details, see:
#' \url{http://www.census.gov/geo/reference/codes/cou.html}
#' 
#' @docType data
#' @keywords datasets
#' @name counties
#' @usage data(counties)
#' @format A data frame with 3235 rows and 6 variables
NULL


#' Data for U.S. Cities by Zip Code
#'
#' This data set considers each zip code throughout the U.S. and provides
#' additional information, including the city and state, latitude and longitude,
#' and the FIPS code for the corresponding county.
#'
#' The ZIP code data was obtained from version 1.0 of the \code{\link[zipcode]{zipcode}}
#' package on CRAN. The county FIPS codes were obtained by querying the FIPS
#' code from each zip's latitude and longitude via the FCC's Census Block
#' Conversions API. For details regarding the API, see
#' \url{http://www.fcc.gov/developers/census-block-conversions-api}.
#'
#' \itemize{
#'   \item zip. U.S. ZIP (postal) code
#'   \item city. ZIP code's city
#'   \item state. ZIP code's state
#'   \item latitude. ZIP code's latitude
#'   \item longitude. ZIP code's longitude
#'   \item fips. County FIPS Code
#' }
#'
#' The FIPS codes are useful for mapping ZIP codes and cities to counties in the
#' \code{\link[noncensus]{counties}} data set.
#'
#' Fun fact: ZIP is an acronym for "Zone Improvement Plan."
#' 
#' @docType data
#' @keywords datasets
#' @name zip_codes
#' @usage data(zip_codes)
#' @format A data frame with 43524 rows and 6 variables
NULL


#' Combined Statistical Areas (CSAs)
#'
#' The U.S. Census Bureau groups counties into CSAs primarily based on county
#' population. NOTE: Not all counties are members of a CSA. For a detailed
#' description, Wikipedia has an excellent discussion:
#' \url{http://en.wikipedia.org/wiki/Combined_Statistical_Area}.  Also, the
#' following map from Wikipedia is excellent to visualize the areas:
#' \url{http://upload.wikimedia.org/wikipedia/commons/7/7b/Combined_statistical_areas_of_the_United_States_and_Puerto_Rico.gif}
#'
#' @docType data
#' @keywords datasets
#' @name combined_areas
#' @usage data(combined_areas)
#' @format A data frame with 166 rows and 2 variables.
NULL


#' Core-based Statistical Area (CBSAs)
#'
#' The U.S. Census Bureau groups counties into CBSAs primarily based on county
#' population. NOTE: Not all counties are members of a CBSA. For a detailed
#' description, Wikipedia has an excellent discussion:
#' \url{http://en.wikipedia.org/wiki/Core_Based_Statistical_Area}.  Also, the
#' following map from Wikipedia is excellent to visualize the areas:
#' \url{http://upload.wikimedia.org/wikipedia/commons/7/7b/Combined_statistical_areas_of_the_United_States_and_Puerto_Rico.gif}
#'
#' @docType data
#' @keywords datasets
#' @name corebased_areas
#' @usage data(corebased_areas)
#' @format A data frame with 917 rows and 4 variables.
NULL
