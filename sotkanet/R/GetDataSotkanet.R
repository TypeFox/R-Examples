#' Retrieves Sotkanet data according to the query arguments and combines the
#' indicator, region, and overall data into a single table
#'
#' Arguments:
#'   @param indicators Dataset identifier(s)
#'   @param years vector of years c(2010, 2012, ... )
#'   @param genders vector of genders ('male' | 'female' | 'total')
#'   @param regions pick selected regions only (default: all regions)
#'   @param region.category return selected regions category (for options, see:
#'    	    unique(SotkanetRegions(type = "table")$region.category)); 
#'	    "ALUEHALLINTOVIRASTO, "ERVA", "EURALUEET", "EUROOPPA", "KUNTA", 
#'	    "MAA", "MAAKUNTA", "NUTS1", "POHJOISMAAT", "SAIRAANHOITOPIIRI", 
#'	    "SEUTUKUNTA", "SUURALUE"   
#'   @param verbose verbose
#'
#' Returns:
#'   @return data.frame
#'
#' @export
#' @references
#' See citation("sotkanet") 
#' @author Einari Happonen. Maintainer: Louhos \email{louhos@@googlegroups.com}
#' @examples # dat <- GetDataSotkanet(indicators = 10013)
#' @keywords utilities

GetDataSotkanet <- function (indicators, years = 1991:2015, genders = c("total"), regions = NULL, region.category = NULL, verbose = TRUE) {

  # List all indicators in Sotkanet database
  sotkanet.indicators <- SotkanetIndicators(type = "table")

  dats <- list()
  for (indicator in indicators) { 
    if (verbose) {message(paste("Retrieving indicator", indicator))}
    dats[[as.character(indicator)]] <- GetDataSotkanetSingleIndicator(indicator, years = years, 
    genders = genders, regions = regions, region.category = region.category) 
  }

  # Merge all data from the different indicators in a single table
  combined.data <- do.call("rbind", dats)

  # Add indicator information
  combined.data$indicator.organization.title.fi <- sotkanet.indicators[match(combined.data$indicator, 
  	sotkanet.indicators$indicator), "indicator.organization.title.fi"]
  
  combined.data

}




#' Description:
#' GetDataSotkanetSingleIndicator retrieves Sotkanet data 
#' for given indicator according to the query arguments and combines
#' indicator, region, and overall data into one table
#'
#' Arguments:
#'   @param indicator Dataset identifier
#'   @param years vector of years c(2010, 2012, ... )
#'   @param genders vector of genders ('male' | 'female' | 'total')
#'   @param regions return selected regions only
#'   @param region.category return selected regions category (for options, see:
#'          unique(SotkanetRegions(type = "table")$region.category)); 
#'	    "ALUEHALLINTOVIRASTO, "ERVA", "EURALUEET", "EUROOPPA", "KUNTA", 
#'	    "MAA", "MAAKUNTA", "NUTS1", "POHJOISMAAT", "SAIRAANHOITOPIIRI", 
#'	    "SEUTUKUNTA", "SUURALUE"   
#'
#' Returns:
#'   @return sotkanet data table
#'
#' @references
#' See citation("sotkanet") 
#' @author Einari Happonen. Maintainer: Louhos/Opasnet \email{louhos@@googlegroups.com}
#' @examples # 
#' @keywords utilities
GetDataSotkanetSingleIndicator <- function (indicator, years = 1990:2000, genders = "total", regions = NULL, region.category = NULL) {

  # FIXME: is it possible to specify already in query which regions we select

  dat <- SotkanetData(indicator = indicator, years = years, genders = genders)

  # Pick corresponding indicator 
  indicator.data <- SotkanetIndicators(indicator)[, c("indicator", "indicator.title.fi")]
  dat <- merge(indicator.data, dat)

  # Pick corresponding region
  #message(paste("Picking region"))
  region.data <- SotkanetRegions()[, c("region", "region.title.fi", "region.code", "region.category")]
  dat <- merge(region.data, dat)

  # Replace comma by point as decimal separator
  #message(paste("Polishing"))
  dat$primary.value <- as.numeric(gsub("\\,", "\\.", as.character(dat$primary.value)))
  dat$absolute.value <- as.numeric(gsub("\\,", "\\.", as.character(dat$absolute.value)))

  # Remove unnecessary column
  #if (all(is.na(dat$absolute.value))) {dat$absolute.value <- NULL}

  # Pick only the selected regions
  if (!is.null(region.category)) {
    dat <- dat[dat$region.category %in% region.category, ]
  }

  if (!is.null(regions)) {
    dat <- dat[dat$region.title.fi %in% regions, ]
  }

  dat

}

