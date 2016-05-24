#' Description:
#' SotkanetData retrieves Sotkanet data 
#' according to the query arguments.
#'
#' Arguments:
#'   @param indicator Dataset identifier
#'   @param years vector of years c(2010, 2012, ... )
#'   @param genders vector of genders ('male' | 'female' | 'total')
#'
#' Returns:
#'   @return sotkanet json query
#'
#' @references
#' See citation("sotkanet") 
#' @author Einari Happonen / Opasnet / Louhos. Maintainer: Louhos/Opasnet \email{louhos@@googlegroups.com}
#' @examples # 
#' @keywords utilities
SotkanetData <- function(indicator, years, genders)
{

  #base.url <- base_url()
  #url <- paste(base.url, 'data/csv?', sep = "")
  # Here the older url is in use for some reason:
  url <- 'http://www.sotkanet.fi/rest/1.0/data/csv?'
  #url <- "https://www.sotkanet.fi/sotkanet/fi/taulukko/?"

  url <- paste(url, 'indicator=',indicator, '&years=', 
      	 	    paste(years, collapse='&years='), 
		    '&genders=', paste(genders, 
		    collapse='&genders='), sep='')

  #res <- sotkanet.json_query(url)
  res <- sotkanet.csv_query(url)

  return(res)

}

