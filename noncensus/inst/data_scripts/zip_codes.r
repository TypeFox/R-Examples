library(dplyr)
library(zipcode)
data(zipcode)

#' Try/catch with exponential backoff
#'
#' Attempts the expression in \code{expr} up to the number of tries specified in
#' \code{max_attempts}. Each time a failure results, the functions sleeps for a
#' random amount of time before re-attempting the expression. The upper bound of
#' the backoff increases exponentially after each failure.
#'
#' For details on exponential backoff, see:
#' \url{http://en.wikipedia.org/wiki/Exponential_backoff}
#'
#' @param expr an R expression to try.
#' @param silent logical: should the report of error messages be suppressed?
#' @param max_tries the maximum number of times to attempt the expression
#' \code{expr}
#' @param verbose logical: Should detailed messages be reported regarding each
#' attempt? Default: no.
#' @return the value of the expression in \code{expr}. If the final attempt was
#' a failure, the objected returned will be of class try-error".
try_backoff <- function(expr, silent=FALSE, max_attempts=10, verbose=FALSE) {
  for (attempt_i in seq_len(max_attempts)) {
    results <- try(expr=expr, silent=silent)
    if (class(results) == "try-error") {
      backoff <- runif(n=1, min=0, max=2^attempt_i - 1)
      if (verbose) {
        message("Backing off for ", backoff, " seconds.")
      }
      Sys.sleep(backoff)
    } else {
       if (verbose) {
         message("Succeeded after ", attempt_i, " attempts.")
       }
      break
    }
 }
 results
}


# FCC's Census Block Conversions API
# http://www.fcc.gov/developers/census-block-conversions-api
latlong2fips <- function(latitude, longitude) {
  out <- try({
    url <- "http://data.fcc.gov/api/block/find?format=json&latitude=%f&longitude=%f"
    url <- sprintf(url, latitude, longitude)
    json <- try_backoff(RCurl::getURL(url), verbose=TRUE)
    json <- RJSONIO::fromJSON(json)
    as.character(json$County['FIPS'])
  })
  out
}

# Gets unique latitude/longitude pairs
lat_long <- zipcode[c("latitude", "longitude")]
lat_long <- lat_long[!duplicated(lat_long), ]

# Gets FIPS codes for each latitude/longitude pair
lat_long <- within(lat_long, {
  fips <- mapply(latlong2fips, latitude, longitude)
})
lat_long$fips <- factor(lat_long$fips)

# Merge FIPS with zipcode to get mapping from zip to FIPS
zip_codes <- inner_join(zipcode, lat_long)

zip_codes$zip <- factor(zip_codes$zip)
zip_codes$city <- factor(zip_codes$city)
zip_codes$state <- factor(zip_codes$state)

rownames(zip_codes) <- NULL

save(zip_codes, file="../../data/zip_codes.RData", compress='xz')
