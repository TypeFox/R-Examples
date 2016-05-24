#'@include pollstr-package.R
NULL

# Create URL for the charts API method
pollstr_charts_url <- function(page, topic, state, showall) {
  query <- list()
  if (!is.null(page)) {
    query[["page"]] <- as.character(page)[1]
  }  
  if (!is.null(topic)) {
    query[["topic"]] <- as.character(topic)[1]
  }
  if (!is.null(state)) {
    query[["state"]] <- as.character(state)[1]
  }
  if (!is.null(showall)) {
    query[["showall"]] <- if (showall) "true" else "false"
  }
  if (!length(query)) {
    query <- NULL
  }
  modify_url(paste(.POLLSTR_API_URL, "charts", sep = "/"), query = query)
}

# clean up the objects returned by the API
charts2df <- function(.data) {
  charts <- ldply(.data, function(x) {
    x[["estimates"]] <- NULL
    if (is.null(x[["topic"]])) {
      x[["topic"]] <- ""
    }
    x[["election_date"]] <- electiondate2date(x[["election_date"]])
    convert_df(x)
  })
  # Convert
  charts[["last_updated"]] <-
    as.POSIXct(charts[["last_updated"]],
               format = "%Y-%m-%dT%H:%M:%OSZ",
               tz = "GMT")
  
  estimates <- ldply(.data,
                     function(x) {
                       if (length(x[["estimates"]])) {
                         y <- ldply(x[["estimates"]], convert_df)
                         y[["slug"]] <- x[["slug"]]
                         y
                       }
                     })
  structure(list(charts = charts, estimates = estimates),
            class = c("pollstr_charts"))
}

get_charts_page <- function(page, topic, state, showall, as = "parsed") {
  url <- pollstr_charts_url(page, topic, state, showall)
  get_url(url, as = as)
}

#' Get list of available charts
#'
#' @param page Return page number
#' @param state Only include charts from a single state. Use 2-letter state abbreviations. "US" will return all national charts.
#' @param topic Only include charts related to a specific topic. See \url{http://elections.huffingtonpost.com/pollster/api} for examples.
#' @param showall logical Include charts for races that were once possible but didn't happen (e.g. Gingrich vs. Obama 2012)
#' @param convert Rearrange the data returned by the API into easier to use data frames.
#' @param max_pages Maximum number of pages to get.
#' 
#' @references \url{http://elections.huffingtonpost.com/pollster/api}
#' @return If \code{convert=TRUE}, a \code{"pollstr_charts"} object with elements
#' \describe{
#'   \item{\code{charts}}{Data frame with data on charts.}
#'   \item{\code{estimates}}{Data frame with current estimates from each chart. The column \code{slug} matches this data frame to \code{charts}}
#' }
#' Otherwise, a \code{"list"} in the original structure of the json returned by the API.
#' @examples
#' \dontrun{
#'  # Get charts related to Washington
#'  wa <- pollstr_charts(state='WA')
#'  # Get national charts
#'  us_charts <- pollstr_charts(state='US')
#'  # Get charts in the topic '2016-president'
#'  gov <- pollstr_charts(topic='2016-president')
#'  # Get all charts
#'  allcharts <- pollstr_charts()
#' }
#' @export
pollstr_charts <- function(page = 1, topic = NULL, state = NULL, showall = NULL,
                           convert = TRUE, max_pages = 1) {
  
  .data <- list()
  i <- 0L
  while (i < max_pages) {
    newdata <- get_charts_page(page + i, topic, state, showall)
    if (length(newdata)) {
      .data <- append(.data, newdata)
    } else {
      break
    }
    i <- i + 1L
  }
  if (convert) .data <- charts2df(.data)
  .data
}


#' @export
print.pollstr_charts <- function(x, ...) {
  print(x$charts[,c('title','slug','state','poll_count','last_updated')])
  return(invisible(x))
}
