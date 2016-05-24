#' Returns the URI constructed from the parameter settings. This also
#' URI-encodes all the values in each query parameter.
#'
#' @param query.builder Name of the Object of the Query Builder Class
#'
#' @param token Token object containing the OAuth2.0 Authentication details
#'
#'
#'
#' @importFrom utils URLencode
#'
#' @return
#'   A full URI that can be used with the Google Analytics API.
ToBody <- function(query.builder,token) {

  query <- c("end.date"  = query.builder$end.date(),
             "metrics"     = query.builder$metrics(),
             "start.date"    = query.builder$start.date(),
             "title" = query.builder$title(),
             "dimensions"  = query.builder$dimensions(),
             "filters"     = query.builder$filters(),
             "segment"     = query.builder$segment())

  postbody <- "{"

  for (name in names(query)) {
    postbody.name <- switch(name,
                       end.date    = "end-date",
                       metrics     = "metrics",
                       start.date  = "start-date",
                       title    = "title",
                       dimensions  = "dimensions",
                       filters     = "filters",
                       segment     = "segment")


    if (!is.null(postbody.name)) {
      postbody <- paste(postbody,'"',
                   URLencode(postbody.name, reserved = TRUE),'"',
                   ":",'"',
                   URLencode(query[[name]], reserved = FALSE),'"',
                   ",",
                   sep = "",
                   collapse = "")
    }
  }
  # remove the last '&' that joins the query parameters together.
  postbody <- sub(",$", "",postbody)
  # remove any spaces that got added in from bad input.
  postbody <- gsub("\\s", "", postbody)
  postbody <- paste(postbody,"}",sep = "")
  return(postbody)
}
