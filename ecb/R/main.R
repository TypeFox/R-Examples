#' Retrieve data frame of all datasets in the ECB Statistical Data Warehouse
#'
#' @return A dataframe
#' @export
#'
#' @examples
#' df <- get_dataflows()
#' head(df)
get_dataflows <- function() {

  query_url <- "https://sdw-wsrest.ecb.europa.eu/service/dataflow"

  req <- make_request(query_url, "metadata")

  res <- xml2::read_xml(httr::content(req, "text"), verbose = TRUE)

  ecb_ns <- xml2::xml_ns(res) # xml namespace

  data_flows_nodes <- xml2::xml_find_all(res, "//str:Dataflow", ecb_ns)
  name_nodes <- xml2::xml_find_all(xml2::xml_children(data_flows_nodes),
                                   "//com:Name", ecb_ns)

  flow_ref <- xml2::xml_attr(data_flows_nodes, "id")
  title <- xml2::xml_text(name_nodes)

  df <- data.frame(flow_ref, title, stringsAsFactors = FALSE)
  structure(df, class = c("tbl_df", "tbl", "data.frame"))
}

#' Retrieve data from the ECB Statistical Data Warehouse API
#'
#' @param key A character string identifying the series to be retrieved
#' @param filter A named list with additional filters (see \code{details})
#'
#' @details
#' The \code{filter} option of \code{get_data()} takes a named list of key-value pairs.
#' If left blank, it returns all data for the current version.
#'
#' Available filter parameters:
#'
#' \itemize{
#' \item \code{startPeriod} & \code{endPeriod}
#'  \itemize{
#'    \item \code{YYYY} for annual data (e.g.: 2013)
#'    \item \code{YYYY-S[1-2]} for semi-annual data (e.g.: 2013-S1)
#'    \item \code{YYYY-Q[1-4]} for quarterly data (e.g.: 2013-Q1)
#'    \item \code{YYYY-MM} for monthly data (e.g.: 2013-01)
#'    \item \code{YYYY-W[01-53]} for weekly data (e.g.: 2013-W01)
#'    \item \code{YYYY-MM-DD} for daily data (e.g.: 2013-01-01)
#'    }
#' \item \code{updatedAfter}
#'  \itemize{
#'    \item A timestamp to retrieve the latest version of changed values in the database since a certain point in time
#'    \item Example: \code{filter = list(updatedAfter = 2009-05-15T14:15:00+01:00)}
#'    }
#' \item \code{firstNObservations} & \code{lastNObservations}
#'  \itemize{
#'    \item Example: \code{filter = list(firstNObservations = 12)} retrieves the first 12 observations of all specified series
#'    }
#' \item \code{detail}
#'  \itemize{
#'    \item Possible options: \code{full/dataonly/serieskeysonly/nodata}
#'    \item \code{dataonly} is the default
#'    \item Use \code{serieskeysonly} or \code{nodata} to list series that match a certain query, without returning the actual data
#'    \item An alternative to using \code{serieskeys/nodata} is the convenience function \code{get_dimensions()}, which returns a list of dataframes with dimensions and explanations (see extended example below).
#'    \item \code{full} returns both the series values and all metadata. This entails retrieving much more data than with the `dataonly` option.
#'    }
#' \item \code{includeHistory} (not currently implemented)
#'  \itemize{
#'    \item \code{false} (default) returns only version currently in production
#'    \item \code{true} returns version currently in production, as well as all previous versions
#'    }
#' }
#' See the \href{https://sdw-wsrest.ecb.europa.eu/}{SDW API} for more details.
#'
#' @return A data frame
#' @export
#'
#' @examples
#' # Get monthly data on annualized euro area headline HICP
#' hicp <- get_data("ICP.M.U2.N.000000.4.ANR")
#' head(hicp)
get_data <- function(key, filter = NULL) {

  if(!"detail" %in% names(filter)) {
    filter <- c(filter, "detail" = "dataonly")
  }

  if(!filter[["detail"]] %in% c("full", "dataonly")) {
    return(get_dimensions(key))
  }

  query_url <- create_query_url(key, filter = filter)

  req <- make_request(query_url, "data")

  tmp <- tempfile()
  writeLines(httr::content(req, "text", encoding = "utf-8"), tmp)

  result <- rsdmx::readSDMX(tmp, FALSE)

  unlink(tmp)

  df <- as.data.frame(result)
  df <- structure(df,
                  class = c("tbl_df", "tbl", "data.frame"),
                  names = tolower(names(df)))
  df
}

#' Retrieve dimensions of series in the ECB's SDW
#'
#' @param key A character string identifying the series to be retrieved
#'
#' @return A list of data frames, one for each series retrieved
#' @export
#'
#' @examples
#' hicp_dims <- get_dimensions("ICP.M.U2.N.000000.4.ANR")
#' hicp_dims[[1]]
get_dimensions <- function(key) {

  query_url <- create_query_url(key, filter = list("detail" = "nodata"))

  # Used in creating names (series_names) below
  flow_ref <- regmatches(key, regexpr("^[[:alnum:]]+", key))

  req <- make_request(query_url, "metadata")

  skeys <- xml2::read_xml(httr::content(req, "text", encoding = "utf-8"),
                          verbose = TRUE)

  skeys_ns <- xml2::xml_ns(skeys) # xml namespace

  series <- xml2::xml_find_all(skeys, "//generic:Series", skeys_ns)

  series_list <- lapply(series, xml2::xml_children)

  # Concatenate dimensions to recreate series code
  series_names <- vapply(series_list, function(x) {

      attrs <- xml2::xml_attr(xml2::xml_children(x[1]), "value")
      name <- paste0(attrs, collapse = ".")
      paste(flow_ref, name, sep = ".")

    }, character(1))

  # Return list of dataframes, one for each series, with dimension-value pairs
  df_dim <- lapply(series_list, function(nodeset) {

    data.frame(dim = xml2::xml_attr(xml2::xml_children(nodeset), "id"),
               value = xml2::xml_attr(xml2::xml_children(nodeset), "value"),
               stringsAsFactors = FALSE)
  })

  names(df_dim) <- series_names
  df_dim
}

#' Get full, human-readable description of a series
#'
#' @param key A character string identifying the series to be retrieved
#'
#' @return A character vector
#' @export
#'
#' @examples
#' get_description("ICP.M.DE.N.000000+XEF000.4.ANR")
get_description <- function(key) {
  vapply(get_dimensions(key), function(x) x$value[x$dim == "TITLE_COMPL"],
         character(1))
}

#' Format date variable retrieved from the SDW into a proper date variable
#'
#' @param x A vector of dates
#'
#' @return A date-formatted vector
#' @export
#'
#' @examples
#' hicp <- get_data("ICP.M.U2.N.000000.4.ANR")
#' hicp$obstime <- convert_dates(hicp$obstime)
#' str(hicp)
convert_dates <- function(x) {

  if(grepl("^[0-9]{4}$", x[1])) {
    return(as.Date(paste0(x, "-01-01"), "%Y-%m-%d"))
  }

  if(grepl("^[0-9]{4}-[0-9]{2}$", x[1])) {
    return(as.Date(paste0(x, "-01"), "%Y-%m-%d"))
  }

  if(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x[1])) {
    return(as.Date(x, "%Y-%m-%d"))
  }
}

create_query_url <- function(key, filter = NULL) {

  url <- "https://sdw-wsrest.ecb.europa.eu/service/data"

  # Get flow reference (= dataset abbreviation, e.g. ICP or BOP)
  flow_ref <- regmatches(key, regexpr("^[[:alnum:]]+", key))
  key_q <- regmatches(key, regexpr("^[[:alnum:]]+\\.", key),
                      invert = TRUE)[[1]][2]

  if(any(names(filter) == "")) {
    stop("All filter parameters must be named!")
  }

  if("updatedAfter" %in% names(filter)) {
    filter$updatedAfter <- curl::curl_escape(filter$updatedAfter)
  }

  # Create parameter part of query string
  names <- curl::curl_escape(names(filter))
  values <- curl::curl_escape(as.character(filter))
  query <- paste0(names, "=", values, collapse = "&")
  query <- paste0("?", query)

  query_url <- paste(url, flow_ref, key_q, query, sep = "/")
  query_url
}

check_status <- function(req) {
  if(req$status_code >= 400)
    stop("HTTP failure: ", req$status_code, "\n", httr::content(req, "text"))
}

make_request <- function(query_url, header_type) {

  accept_headers <-
    c("metadata" = "application/vnd.sdmx.genericdata+xml;version=2.1",
      "data" = "application/vnd.sdmx.structurespecificdata+xml;version=2.1")

  req <- httr::GET(query_url, httr::add_headers(
    "Accept" = accept_headers[header_type],
    "Accept-Encoding" = "gzip, deflate"))

  check_status(req)
  req
}
