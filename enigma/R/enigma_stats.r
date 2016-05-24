#' Get statistics on columns of a dataset from Enigma.
#'
#' @importFrom httr GET content write_disk parse_url
#' @importFrom jsonlite fromJSON
#' @importFrom plyr rbind.fill
#'
#' @export
#'
#' @param dataset Dataset name. Required.
#' @param select (character) Column to get statistics on. Required.
#' @param conjunction one of "and" or "or". Only applicable when more than one \code{search}
#' or \code{where} parameter is provided. Default: "and"
#' @param operation (character) Operation to run on a given column. For a numerical column, valid operations
#' are sum, avg, stddev, variance, max, min and frequency. For a date column, valid operations are
#' max, min and frequency. For all other columns, the only valid operation is frequency. Defaults
#' to all available operations based on the column's type.
#' @param by (character) Compound operation to run on a given pair of columns. Valid compound
#' operations are sum and avg. When running a compound operation query, the \code{of} parameter is
#' required (see below).
#' @param of (character) Numerical column to compare against when running a compound operation.
#' Required when using the \code{by} parameter. Must be a numerical column.
#' @param limit (numeric) Limit the number of frequency, compound sum, or compound average results
#' returned. Max: 500; Default: 500.
#' @param search (character) Filter results by only returning rows that match a search query. By
#' default this searches the entire table for matching text. To search particular fields only, use
#' the query format "@@fieldname query". To match multiple queries, the | (or) operator can be used
#' eg. "query1|query2".
#' @param where (character) Filter results with a SQL-style "where" clause. Only applies to
#' numerical columns - use the \code{search} parameter for strings. Valid operators are >, < and =.
#' Only one \code{where} clause per request is currently supported.
#' @param sort (character) Sort frequency, compound sum, or compound average results in a given
#' direction. + denotes ascending order, - denotes descending
#' @param page (numeric) Paginate frequency, compound sum, or compound average results and return
#' the nth page of results. Pages are calculated based on the current limit, which defaults to 500.
#' @param key (character) Required. An Enigma API key. Supply in the function call, or store in
#' your \code{.Rprofile} file, or do \code{options(enigmaKey = "<your key>")}. Obtain an API key
#' by creating an account with Enigma at \url{http://enigma.io}, then obtain an API key from
#' your account page.
#' @param ... Named options passed on to \code{\link[httr]{GET}}
#' @examples \dontrun{
#' # After obtaining an API key from Enigma's website, pass in your key to the function call
#' # or set in your options (see above instructions for the key parameter)
#' # If you pass in your key to the function call use the key parameter
#'
#' # stats on a varchar column
#' cbase <- 'com.crunchbase.info.companies.acquisition'
#' enigma_stats(dataset=cbase, select='acquired_month')
#'
#' # stats on a numeric column
#' enigma_stats(dataset=cbase, select='price_amount')
#'
#' # stats on a date column
#' pakistan <- 'gov.pk.secp.business-registry.all-entities'
#' enigma_metadata(dataset=pakistan)
#' enigma_stats(dataset=pakistan, select='registration_date')
#'
#' # stats on a date column, by the average of a numeric column
#' aust <- 'gov.au.government-spending.federal-contracts'
#' enigma_metadata(dataset=aust)
#' enigma_stats(dataset=aust, select='contractstart', by='avg', of='value')
#'
#' # Get frequency of distances traveled, and plot
#' ## get columns for the air carrier dataset
#' dset <- 'us.gov.dot.rita.trans-stats.air-carrier-statistics.t100d-market-all-carrier'
#' enigma_metadata(dset)$columns$table[,c(1:4)]
#' out <- enigma_stats(dset, select='distance')
#' library("ggplot2")
#' library("ggthemes")
#' df <- out$result$frequency
#' df <- data.frame(distance=as.numeric(df$distance), count=as.numeric(df$count))
#' ggplot(df, aes(distance, count)) +
#'  geom_bar(stat="identity") +
#'  geom_point() +
#'  theme_grey(base_size = 18) +
#'  labs(y="flights", x="distance (miles)")
#'  
#' # conjunction parmeter, compare these two queries
#' enigma_stats(dataset = 'us.gov.dol.ogesdw.msha.msha-accident', select = "no_injuries",
#'    where = c('degree_injury_cd > 2', 'no_injuries > 1'), conjunction = "and")
#' enigma_stats(dataset = 'us.gov.dol.ogesdw.msha.msha-accident', select = "no_injuries",
#'    where = c('degree_injury_cd > 2', 'no_injuries > 1'), conjunction = "or")
#' }

enigma_stats <- function(dataset=NULL, select, conjunction = NULL, operation=NULL, 
  by=NULL, of=NULL, limit=500, search=NULL, where=NULL, sort=NULL, page=NULL, key=NULL, ...) {

  key <- check_key(key)
  check_dataset(dataset)

  url <- sprintf('%s/stats/%s/%s/select/%s', en_base(), key, dataset, select)
  sw <- proc_search_where(search, where)
  args <- list(operation = operation, conjunction = conjunction, by = by, of = of, 
               limit = limit, sort = sort, page = page)
  args <- as.list(unlist(ec(c(sw, args))))
  json <- enigma_GET(url, args, ...)
  sum_stats <- enigma_stats_dat_parser(json)
  structure(list(success = json$success, datapath = json$datapath, info = json$info, result = sum_stats), class = "enigma_stats")
}

enigma_stats_dat_parser <- function(x) {
  nn <- names(x$result)
  res <- lapply(nn, function(z){
    tmp <- x$result[[z]]
    if (length(tmp) > 1) {
      do.call(rbind.fill, lapply(tmp, function(w){
        b <- as.list(w)
        b[sapply(b, is.null)] <- "null"
        data.frame(b, stringsAsFactors = FALSE)
      }))
    } else { 
      tmp 
    }
  })
  names(res) <- nn
  res
}
