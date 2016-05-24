#' Fetch a dataset from Enigma.
#'
#' @export
#'
#' @param dataset Dataset name. Required.
#' @param limit (numeric) Number of rows of the dataset to return. Default (and max): 500
#' @param select (character) Vector of columns to be returned with each row. Default is to return
#' all columns.
#' @param conjunction one of "and" or "or". Only applicable when more than one \code{search}
#' or \code{where} parameter is provided. Default: "and"
#' @param sort (character) Sort rows by a particular column in a given direction. + denotes
#' ascending order, - denotes descending. See examples.
#' @param page (numeric) Paginate row results and return the nth page of results. Pages are
#' calculated based on the current limit, which defaults to 50.
#' @param where (character) Filter results with a SQL-style "where" clause. Only applies to
#' numerical columns - use the \code{search} parameter for strings. Valid operators are >, < and =.
#' Only one \code{where} clause per request is currently supported.
#' @param search (character) Filter results by only returning rows that match a search query. By
#' default this searches the entire table for matching text. To search particular fields only, use
#' the query format "@@fieldname query". To match multiple queries, the | (or) operator can be used
#' eg. "query1|query2".
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
#' # White house visitor list
#' enigma_data(dataset='us.gov.whitehouse.visitor-list')
#'
#' # White house visitor list - selecting three columns
#' enigma_data(dataset='us.gov.whitehouse.visitor-list', select=c('namelast','visitee_namelast',
#' 'last_updatedby'))
#'
#' # White house visitor list - sort by last name
#' enigma_data(dataset='us.gov.whitehouse.visitor-list', sort='+namelast')
#'
#' # White house visitor list - get rows where total_people > 5
#' enigma_data(dataset='us.gov.whitehouse.visitor-list', where='total_people > 5')
#'
#' # White house visitor list - search for Vitale in full name field
#' ## remove the 2nd at symbol before running
#' # enigma_data(dataset='us.gov.whitehouse.visitor-list', search='@@namefull=Vitale')
#'
#' # White house visitor list - search for SOPHIA in first name field
#' ## remove the 2nd at symbol before running
#' # enigma_data(dataset='us.gov.whitehouse.visitor-list', search='@@namefirst=SOPHIA')
#'
#' # Domestic Market Flight Statistics (Final Destination)
#' dataset='us.gov.dot.rita.trans-stats.air-carrier-statistics.t100d-market-all-carrier'
#' enigma_data(dataset=dataset)
#'
#' # Search for 'apple' in the Crunchbase dataset, and get two columns back
#' enigma_data(dataset='com.crunchbase.info.companies.acquisition', search='apple',
#' select=c('acquisition','price_amount'))
#'
#' # conjunction parmeter, compare these two queries
#' enigma_data(dataset = 'us.gov.dol.ogesdw.msha.msha-accident',
#'    where = c('degree_injury_cd > 2', 'no_injuries > 1'), conjunction = "and")
#' enigma_data(dataset = 'us.gov.dol.ogesdw.msha.msha-accident',
#'    where = c('degree_injury_cd > 2', 'no_injuries > 1'), conjunction = "or")
#' }

enigma_data <- function(dataset=NULL, limit=500, select=NULL, conjunction=NULL,
                        sort=NULL, page=NULL, where=NULL, search=NULL, key=NULL, ...) {

  key <- check_key(key)
  check_dataset(dataset)
  if (!is.null(select)) select <- paste(select, collapse = ",")

  url <- sprintf('%s/data/%s/%s', en_base(), key, dataset)
  sw <- proc_search_where(search, where)
  args <- list(limit = limit, select = select, conjunction = conjunction,
                               sort = sort, page = page)
  args <- as.list(unlist(ec(c(sw, args))))
  json <- enigma_GET(url, args, ...)
  meta <- json$info
  json$result <- lapply(json$result, as.list)
  dat2 <- do.call(rbind.fill,
                  lapply(json$result, function(x){
                    x[sapply(x, is.null)] <- NA; data.frame(x, stringsAsFactors = FALSE)
                  }))
  structure(list(success = json$success, datapath = json$datapath, info = meta, result = dat2),
            class = "enigma", dataset = dataset)
}

#' @export
print.enigma <- function(x, ..., n = 10) {
  cat("<<enigma data>>", sep = "\n")
  cat(paste0("  Dataset: ", attr(x, "dataset")), sep = "\n")
  cat(paste0("  Found/returned: ", sprintf("[%s/%s]", x$info$total_results, NROW(x$result))), "\n", sep = "\n")
  trunc_mat(x$result, n = n)
}

proc_search_where <- function(x, y) {
  if (!is.null(x)) {
    if (length(x) > 1) {
      x <- lapply(x, function(z) list(search = z))
    } else {
      x <- list(search = x)
    }
  }
  if (!is.null(y)) {
    if (length(y) > 1) {
      y <- lapply(y, function(z) list(where = z))
    } else {
      y <- list(where = y)
    }
  }
  c(x, y)
}
