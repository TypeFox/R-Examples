#' Search CrossRef types
#'
#' @export
#'
#' @param types (character) Type identifier, e.g., journal
#' @template args
#' @template moreargs
#' @template cursor_args
#' @param facet (logical) Include facet results. Default: \code{FALSE}
#' @param works (logical) If TRUE, works returned as well, if not then not.
#' @param parse (logical) Whether to output json \code{FALSE} or parse to 
#' list \code{TRUE}. Default: \code{FALSE}
#' 
#' @details BEWARE: The API will only work for CrossRef DOIs.
#'
#' @references \url{https://github.com/CrossRef/rest-api-doc/blob/master/rest_api.md}
#' @examples \dontrun{
#' cr_types()
#' cr_types("monograph")
#' cr_types("monograph", works=TRUE)
#' cr_types(c('monograph', 'book-set', 'book', 'book-track'))
#' cr_types(c('monograph', 'book-set'), works=TRUE)
#' 
#' ## get facets back
#' cr_types("journal-article", works=TRUE, facet=TRUE)$facets
#' cr_types(c('monograph', 'book-set'), works=TRUE, facet=TRUE)
#' 
#' # Use the cursor for deep paging
#' cr_types("journal-article", works = TRUE, cursor = "*", cursor_max = 500, limit = 100)
#' cr_types(c('monograph', 'book-set'), works = TRUE, cursor = "*",
#'    cursor_max = 300, limit = 100)
#'    
#' # query doesn't work unless using works=TRUE
#' ### you get results, but query param is silently dropped
#' cr_types(query = "ecology")
#' 
#' # print progress - only works when passing more than one type
#' cr_types(c('monograph', 'book-set'), works=TRUE, .progress='text')
#' 
#' # Low level function - does no parsing to data.frame, get json or a list
#' cr_types_('monograph')
#' cr_types_('monograph', parse = TRUE)
#' cr_types_("journal-article", works = TRUE, cursor = "*", 
#'    cursor_max = 300, limit = 100)
#' cr_types_("journal-article", works = TRUE, cursor = "*", 
#'    cursor_max = 300, limit = 100, parse = TRUE)
#' }

`cr_types` <- function(types = NULL, query = NULL, filter = NULL, offset = NULL, limit = NULL, 
  sample = NULL, sort = NULL, order = NULL, facet = FALSE, works = FALSE, 
  cursor = NULL, cursor_max = 5000, .progress="none", ...) {
  
  args <- prep_args(query, filter, offset, limit, sample, sort, order, facet, cursor)
  if (length(types) > 1) {
    res <- llply(types, types_GET, args = args, works = works, 
                 cursor = cursor, cursor_max = cursor_max, ..., .progress = .progress)
    if (!is.null(cursor)) {
      out <- lapply(res, "[[", "data")
      bind_rows(out)
    } else {
      out <- lapply(res, "[[", "message")
      out <- if (works) do.call(c, lapply(out, function(x) lapply(x$items, parse_works))) else lapply(out, DataFrame)
      df <- rbind_all(out)
      meta <- if (works) data.frame(types = types, do.call(rbind, lapply(res, parse_meta)), stringsAsFactors = FALSE) else NULL
      if (facet) { 
        ft <- Map(function(x, y){ 
          rr <- list(parse_facets(x$message$facets)); names(rr) <- y; rr 
        }, res, types) 
      } else {
        ft <- list() 
      }
      list(meta = meta, data = df, facets = ft)
    }
  } else {
    res <- types_GET(types, args, works = works, cursor = cursor, 
                     cursor_max = cursor_max, ...)
    if (!is.null(cursor)) {
      res
    } else {
      if (all(is.na(res))) {
        list(meta = NA, data = NA)
      } else {
        if (is.null(types)) {
          list(meta = list(count = res$message$`total-results`), 
               data = bind_rows(lapply(res$message$items, parse_type)))
        } else {
          if (works) {
            wout <- rbind_all(lapply(res$message$items, parse_works))
            meta <- parse_meta(res)
          } else {
            wout <- parse_type(res$message)
            meta <- NULL
          }
          list(meta = meta, data = wout, facets = parse_facets(res$message$facets))
        }
      }
    }
  }
}

#' @export
#' @rdname cr_types
`cr_types_` <- function(types = NULL, query = NULL, filter = NULL, offset = NULL, limit = NULL, 
  sample = NULL, sort = NULL, order = NULL, facet=FALSE, works = FALSE, 
  cursor = NULL, cursor_max = 5000, .progress="none", parse=FALSE, ...) {
  
  args <- prep_args(query, filter, offset, limit, sample, sort, order, facet, cursor)
  if (length(types) > 1) {
    llply(types, types_GET_, args = args, works = works, 
          cursor = cursor, cursor_max = cursor_max, parse = parse, ..., .progress = .progress)
  } else {
    types_GET_(types, args, works = works, cursor = cursor, 
              cursor_max = cursor_max, parse = parse, ...)
  }
}

types_GET <- function(x, args, works, cursor = NULL, cursor_max = NULL, ...){
  path <- if (!is.null(x)) {
    if (works) sprintf("types/%s/works", x) else sprintf("types/%s", x)
  } else { 
    "types" 
  }
  
  if (!is.null(cursor) && works) {
    rr <- Requestor$new(path = path, args = args, cursor_max = cursor_max,
                        should_parse = TRUE, ...)
    rr$GETcursor()
    rr$parse()
  } else {
    cr_GET(path, args, todf = FALSE, ...)
  }
}

types_GET_ <- function(x, args, works, cursor = NULL, cursor_max = NULL, parse, ...){
  path <- if (!is.null(x)) {
    if (works) sprintf("types/%s/works", x) else sprintf("types/%s", x)
  } else { 
    "types" 
  }
  
  if (!is.null(cursor) && works) {
    rr <- Requestor$new(path = path, args = args, cursor_max = cursor_max,
                        should_parse = parse, ...)
    rr$GETcursor()
    rr$cursor_out
  } else {
    cr_GET(path, args, todf = FALSE, parse = parse, ...)
  }
}

parse_type <- function(x){
  x[vapply(x, length, 1) == 0] <- NA
  data.frame(x, stringsAsFactors = FALSE)
}
