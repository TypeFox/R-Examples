#' Paging helpers
#'
#' @export
#' @param .data Result of a call to \code{api}
#' @param limit	Maximum number results desired.
#' @param limit_max Maximum number results allowed in each request.
#' @param offset Record to start at
#' @param by Chunk size, if chunking desired. Default:
#' @family dsl
#' @examples
#' quer <- api('https://api.github.com/') %>%
#'   api_path(repos, ropensci, rgbif, issues) %>%
#'   api_query(state = open)
#'
#' # per_page & page, w/ known max_limit
#' api('https://api.github.com/') %>%
#'   api_paging(limit = 220, limit_max = 100)
#'
#' ##### Not working yet
#' # per_page & page
#' # quer %>%
#' #   api_paging(per_page = 10, page = 2)
#'
#' # limit & offset
#' # quer %>%
#' #   api_paging(limit = 10, offset = 20)
#'
#' # rows & start
#' # quer %>%
#' #   api_paging(rows = 10, start = 5)
api_paging <- function(.data, limit, limit_max, offset = 0, by = NULL) {
  .data <- as.req(.data)
  stopifnot(is.numeric(limit), is.numeric(limit_max), is.numeric(offset))
  by <- get_by(by, limit, limit_max)
  args <- list(limit = limit, limit_max = limit_max,
               offset = offset, by = by)
  modifyList(.data, list(paging = args))
}

get_by <- function(by, limit, limit_max) {
  if (!is.null(by)) {
    stopifnot(is.numeric(by))
    stopifnot(by < limit_max)
    return(by)
  } else {
    if (limit > limit_max) {
      return(limit_max)
    } else {
      return(limit)
    }
  }
}

# ## four headers
# x=HEAD("https://api.github.com/repos/ropensci/taxize/issues?state=open&per_page=5&page=4")
# get_links(x$headers)
# ## two headers
# x=HEAD("https://api.github.com/repos/ropensci/taxize/issues")
# get_links(x$headers)
# ## no headers
# x=HEAD("https://api.github.com/repos/ropensci/pangaear/issues")
# get_links(x$headers)
get_links <- function(w) {
  lk <- w$link
  if (is.null(lk)) {
    NULL
  } else {
    if (is(lk, "character")) {
      links <- strtrim(strsplit(lk, ",")[[1]])
      lapply(links, each_link)
    } else {
      nms <- sapply(w, "[[", "name")
      tmp <- unlist(w[nms %in% "next"])
      grep("http", tmp, value = TRUE)
    }
  }
}
# get_links <- function(w) {
#   lk <- w$link
#   urls <- comp(sapply(w, "[[", "url"))
#   if (is.null(lk) && length(urls) == 0) {
#     NULL
#   } else {
#     if (is(w, "character")) {
#       links <- strtrim(strsplit(lk, ",")[[1]])
#       lapply(links, each_link)
#     } else {
#       nms <- sapply(w, "[[", "name")
#       tmp <- unlist(w[nms %in% "next"])
#       grep("http", tmp, value = TRUE)
#     }
#   }
# }

each_link <- function(z) {
  tmp <- strtrim(strsplit(z, ";")[[1]])
  nm <- gsub("\"|(rel)|=", "", tmp[2])
  url <- gsub("^<|>$", "", tmp[1])
  list(name = nm, url = url)
}
