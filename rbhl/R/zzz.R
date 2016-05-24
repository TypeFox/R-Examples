query2message <- function(url, x) {
  mylist <- list()
  for (i in 1:length(x)) {
    mylist[i] <- paste(names(x[i]), '=', x[i][[1]], sep = '')
  }
  paste('API call:\n', paste(url, "?", paste(mylist, collapse = "&"), sep = ''))
}

return_results <- function(x, y, z) {
  if (y == 'raw') {
    return( x )
  } else if (y == 'list') {
    return( jsonlite::fromJSON(I(x)) )
  } else {
    if (z == "json") {
      return(jsonlite::fromJSON(I(x)))
    } else {
      return(XML::xmlTreeParse(I(x)))
    }
  }
}

check_key <- function(x) {
  tmp <- if (is.null(x)) Sys.getenv("BHL_KEY", "") else x
  if (tmp == "") getOption("bhl_key", stop("need an API key for BHL")) else tmp
}

bhl_GET <- function(as, args, ...){
  out <- GET(bhl_url(), query = args, ...)
  stop_for_status(out)
  res <- switch(as,
         xml = xmlSize(xpathSApply(content(out), "//Result")[[1]]),
         json = length(content(out)$Result),
         list = length(content(out)$Result),
         table = length(content(out)$Result))
  if (is.null(res) || res == 0) {
    stop("No results found", call. = FALSE)
  }
  tt <- content(out, as = "text")
  switch(as, json = tt, xml = tt, list = fjson(tt), table = todf(tt))
}

todf <- function(x){
  temp <- jsonlite::fromJSON(I(x), TRUE)$Result
  if (is.character(temp)) {
    temp
  } else {
    tmp <- if (!is.null(names(temp))) {
      data.frame(bhlc(temp), stringsAsFactors = FALSE)
    } else {
      do.call(rbind.fill, lapply(bhlc(temp), data.frame))
    }
    structure(list(data = tmp), class = "bhldf")
  }
}

fjson <- function(x) jsonlite::fromJSON(x, FALSE)

bhl_url <- function() "http://www.biodiversitylibrary.org/api2/httpquery.ashx"

bhlc <- function(l) Filter(Negate(is.null), l)

as_f <- function(x) {
  as <- match.arg(x, c("table","list","json","xml"))
  if (as %in% c('list','table','json')) 'json' else 'xml'
}

#' @export
print.bhldf <- function(x, ..., n = 10) {
  cat(sprintf("<bhl data> [%d, %d]", NROW(x$data), NCOL(x$data)), sep = "\n")
  trunc_mat(x$data, n = n)
}
