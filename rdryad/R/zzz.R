strextract <- function(str, pattern) regmatches(str, regexpr(pattern, str))

dGET <- function(x, ...) {
  res <- GET(x, ...)
  stop_for_status(res)
  content(res, "text")
}

dr_base_oai <- function() "http://www.datadryad.org/oai/request"
