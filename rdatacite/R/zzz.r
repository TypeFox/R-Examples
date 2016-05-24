dc_compact <- function(l) Filter(Negate(is.null), l)

dc_base <- function() "http://search.datacite.org/api"

dc_oai_base <- function() "http://oai.datacite.org/oai"

pluck <- function(x, name, type) {
  if (missing(type)) {
    lapply(x, "[[", name)
  } else {
    vapply(x, "[[", name, FUN.VALUE = type)
  }
}

last <- function(x) x[length(x)][[1]]

strextract <- function(str, pattern) regmatches(str, regexpr(pattern, str))

