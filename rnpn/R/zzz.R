base <- function() 'https://www.usanpn.org/npn_portal/'

npnc <- function(l) Filter(Negate(is.null), l)

pop <- function(x, y) {
  x[!names(x) %in% y]
}

ldfply <- function(y){
  res <- lapply(y, function(x){
    x[ sapply(x, is.null) ] <- NA
    data.frame(x, stringsAsFactors = FALSE)
  })
  do.call(rbind.fill, res)
}

npn_GET <- function(url, args, ...){
  tmp <- GET(url, query = args, ...)
  stop_for_status(tmp)
  tt <- content(tmp, as = "text", encoding = "UTF-8")
  jsonlite::fromJSON(tt, FALSE)
}
