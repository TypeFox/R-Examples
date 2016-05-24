#' @importFrom jsonlite fromJSON
#' @importFrom utils URLencode
resolve_resource <- function(url, ..., cache=TRUE){
  url <- utils::URLencode(url)
  
  if (isTRUE(cache)){
    res <- cache_get(url)
    if (!is.null(res)){
      return(res)
    }
  }
  
  message(...," ", url)
  res <- jsonlite::fromJSON(url)$value
  
  if (isTRUE(cache)){
    cache_add(url, res)
  }
  res
}
