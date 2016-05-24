.cache <- new.env(TRUE)

#' clears the cache
#' 
#' @export
cache_clear <- function(){
  rm( list  = ls(envir=.cache)
    , envir = .cache
    )
}

cache_add <- function(url, result){
  .cache[[url]] <- result
  result
}

cache_get <- function(url){
  .cache[[url]]
}