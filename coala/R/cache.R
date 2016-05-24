cache_env <- new.env()

cache_env[["counter"]] <- 0

get_id <- function() {
  cache_env[["counter"]] <- cache_env[["counter"]] + 1
  as.character(cache_env[["counter"]])
}

cache <- function(model, name, variable) {
  if (is.null(cache_env[[model$id]])) cache_env[[model$id]] <- list()
  cache_env[[model$id]][[name]] <- variable
}

read_cache <- function(model, name) {
  cache_env[[model$id]][[name]]
}

reset_cache <- function() {
  rm(list = ls(envir = cache_env), envir = cache_env)
  cache_env[["counter"]] <- 0
}
