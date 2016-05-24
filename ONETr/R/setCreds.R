cacheEnv <- new.env()

setCreds <- function(user,pass){
  assign("creds",list(user,pass), envir=cacheEnv)
  message("API credentials saved. You may now use package functions.")
}