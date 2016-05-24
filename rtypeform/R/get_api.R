get_api = function(api){
  if(is.null(api))
    api = Sys.getenv("typeform_api")
  if(nchar(api) == 0)
    stop("Invalid api key.")
  api
}