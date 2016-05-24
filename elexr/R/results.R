#' Fetch election results
#'
#' @param election_date Election date as a string in MM-DD-YYYY format.
#' @export
results <- function(election_date) {
  errcode <- system2("which", "elex", stdout = FALSE, stderr = FALSE)
  if (errcode > 0) {
    stop(paste("The elex command does not appear to be installed.  You'll need to ",
               "install it following the instructions at ",
               "http://elex.readthedocs.org/ before you can use this ",
               "function", sep = ""))
  }
  
  api_key <- Sys.getenv("AP_API_KEY")
  if (api_key == "") {
    stop(paste("elex needs an API key specified in the AP_API_KEY ",
               "environment variable.  Set it using Sys.setenv() before ",
               "calling this function", sep = ""))
  }
  
  args <- paste("results", election_date)
  out <- system2("elex",  args = args, stdout = TRUE, stderr = FALSE)
  con <- textConnection(out)
  return(utils::read.csv(con, header = TRUE))
}