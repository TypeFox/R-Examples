#' Return number of cpus (or a default on failure)
#' 
#' @param default Number of cores to assume if detectCores fails
#' @return Integer number of cores
#' @author jefferis
#' @export
#' @seealso \code{\link{detectCores}}
#' @return integer number of cores always >=1 for default values
#' @examples 
#' ncpus()
ncpus<-function(default=1L){
  # nb now in base
  cores=parallel::detectCores()
  if(is.na(cores)) {
    cores=parallel::detectCores(TRUE)
    if(is.na(cores)){
      warning("I can't identify the number of cpus. Defaulting to",default)
      return(as.integer(default))
    }
  }
  cores
}
