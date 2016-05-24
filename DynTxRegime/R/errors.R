UserError <- function(reason, msg){
  if( tolower(reason) == 'input' ) {
    e <- simpleError(msg)
    stop(e)
  }
}

DeveloperError <- function(msg, place){

  cat("You have discovered a bug in DynTxRegime.\n")
  cat("Please forward the following information to sthollow@ncsu.edu.\n")
  cat(msg, " received from ", place, "\n")
  e <- simpleError("End Error Report")
  stop(e)

}
