
#' @encoding UTF-8
#' @title Pause
#' @description A replication of MatLab pause function.
#' @param x is optional. If x>0 a call is made to \code{\link{Sys.sleep}}. Else, execution pauses until a key is entered.
#' @export
`pause` <-
  function (x=0) {
    if(x > 0){
      Sys.sleep(x)
    }else{
      cat("Hit <enter> to continue...","green")
      readline()
      invisible()
    }
  }
NULL
