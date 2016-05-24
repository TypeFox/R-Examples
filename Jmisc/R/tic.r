#' Start/clock to measure performance. Same as tic and toc in matlab
#' @name tic
#' @aliases tic toc
#' @rdname timer
#' @title Start Stop clock to measure performance
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @param name Name of the temporary time variable
#' @param envir environment of the temporary time variable
#' @export
#' @examples
#' \dontrun{   
#' tic()
#' Sys.sleep(1)
#' toc
#' }
tic <- function(name = ".time_Jmisc", envir = .GlobalEnv){
  assign(name, Sys.time(), envir= envir)
}
