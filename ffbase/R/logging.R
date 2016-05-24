emptyLogger <- function(...) invisible()

Log <- new.env()
Log$info <- if (interactive()) cat else emptyLogger

Log$chunk <- function(i){
  if (is.na(i[3])){
    Log$info("\r< Processing chunk:",i," >")    
  } else {
    if (i[1]==1) Log$info("\n")
    Log$info("\r< Processing :",round(100*(i[2])/i[3]), "% >" , sep="")
    if (i[2] == i[3]){
      Log$info("\r")
    }
  } 
}

#' sets the logging of ffbase
#' @param level logging level: info/debug
#' @param logger function to be called for logging statements, by default this is \code{cat}
#' @export
set_ffbase_logging <- function(level = c("info"), logger=if (interactive()) cat){
  
  if (isTRUE(logger)){
    logger <- cat
  } else if (!is.function(logger)){
    message("switching logging of. Logging can be restored using
            'set_ffbase_logging(logger=TRUE)'")
    logger <- emptyLogger
  }
  
  assign(level, logger, Log)
}
