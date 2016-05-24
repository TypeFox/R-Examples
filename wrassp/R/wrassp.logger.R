##' Designated logger for the wrassp signal processing functions
##'
##' Function logs the call to a signal processing function (spf) of wrassp. 
##' It is called by default if the forceToLog option of the spf is not set to 
##' false. I tries to format the output in an easily readable fashion. 
##' 
##' @title wrassp.logger
##' @param fName the name of the function calling the logger
##' @param fOpts are the function options given by the user aquired by match.call 
##' @param optLogFilePath path to option log file 
##' @param listOfFiles vector of file paths that the spf calling the logger processed
##' @author Raphael Winkelmann
##' @seealso \code{\link{match.call}}
wrassp.logger <- function(fName, fOpts,
                          optLogFilePath, listOfFiles){
  
  fid = file(optLogFilePath, open="a")

  cat("\n##################################\n", file = fid, append = T)
  cat("##################################\n", file = fid, append = T)
  cat(paste("########", fName, "performed ########\n"), file = fid, append = T)
  
  cat("Timestamp: ", paste(Sys.time()), '\n', file = fid, append = T)
  for (opt in names(fOpts)){
    if(opt != "listOfFiles" && opt != "optLogFilePath"){
      cat(paste(opt, ":", fOpts[opt][[1]], "\n"), file = fid, append = T)      
    }
  }
  cat(" => on files:\n\t", file = fid, append = T)
  cat(paste(listOfFiles, collapse="\n\t"), file = fid, append = T)

  close(fid)
}
