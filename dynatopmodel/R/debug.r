# write a log entry console and/or file
LogEvent <- function(msg, tm=NULL, # time can be actual run time or simulation time
                     warn=F)  #, logout=stderr())  #System.getenv("log"))  # destination
{

     # browser()
		if(is.null(tm)){tm=Sys.time()}
    msg <- paste(as.character(tm), ": ", msg)
    if(warn)
    {
      warning(msg, immediate.=T)
    }
    else
    {
      message(msg)
    }

}

