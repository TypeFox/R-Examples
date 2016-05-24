# Copyright 2014 by Dmitry Pavlyuk <Dmitry.V.Pavlyuk@gmail.com>

#
# Routines for logging
#

# Available logging levels
.loggingLevels <- c("quiet","warn","info","debug")

#
# Logs a message, an (optional) object, and caller function's name to console
#
logging <- function(message, obj = NULL, level = "debug", caller = sys.call(-1)){
        lvl <- envirGet("logging")
        if (!is.null(lvl) 
                && ((which(.loggingLevels == lvl) >=    which(.loggingLevels == level)))){
                message(toupper(level)," ", toString(caller[[1]])," (pid=",Sys.getpid(),") ", message)
                if (!is.null(obj)){
                        print(obj)
                }
                flush.console()
        }
}