# AB: 21/02/2010
# source all the files of a directory
sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[R]$")) {
       if(trace) cat("source:", nm)           
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
 }
