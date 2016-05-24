
# Name:        go(..., add=FALSE,timer=FALSE)
# Description: Like source() but recalls the last source file names by default. Multiple source files can be specified.
# Parameters:  ... = list of filenames as character strings;
#              add = add these names to the current list? if replace, then FALSE
# Note:        does not pass parameters to source()
# Example:     go('myprog')  # will run source('myprog.r')
#              go()      # will run source('myprog.r') again
#              go('somelib',add=TRUE) # will run source('myprog.r') and source('somelib.r')
#              go('myprog','somelib')  # same as above
#              go('mytest') # will run source('mytest') only
#              go()      # runs source('mytest') again
# Reference:   jouni@kerman.com, kerman@stat.columbia.edu
# Modified:    2004-06-22
#

go <- function(..., add=FALSE, timer=FALSE)
{
  last.sources <- getOption(".Last.Source")
  sources <- unlist(list(...))
  if (length(sources)<1) {
    sources <- last.sources
  } else if (add) {
    sources <- c(last.sources,sources)
  }
  if (length(sources)<1) {
    return(cat("Usage: go('sourcefile', 'sourcefile2', ..., add=?, timer=?)\n"))
  }
  options(".Last.Source"=sources)
  cat("Source file(s): ",sources,"\n")
  yy <- NULL
  for (src in sources) {
    if (is.na(src)) {
      next
    }
    if (!file.exists(src)) {
      src2 <- paste(src, ".R", sep="")
      if (file.exists(src2))
        src <- src2
      else {
        cat("source('",src,"') : file does not exist.\n",sep='')
        next
      }
    }
    cat("source('",src,"')\n",sep="")
    if (timer)
      cat("source('",src,"') : ",max(na.omit(system.time(source(src)))),
        " seconds elapsed.\n", sep='')
    else 
      yy[[src]] <- source(src)
  }
  invisible(yy)
}


# By entering "G" on the console, go() is run. This is faster than typing "go()"...
print.GO <- function(x,...) {go()}
G <- structure(NA, class="GO")
#class(G) <- "GO"

# end of go.R
