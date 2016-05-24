runscript <- function(x, method=c('run', 'copy', 'view', 'show', 'dir'),
                      ask = TRUE, fmt="ch%02d.R", package="FinTS",
                      subdir="scripts", lib.loc=NULL){
##
## 1.  Set up
##
  method <- match.arg(method)
    # chnames <- c(...) # chapter names
    # if (missing(x)) x <- match(select.list(chnames), chnames)
#    s <- system.file("scripts", package = "FinTS")
  s <- system.file(subdir, package = package, lib.loc=lib.loc)
  {
    if(missing(x)){
      Ch0 <- dir(s, full.names=TRUE)
      Ch.info <- file.info(Ch0)
      Chs <- Ch0[!Ch.info$isdir]
      chs <- dir(s)[!Ch.info$isdir]
      ns <- length(chs)
      if(ns<1){
        cat("no files found in directory", s, "\n")
        return()
      }
      firstLine <- chs
      for(i in seq(1, length=ns)){
        fL <- try(readLines(Chs[i], 1))
        if(class(fL) != "try-error")
          firstLine[i] <- paste(chs[i], fL, sep=" - ")
      }
      fL. <- (select.list(firstLine) == firstLine)
      ch <- chs[fL.]
      Ch <- Chs[fL.]
    }
    else {
#    ch <- sprintf("ch%02d.R", x)
      x <- as.numeric(x)
      ch <- sprintf(fmt, x)
      Ch <- paste(s, ch, sep="/")
    }
  }
##
## 2.  method == 'dir'
##
  if(method=='dir')return(Ch)
##
## 3.  method == 'run'
##
  if(method=='run'){
    if(ask){
      op <- par(ask=TRUE)
      on.exit(par(op))
    }
    source(Ch, echo = TRUE)
    return(invisible(Ch))
  }
##
## 4.  method == 'view'
##
  ch. <- readLines(Ch)
#
  if(method=='view'){
    cat(ch., sep="\n")
    return(invisible(Ch))
  }
##
## 4a.  method == 'show'
##
  if(method=='show'){
    file.show(Ch)
    return(invisible(Ch))
  }
##
## 5.  method == 'copy'
##
  writeLines(ch., ch)
  return(invisible(Ch))
}
