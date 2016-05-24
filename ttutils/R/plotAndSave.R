.saveDevice <- function(plot.func, plot.name, options, ext, ...) {
  fileName <- paste(plot.name, ext, sep=".")
  options <- c(options, list(file=fileName))
  device <- switch(ext,
                   "eps"  = "postscript",
                   "ps"   = "postscript",
                   "pdf"  = "pdf",
                   "jpg"  = "jpeg",
                   "jpeg"  = "jpeg",
                   "png"  = "png",
                   "bmp"  = "bmp",
                   "tiff" = "tiff",
                   "emf"  = "win.metafile",
                   "wmf"  = "win.metafile")
  do.call(device, options)
  plot.func(...)
  dev.off()
}
  
plotAndSave <- function(plot.func, plot.name, ..., folder=getwd(),
                        format=c("eps", "pdf"),
                        options=list(eps = list(onefile=TRUE, horizontal=FALSE,
                                       paper="special", width=7, height=7),
                          ps = list(onefile=TRUE, horizontal=FALSE,
                            paper="special", width=7, height=7),
                          pdf=list(onefile=TRUE)),
                        do.plot=TRUE, do.return=do.plot) { 
  n <- nchar(folder)
  ch <- substr(folder, n, n)
  if (ch != .Platform$file.sep) {
    folder <- paste(folder, .Platform$file.sep, sep="")
  }
  folder <- path.expand(folder)
  if (do.return && !do.plot) {
    warning("Plot will not be displayed, hence there will be no return value!")
  }
  fileName <- paste(folder, plot.name, sep="")
  chosenFormats <- unique(match.arg(tolower(format),
                                    c("eps", "pdf", "ps", "jpg", "png", "bmp",
                                      "tiff", "emf", "wmf"), several.ok=T))
  
  winFormats <- chosenFormats %in% c("emf", "wmf") 
  if ((.Platform$OS.type != "windows") && any(winFormats)) {
    warning("Windows metafiles can be generated only under Windows!")
    chosenFormats <- chosenFormats[!winFormats]
  }
  plotFunction <- match.fun(plot.func)
  lapply(chosenFormats, function(ext)
         .saveDevice(plotFunction, fileName, options[[ext]], ext, ...))
  if (do.plot) {
    if (do.return) {
      return(plot.func(...))
    } else {
      invisible(plot.func(...))
    }
  }
}
