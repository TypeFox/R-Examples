saveOldFile <- function(filename) {
  preexists <- file.info(filename)[1,1]
  if ( ! is.na(preexists)) { ## if file preexists on disk, we save it under another name
    # first generate the name
    namesplit <- strsplit(filename,split=".",fixed=TRUE)[[1]] ## e.g. "Rplots_1" "eps"
    len <- length(namesplit)
    begname <- paste(paste(namesplit[-len],collapse="."),".old_",sep="") ## "Rplots_1.old_"
    endname <- paste(".",namesplit[len],sep="") ## ".eps"
    copyname <- generateFileName(begname,endname) ## spaMM::generateFileName => "Rplots_1.old_<#>.eps"
    #
    unlink(copyname)
    success <- file.copy(filename,copyname)
    if (success) {
      return(copyname)
    } else return(FALSE)
  } else return("")
} ## returns -1 if no old file, FALSE is failed to copy an old file, TRUE if successfully copied such a file

## Ouvertures de fichiers graphiques doivent tous passer par cette function: except one call in preprocess
## Drawback of this tracking mechanism: users should not explicitly close (by dev.off or graphics.off) any graphic file
## because $plotFiles will not be updated
providePlotFile <- function(filename, verbose=FALSE) { ## to open, keep track of, and reopen files
  if(verbose) message.redef(paste("providePlotFile(", filename, ") called."))
  plotFiles <- blackbox.getOption("plotFiles")
  newFile <- is.null(plotFiles[[filename]])
  if ( ! newFile ) { ## file already created in session
    dev.set(plotFiles[[filename]]) ## set it as current output file
    if (verbose) message.redef(paste("dev.set(...) called for", filename, "."))
  } else { ## create new file
    if (length(dev.list())==62) {
      message.redef("(!) Maximum number of graphic devices reached. Closing the first in the list...")
      dev.off(dev.list()[1])
    }
    abyss <- saveOldFile(filename)
    eval(call(blackbox.getOption("graphicsFormat"), file=filename))
    .blackbox.data$options$plotFiles[[filename]] <- dev.cur() ## .blackbox required here
    ## calling fn provideDevice will set par() argsfor the device
    if (verbose) message.redef(paste(filename, " created."))  ## 'zut[[filename]]' ou simplement le filename ??
  }
  invisible(newFile)
}
