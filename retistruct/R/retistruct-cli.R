##' This calls \code{\link{retistruct.cli.process}} with a time limit
##' specified by \code{cpu.time.limit}.
##' 
##' @title Process a dataset with a time limit
##' @param dataset Path to dataset to process 
##' @param cpu.time.limit Time limit in seconds
##' @param outputdir Directory in which to save any figures
##' @param device String representing device to print figures to
##' @param ... Other arguments to pass to \code{\link{retistruct.cli.process}}
##' @return A list comprising
##' \item{\code{status}}{0 for success, 1 for reaching
##' \code{cpu.time.limit} and 2 for an unknown error}
##' \item{\code{time}}{The time take in seconds}
##' \item{\code{mess}}{Any error message}
##' @author David Sterratt
##' @export
retistruct.cli <- function(dataset, cpu.time.limit=Inf, outputdir=NA,
                           device="pdf", ...) {
  ## Return code
  status <- 0
  setTimeLimit(cpu=cpu.time.limit)
  syst <- system.time(out <- tryCatch(retistruct.cli.process(dataset,
                                                             outputdir=outputdir, device=device, ...),
                                      error=function(e) {return(e)}))
  mess <- "Success"
  if (inherits(out, "error")) {
    mess <- as.character(out)
    if (grepl("reached CPU time limit", mess)) {
      status <- 1
    } else {
      ## Unknown error
      status <- 2
    }
  }
  ## Success
  return(list(status=status, time=syst["user.self"], mess=mess))
}

##' This function processes a \code{dataset}, saving the
##' reconstruction data and matlab export data to the \code{dataset}
##' directory and printing figures to \code{outputdir}.
##'
##' @title Process a dataset, saving results to disk
##' @param dataset Path to dataset to process 
##' @param outputdir Directory in which to save any figures
##' @param device String representing device to print figures to
##' @param titrate If \code{TRUE} add output of
##' \code{\link{titrate.reconstructedOutline}}  to object saved.
##' @author David Sterratt
##' @export
retistruct.cli.process <- function(dataset, outputdir=NA, device="pdf",
                                   titrate=FALSE) {
  ## Processing
  warn.opt <- getOption("warn")
  options(warn=1)
  r <- retistruct.read.dataset(dataset)
  r <- retistruct.read.markup(r)
  r <- retistruct.reconstruct(r)
  if (titrate) {
    r$titration <- titrate.reconstructedOutline(r)
  }
  ## Output
  retistruct.save.recdata(r)
  
  if (!is.na(outputdir)) {
    message("Producing figures")
    retistruct.cli.figure(dataset, outputdir, device=device)
  }

  ## Export to matlab
  message("Exporting to matlab")
  retistruct.export.matlab(r)
  options(warn=warn.opt)
}

## retistruct.cli.basepath - generate a path based on the elided directory name
##
retistruct.cli.basepath <- function(dataset) {
  basepath <- gsub("\\./", "", dataset)
  basepath <- gsub("/", "_", basepath)
  basepath <- gsub(" ", "_", basepath)
  return(basepath)
}

##' @title Print a figure to file
##' @param dataset Path to dataset to process 
##' @param outputdir Directory in which to save any figures
##' @param device String representing device to print figures to
##' @param width Width of figures in inches
##' @param height Height of figures in inches
##' @param res Resolution of figures in dpi (only applies to bitmap
##' devices)
##' @export
##' @author David Sterratt
retistruct.cli.figure <- function(dataset,
                                  outputdir, device="pdf", width=6, height=6,
                                  res=100) {
  suppressMessages(r <- retistruct.read.recdata(list(dataset=dataset), check=FALSE))
  units <- NULL
  if (device!="pdf") {
    height <- height*res
    width  <- width*res
  }
  suffix <- paste(".", device, sep="")
  dev <- switch(device,
                pdf=pdf,
                png=png,
                jpeg=jpeg,
                tiff=tiff)
  if (is.null(dev)) {
    stop(paste("Device", device, "is not supported"))
  }
  if (!is.null(r)) {
    ## Determine the name of a figure
    basepath <- retistruct.cli.basepath(dataset)
    
    ## Flat plot
    dev(file=file.path(outputdir, paste(basepath, "-flat", suffix, sep="")),
           width=width, height=height)
    par(mar=c(1, 1, 1, 1))
    flatplot(r, axt="n",
              datapoints=TRUE,
              landmarks=TRUE,
              markup=FALSE,
              stitch=TRUE,
              grid=TRUE,
              mesh=FALSE,
              strain=FALSE)
    title(dataset)
    dev.off()

    ## Polar plot with KDE contours
    dev(file=file.path(outputdir, paste(basepath, "-polar-kde", suffix, sep="")),
           width=width, height=height)
    par(mar=c(2, 2, 2, 2))
    projection(r, datapoint.contours=TRUE, grouped.contours=FALSE)
    title(paste("KDE:", dataset))
    if (!is.null(r$EOD)) {
      polartext(paste("OD displacement:", format(r$EOD, digits=3, nsmall=2), "deg"))
    }
    dev.off()

    ## Polar plot with KR contours
    dev(file=file.path(outputdir, paste(basepath, "-polar-kr", suffix, sep="")),
           width=width, height=height)
    par(mar=c(2, 2, 2, 2))
    projection(r, datapoint.contours=FALSE, grouped.contours=TRUE)
    title(paste("KR:", dataset))
    if (!is.null(r$EOD)) {
      polartext(paste("OD displacement:", format(r$EOD, digits=3, nsmall=2), "deg"))
    }
    dev.off()

    
    ## Strain plot
    dev(file=file.path(outputdir, paste(basepath, "-strain", suffix, sep="")),
           width=width, height=height)
    par(mar=c(1, 1, 1, 1))
    flatplot(r, axt="n",
              datapoints=FALSE,
              landmarks=FALSE,
              markup=FALSE,
              stitch=FALSE,
              grid=FALSE,
              mesh=FALSE,
              strain=TRUE)
    title(dataset)
    dev.off()

    ## l.vs.L plot
    dev(file=file.path(outputdir, paste(basepath, "-strain-lvsL", suffix, sep="")),
           width=width, height=height)
    par(mar=c(3.0, 3.0, 1.5, 0.5))
    par(mgp=c(1.5, 0.5, 0))
    par(tcl=-0.3)
    
    lvsLplot(r)
    title(dataset)
    dev.off()
  }
}
