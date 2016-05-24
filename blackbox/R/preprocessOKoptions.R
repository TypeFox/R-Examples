preprocessbboptions <- function(optionList) {
  miscOptions <- optionList$miscOptions
  if("massBarycenter" %in% miscOptions) {## option added (03/2014)
    optionList$barycenterFn = "massBarycenter"
  } else optionList$barycenterFn <- "verticesBarycenter"
  # GRAPHICS
  # graphicPars -> for par() AND for postscript's width and height AND anything else like xmaxticks...
  graphicPars <- optionList$graphicPars
  if (is.null(graphicPars$cex.main)) graphicPars$cex.main <- 1.5 ## larger default font for main title
  if (is.null(graphicPars$cex.lab)) graphicPars$cex.lab <- 1.5 ## larger default font for axis title
  if (is.null(graphicPars$cex.axis)) graphicPars$cex.axis <- 1.4 ## larger default font for tick labels
  if (is.null(graphicPars$labcex)) graphicPars$labcex <- 0.9 ## for contour; default 0.6.
  if (is.null(graphicPars$tcl)) graphicPars$tcl <- NA ## smaller ticks
  if (is.null(graphicPars$mgp)) graphicPars$mgp <- c(3.5, 0.5, 0) ## tick labels closer to axis than default, axis label more distant
  if (is.null(graphicPars$mar)) graphicPars$mar <- c(5, 5.1, 4, 2)+0.1 ## bltr axis label positions on each side; default c(5, 4, 4, 2) + 0.1, starting from 0
  if (is.null(graphicPars$xmaxticks)) graphicPars$xmaxticks <- 9 ## for maketicks -> nice, not for par
  if (is.null(graphicPars$ymaxticks)) graphicPars$ymaxticks <- 10 ## for maketicks -> nice, not for par
  epsArgsNames <- intersect(names(graphicPars), c("height", "width")) ## not even in par() for ps file ...
  optionList$epsArgs <- graphicPars[epsArgsNames] ## only  'eps' arguments
  optionList$graphicPars <- graphicPars
  ##Let us call par() to get names(graphics:::.Pars)
  ParArgsNames <- intersect(names(graphicPars), names(par())) ## using names of 'par' arguments, see ?par
  optionList$ParArgs <- graphicPars[ParArgsNames] ## only  'par' arguments
  ## par() opens a device (not even documented...),
  ##    typically Rplots.pdf in batch session (visible when/if dev.off() is called)
  ## so either we call par() and have an empty Rplots.pdf (overwriting any previous one),
  ## or we would already called something like providePlotFile here (but not using the $options at this point... heavy).
  ## A C++ kind of file control system is missing.
  dev.off() ## close the parasitic device
  unlink("Rplots.pdf") ## removes the parasitic file
  # Set device type and basic filename, open file
  if (is.null(optionList$graphicsFormat)) { ## then gets the default (eps)
    optionList$graphicsFormat <- blackbox.getOption("graphicsFormat")
  } ## else it should be a valid device _function (name)_ such as eps or pdf
  interactiveGraphics <- optionList$interactiveGraphics ## check user's non-default
  if (is.null(interactiveGraphics)) { ## then uses default set by .onLoad
    interactiveGraphics <- blackbox.getOption("interactiveGraphics")
  }
  if ( interactiveGraphics ) {
    optionList["basicRplotsfile"] <- list(NULL) ## probably already is, but ... ?
  } else {
    if (identical(tolower(optionList$graphicsFormat),"postscript")) {
      graphicsExt <- ".ps"
    } else { ## here Migraine's default 'eps' which calls OKsmooth's eps function with 'good' height/width
      graphicsExt <- paste(".", tolower(optionList$graphicsFormat), sep="")
    }
    optionList$basicRplotsfile <- paste("Rplots_", optionList$jobSampleNbr, graphicsExt, sep="");
  }
  #
  if ("NLOPT_LD_MMA_for_CI" %in% optionList$miscOptions) optionList$optimizers <- c(optionList$optimizers,"NLOPT_LD_MMA_for_CI")
  optionList$redundant.mode <- intersect(optionList$miscOptions, c("alwaysRational", "noElim", "alwaysDouble", "defaultPrecision"))
  if (length(optionList$redundant.mode)==0) {
    optionList$redundant.mode <- "defaultPrecision"
  } else {
    optionList$redundant.mode <- optionList$redundant.mode[length(optionList$redundant.mode)] ## get the last user's choice
  }
  #
  optionList["respCols"] <- optionList["fittedLoci"] # if $fittedLoci is NULL, $respCols is filled by buildPointls.
  optionList["fittedLoci"] <- list(NULL) ## explicit element with NULL value => to over-write any preexisting option value
  #
  optionList$GCVrangeFactors <- c(optionList$GCVlowerFactor, optionList$GCVupperFactor) # FR->FR but not yet used !!
  optionList["GCVlowerFactor"] <- list(NULL) ## explicit element with NULL value => to over-write any preexisting option value
  optionList["GCVupperFactor"] <- list(NULL) ## explicit element with NULL value => to over-write any preexisting option value
  # Migraine  defaults (will overwrite OKsmooth defaults)
  if (optionList$covFamily %in% c("Matern", "GneitingMatern")) {
    if (is.null(optionList$initSmoothness)) optionList$initSmoothness <- 3.99 ## but CKrigcoefs uses another value
    if (is.null(optionList$minSmoothness)) optionList$minSmoothness <- 2 ## this one used by CKrigcoefs
    if (is.null(optionList$maxSmoothness)) optionList$maxSmoothness <- 4
  } else {
    if (is.null(optionList$initSmoothness)) optionList$initSmoothness <- 1
    if (is.null(optionList$minSmoothness)) optionList$minSmoothness <- 0.5
    if (is.null(optionList$maxSmoothness)) optionList$maxSmoothness <- 2
  }
  #
  if (optionList$D2IBDbool) {
    optionList$DemographicModel <- c("IBD", "2D")
  } else if (optionList$D1IBDbool) {
    optionList$DemographicModel <- c("IBD", "1D")
    locstring <- paste("Geographic bin width= ", optionList$Nbfactor, sep="")
    print(locstring)
    write(locstring, file=optionList$cleanResu)
  }
  invisible(optionList)
}
