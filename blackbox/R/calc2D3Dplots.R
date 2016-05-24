calc2D3Dplots <- function(plotFile=NULL,pairlist=list()) { ## plain or slice plots of predicted likelihood
  ## NOT PROFILES!!!
  plotobject <- blackbox.getOption("fitobject")
  fittedparamnbr <- blackbox.getOption("fittedparamnbr")
  fittedNames <- blackbox.getOption("fittedNames")
  # (1) determine number of plots
  spec2DProfiles <- blackbox.getOption("spec2DProfiles")
  if (fittedparamnbr<3L) {
    nplots <- 1L ## maybe not quite so, but only test < 10 below
  } else {
    if (length(pairlist)==0L) {
      if (length(spec2DProfiles)>0L) {
        for (locit in seq_len(length(spec2DProfiles))) {
          locstring <- spec2DProfiles[locit]
          locstring <- substr(locstring,2,(nchar(locstring)-1L)) ## removes "(" ")"
          vars <- strsplit(locstring,",")[[1]]
          pairlist <-  c(pairlist, list(c(vars[1], vars[2])))
        }
      } else { ## does not test all2DProfiles (not profiles here)
        varNames <- blackbox.getOption("fittedNames")
        for (kk in 1:fittedparamnbr) for (ll in 1:fittedparamnbr) {
          if (kk<ll) pairlist <- c(pairlist, list(c(varNames[kk], varNames[ll])))
        }
      }
    }
    nplots <- length(pairlist)
  }
  # (2) determine plotFile
  if (is.null(plotFile)) {
    plotFile <- ! (blackbox.getOption("interactiveGraphics") && nplots<10) ## no screen output if too many plots
    if (plotFile) plotFile <- blackbox.getOption("basicRplotsfile")
  }
  # (3) actual plots
  if (nplots>0L) {
    if ( is.character(plotFile)) providePlotFile(plotFile) ## exclude NULL,FALSE...
    if (fittedparamnbr==1L) {
      message.redef("*** A plot of the predicted likelihood is produced. ***")
      rawProfiles(plottype="prediction") ## NOT a 'raw' Profile !
    } else if (fittedparamnbr==2) { ## plotSlice really not appropriate for this case
      message.redef("*** Plots of the predicted likelihood are produced. ***")
      xvar <- fittedNames[1];yvar <- fittedNames[2]
      ## 3-D (persp) and 2-D (contour) plots of the likelihood surface
      makeplottypes(xvar, yvar, types=c("C", "p"),
                    main="",
                    grid.list= list( v1=gridfn(varname=xvar, varnameS=fittedNames), v2=gridfn(yvar, varnameS=fittedNames)))
      if ("BWPlots" %innc% blackbox.getOption("plotOptions")) {
        ## Black & White 2D contour plot
        makeplottypes(xvar, yvar, types=c("C", "p"),
                      main="",
                      grid.list= list( v1=gridfn(varname=xvar, varnameS=fittedNames), v2=gridfn(yvar, varnameS=fittedNames)),
                      bw=TRUE) ## bw is surfaceKrig argument
      }
    } else if (fittedparamnbr>2) { #FR->FR BW plots missing here
      message.redef("*** 'Slice plots' of the predicted likelihood are produced. ***")
      ##    par(mfrow=c(1, 1))
      for (pairit in seq_len(nplots)) {
        varParams <- pairlist[[pairit]]
        FixedParams <- (fittedNames %w/o% varParams)
        FixedValue <- blackbox.getOption("rosglobal")$par[FixedParams]
        plotSlice(varParams[1], varParams[2], as.list(FixedValue))
      }
    } ## endif fittedparamnbr==2 else...
  }
  invisible(NULL)
}
