## possible usage
# rawProfiles(plottype="showIters",unselectedKgpt=pointls,iterThresholds = 100)
# plots in red the points >= 100th. No FONKgpointls argument.
rawProfiles <- function(unselectedKgpt, ## not needed if plottype="prediction"
                      FONKgpointls, ## not needed if plottype="prediction" or "showIters"
                      plottype,
                      FONKgLow= blackbox.getOption("FONKgLow"),
                      FONKgUp= blackbox.getOption("FONKgUp"),
                      style="compact",
                      iterThresholds=NA,
                      addPrediction=F) {
  ## addPrediction currently works only for OnePop (one var) case
  ## usage: rawProfiles("both") or rawProfiles("unselected") or rawProfiles("showIters", "", 531)
  ## very raw profiling version (except if L is computed without error)
  ## Plots of the points values computed by Migraine before kriging
  ## X scale may not be exact for the selected pojnts with high likelihood
  ## default style 'compact' attemprs to put all panels in a single plot
  fittedLoci <- blackbox.getOption("respCols")
  fittedNames <- blackbox.getOption("fittedNames")
  fittedparamnbr <- blackbox.getOption("fittedparamnbr")
  GCVrownames <- blackbox.getOption("GCVrownames")
  axesBool <- FALSE
  provideDevice(bbdefaultPars=TRUE) ## provides a file or not depending on its default file argument
  intsqrt <- as.integer(sqrt(fittedparamnbr))
  op <- par(no.readonly = TRUE) # the whole list of settable par's to be reset afterwards
  if (style=="compact") {
    if (intsqrt>1) {loccex.axis <- blackbox.getOption("graphicPars")$cex.axis*0.6} else {loccex.axis <- blackbox.getOption("graphicPars")$cex.axis}
    par(mfrow=c(ceiling(fittedparamnbr/intsqrt), intsqrt), cex.axis=loccex.axis)
  }
  if (length(fittedLoci)==0) {locfittedLoci <- blackbox.getOption("ycolname")} else {locfittedLoci <- fittedLoci}
  needFONKgpointls <- ! all( (plottype %in%  c("prediction","showIters" )))
  if ( ! all( (plottype ==  "prediction" )) ) { 
    unsyvalues <- - apply(unselectedKgpt[, locfittedLoci, drop=FALSE], 1, sum)
    if ( needFONKgpointls) {
      syvalues <- - apply(FONKgpointls[, locfittedLoci, drop=FALSE], 1, sum)
      gcvKgpointls <- FONKgpointls[GCVrownames, ]
      gcvyvalues <- - apply(gcvKgpointls[, locfittedLoci, drop=FALSE], 1, sum)
    }
  }
  for(st in fittedNames) {
    if ( ! all( (plottype ==  "prediction" ))) { ## rajout test 01/2015
      rown <- rownames(unselectedKgpt)
      ## everything goes in canonical scale; the log plot is controlled through the log=logv argument
      if (islogscale(st)) {
        unsxx <- exp(unselectedKgpt[, st])
        if (needFONKgpointls) {
          xx <- exp(FONKgpointls[, st])
          gcvxx <- exp(gcvKgpointls[, st])
        }
        locstring <- " on a log scale"
        logv <- "x"
      } else {
        unsxx <- unselectedKgpt[, st]
        if (needFONKgpointls) {
          xx <- FONKgpointls[, st]
          gcvxx <- gcvKgpointls[, st]
        }
        locstring <- ""
        logv <- ""
      }
      if (st =="latt2Ns2") {
        Nbfactor <- blackbox.getOption("Nbfactor")
        unsxx <- unsxx*Nbfactor
        if (needFONKgpointls) {
          xx <- xx*Nbfactor
          gcvxx <- gcvxx*Nbfactor
        }
      }
      if (st =="condS2") {
        S2factor <- blackbox.getOption("S2factor")
        unsxx <- unsxx*S2factor
        if (needFONKgpointls) {
          xx <- xx*S2factor
          gcvxx <- gcvxx*S2factor
        }
      }
      names(unsxx) <- rown
      if (plottype=="both") {
        if (nrow(unselectedKgpt)!=nrow(FONKgpointls)) {
          plot(unsxx, unsyvalues,
               cex=0.5, cex.axis=loccex.axis, col=grey(0.7), xlab="", ylab="", log=logv)
          par(new=TRUE) ## will ignore previous output in the same plot (eg for axis(3), (4))
          axesBool <- TRUE ## already axes on an existing plot
        }
      } else if (plottype=="showIters") { ## plots only unselectedKgpt
        it2rownames <- as.character(seq(iterThresholds+1, nrow(unselectedKgpt))) ## iterThresholds is of length 1 despite its name
        plot(unsxx, unsyvalues,
             pch=".", cex=2, cex.axis=loccex.axis, xlab=userunit(st, locstring),
             ylab="ln(L)", log=logv) ## all points in black to have correct plot size
        points(unsxx[it2rownames], unsyvalues[it2rownames],
               pch=".", cex=2, col="red") ## red over black
      } else if (plottype=="unselected") {
        plot(unsxx, unsyvalues,
             pch=".", cex=2, cex.axis=loccex.axis, xlab=userunit(st, locstring),
             ylab="ln(L)", log=logv)
        axesBool <- TRUE ## already axes on an existing plot
      }
      if ( ! plottype=="unselected" ) {
        plot(xx, syvalues, pch=".", cex=2, cex.axis=loccex.axis, ##pch=15,
             axes=(!axesBool), xlab=userunit(st, locstring),
             ylab="ln(L)", log=logv)
        #       if(DemographicModel=="OnePopVarSize") {
        #            abline(h=max(syvalues)- qchisq(0.95, fittedparamnbr)/2, col=2)
        #      abline(h=max(syvalues)- qchisq(0.99, fittedparamnbr)/2, col=3)
        #          abline(h=max(syvalues)- qchisq(0.999, fittedparamnbr)/2, col=4)
        #	        abline(h=max(syvalues)- qchisq(0.9999, fittedparamnbr)/2, col=5)
        #        }
        if ( plottype %in% c("selected", "both") ) points(gcvxx, gcvyvalues, col="red")
      }
    }
    if ( plottype=="prediction" || addPrediction) {
      ## here problem is because we have used the log option, contrary to the other plots
      if (islogscale(st)) {
        xmin <- exp(FONKgLow[st])
        xmax <- exp(FONKgUp[st])
        plotfn <- function(z) {z <- matrix(log(z), ncol=1);colnames(z) <- st;purefn(z)}
      } else {
        xmin <- FONKgLow[st]
        xmax <- FONKgUp[st]
        plotfn <- function(z) {z <- matrix(z, ncol=1);colnames(z) <- st;purefn(z)}
      }
      plot(plotfn, xmin, xmax, add=(dev.cur()>1),xlab=st,ylab="ln(L)",cex.axis=loccex.axis)
    }
    if (plottype=="both") {
      axis(3,cex.axis=loccex.axis);axis(4,cex.axis=loccex.axis)
    }
  }## end loop over fittedNames
  par(mfrow=c(1, 1))
  par(op)
}
