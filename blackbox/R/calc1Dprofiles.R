calc1Dprofiles <- function(varNames=blackbox.getOption("spec1DProfiles")) {
  # **** 1D Profiling Computations: ****
  ##Grid (RelLik) of relative profile likelihoods:
  ##arguments to control grid: gridsteps, xGrid, yGrid
  ##have defaults (see calcGridRelProfile definition) which can be modified.
  #fixedvars <- c()
  INFO <- list(DemographicModel=blackbox.getOption("DemographicModel"),
               fittedparamnbr=blackbox.getOption("fittedparamnbr"),
               #rosglobal=blackbox.getOption("rosglobal"), ## NOTE this cannot be here
               #   because rosglobal may be modified by this function
               plotOptions=blackbox.getOption("plotOptions"),
               ParameterNames=blackbox.getOption("ParameterNames"),
               constantNames=blackbox.getOption("constantNames"),
               interactiveGraphics=blackbox.getOption("interactiveGraphics"),
               graphicPars=blackbox.getOption("graphicPars"))
  if (length(varNames)==0L) {
    if("all1DProfiles" %innc% INFO$plotOptions) {
      varNames <- INFO$ParameterNames ## original conception, does not use fittedNames
      if ("IBD" %in% INFO$DemographicModel) varNames <- c(varNames, "latt2Ns2")
      if((length(intersect(INFO$DemographicModel, c("OnePopVarSize", "IM")))) && ("PopSizeRatioProf" %innc% INFO$plotOptions)) varNames <- c(varNames, "Nratio")
      else if("OnePopFounderFlush" %in% INFO$DemographicModel && ("PopSizeRatioProf" %innc% INFO$plotOptions)) varNames <- c(varNames, "Nratio", "NactNfounderratio", "NfounderNancratio")
    } else varNames <- c()
  }
  varNames <- varNames %w/o% INFO$constantNames
  varNames <- unique(varNames)
  nprofs <- length(varNames)
  It <- 0
  intsqrt <- as.integer(sqrt(nprofs))
  if (intsqrt>0) {
    message.redef(paste("*** Computing ", nprofs, " one-dimensional profile confidence plots *** (may be slow)", sep=""))
    if (INFO$interactiveGraphics) provideDevice(bbdefaultPars=TRUE)
    op <- par(no.readonly = TRUE) # the whole list of settable par's to be reset afterwards
    if (intsqrt>1) {loccex.axis <- INFO$graphicPars$cex.axis*0.6} else {loccex.axis <- INFO$graphicPars$cex.axis}
    if(nprofs<=4) {par(mfrow=c(ceiling(nprofs/intsqrt), intsqrt), cex.axis=loccex.axis); It2=intsqrt*intsqrt;} else {par(mfrow=c(2, 2), cex.axis=loccex.axis); It2=4}
    prevmsglength <- 0
    for (st in varNames) {
      if(It==It2) {
        title("One-parameter likelihood ratio profiles", outer = TRUE, line=-1)
        if (INFO$interactiveGraphics) provideDevice(bbdefaultPars=FALSE)
      }
      It <- It+1
      if (interactive()) {cat(paste("Computing profile grid for ", formatName(st, format="ASCII")), "\n")}
      loclist <- calcGridRelProfile(fixed=st)
      xGrid <- loclist$xGrid;yGrid <- loclist$yGrid;RelLik <- loclist$RelLik;inKrigSpace <- loclist$inKrigSpace
      xy <- cbind(xGrid, RelLik)
      rosglobal <- blackbox.getOption("rosglobal") ## get the very latest value
      if (st %in% names(rosglobal$canonVP)) {
        xMLcoord <- rosglobal$canonVP[st]
      } else if (st=="Nratio") {
        xMLcoord <- rosglobal$Nratio
      } else if (st=="Nancratio") {
        xMLcoord <- rosglobal$Nratio  ## Nratio et non Nancratio OK
      } else if (st=="NactNfounderratio") {
        xMLcoord <- rosglobal$NactNfounderratio
      } else if (st=="NfounderNancratio") {
        xMLcoord <- rosglobal$NfounderNancratio
      } else {
        if (st=="latt2Ns2") xMLcoord <- rosglobal$latt2Ns2 ## internal
      }
      if (islogscale(st)) xMLcoord <- log(xMLcoord)
      xy <- rbind(xy, c(xMLcoord, 1))
      xy <- xy[ ! is.na(xy[, 2]), , drop=FALSE]
      # 10/2015: sets plot()'s ylim argument but does not select points
      if (nrow(xy)>1) xy <- xy[order(xy[, 2], decreasing=TRUE), ]
      lowlogLR <- qchisq(0.99, 1)*0.525 ## so that exp(-lowlogLR) is slightly below the level for 99% coverage
      ylimmin <- min(exp(-lowlogLR),xy[min(5, nrow(xy)), 2]/2 ) ## to see the 0.99 level, and at least 5 points
      ylimmin <- max(1e-100,ylimmin) ## ylimmin final must be >0
      xy <- xy[order(xy[, 1], decreasing=FALSE), ] ## x order important for plot joined...
      xargs <- maketicks(st, labels=xy[, 1], axis=1, maxticks=INFO$graphicPars$xmaxticks)
      xat <- xargs$at;xlabels <- xargs$labels;xlegend <- xargs$legend;xph <- xargs$phantomat
      plot(x=xy[, 1], y=xy[, 2], xlab=xlegend, ylab="Likelihood ratio",
           xaxt="n", main="", xlim=c(min(xat), max(xat)), ylim=c(ylimmin, 1.05), log="y", #log=loclog
           cex.axis=loccex.axis, pch=3, type="b", cex=0.5
      )
      points(xMLcoord, 1, pch=1)
      axis(1, at=xat, labels=mantissExp(xlabels), cex.axis=1.4*loccex.axis)
      if (!is.null(xph)) axis(1, at=xph, labels=rep("", length(xph)), tcl=-0.1)
      abline(h=exp(-qchisq(0.95, df=1)/2), col=2)
      text(xat[1], pos=4, exp(-qchisq(0.95, df=1)/2)*1.2, "0.95", col=2, cex=loccex.axis, offset=-0.3)
      abline(h=exp(-qchisq(0.99, df=1)/2), col=3)
      text(xat[1], pos=4, exp(-qchisq(0.99, df=1)/2)*1.2, "0.99", col=3, cex=loccex.axis, offset=-0.3)
    }
    title("One-parameter likelihood ratio profiles", outer = TRUE, line=-1)
    par(mfrow=c(1, 1))
    par(op)
  }
  # **** end 1D Profiling Computations: ****
  invisible(NULL)
}

plot1DprofFrom2D <- function(varNames=blackbox.getOption("spec1DProfiles")) {
  INFO <- list(margProfsInfo=blackbox.getOption("margProfsInfo"),
               DemographicModel=blackbox.getOption("DemographicModel"),
               fittedparamnbr=blackbox.getOption("fittedparamnbr"),
               plotOptions=blackbox.getOption("plotOptions"),
               rosglobal=blackbox.getOption("rosglobal"), ## should not vary within the function
               ParameterNames=blackbox.getOption("ParameterNames"),
               constantNames=blackbox.getOption("constantNames"),
               interactiveGraphics=blackbox.getOption("interactiveGraphics"),
               graphicPars=blackbox.getOption("graphicPars"))
  margProfsInfo <- INFO$margProfsInfo
  if (length(varNames)==0L) {
    if("all1DProfiles" %innc% INFO$plotOptions) {
      varNames <- INFO$ParameterNames ## original conception, does not use fittedNames
      if ("IBD" %in% INFO$DemographicModel) varNames <- c(varNames, "latt2Ns2")
      if((length(intersect(INFO$DemographicModel, c("OnePopVarSize", "IM")))) && ("PopSizeRatioProf" %innc% INFO$plotOptions)) varNames <- c(varNames, "Nratio")
      if("OnePopFounderFlush" %in% INFO$DemographicModel && ("PopSizeRatioProf" %innc% INFO$plotOptions)) varNames <- c(varNames, "Nratio", "NactNfounderratio", "NfounderNancratio")
    } #else varNames <- c()
  }
  varNames <- varNames %w/o% INFO$constantNames
  varNames <- unique(varNames)
  if (length(varNames)>0L) margProfsInfo <- margProfsInfo[intersect(names(margProfsInfo),varNames)]
  nprofs <- length(margProfsInfo)
  It <- 0
  intsqrt <- as.integer(sqrt(nprofs))
  if (intsqrt>0) {
    if (interactive()) {cat("Drawing 1D profiles after 2D profiles computation\n")}
    if (INFO$interactiveGraphics) provideDevice(bbdefaultPars=TRUE)
    op <- par(no.readonly = TRUE) # the whole list of settable par's to be reset afterwards
    if (intsqrt>1) {loccex.axis <- INFO$graphicPars$cex.axis*0.6} else {loccex.axis <- INFO$graphicPars$cex.axis}
    if(nprofs<=4) {par(mfrow=c(ceiling(nprofs/intsqrt), intsqrt), cex.axis=loccex.axis); It2=intsqrt*intsqrt;} else {par(mfrow=c(2, 2), cex.axis=loccex.axis); It2=4}    
    prevmsglength <- 0
    for (st in names(margProfsInfo)) {
      if(It==It2) {
        title("One-parameter likelihood ratio profiles", outer = TRUE, line=-1)
        if (INFO$interactiveGraphics) provideDevice(bbdefaultPars=FALSE)
      }
      It <- It+1
      #if (interactive()) {cat(paste("Computing profile grid for ", formatName(st, format="ASCII")), "\n")}
      loclist <- margProfsInfo[[st]]
      xGrid <- loclist$xGrid
      #yGrid <- loclist$yGrid
      RelLik <- loclist$RelLik
      inKrigSpace <- loclist$inKrigSpace
      xy <- cbind(xGrid, RelLik)
      if (st %in% names(INFO$rosglobal$canonVP)) {
        xMLcoord <- INFO$rosglobal$canonVP[st]
      } else if (st=="Nratio") {
        xMLcoord <- INFO$rosglobal$Nratio
      } else if (st=="Nancratio") {
        xMLcoord <- INFO$rosglobal$Nratio  ## Nratio et non Nancratio OK
      } else if (st=="NactNfounderratio") {
        xMLcoord <- INFO$rosglobal$NactNfounderratio
      } else if (st=="NfounderNancratio") {
        xMLcoord <- INFO$rosglobal$NfounderNancratio
      } else {
        if (st=="latt2Ns2") xMLcoord <- INFO$rosglobal$latt2Ns2 ## internal
      }
      if (islogscale(st)) xMLcoord <- log(xMLcoord)
      xy <- rbind(xy, c(xMLcoord, 1))
      xy <- xy[ ! is.na(xy[, 2]), , drop=FALSE]
      # 10/2015: sets plot()'s ylim argument but does not select points
      if (nrow(xy)>1) xy <- xy[order(xy[, 2], decreasing=TRUE), ]
      lowlogLR <- qchisq(0.99, 1)*0.525 ## so that exp(-lowlogLR) is slightly below the level for 99% coverage
      ylimmin <- min(exp(-lowlogLR),xy[min(5, nrow(xy)), 2]/2 ) ## to see the 0.99 level, and at least 5 points
      ylimmin <- max(1e-100,ylimmin) ## ylimmin final must be >0
      xy <- xy[order(xy[, 1], decreasing=FALSE), ] ## x order important for plot joined...
      xargs <- maketicks(st, labels=xy[, 1], axis=1, maxticks=INFO$graphicPars$xmaxticks)
      xat <- xargs$at;xlabels <- xargs$labels;xlegend <- xargs$legend;xph <- xargs$phantomat
      plot(x=xy[, 1], y=xy[, 2], xlab=xlegend, ylab="Likelihood ratio",
           xaxt="n", main="", xlim=c(min(xat), max(xat)), ylim=c(ylimmin, 1.05), log="y", #log=loclog
           cex.axis=loccex.axis, pch=3, type="b", cex=0.5
      )
      points(xMLcoord, 1, pch=1)
      axis(1, at=xat, labels=mantissExp(xlabels), cex.axis=1.4*loccex.axis)
      if (!is.null(xph)) axis(1, at=xph, labels=rep("", length(xph)), tcl=-0.1)
      abline(h=exp(-qchisq(0.95, df=1)/2), col=2)
      text(xat[1], pos=4, exp(-qchisq(0.95, df=1)/2)*1.2, "0.95", col=2, cex=loccex.axis, offset=-0.3)
      abline(h=exp(-qchisq(0.99, df=1)/2), col=3)
      text(xat[1], pos=4, exp(-qchisq(0.99, df=1)/2)*1.2, "0.99", col=3, cex=loccex.axis, offset=-0.3)
    }
    title("One-parameter likelihood ratio profiles", outer = TRUE, line=-1)
    par(mfrow=c(1, 1))
    par(op)
  }
  invisible(NULL)
}
