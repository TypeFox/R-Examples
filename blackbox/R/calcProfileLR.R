## unique procedure for all optimization methods
calcProfileLR <- function(varNames=blackbox.getOption("fittedNames"),
                          pairlist=list(),
                          cleanResu="") {
  INFO <- list(DemographicModel=blackbox.getOption("DemographicModel"),
              spec2DProfiles=blackbox.getOption("spec2DProfiles"),
              plotOptions=blackbox.getOption("plotOptions"),
              #rosglobal=blackbox.getOption("rosglobal"), ## NOTE this cannot be here
              #   because rosglobal may be modified by this function
              fittedparamnbr=blackbox.getOption("fittedparamnbr"),
              graphicPars=blackbox.getOption("graphicPars"),
              interactiveGraphics=blackbox.getOption("interactiveGraphics")
              )
  ## first construct a list of all desired pairs of variables
  if (length(pairlist)==0L) {
    if (INFO$fittedparamnbr>2L) { ## if profiles required
      if (length(INFO$spec2DProfiles)>0L) { ## if speficic profiles required
        for (locit in seq_len(length(INFO$spec2DProfiles))) {
          locstring <- INFO$spec2DProfiles[locit]
          locstring <- substr(locstring,2,(nchar(locstring)-1L)) ## removes "(" ")"
          vars <- strsplit(locstring,",")[[1]]
          pairlist <-  c(pairlist, list(c(vars[1], vars[2])))
        }
      } else if ("All2DProfiles" %innc% INFO$plotOptions ) {
        for (kk in 1:INFO$fittedparamnbr) for (ll in 1:INFO$fittedparamnbr) {
          if (kk<ll) pairlist <- c(pairlist, list(c(varNames[kk], varNames[ll])))
        }
      } # else pairlist <- c()
      ## adds 2Nmu, 2Nm if IBD and not already there
      if (("IBD" %in% INFO$DemographicModel) &&
          ! any(as.logical(lapply(pairlist, function(v){all(c("twoNmu", "twoNm")==v)})))) {
        pairlist <- c(pairlist, list(c("twoNmu", "twoNm")))
      }
      ## adds 2Nm, g if IBD and not already there ## FR->FR voir a un meilleur controle utilisateur...
      if (("IBD" %in% INFO$DemographicModel) && ! any(as.logical(lapply(pairlist, function(v){all(c("twoNm", "g")==v)})))) {
        pairlist <- c(pairlist, list(c("twoNm", "g")))
      }
      ## adds 2Nmu, Nb if IBD and not already there
      if (("IBD" %in% INFO$DemographicModel) && ! any(as.logical(lapply(pairlist, function(v){all(c("twoNmu", "latt2Ns2")==v)})))) {
        pairlist <- c(pairlist, list(c("twoNmu", "latt2Ns2")))
      }
      ## adds 2Nm, condS2 if IBD, condS2 in varNames and not already there
      if (("IBD" %in% INFO$DemographicModel) && "condS2" %in% varNames && ! any(as.logical(lapply(pairlist, function(v){all(c("twoNm", "condS2")==v)})))) {
        pairlist <- c(pairlist, list(c("twoNm", "condS2")))
      }
    } else if (INFO$fittedparamnbr==2L) {pairlist <- list(varNames)}
  }
  if (INFO$fittedparamnbr>2L) {
    localmain <- "Profile likelihood ratio"
  } else localmain <- "Likelihood ratio";
  prevmsglength <- 0
  npaires <- length(pairlist);paireIt <- 0
  ## a check of feasible pairs: same algo as in calcGridRelProfile
  checklist <- list()
  for(paire in pairlist) {
    testvars <- (paire %w/o% varNames) ## should be empty... except in case adhocked below
    ## ad hoc fix for Nm bound /Nb profile
    if ("latt2Ns2" %in% paire && "twoNm" %in% varNames) {testvars <- (testvars %w/o% "latt2Ns2")}
    ## ad hoc fix for Nb bound/Nm profile
    if ("twoNm" %in% paire && "latt2Ns2" %in% varNames) {testvars <- (testvars %w/o% "twoNm")}
    ## ad hoc fix...
    if ("Nratio" %in% paire && "twoNmu" %in% varNames) {testvars <- (testvars %w/o% "Nratio")}
    if ("Nancratio" %in% paire && "twoNmu" %in% varNames) {testvars <- (testvars %w/o% "Nancratio")}
    if ("NactNfounderratio" %in% paire && "twoNmu" %in% varNames) {testvars <- (testvars %w/o% "NactNfounderratio")}
    if ("NfounderNancratio" %in% paire && "twoNancmu" %in% varNames) {testvars <- (testvars %w/o% "NfounderNancratio")}
    if(length(testvars)==0) checklist=c(checklist, list(paire))
  }
  pairlist <- checklist ## now contains only pairs valid for calcGridRelProfile
  if (length(pairlist)>0) {
    message.redef(paste("*** Computing ", length(pairlist), " two-dimensional profile confidence plots *** (may be slow)", sep=""))
    write("\n", file=cleanResu)
  }
  loccex.axis=INFO$graphicPars$cex.axis
  ## main loop
  for (paire in pairlist) {
    paireIt <- paireIt+1
    stk <- paire[1];stl <- paire[2]
    locxstring <- "";locystring <- "";
    if (interactive()) {
      tmp <- sapply(paire, formatName, format="ASCII")
      locmess <- paste("Computing profile grid for (", paste(tmp, collapse=", "), ")")
      locmess <- paste(locmess, " (profile #", paireIt, " out of ", npaires, ")", sep="")
      cat(locmess, "\n")
    }
    loclist <- do.call("calcGridRelProfile",list(fixed=paire))
    xGrid <- loclist$xGrid;yGrid <- loclist$yGrid;RelLik <- loclist$RelLik;inKrigSpace <- loclist$inKrigSpace
    if ( any( ! is.na(RelLik) )) {
      rosglobal <- blackbox.getOption("rosglobal") ## latest one
      if (min(length(xGrid), length(yGrid))<2) next; ## no variation for one variable
      ## we want pretty labels so we must determine xat from pretty xlabs values, not xlabs from pretty xat values...
      xargs <- maketicks(stk, labels=xGrid, axis=1, maxticks=INFO$graphicPars$xmaxticks)
      xat <- xargs$at;xlabs <- xargs$labels;xlegend <- xargs$legend;xph <- xargs$phantomat
      if (stk=="latt2Ns2") {
        xMLcoord <- rosglobal$latt2Ns2 ## not times GV$Nbfactor since internal coord of the plot are (log) latt2Ns2
      } else xMLcoord <- rosglobal$canonVP[stk]
      if(islogscale(stk)) {xMLcoord <- log(xMLcoord)}
      ## same for y
      yargs <- maketicks(stl, labels=yGrid, axis=2, maxticks=INFO$graphicPars$ymaxticks)
      yat <- yargs$at;ylabs <- yargs$labels;ylegend <- yargs$legend;yph <- yargs$phantomat
      if (stl=="latt2Ns2") {
        yMLcoord <- rosglobal$latt2Ns2 ## not times GV$Nbfactor since internal coord of the plot are (log) latt2Ns2
      } else yMLcoord <- rosglobal$canonVP[stl]
      if(islogscale(stl)) {yMLcoord <- log(yMLcoord)}
      if (INFO$interactiveGraphics) provideDevice(bbdefaultPars=TRUE) ## one window for _each_ graphic
      RelLik[RelLik>1] <- 1
      BadExtrapol <- which( ! inKrigSpace, arr.ind=T) ## inKrigSpace is another (no longer global) array created by calcGridRelProfile
      cautiousRelLik <- RelLik;cautiousRelLik[BadExtrapol] <- NA
      if (nrow(BadExtrapol)>1) {
        maxBad <- max(RelLik[BadExtrapol])
        usernames <- sapply(paire, formatName, format="ASCII")
        if (maxBad==1) {
          warningue <- paste("(!)    The highest profile likelihood is extrapolated in the ", paste(usernames, collapse=", "), " profile", sep="")
        } else if (maxBad>0.05) {
          warningue <- paste("    Some high profile likelihoods are extrapolated in the ", paste(usernames, collapse=", "), " profile", sep="")
        } else {
          warningue <- paste("    Some profile likelihoods are extrapolated in the ", paste(usernames, collapse=", "), " profile", sep="")
        }
        message.redef(warningue)
        write(warningue, file=cleanResu)
      }
      if ("cautious" %innc% INFO$plotOptions) {
        RelLik <- cautiousRelLik
      } else if ("Incautious" %innc% INFO$plotOptions) {
        cautiousRelLik <- RelLik
      } else { ## default
        ## contour lines show the uncautious and shading shows the cautious. Analogous, though not exactly, to slice plots
      }
      ## each lik is NA either because (x,y) vars are not in krigSpace or that full (x,y,...) point is not in KrigSpace
      if (all(is.na(cautiousRelLik))) {
        warningue <- paste("    All profile likelihoods are extrapolated in the ", paste(usernames, collapse=", "), " profile", sep="")
        cautiousRelLik <- matrix(-1,nrow=nrow(cautiousRelLik),ncol=ncol(cautiousRelLik)) ## FR->FR quick patch for plot fns
      }
      niveaux <- c(0.001, 0.01, 0.05, seq(1L, 9L, 1L)/10) ## one contour() call => non_overlapping labels
      filled.contour(x=xGrid, y=yGrid, cautiousRelLik, xlab=xlegend,
                     ylab=ylegend, main=localmain, color.palette=spaMM.colors, #topo.colors, ## col and color have different syntaxes
                     plot.axes={contour(xGrid, yGrid, RelLik, add=T, nlevels=1, levels=niveaux, labcex=INFO$graphicPars$labcex);
                                if (!is.null(xph)) axis(1, at=xph, labels=rep("", length(xph)), tcl=-0.1);
                                if (!is.null(yph)) axis(2, at=yph, labels=rep("", length(yph)), tcl=-0.1);
                                #      	points(xGrid[BadExtrapol[, 1]], yGrid[BadExtrapol[, 2]], pch=19); ## (NULL OK on points since limits of plot are known)
                                points(xMLcoord, yMLcoord, pch="+");
                                axis(1, at=xat, labels=mantissExp(xlabs), cex.axis=loccex.axis);
                                axis(2, at=yat, labels=mantissExp(ylabs), cex.axis=loccex.axis)})
      ## for grey scale figure one needs to have package gplots installed;
      if (("BWPlots" %innc% INFO$plotOptions) && length(yGrid)>1) {
        if (INFO$interactiveGraphics) provideDevice(bbdefaultPars=TRUE) ## one window for _each_ graphic
        filled.contour(x=xGrid, y=yGrid, cautiousRelLik, xlab=xlegend,
                       ylab=ylegend, main=localmain,
                       col=colorpanel(20, "grey75", "white"),
                       plot.axes={contour(xGrid, yGrid, RelLik, add=T, nlevels=1, levels=niveaux, labcex=INFO$graphicPars$labcex);
                                  #				   points(xGrid[BadExtrapol[, 1]], yGrid[BadExtrapol[, 2]], pch=19);
                                  points(xMLcoord, yMLcoord, pch="+");
                                  axis(1, at=xat, labels=mantissExp(xlabs), cex.axis=loccex.axis);
                                  axis(2, at=yat, labels=mantissExp(ylabs), cex.axis=loccex.axis)
                       }
        )
      }
    } else {
      warningue <- paste("(!)    No grid point for the ", paste(paire), " profile was within convex hull of Kriged points", sep="")
      message.redef(warningue)
      write(warningue, file=cleanResu)
      warningue <- "Perhaps increase 'gridSteps' so that some grid points fall within the hull." ## dLnLthreshold currently cannot be changed....
      message.redef(warningue)
      write(warningue, file=cleanResu)
    }
  } ## end for paire...
  plot1DprofFrom2D()
}
## for labels on contours: the default method="flattest"
##    does not always show anything, and the others may be ugly
## for contours only without shading: execute the contour(...) above
##    only and without the add=T statement.
