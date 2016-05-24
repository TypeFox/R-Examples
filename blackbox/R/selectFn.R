selectFn <- function(ptls,
                    topmode=FALSE,
                    dLnLthreshold=NA,
                    designRetain=1,
                    minPtNbr=30,
                    maxPtNbr=NA,
                    replicatesInfo=NA,
                    ycolname=blackbox.getOption("ycolname"),
                    fittedNames =blackbox.getOption("fittedNames"),
                    fittedLoci= blackbox.getOption("respCols"),
                    verbose=F) {
  fittedparamnbr <- length(fittedNames)
  designRetain <- min(1,designRetain)
  if ( ! is.numeric(maxPtNbr)) {
    maxPtNbr <- as.integer(nrow(ptls)*designRetain) ## default is some cover design
    if (verbose) {print(paste("non-numeric maxPtNbr set to ", maxPtNbr, sep=""), quote=F)}
  }
  # three options for topmode:
  # topmode is "dLnL" => cover design dLnL highest points
  # topmode is FALSE => cover design in all points
  # topmode is number=> take topmode highest points
  if(!(topmode==FALSE || topmode=="dLnL" || is.numeric(topmode))) {
    message.redef("(!) Incorrect topmode specification in selectFn")
  }
  # dLnLthreshold is number => use non default value for "dLnL" likelihood threshold
  if(!(is.na(dLnLthreshold) || is.numeric(dLnLthreshold))) {
    message.redef("(!) Incorrect dLnLthreshold specification in selectFn")
  }
  if (minPtNbr>nrow(ptls)) {
    locstring <- paste("(!) minPtNbr>nrow(ptls) (", minPtNbr, " vs. ", nrow(ptls), ").", sep="")
    message.redef(locstring)
  }
  if (is.numeric(topmode) && topmode>nrow(ptls)) {
    locstring <- paste("(!) topmode>nrow(ptls) (", topmode, " vs. ", nrow(ptls), "). Computation will likely fail (?).", sep="")
    message.redef(locstring)
  }
  if (is.numeric(topmode) && is.numeric(maxPtNbr)) {
    cat("topmode argument overrides maxPtNbr argument", "\n")
    maxPtNbr <- topmode
  }
  if(!("replicatesInfoClass" %in% class(replicatesInfo))) {
    replicatesInfo <- findReplicates(ptls)
    if (verbose) cat("replicatesInfo recomputed within selectFn", "\n")
  }
  warningue <- NA
  singletonsnames <- replicatesInfo$singletonsnames
  doublonsFirstsnames <- replicatesInfo$doublonsFirstsnames
  doublonsSecondsnames <- replicatesInfo$doublonsSecondsnames
  if (is.null(doublonsFirstsnames) || length(doublonsFirstsnames)==0) {
    stop.redef("(!!) No replicates found among input points. I exit.")
  } else {
    RMSy <- replicatesInfo$pureRMSE
  }
  ## below are *typically tables* of y values (par locus and total)
  ## one reason for keeping per-locus values may be the 'modified' case below
  ## sortandxy will contain both singletons and mean values of replicates
  sortandxy <- ptls[singletonsnames, c(fittedNames, fittedLoci, ycolname)]
  if (verbose) {
    locstring <- paste("sortandxy has ", nrow(sortandxy), " singlets", ep="")
    print(locstring, quote=F)
  }
  tmp <- cbind(ptls[doublonsFirstsnames, fittedNames, drop=FALSE], replicatesInfo$doublonsymeans) ## both in the same row order !
  sortandxy <- rbind(sortandxy, tmp)
  if (verbose) {
    locstring <- paste("sortandxy now has ", nrow(sortandxy), " rows from singlets and doublets", sep="")
    print(locstring, quote=F)
  }
  if(nrow(sortandxy)<minPtNbr/designRetain) {
    message.redef("    Not enough points satisfy the minPtNbr constraint")
    message.redef("          (as controlled by the minKrigPtNbr setting).")
  }
  #
  if (topmode=="dLnL") {
    if (verbose) {cat("Performing topmode selection", "\n")}
    # evaluates dLnLthreshold if no explicit value taken
    if (is.na(dLnLthreshold) || dLnLthreshold<=0) { #default
      dlrq <- qchisq(0.9999, fittedparamnbr)/2 ## log L difference
      #sqrt(2*ln(n)) is the likely extreme value of n Normal(0, 1) random deviates (e.g. Durrett05, p.83)
      extremev <- RMSy*sqrt(2*log(nrow(ptls)))
      #dlr <- dlrq+extremev ## pb if extremev> dlrq: variation of signal < noise. Try
      if (dlrq> 2*extremev) {
        dlr <- dlrq+extremev
        message.redef("default dLnL threshold value (Chi-square threshold + probable extreme deviate): ")
        locstring <- paste(prettynum(dlr), " ( = ", prettynum(dlrq), " + ", prettynum(extremev), " )", sep="")
      } else {
        dlr <- 3*extremev ## facteur 3 choisi empiriquement (2 n'est pas assez)
        message.redef("default dLnL threshold value, given large variance of replicates: ")
        locstring <- paste(prettynum(dlr), " ( = 3 * probable extreme deviate )", sep="")

      }
      message.redef(locstring)
    } else {
      dlr <- dLnLthreshold ## a positive value has been provided
      locstring <- paste("dLnL threshold value ", dlr, " as given in function call.", sep="")
      message.redef(locstring)
    }
  } else if (is.numeric(topmode)) { ## translate topmode in dlr for later simplicity of code
    message.redef("(!) From selectFn(): numeric(topmode) is probably a bad idea !!!!! ")
    tmp <- sort(ptls[, ycolname])
    if (topmode<nrow(tmp)) {
      dlr <- (tmp[topmode]+tmp[topmode+1])/2-tmp[1]
      if (verbose) {
        locstring <- paste("dLnL threshold value ", dlr, " as implied by numeric topmode value.", sep="")
        print(locstring, quote=F)
      }
    } else {
      dlr <- (tmp[topmode]+1)/2-tmp[1]
      if (verbose) {
        locstring <- paste("dLnL threshold value ", dlr, " as implied by numeric topmode value>nrow(tmp).", sep="")
        print(locstring, quote=F)
      }
    }
  } else dlr <- NA
  # then selects points according to this threshold
  if ( ! is.na(dlr)) {
    minMLL <- min(ptls[, ycolname])
    uppersurf <- which(sortandxy[, ycolname]<(minMLL+dlr)) ## indices
    if (verbose) {
      locstring <- paste(length(uppersurf), " rows selected according to dLnL threshold", sep="")
      print(locstring, quote=F)
    }
  } else { ## is.na(dlr)
    uppersurf <- 1:(nrow(sortandxy)) ## not dLnL/topmode selection, all points are still there
    if (verbose) {
      cat("All sortandxy rows retained as no dLnL threshold is available", "\n")
    }
  }
  # now looks whether there are enough coordinates for cover.design selection
  if(length(uppersurf)>=minPtNbr/designRetain) { ## ie if enough points
    ## then keep all of them for later selection
    selectedxy <- sortandxy[uppersurf, ] ## no need to sort according to y values
    if (verbose) {
      locstring <- paste("All ", length(uppersurf), " rows retained as length(uppersurf)>=minPtNbr/designRetain", sep="")
      print(locstring, quote=F)
    }
  } else { ## won't be able to retain at least minPtNbr points
    ## => work with the original sortandxy, not with the selection sortandxy[uppersurf, ]
    ## The maximum number that can be retained in sortandxy:
    selectableNbr <- min(as.integer(minPtNbr/designRetain+1), nrow(sortandxy))
    if (verbose) {
      locstring <- paste("only ", selectableNbr, " rows retained as length(uppersurf)<minPtNbr/designRetain", sep="")
      print(locstring, quote=F)
    }
    locstring <- paste("(!) Few points in upper ", prettynum(dlr), " [ln(L) units] range:", sep="")
    message.redef(locstring)
    warningue <- locstring
    locstring <- paste("	   only ", length(uppersurf), " points in this range. ", sep="") ##
    message.redef(locstring)
    warningue <- paste(warningue, locstring, sep="\n")
    locstring <- paste("	   Pre-selecting ", selectableNbr, " highest mean likelihoods instead ", sep="")
    message.redef(locstring)
    message <- paste("     (this can be controlled by the minKrigPtNbr setting).", sep="")
    message.redef(locstring)
    if (length(fittedLoci)==0) {locfittedLoci <- ycolname} else {locfittedLoci <- fittedLoci}
    yvalues <- apply(sortandxy[, locfittedLoci, drop=FALSE], 1, sum) ##
    selectedxy <- sortandxy[(order(yvalues, decreasing=F))[1:selectableNbr], ] ## lowest *-*LnL first!
  }
  # distinguishes singletons and doublons in the above result
  remainingsingletons <- (rownames(selectedxy) %w/o% doublonsFirstsnames) %w/o% doublonsSecondsnames
  remainingdoublonsFirst <- (rownames(selectedxy) %w/o% remainingsingletons)
  locstring <- paste("    ", length(remainingsingletons)+2*length(remainingdoublonsFirst),
                   " points pre-selected, (", length(remainingsingletons), " singlets plus ",
                   length(remainingdoublonsFirst), " doublets)", sep="")
  message.redef(locstring)
  CovFnParam <- blackbox.getOption("CovFnParam")
  if (length(CovFnParam)>0 && all(CovFnParam[fittedNames]>0)) scales <- CovFnParam[fittedNames] else {
    scales <- numeric(fittedparamnbr)
    for (ff in 1:fittedparamnbr) {scales[ff] <- max(selectedxy[, ff]-min(selectedxy[, ff]))}
    ## keep in mind selectedxy has only fittedparamnbr 'x' columns
  }
  if (any(is.na(scales))) {
    stop.redef("(!) From selectFn(): 'scales' contains NA; check contents of selectedxy")
  }
  ## the following block provides remainingxynames containing singletonsnames (possibly none) and first names of some doublons
  if (length(remainingdoublonsFirst)==maxPtNbr) {
    ## keep them and only them, as they are well separated from each other
    remainingxynames <- remainingdoublonsFirst
    if (verbose) cat(" # of doublets exactly equal to maxPtNbr", "\n")
  } else if (length(remainingdoublonsFirst)>maxPtNbr) {
    ## keep as much of them as possible, and no singleton
    ## Not tested for length(fittedNames)==1
    if (verbose) cat(" Selecting maxPtNbr points out of more doublets", "\n")
    remainingxynames <- stripclosestpairsWrapper(selectedxy[remainingdoublonsFirst, fittedNames, drop=FALSE],
                                               finallength=maxPtNbr,
                                               scales=scales,
                                               0) ############################################################################# 1 is just a trick to go forward. We need *some* fixed point.
  } else {
    ## keep all doublons, but singletons are needed
    remainingxy <- selectedxy[c(remainingsingletons, remainingdoublonsFirst), fittedNames, drop=FALSE] ## doublons last for stripclosestpairs!
    ndsafe <- min(maxPtNbr, as.integer(nrow(remainingxy)*designRetain))
    # so that final length is min(maxPtNbr, as.integer(nrow(remainingxy)*designRetain)
    if (ndsafe <= length(remainingdoublonsFirst)) {
      remainingxynames <- remainingdoublonsFirst ## FR 18/12/2014
    } else if (length(remainingdoublonsFirst)==0) {
      stop.redef("(!) No doublets remaining!")
    } else {
      if (verbose) cat(" retaining all doublets and selecting additional singlets", "\n")
      remainingxynames <- stripclosestpairsWrapper(remainingxy,
                                                 finallength=ndsafe,
                                                 scales=scales,
                                                 fixedNbr=length(remainingdoublonsFirst))
    }
  }
  ## then find matching second names of doublons first namesfrom positions (not names) in each list
  posfirsts <- which(doublonsFirstsnames %in% remainingxynames)
  allnames <- c(remainingxynames, doublonsSecondsnames[posfirsts]) ## assumes nothing about rownames (indices) in a doublon
  ## but assumes that pair mates are at the same positions in the two doublons names lists
  ptls <- ptls[allnames, ]
  ## ne pas trier sur les names car incorrect si plusieurs iterations de Migraine dans le pointls.
  ptls <- ptls[ do.call(order, ptls) , ] ## important sinon probleme avec CKrigcoefs
  locstring <- paste("    ", length(remainingxynames)+length(remainingdoublonsFirst),
                   " points selected, (", length(remainingxynames)-length(remainingdoublonsFirst), " singlets plus ",
                   length(remainingdoublonsFirst), " doublets)", sep="")
  message.redef(locstring)
  return(list(ptls=ptls, nuniquerows=length(remainingxynames), RMSy=RMSy, warningue=warningue))
}  ## end def selectFn
