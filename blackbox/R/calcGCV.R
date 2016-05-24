## data would need attributes fittedparamnbr, CovFnParam,fittedNames,GCVdesignRetain
calcGCV <- function(sorted_data=data, ## assumes sorted data
                    data, ## back compatibility
                    CovFnParam=NULL,
                    GCVptnbr=Inf,
                    topmode=FALSE,
                    verbose=FALSE,
                    cleanResu="",
                    force=FALSE,
                    decreasing=FALSE, ## important for call to selectFn
                    verbosity=blackbox.getOption("verbosity"),
                    optimizers=blackbox.getOption("optimizers") ## change this manually to NLOPT_LN_BOBYQA if this is the desired optimizer
) {
  if ("Migraine" %in% blackbox.getOption("usedBy")) {
    fittedNames <- blackbox.getOption("fittedNames")
    ycolname <- blackbox.getOption("ycolname")
  } else {
    fittedNames <- attr(sorted_data,"fittedNames")
    ycolname <- attr(sorted_data,"ycolname")
  }
  fittedparamnbr <- length(fittedNames)
  if (is.null(GCVptnbr) || GCVptnbr == -1) {
    GCVptnbr <- floor(250**((fittedparamnbr/3)**(1/2))) ## default minimum point nbr for cross validation
    if (fittedparamnbr==1) GCVptnbr <- 90 ## default minimum point nbr for cross validation for OnePop 90 = floor(250**((2/3)**(1/2)))
  }
  ## add missing names to params given by users
  if ( length(CovFnParam) > 0) {
    if(length(CovFnParam)==fittedparamnbr) {
      names(CovFnParam) <- fittedNames ## selectFn can use this before any CovFnParam estimation
    } else {
      stop.redef("The length of 'CovFnParam' should be that of blackbox.getOption('fittedNames')")
    }
    if (verbose) {
      locstring <- paste("\n*** Estimation of lambda parameter by generalized cross-validation ***", sep="")
      message.redef(locstring)
  }
  } else {
  if (verbose) {
    locstring <- paste("\n*** Estimation of Kriging parameters by generalized cross-validation ***", sep="")
    message.redef(locstring)
    }
  }
    locstring <- paste("    Selection of points for GCV:", sep="")
    message.redef(locstring)
  if ("Migraine" %in% blackbox.getOption("usedBy")) blackbox.options(CovFnParam=CovFnParam)
  ##
  ## lambdaEst must always be estimated even if the CovFnParams are given.
  replicatesInfo <- findReplicates(sorted_data) ##
  if ( length(replicatesInfo$doublonsFirstsnames)==0L) {
    if (force) {
      resp_order <- order(sorted_data[,ycolname],decreasing=decreasing)
      bestvals_idxs <- resp_order[seq(min(19L,nrow(sorted_data)))] ## 19 best values if possible
      data <- rbind(sorted_data,sorted_data[bestvals_idxs,,drop=FALSE])
    } else {
      message.redef("No replicates found among input points.")
      message.redef("Use force=TRUE if this is really what you want,")
      message.redef("otherwise check input.")
      stop.redef("Exiting because no replicates among input points.")
    }
  }
  if (verbose) cat(paste("Pure RMSE from all smoothed duplicated points= ", prettynum(replicatesInfo$pureRMSE), sep=""), "\n")
  if (GCVptnbr<nrow(sorted_data)) {
    blob <- selectFn(sorted_data,
                     topmode=topmode, ## default is FALSE => selection only by stripClosestPairs with distinct handling of replicates
                     , #dLnLthreshold default
                     designRetain=blackbox.getOption("GCVdesignRetain"), ## default=1 since close points are useful for GCV
                     minPtNbr=GCVptnbr,
                     maxPtNbr=GCVptnbr,
                     replicatesInfo=replicatesInfo)
    gcvKgpointls <- blob$ptls
    gcvnuniquerows <- blob$nuniquerows
  } else {
    gcvnuniquerows <- length(replicatesInfo$doublonsFirstsnames)+length(replicatesInfo$singletonsnames)
    gcvKgpointls <- sorted_data
  }
  if ("Migraine" %in% blackbox.getOption("usedBy")) blackbox.options(GCVrownames=rownames(gcvKgpointls)) ## stor only the row names for future reference
  if (verbose) {
    locstring <- paste("  (see GCVptnbr setting(Migraine)/argument(R) for controlling this).", sep="")
    message.redef(locstring)
  }
  if( ! is.null(cleanResu) && gcvnuniquerows<floor(80**((fittedparamnbr/3)**(1/3)))) {## threshold for warning in cleanResu: 20  45  80 124 180 249
    write(paste("    Only ", gcvnuniquerows, " unique coordinates selected for cross-validation.", sep=""), file=cleanResu)
  }
  initCovFnParam <- NULL
  if ("Migraine" %in% blackbox.getOption("usedBy")) { ## FR->FR try with bboptim next time...
    if ( ! is.null(nextPointsf <- blackbox.getOption("nextPointsf"))) { ## name provided by *Migraine* default settings
      tmp <- try(suppressWarnings(readLines(con = nextPointsf)),silent=TRUE)
      if( ! inherits(tmp,"try-error") ) {
        tmp <- tmp[grep("initCovFnParam",x=tmp)]
        if (length(tmp)>0L) {
          if (blackbox.getOption("verbosity")) cat(paste("Found",tmp,"in 'nextpoints' file\n"))
          tmp <- strsplit(tmp,"= ")[[1]][2]
          tmp <- strsplit(tmp," ")
          initCovFnParam <- as.numeric(unlist(tmp))
        }
      }
    }
  }
  ## estimates CovFnParam vector by GCV
  loc <- provideCovFnParams(gcvKgpointls=gcvKgpointls, gcvnuniquerows=gcvnuniquerows,
                            fittedNames=fittedNames,
                            ycolname=ycolname,
                            initCovFnParam=initCovFnParam,
                            cleanResu=cleanResu,
                            verbosity=verbosity, optimizers=optimizers)
  resu <- c(loc[c("CovFnParam", "lambdaEst", "hglmLambdaEst", "hglmPhiEst")], pureRMSE=replicatesInfo$pureRMSE)
  resu[unlist(lapply(resu,is.null))] <- NULL
  if ("Migraine" %in% blackbox.getOption("usedBy")) do.call(blackbox.options, resu)
  resu$GCVmethod <- loc$method
  invisible(resu)
}
