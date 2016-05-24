calcPredictorOK <- function(FONKgpointls,
                            minKrigPtNbr=blackbox.getOption("minKrigPtNbr"),
                            krigmax=NULL,
                            topmode=FALSE,
                            rawPlots=TRUE,
                            cleanResu=""
) {
  fittedNames <- blackbox.getOption("fittedNames")
  FONKgNames <- blackbox.getOption("FONKgNames")
  fittedparamnbr <- blackbox.getOption("fittedparamnbr")
  if (is.null(minKrigPtNbr) || minKrigPtNbr==-1) {
    minKrigPtNbr <- floor(500**((fittedparamnbr/3)**(1/2))) ## default minimum point nbr for Kriging 36  159  500 1307 3050 6560
  } # else if (minKrigPtNbr<=GCVptnbr) {minKrigPtNbr <- GCVptnbr} ## safe, but restricts experimentation
  if (fittedparamnbr==1) minKrigPtNbr <- 90 ## default minimum point nbr for Kriging for OnePop 90 = floor(250**((2/3)**(1/2)))
  locstring <- paste("  Selection of points for Kriging:", sep="")
  message.redef(locstring)
  unselectedKgpt <- FONKgpointls ## for figures and second selectFn() call
  ## selectFn -> findReplicates -> (in case incorrect imput is detected) -> canonize
  ##  -> uses getOption("FONKgLow") to build message
  ## => FONKgLow must have been set previously by buildFONKgpointls (and will be reset below)
  blob <- selectFn(FONKgpointls,
                   topmode=topmode,
                   #, #dLnLthreshold default
                   designRetain=blackbox.getOption("designRetain"),
                   minPtNbr=minKrigPtNbr, ## il faut bien utiliser la valeur locale !
                   maxPtNbr=blackbox.getOption("maxKrigPtNbr"),
                   replicatesInfo=NA)
  FONKgpointls <- blob$ptls
  FONKgLow <- apply(FONKgpointls[,FONKgNames, drop=FALSE], 2, min)
  FONKgUp <- apply(FONKgpointls[, FONKgNames, drop=FALSE], 2, max)
  ## updating of previous existing vars: (cf buildFONKgpointls or prepareData)
  do.call(blackbox.options, list(FONKgpointls=FONKgpointls,FONKgLow=FONKgLow, FONKgUp=FONKgUp))
  nuniquerows <- blob$nuniquerows
  blackbox.options(RMSy = blob$RMSy)
  if (! is.na(blob$warningue)) write(blob$warningue, file=cleanResu) ## writing this at this step prevents writing it at each selectFn call
  if(nuniquerows<floor(250**((fittedparamnbr/3)**(1/3)))) {## threshold for warning in cleanResu: 24    90   250   587  1246  2461
    write(paste("    Only ", nuniquerows, " unique coordinates selected for Kriging.\n", sep=""), file=cleanResu)
  }
  #
  blackbox.options(hulls=NULL) ## important to make sure that GV$hulls are computed once and only once for a given FONKgpointls
  locstring <- paste("  (see min/maxKrigPtNbr setting(Migraine)/argument(R) for controlling this).", sep="")
  message.redef(locstring)
  if (rawPlots) {
    if (length(intersect(blackbox.getOption("DemographicModel"), c("IBD", "OnePop" )))>0) {
      rawProfiles(unselectedKgpt=unselectedKgpt, FONKgpointls=FONKgpointls, "both")
    } else {
      rawProfiles(unselectedKgpt=unselectedKgpt, FONKgpointls=FONKgpointls, "unselected")
      rawProfiles(unselectedKgpt=unselectedKgpt, FONKgpointls=FONKgpointls, "selected")
    }
  }
  if ( ! ("extrapolateOutOfHull" %innc% blackbox.getOption("miscOptions"))) {
    providefullhull(fittedNames) ## sets .blackbox.data$options$hulls$Kgtotal
  }
  krigit <- 0
  if (is.null(krigmax) || krigmax<0) { ## no explicit value: apply defaults
    if (.Platform$r_arch=="i386") {
      #do.call(blackbox.options, list(krigmax=2198, kriglength=2000, krigoverlap=500))
      krigmax <- 2198
      kriglength <- 2000
      krigoverlap <- 500
    } else { ## x64
      #do.call(blackbox.options, list(krigmax=10000, kriglength=6000, krigoverlap=1000))
      krigmax <- 10000
      kriglength <- 6000
      krigoverlap <- 1000
    }
  } else {
    kriglength <- blackbox.getOption("kriglength")
    krigoverlap <- blackbox.getOption("krigoverlap")
  }
  message.redef("\n*** Kriging ***")
  ycolname <- blackbox.getOption("ycolname")
  flushCSmoothTable() ## although calcPredictorOK() is called only in migraine use so flushing the pointer table is unlikely to be required
  if (nrow(FONKgpointls)<=krigmax) {
    fitobject <- generatePredictor(,first=1L, last=nrow(FONKgpointls)) ## -> OKrig -> CKrigcoefs
    #fitobject$shat.GCV <- NA ## otherwise error on 'summarizing' fit ## (old comment pre 2015)
  } else {
    fitobject <- list() ## will contain $x, $y, $fitted.values, $blocmin $blocmax, $Kriglist
    class(fitobject) <- c("OKriglistplus", class(fitobject))
    krigit <- (nrow(FONKgpointls)-krigoverlap)%/%(kriglength-krigoverlap) ## for krigit+1 blocks
    fitobject$blocmin <- numeric(krigit)
    fitobject$blocmax <- numeric(krigit) ## length=krigit sufficient for krigit+1 blocks
    fitlist <- list()
    class(fitlist) <- c("OKriglist", class(fitlist))
    for (bloc in 0:krigit) {
      fit <- generatePredictor(,first=bloc*(kriglength-krigoverlap)+1, last=min(nrow(FONKgpointls), (bloc+1)*kriglength-bloc*krigoverlap))
      #fit$shat.GCV <- NA ## otherwise error on 'summarizing' fit
      fitobject$blocmin[bloc+1] <- min(fit$x[, 1]);fitobject$blocmax[bloc+1] <- max(fit$x[, 1])
      fitlist[[bloc+1]] <- fit
    }
    #fitlist[[krigit+1]] <- FONKgpointls ## could be convenient for predict(fitobject)
    fitobject$Kriglist <- fitlist
    fitobject$x <- FONKgpointls[, fittedNames, drop=FALSE]
    fitobject$y <- - FONKgpointls[, ycolname]
  }
  ## A single Krig object or a list of such object is stored in the global 'fitobject' above.
  ## This fitobject should serve as default argument for purefn.
  ## Some diagnostics...
  ## next code should be equiv to
  ##  fitobject$fitted.values <- as.matrix(predict(fitobject), ncol=1)
  ## tmp <- fitobject$fitted.values-fit$y
  ## when there is a single block
  ## differences between predicted and observed values
  tmp <- predict(fitobject, testhull=FALSE)
  fitobject$fitted.values <- tmp ## prediction of kriged points ## removed t(t()) 03/2015
  tmp <- tmp+as.vector(FONKgpointls[, ycolname]) #+ car y/-y
  blackbox.options(RMSpred=sqrt(mean((tmp)^2)))  #note replicated coordinates count twice (as they should)
  llocalst <- paste("RMS residual error from predictions of smoothed points: ", prettynum(blackbox.getOption("RMSpred")),
                    "\n (RMS of variance from all pairs, smoothed or not: ", prettynum(blackbox.getOption("RMSy")),")")
  cat(llocalst, "\n") ## FR->FR for call through migraine, cat() goes in Rout_1 but not on screen. llocalst could deserve to go on screen
  blackbox.options(fitobject=fitobject)
  invisible(fitobject)
}
