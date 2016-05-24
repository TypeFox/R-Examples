# ptls must be ordered...
findReplicates <- function(ptls) {
  if ("Migraine" %in% blackbox.getOption("usedBy")) {
    ycolname <- blackbox.getOption("ycolname")
    fittedLoci <- blackbox.getOption("respCols")
    fittedNames <- blackbox.getOption("fittedNames")
  } else {
    ycolname <- attr(ptls,"ycolname")
    fittedLoci <- NULL
    fittedNames <- attr(ptls,"fittedNames")
  }
  if("automatedCleaning" %innc% blackbox.getOption("miscOptions")) {
    ## then we should not have suspicious replicates here
    ## we can check that automated cleaning worked:
    absD <- (apply(abs(diff(t(t(ptls)), lag=1)), 1, max)) ## SSD no longer has rownames; t(t(.)) : conversion of data.frame to matrix...
    nullabsD <- (absD==0)
    if (length(which(nullabsD))>0) {
      message.redef("(!) Some likelihood estimates  from independent replicates appear identical. ")
      stop.redef("    this should not occur in selectFn() after automated cleaning. I exit. ")
    }
  } else {
    ## there is no point in checking again as a warning was already issued before
  }
  ## the check on replicates x coordinates
  ptlsx <- ptls[, fittedNames, drop=FALSE] ## selection of relevant x columns
  doublonsv <- numeric();
  nr <- nrow(ptlsx)
  #  pif <- ptlsx[-nrow(ptlsx), , drop=FALSE]
  #	paf <- ptlsx[-1, , drop=FALSE]
  #    absD <- (apply(abs(pif-paf), 1, max)) ## absD no longer has rownames
  absD <- (apply(abs(diff(t(t(ptlsx)), lag=1)), 1, max)) ## SSD no longer has rownames; t(t(.)) : conversion of data.frame to matrix...
  nullabsD <- (absD==0)
  ptlnames <- rownames(ptlsx)
  pifnames <- ptlnames[-nrow(ptlsx)]
  pafnames <- ptlnames[-1]
  pifreplnames <- pifnames[nullabsD]
  pafreplnames <- pafnames[nullabsD]
  errorcheck <- intersect(pifreplnames, pafreplnames)
  ## no automatedCleaning below...
  if (length(errorcheck)>0) { ## check on x coordinates only
    message.redef("(!) Some parameter point(s) occur more than twice.")
    message.redef("    This typically results from appending twice or more a single likelihood estimate to the pointls file.")
    message.redef("    Look in particular for replicates of the following cordinates in the pointls file:")
    apply(ptlsx[errorcheck, , drop=FALSE], 1, function(v) {message.redef(canonize(v)$canonVP)})
    stop.redef("Exiting as a result of some parameter point(s) occurring more than twice.")
  }
  sorteddoublonsnames <- as.character(sort(as.numeric(c(pifreplnames, pafreplnames)))) ## length is twice the nbr of repls
  if (length(sorteddoublonsnames)>0) {
    doublonsFirstsnames <- sorteddoublonsnames[seq(1, length(sorteddoublonsnames), 2)] #rownames
    doublonsSecondsnames <- sorteddoublonsnames[seq(2, length(sorteddoublonsnames), 2)] #rownames
    replyFirsts <- ptls[doublonsFirstsnames, c(fittedLoci, ycolname), drop=FALSE]
    replySeconds <- ptls[doublonsSecondsnames, c(fittedLoci, ycolname), drop=FALSE]
    pureRMSE <- sqrt(mean((replyFirsts[, ycolname]-replySeconds[, ycolname])^2)/2)
    ## gets the doublonsFirsts rownames and order:
    doublonsymeans <- (replyFirsts+replySeconds)/2 ## columns c(fittedLoci, ycolname) ## nrows is the number of replicates
    singletonsnames <- ptlnames %w/o% sorteddoublonsnames
  } else {
    doublonsFirstsnames <- character(0)
    doublonsSecondsnames <- character(0)
    doublonsymeans <- numeric(0) ## columns c(fittedLoci, ycolname) ## nrows is the number of replicates
    singletonsnames <- ptlnames
    pureRMSE <- NA
  }
  tmp <- list(ptlsx=ptlsx, doublonsFirstsnames=doublonsFirstsnames, doublonsSecondsnames=doublonsSecondsnames,
            singletonsnames=singletonsnames, pureRMSE=pureRMSE, doublonsymeans=doublonsymeans)
  class(tmp) <- c(class(tmp), "replicatesInfoClass")
  return(tmp)
}
