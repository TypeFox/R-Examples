cleanResults <- function(rawTMLE, verbose=FALSE) {
  ## - - - - - - - - - - - 
  ## CLEANING THE RESULTS
  ## - - - - - - - - - - -
  TMLE <- rawTMLE

  ## - - - - - - - - - - - - - - - -
  ## failures due to 'NA' in 'obs'
  ## - - - - - - - - - - - - - - - -

  msg <- "contains at least one 'NA'"
  idx <- sapply(TMLE,
                function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
  if (sum(idx)) {
    verbose && enter(verbose, "at least one 'NA': ", sum(idx))
    ## get names of an instance of guilty genes
    verbose && enter(verbose, paste("example: ", names(TMLE[min(which(idx))])))
    verbose && exit(verbose)
  }
  TMLE <- TMLE[!idx]

  ## - - - - - - - -
  ## other failures
  ## - - - - - - - -

  ##
  ## 'simpleError' 
  ##

  msg <- "simpleError"
  idx <- sapply(TMLE, inherits, msg)
  if (FALSE) {
    verbose && enter(verbose, paste("simpleError: ", sum(idx)))
    ## get names of an instance of guilty genes
    verbose && cat(verbose, "example: ", names(TMLE[min(which(idx))]))
    verbose && exit(verbose)
  }
  err <- sapply(TMLE[idx], as.character) 

  ##
  ## compare to known error messages
  ##

  knownMsgs <- c("Error in findInterval",
                 "impossible",
                 "Parsimonious conditional simulation of X given W failed",
                 "sLeftA",
                 "Range of argument",
                 "Error in if \\(ps",
                 "must be positive",
                 "Missing element in argument 'lib'")
  msgs <- lapply(knownMsgs, grep, err)
  names(msgs) <- knownMsgs
  
  if (sum(sapply(msgs, length)) < length(err)) {
    stop("some messages are not in 'knownMsgs'")
  }

  ##
  ##  Missing element in argument 'lib'
  ##

  msg <- "Missing element in argument 'lib'"
  idx <- sapply(TMLE,
                function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
  if (sum(idx)) {
    verbose && enter(verbose, paste("Missing element in argument 'lib': ", sum(idx)))
    ## get names of an instance of guilty genes
    verbose && cat(verbose, "example: ", names(TMLE[min(which(idx))]))
    verbose && exit(verbose)
  }
  TMLE <- TMLE[!idx]

  ##
  ##  impossible... 
  ##

  msg <- "impossible"
  idx <- sapply(TMLE,
                function(xx){length(xx)==2 && length(grep(msg, as.character(xx)))})
  if (sum(idx)) {
    verbose && enter(verbose, paste("impossible: ", sum(idx)))
    ## get names of an instance of guilty genes
    verbose && cat(verbose, "example: ", names(TMLE[min(which(idx))]))
    verbose && exit(verbose)
  }
  TMLE <- TMLE[!idx]

  ##
  ##  Error in if (ps[3]==0)
  ##

  msg <- "Error in if \\("
  idx <- sapply(TMLE,
                function(xx){length(xx)==2 && length(grep(msg, as.character(xx)))})
  if (sum(idx)) {
    verbose && enter(verbose, paste("Error in if \\(: ", sum(idx)))
    ## get names of an instance of guilty genes
    verbose && cat(verbose, "example: ", names(TMLE[min(which(idx))]))
    verbose && exit(verbose)
  }
  TMLE <- TMLE[!idx]

  ##
  ##  Parameter 'sigma2' must be positive...
  ##

  msg <- "must be positive"
  idx <- sapply(TMLE,
                function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
  if (sum(idx)) {
    verbose && enter(verbose, paste("must be positive: ", sum(idx)))
    ## get names of an instance of guilty genes
    verbose && cat(verbose, "example: ", names(TMLE[min(which(idx))]))
    verbose && exit(verbose)
  }
  TMLE <- TMLE[!idx]

  ##
  ## 'Exception: Range of argument 'weightsW' is out of range 
  ##

  msg <- "Range of argument 'weightsW'"
  idx <- sapply(TMLE,
                function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
  if (sum(idx)) {
    verbose && enter(verbose, paste("Range of argument 'weightsW': ", sum(idx)))
    ## get names of an instance of guilty genes
    verbose && cat(verbose, "example: ", names(TMLE[min(which(idx))]))
    verbose && exit(verbose)
  }
  TMLE <- TMLE[!idx]

  ##
  ## <simpleError in findInterval(V[xx], cumsum(ps)): 'vec' contient des valeurs manquantes (NAs)>
  ##

  msg <- "Error in findInterval"
  idx <- sapply(TMLE,
                function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
  if (sum(idx)) {
    verbose && enter(verbose, paste("Error in findInterval: ", sum(idx)))
    ## get names of an instance of guilty genes
    verbose && cat(verbose, "example: ", names(TMLE[min(which(idx))]))
    verbose && exit(verbose)
  }
  TMLE <- TMLE[!idx]

  ##
  ##  Parsimonious conditional simulation of X given W failed...
  ##

  msg <- "Parsimonious conditional simulation of X given W failed"
  idx <- sapply(TMLE,
                function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
  if (sum(idx)) {
    verbose && enter(verbose, paste("Parsimonious conditional simulation of X given W failed: ", sum(idx)))
    ## get names of an instance of guilty genes
    verbose && cat(verbose, "example: ", names(TMLE[min(which(idx))]))
    verbose && exit(verbose)
  }
  TMLE <- TMLE[!idx]

  ##
  ## <simpleError in sLeftA[bases[, 1]]: type 'list' d'indice incorrect>
  ##

  msg <- "Error in sLeftA\\[bases"
  idx <- sapply(TMLE,
                function(xx){length(xx)==2 && length(grep(msg, as.character(xx)))})
  if (sum(idx)) {
    verbose && enter(verbose, paste("Error in sLeftA\\[bases: ", sum(idx)))
    ## get names of an instance of guilty genes
    verbose && cat(verbose, "example: ", names(TMLE[min(which(idx))]))
    verbose && exit(verbose)
  }  
  TMLE <- TMLE[!idx]
  TMLE <- TMLE[!sapply(TMLE, is.null)]
  return(TMLE)
}




enrich <- function ## calculates and plots the degree of over-résentation
(p.values, chr, thr=1e-5, ...) {
    p <- length(p.values)
    rk <- rank(1-p.values)
    rkChr <- tapply(rk, chr, FUN=mean)
    isRej <- (rk>length(rk)-nbRej)
    rejByChr <- tapply(isRej, chr, sum)
    nbByChr <- table(chr)
    phg <- 1-phyper(rejByChr, nbByChr, p-nbByChr, sum(isRej))
    lphg <- -log10(phg)
    ww <- which(is.finite(lphg))
    ymax <- max(lphg[ww], na.rm=TRUE)
    ylim <- c(0, ymax+1)
    lphg[-ww] <- ymax
    cols <- rep("purple", length(lphg))
    cols[-ww] <- "pink"
    bp <- barplot(lphg, ..., col=cols, ylim=ylim)
    sa <- sapply(bp[-ww], points, ymax, pch=2)
}
