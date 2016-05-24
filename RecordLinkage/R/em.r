# moved to em-methods.r:
#
# emWeights <- function (rpairs, cutoff=0.95,...)
#
#emClassify <- function (rpairs,threshold.upper=Inf, 
#                        threshold.lower=threshold.upper,my=Inf, ny=Inf)


if(getRversion() >= "2.15.1")  utils::globalVariables(c("is_match", "W",
  "nMatch", "nAll", "nNonMatch", "Wdata"))

setGeneric(
  name = "optimalThreshold",
  def = function(rpairs, my=NaN, ny=NaN) standardGeneric("optimalThreshold")
)

setMethod(
  f = "optimalThreshold",
  signature = "RecLinkData",
  definition = function (rpairs, my=NaN, ny=NaN)
  {
    if (!("RecLinkData" %in% class(rpairs) || "RecLinkResult" %in% class(rpairs)))
      stop(sprintf("Wrong class for rpairs: %s", class(rpairs)))

    if (nrow(rpairs$pairs) == 0)
      stop("No record pairs!")

    if (is.null(rpairs$Wdata))
      stop("No weights in rpairs!")

    if (!is.numeric(my))
      stop(sprintf("Illegal type for my: %s", class(my)))
    if (!missing(my) && (my < 0 || my > 1))
      stop(sprintf("Illegal value for my: %g", my))

    if (!is.numeric(ny))
      stop(sprintf("Illegal type for ny: %s", class(ny)))
    if (!missing(ny) && (ny < 0 || ny > 1))
      stop(sprintf("Illegal value for ny: %g", ny))

    # remove missing values for matching status
    indMissing <- which(is.na(rpairs$pairs$is_match))
    if(length(indMissing)==nrow(rpairs$pairs))
      stop("Only pairs with unknown status in rpairs!")

    if (length(indMissing >0)) rpairs <- rpairs[-indMissing]

    # grouped by weight, count number of non-matches and matches
    errStatFun <- function(x)
    {
        nMatch <- sum(x)
        nNonMatch <- length(x) - sum(x)
        list(nMatch = nMatch, nNonMatch = nNonMatch)
    }
    dt <- data.table(is_match=rpairs$pairs$is_match,
      Wdata=match(rpairs$Wdata, sort(rpairs$W)), key="Wdata")
    errStat <- dt[,errStatFun(is_match), by=Wdata]

    # Count false negatives and positives per weight. The vectors correspond to
    # unique(sort(rpairsWdata)). Each position holds
    #   for FN: the number of false negatives if all record pairs with
    #           smaller or equal weight are classified as non-links
    #   for FP: the number of false positives if all record pairs with
    #           greater or equal weight are classified as links
    FN <- cumsum(errStat$nMatch)
    FP <- rev(cumsum(rev(errStat$nNonMatch)))

    # Construct error measures therefrom. Because thresholds define
    # right-closed and left-open intervals, the range of weights has to be extended
    # by a value greater than all exisiting weights (the case that no pairs are
    # classified as links at all). The error rates are extended by zeros.
    alphaErr <- c(0, FN/sum(errStat$nMatch))
    betaErr <- c(FP/sum(errStat$nNonMatch), 0)
    accuracy <- (nrow(rpairs$pairs) - (FP + FN)) / nrow(rpairs$pairs)

    classWeights <- c(sort(unique(rpairs$Wdata)), Inf)

    # set thresholds

    # no error bounds given: maximize accuracy
    if (missing(my) && missing(ny))
      return(as.numeric(classWeights[which.max(accuracy)]))

    # only bound for alpha error given: minimize beta error under constraint that
    # bound for alpha error is met
    if (!missing(ny))
    {
      min_ind <- which.min(betaErr[alphaErr<=ny])
      return(as.numeric(classWeights[alphaErr<=ny][min_ind]))
    }

    # only bound for beta error given: minimize alpha error under constraint that
    # bound for beta error is met
    if (!missing(my))
    {
      min_ind <- which.min(alphaErr[betaErr<=my])
      return(as.numeric(classWeights[betaErr<=my][min_ind]))
    }
  }
)


setMethod(
  f = "optimalThreshold",
  signature = "RLBigData",
  definition = function (rpairs, my=NaN, ny=NaN)
  {
    if (nrow(rpairs@pairs) == 0)
      stop("No record pairs!")

    if (!hasWeights(rpairs))
      stop("No weights in rpairs!")

    if (!is.numeric(my))
      stop(sprintf("Illegal type for my: %s", class(my)))
    if (!missing(my) && (my < 0 || my > 1))
      stop(sprintf("Illegal value for my: %g", my))

    if (!is.numeric(ny))
      stop(sprintf("Illegal type for ny: %s", class(ny)))
    if (!missing(ny) && (ny < 0 || ny > 1))
      stop(sprintf("Illegal value for ny: %g", ny))


    # remove missing values for matching status
    # TODO
#    indMissing <- which(is.na(rpairs$pairs$is_match))
#    if(length(indMissing)==nrow(rpairs$pairs))
#      stop("Only pairs with unknown status in rpairs!")

#    if (length(indMissing >0)) rpairs <- rpairs[-indMissing]

    wMatch <- ffdf(W = rpairs@Wdata, is_match = rpairs@pairs$is_match)

    pgb <- txtProgressBar(0, nrow(wMatch))
    summaryTable <- ffrowapply(
      {
        setTxtProgressBar(pgb, i2)
        slice <- data.table(wMatch[i1:i2,])
        # format weight as factor with as many digits as reasonable
        slice$W <- factor(format(slice$W, digits = 20))
        # count number of pairs and matches, ignore pairs with match status NA
        slice[,list(nMatch=sum(is_match, na.rm=TRUE),
          nAll=sum(!is.na(is_match))), by=W]
      }, X = wMatch, RETURN = TRUE, CFUN = "rbind")
    close(pgb)
    setkeyv(summaryTable, "W")
    summaryTable <- summaryTable[,list(nMatch=sum(nMatch), nAll=sum(nAll)), by=W]
    summaryTable$nNonMatch <- summaryTable[,nAll - nMatch]
#browser()
    # sort summaryTable in numeric order (data.table converts to factor and
    # performs sorting by character string

    summaryTable <- summaryTable[order(as.numeric(levels(W)[W])),]

    FN <- summaryTable[,cumsum(nMatch)]
    FP <- summaryTable[,rev(cumsum(rev(nNonMatch)))]


    # Construct error measures therefrom. Because thresholds define
    # right-closed and left-open intervals, the range of weights has to be extended
    # by a value greater than all exisiting weights (the case that no pairs are
    # classified as links at all). The error rates are extended by zeros.
    alphaErr <- c(0, FN/sum(summaryTable$nMatch))
    betaErr <- c(FP/sum(summaryTable$nNonMatch), 0)
    accuracy <- (nrow(rpairs@pairs) - (FP + FN)) / nrow(rpairs@pairs)

    # set thresholds
    # the returned threshold is always the mean of the calculated threshold
    # (equal to a weight in rpairs) and the next lowest weight. This prevents
    # unexpected classification results due to rounding errors

    classWeights <- c(as.numeric(summaryTable[,levels(W)[W]]), Inf)
    # no error bounds given: maximize accuracy
    if (missing(my) && missing(ny))
    {
      max_ind <- which.max(accuracy)
      return(mean(as.numeric(classWeights[c(max_ind - 1, max_ind)])))
    }
    # only bound for alpha error given: minimize beta error under constraint that
    # bound for alpha error is met
    if (!missing(ny))
    {
      min_ind <- which.min(betaErr[alphaErr<=ny])
      return(mean(as.numeric(classWeights[alphaErr<=ny][c(min_ind - 1, min_ind)])))
    }

    # only bound for beta error given: minimize alpha error under constraint that
    # bound for beta error is met
    if (!missing(my))
    {
      min_ind <- which.min(alphaErr[betaErr<=my])
      return(mean(as.numeric(classWeights[betaErr<=my][c(min_ind - 1, min_ind)])))
    }
  }
)

