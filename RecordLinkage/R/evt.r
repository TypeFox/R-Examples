# evt.r: Functions for applying Extreme Value Theory to Record Linkage

# simplified version of mrl.plot in package ismev
mrl <- function(data, umin = min(data), umax = max(data) - 0.1, nint = 
	round(max(data)-min(data))*20)
{
#
# function to produce empirical mean residual life plot
# as function of threshold.
	x <- xu <- xl <- numeric(nint)
	u <- seq(umin, umax, length = nint)
	for(i in 1:nint) {
		data <- data[data > u[i]]
		x[i] <- mean(data - u[i])
	}
	return(list(x=u,y=x))
}

# Estimation of quantile in pareto distribution
gpdEst <- function(Wdata, thresh=-Inf, quantil=0.95)
{
    gpd=fpot(x=Wdata, threshold=thresh,std.err=FALSE)
    n=length(Wdata)
    scale=gpd$estimate[1]
    shape=gpd$estimate[2]
   	k=length(gpd$exceedances) # number of exceedances over thresh
    x_quantil=thresh+scale/shape*((n/k*(1-quantil))^(-shape) -1)
    # adjust to reasonable value
    if (x_quantil > -scale/shape) x_quantil <- (-scale/shape)
    if (x_quantil < thresh) x_quantil <- (median(c(thresh,max(Wdata))))
    return (x_quantil)
} 




# More efficient implementation of MRL function for weight data.
# Computation is done only for input weights and W-epsilon for every weight
# (where W is a weight unless the lowest one, epsilon is a small number
# (.Machine$double.eps). The graph is linear between these points
#
# Args:
#   W: sorted vector of weights
if(getRversion() >= "2.15.1")  utils::globalVariables(".N")
.computeMRL <- function(W)
{
  # convert weights so that a data.table key can be set on them
  Wtable <- data.table(W=factor(W))

  # This sets W as a key through the backdoor. setkey() does not work here
  # (data.table version 1.6.6) because it reorders the levels alphabetical
  # (numerical order matters here)
  attr(Wtable, "sorted") <- "W"

  # count occurences of every weight
  histogram <- Wtable[,.N,by=W]
  histogram$W <- as.numeric(levels(histogram$W)[histogram$W])

  nVal <- nrow(histogram) # number of distinct weights
  mrlX <- mrlVal <- numeric(2 * nrow(histogram) - 1) # holds x and y values of MRL function

  # compute MRL at values of weights (uneven positions in the result)
  for (i in 1:(nVal - 1))
  {
    thisW <- histogram[i, W]
    # calculate weighted distance to weights greater than the current one
    mrlVal[2 * i - 1] <- sum(histogram[(i+1):nVal, .N * (W - thisW)]) / sum(histogram[(i+1):nVal, .N])
    mrlX[2 * i - 1] <- thisW
  }

  # compute MRL at values of weights minus epsilon (even positions in the result)
  epsilon <- sqrt(.Machine$double.eps)
  for (i in 2:nVal) # do not calculate anything before the first weight
  {
    thisW <- histogram[i, W]
    # calculate weighted distance to weights greater than the current one
    mrlVal[2 * (i - 1)] <- sum(histogram[i:nVal, .N * (W - thisW + epsilon)]) / sum(histogram[i:nVal, .N])
    mrlX[2 * (i - 1)] <- thisW - epsilon
  }
  # complete vector of x-values with highest weight
  mrlX[length(mrlX)] <- histogram[nVal,W]

  list(x=mrlX, y=mrlVal)
}

# MRL plot
plotMRL <- function(rpairs,l = .computeMRL(sort(as.ram((rpairs$Wdata)))))
{
  plot(l$x,l$y,type="l",lty="blank",xlab="Threshold",ylab="MRL")
  # Draw grid
  abline(v=pretty(extendrange(l$x),n=40),h=pretty(extendrange(l$y),n=40),col="lightgray")
  abline(v=pretty(extendrange(l$x),n=8),h=pretty(extendrange(l$y),n=8),col="gray")
  box()
  points(l$x,l$y,type="l")
}


setGeneric(
  name = "getParetoThreshold",
  def = function(rpairs, quantil=0.95, interval=NA) standardGeneric("getParetoThreshold")
)

setMethod(
  f = "getParetoThreshold",
  signature = "RecLinkData",
  definition = function(rpairs, quantil=0.95, interval=NA)
  {
    if (!("RecLinkData" %in% class(rpairs) || "RecLinkResult" %in% class(rpairs)))
      stop(sprintf("Wrong class for rpairs: %s", class(rpairs)))

    if (nrow(rpairs$pairs) == 0)
      stop("No record pairs!")

    if (is.null(rpairs$Wdata))
      stop("No weights in rpairs!")

    if (!is.numeric(quantil))
      stop(sprintf("Illegal type for quantil: %s", class(quantil)))

    if (quantil <0 || quantil > 1)
      stop(sprintf("Illegal value for quantil: %g!", quantil))

    if (!missing(interval) && !is.numeric(interval))
      stop(sprintf("Illegal class for interval: %s!", class(interval)))


    # Choose interval in plot
    l=.computeMRL(sort(rpairs$Wdata))
    if (!is.numeric(interval))
    {
      message("Choose interval for pareto distribution")
      flush.console()
      repeat
      {
        plotMRL(NULL,l=l)
  #      title(main=rpairs$description)
        if (existsFunction("bringToTop")) bringToTop()
        indices=sort(identify(l$x,l$y,n=2,labels=signif(l$x,4)))
        interval=l$x[indices]
        if (length(indices)==0)
          stop("At least the left endpoint of the interval must be chosen!")
        if (length(interval)==1)
          interval=c(interval,max(rpairs$Wdata))
        if (any(rpairs$Wdata > interval[1] & rpairs$Wdata <=interval[2])) break
        message("No data in selected range! Choose a larger interval.")
        flush.console()
      }
      if (existsFunction("bringToTop")) bringToTop(-1)
      if (length(indices)==0)
        stop("At least the left endpoint of the interval must be chosen!")
    }
    # If only left limit has been chosen, leave interval open to the right
    fatTail=rpairs$Wdata[rpairs$Wdata <= interval[2]]
    threshold=gpdEst(fatTail,interval[1],quantil)
    return(as.vector(threshold))
  }
)

setMethod(
  f = "getParetoThreshold",
  signature = "RLBigData",
  definition = function(rpairs, quantil=0.95, interval=NA)
  {

    if (nrow(rpairs@pairs) == 0)
      stop("No record pairs!")

    if (is.null(rpairs@Wdata))
      stop("No weights in rpairs!")

    if (!is.numeric(quantil))
      stop(sprintf("Illegal type for quantil: %s", class(quantil)))

    if (quantil <0 || quantil > 1)
      stop(sprintf("Illegal value for quantil: %g!", quantil))

    if (!missing(interval) && !is.numeric(interval))
      stop(sprintf("Illegal class for interval: %s!", class(interval)))


    # Choose interval from plot
    l=.computeMRL(as.ram(ffsort(rpairs@Wdata)))
    if (!is.numeric(interval))
    {
      message("Choose interval for pareto distribution")
      flush.console()
      repeat
      {
        plotMRL(NULL,l=l)
  #      title(main=rpairs$description)
        if (existsFunction("bringToTop")) bringToTop()
        indices=sort(identify(l$x,l$y,n=2,labels=signif(l$x,4)))
        interval=l$x[indices]
        if (length(indices)==0)
          stop("At least the left endpoint of the interval must be chosen!")
        if (length(interval)==1)
          interval=c(interval,max(rpairs$Wdata))
        # determine if there are weights in the selected range
        weightsInRange <- FALSE
        ffvecapply(
          {
            slice <- rpairs@Wdata[i1:i2]
            if (any(slice > interval[1] & slice <= interval[2]))
            {
              weightsInRange <- TRUE
              .break <- TRUE # breaks the ffvecapply loop
            }
          }
          , X = rpairs@Wdata)
        if (weightsInRange) break
        message("No data in selected range! Choose a larger interval.")
        flush.console()
      }
      if (existsFunction("bringToTop")) bringToTop(-1)
      if (length(indices)==0)
        stop("At least the left endpoint of the interval must be chosen!")
    }
    fatTailInd <- ffvecapply(
      {
        slice <- rpairs@Wdata[i1:i2]
        which(slice <= interval[2]) + i1 - 1
      }, X = rpairs@Wdata, RETURN = TRUE, CFUN = "c")
    fatTail=rpairs@Wdata[fatTailInd]
    threshold=gpdEst(fatTail,interval[1],quantil)
    return(as.vector(threshold))
  }
)
