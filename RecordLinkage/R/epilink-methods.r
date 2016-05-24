# epilink-methods.r: Functions for the Epilink matching procedure
# See Continiero et al.: The EpiLink record linkage software, in:
# Methods of Information in Medicine 2005, 44(1):66-71.

setGeneric(
  name = "epiClassify",
  def = function(rpairs, threshold.upper, threshold.lower=threshold.upper, ...)
    standardGeneric("epiClassify")
)

setMethod(
  f = "epiClassify",
  signature = "RecLinkData",
  definition = function (rpairs,threshold.upper, 
                        threshold.lower=threshold.upper)
  {    
  
    if (!("RecLinkData" %in% class(rpairs) || "RecLinkResult" %in% class(rpairs)))
      stop(sprintf("Wrong class for rpairs: %s", class(rpairs)))
  
    if (nrow(rpairs$pairs) == 0)
      stop("No record pairs!")
  
    if (is.null(rpairs$Wdata))
      stop("No weights in rpairs!")
  
    if (!is.numeric(threshold.upper))
      stop(sprintf("Illegal type for threshold.upper: %s", class(threshold.upper)))
  
    if (!is.numeric(threshold.lower))
      stop(sprintf("Illegal type for threshold.lower: %s", class(threshold.lower)))
  
    if (threshold.upper < threshold.lower)
      stop(sprintf("Upper threshold %g lower than lower threshold %g",
      threshold.upper, threshold.lower))
      
    prediction=rep("P",nrow(rpairs$pairs))
    prediction[rpairs$Wdata>=threshold.upper]="L"
    prediction[rpairs$Wdata<threshold.lower]="N"
    
    ret=rpairs # keeps all components of rpairs
    ret$prediction=factor(prediction,levels=c("N","P","L"))
    ret$threshold=threshold.upper
    class(ret)="RecLinkResult"
    return(ret)
  }
) # end of setMethod

setMethod(
  f = "epiClassify",
  signature = "RLBigData",
  definition = function (rpairs,threshold.upper, 
                        threshold.lower=threshold.upper, e=0.01, 
                        f=getFrequencies(rpairs), withProgressBar = (sink.number()==0))
  {    
    if (!is.numeric(threshold.upper))
      stop(sprintf("Illegal type for threshold.upper: %s", class(threshold.upper)))

    if (!is.numeric(threshold.lower))
      stop(sprintf("Illegal type for threshold.lower: %s", class(threshold.lower)))

    if (threshold.upper < threshold.lower)
      stop(sprintf("Upper threshold %g lower than lower threshold %g",
      threshold.upper, threshold.lower))

    .ffWeightClassify(rpairs, threshold.upper, threshold.lower)
  }
) # end of setMethod


setGeneric(
  name = "epiWeights",
  def = function(rpairs, e=0.01, f=getFrequencies(rpairs), ...)
    standardGeneric("epiWeights")
)

setMethod(
  f = "epiWeights",
  signature = "RecLinkData",
  definition = function (rpairs, e=0.01, f=rpairs$frequencies)
  {
    # check for erronous input

    if (!("RecLinkData" %in% class(rpairs) || "RecLinkResult" %in% class(rpairs)))
      stop(sprintf("Wrong class for rpairs: %s", class(rpairs)))

    if (nrow(rpairs$pairs) == 0)
      stop("No record pairs!")

    if (!is.numeric(f))
      stop(sprintf("Illegal type for f: %s", class(f)))

    if (any(f <=0 | f > 1))
      stop(sprintf("Illegal value for f: %s!", paste(f, collapse=", ")))

    if (!is.numeric(e))
      stop(sprintf("Illegal type for e: %s", class(f)))

    if (any(e <0 | e >= 1))
      stop(sprintf("Illegal value for e: %s!", paste(e, collapse=", ")))

    # check condition e <= 1-f, otherwise illegal weights can occur
    if(any(e > 1-f))
      stop("Condition e <= 1-f does not hold, adjust error rate!")

    # leave out ids and matching status
    # weights are computed on transposed matrix of comparison patterns
    # (each column is an observation) to allow vectorization
    pairs=t(as.matrix(rpairs$pairs[,-c(1,2,ncol(rpairs$pairs))]))
    pairs[is.na(pairs)]=0
    w=log((1-e)/f, base=2)

    # recycle w (in order to get the sum of attribute weights right)
    w <- numeric(nrow(pairs)) + w

    # multiply each pattern with the epilink weights
    # the weight for each pattern is the sum of its attribute weights
    S <- colSums(pairs * w)/sum(w)

    if (any(is.na(S) | S < 0 | S > 1))
      warning("Some weights have illegal values. Check error rate and frequencies!")
    rpairs$Wdata=S
    return(rpairs)
  }
)

.epiWeightsBackend <- function(patterns, e, f)
{
  # convert NA to 0
  for (colInd in 1:ncol(patterns))
  {
    patterns[is.na(patterns[,colInd]), colInd] <- 0
  }
  # converting to a matrix speed up the following operations
  patterns <- t(as.matrix(patterns, rownames.force=FALSE))
  w=log((1-e)/f, base=2)

  # recycle w (in order to get the sum of attribute weights right)
  w <- numeric(nrow(patterns)) + w
  S <- colSums(patterns * w) / sum(w)
  if (any(is.na(S) | S < 0 | S > 1))
    warning("Some weights have illegal values. Check error rate and frequencies!")
  S
}

setMethod(
  f = "epiWeights",
  signature = c("RLBigData", "ANY", "ANY"),
  definition = function (rpairs, e=0.01, f=getFrequencies(rpairs),
      withProgressBar = (sink.number()==0))
  {

    rpairs@Wdata <- ffrowapply(.epiWeightsBackend(rpairs@pairs[i1:i2,3:(ncol(rpairs@pairs)-1), drop=FALSE],
    e=e, f=f), X=rpairs@pairs, RETURN = TRUE, RETCOL=NULL, VMODE="double")
    rpairs@WdataInd <- fforder(rpairs@Wdata)
    rpairs
  }
) # end of setMethod
