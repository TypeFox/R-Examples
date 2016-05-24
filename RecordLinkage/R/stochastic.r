# methods for the "classical" stochastic record linkage as in the work of
# Fellegi and Sunter: F
# Fellegi, Ivan; Sunter, Alan (December 1969). "A Theory for Record Linkage". 
# Journal of the American Statistical Association 64 (328): pp. 1183-1210.


# backend function to calculate weights for a matrix of comparison patterns
# returns list with components:
#   M:  m-probabilites of distinct binary patterns
#   U:  u-probabilities of distinct binary patterns
#   W:  weights of distinct binary patterns
#   Wdata: weights of actual data
.fsWeightsBackend <- function(patterns, m, u, cutoff)
{
  # calculate weights for set of distinct binary patterns
  distinctPat <- t(bincombinations(ncol(patterns)))
  logM <- colSums(log(m * distinctPat + (1-m) * (1-distinctPat), base=2))
  logU <- colSums(log(u * distinctPat + (1-u) * (1-distinctPat), base=2))
  W <- logM -logU
  # converting to a matrix speed up the following operations
  patterns <- as.matrix(patterns, rownames.force=FALSE)
  # convert NA to 0, string comparison values to 0 or one
  for (colInd in 1:ncol(patterns))
  {
    patterns[is.na(patterns[,colInd]), colInd] <- 0
    patterns[patterns[,colInd] >= cutoff, colInd] <- 1
    patterns[patterns[,colInd] < cutoff, colInd] <- 0
  }
  # for every actual pattern get its position in the table of distinct patterns
  # (and in the vector of distinct weights)
  # this turned out to be faster than to do the calculation for every
  # row in patterns
  counts <- countpattern(patterns, matching=TRUE)
  Wdata <- W[counts$matching]
  list(M = 2^logM, U = 2^logU, W = W, Wdata = Wdata)
}





setGeneric(
  name = "fsWeights",
  def = function(rpairs, ...) standardGeneric("fsWeights")
)

setMethod(
  f = "fsWeights",
  signature = "RecLinkData",
  definition = function(rpairs, m = 0.95, u = rpairs$frequencies, cutoff = 1)
  {
    if (nrow(rpairs$pairs) == 0)
      stop("No record pairs!")

    if (!is.numeric(m))
      stop(sprintf("Illegal type for m: %s", class(m)))

    if (any(m < 0 | m > 1))
      stop(sprintf("Illegal value for m: %s!", paste(m, collapse=", ")))

    if (!is.numeric(u))
      stop(sprintf("Illegal type for u: %s", class(u)))

    if (any(u < 0 | u > 1))
      stop(sprintf("Illegal value for u: %s!", paste(u, collapse=", ")))

    # check condition u <= m
    if(any(u > m))
      stop("m-probability must be greater than or equal to u-probability for every attribute!")

    pairs <-rpairs$pairs[,-c(1,2,ncol(rpairs$pairs)), drop=FALSE]
    result <- .fsWeightsBackend(pairs, m=m, u=u, cutoff=cutoff)
    rpairs$M <- result$M
    rpairs$U <- result$U
    rpairs$W <- result$W
    rpairs$Wdata <- result$Wdata
    rpairs
  }
)

setMethod(
  f = "fsWeights",
  signature = "RLBigData",
  definition = function (rpairs, m=0.95, u=getFrequencies(rpairs),
    cutoff=1, withProgressBar = (sink.number()==0))
  {
    if (withProgressBar)
      pgb <- txtProgressBar(0, nrow(rpairs@pairs))
    rpairs@Wdata <- ffrowapply(
      {
        if (withProgressBar) setTxtProgressBar(pgb, i2)
        .fsWeightsBackend(rpairs@pairs[i1:i2,3:(ncol(rpairs@pairs)-1)],
          m=m, u=u, cutoff=cutoff)$Wdata
      }, X=rpairs@pairs, RETURN = TRUE, RETCOL=NULL, VMODE="double")
    if (withProgressBar) close(pgb)
    rpairs@WdataInd <- fforder(rpairs@Wdata)
    rpairs
  }
) # end of setMethod



# TODO: my / ny erlauben
fsClassify <- function(...) emClassify(...)


#setGeneric(
#  name = "fsClassify",
#  def = function(rpairs, ...) standardGeneric("fsClassify")
#)
#
#
#
## classification works the same as for EM weights
#setMethod(
#  f = "fsClassify",
#  signature = "RecLinkData",
#  definition = function(rpairs, ...) emClassify(rpairs, ...)
#)
#
#
#
#setMethod(
#  f = "fsClassify",
#  signature = "RLBigData",
#  definition = function (rpairs, threshold.upper,
#                        threshold.lower=threshold.upper, m=0.95,
#                        u=getFrequencies(rpairs), withProgressBar = (sink.number()==0),
#                        cutoff=1)
#  {
#    if (!is.numeric(threshold.upper))
#      stop(sprintf("Illegal type for threshold.upper: %s", class(threshold.upper)))
#
#    if (!is.numeric(threshold.lower))
#      stop(sprintf("Illegal type for threshold.lower: %s", class(threshold.lower)))
#
#    if (threshold.upper < threshold.lower)
#      stop(sprintf("Upper threshold %g lower than lower threshold %g",
#      threshold.upper, threshold.lower))
#
#    if(!isIdCurrent(rpairs@con)) stop(paste("Invalid SQLite connection in rpairs!",
#      "See '?saveRLObject' on how to make persistant copies of such objects."))
#
#    if (dbExistsTable(rpairs@con, "Wdata"))
#    {
#      query <- "select id1, id2 from Wdata where W >= :upper"
#      links <- dbGetPreparedQuery(rpairs@con, query, data.frame(upper = threshold.upper))
#      query <- "select id1, id2 from Wdata where W < :upper and W >= :lower"
#      possibleLinks <- dbGetPreparedQuery(rpairs@con, query,
#        data.frame(upper = threshold.upper, lower = threshold.lower))
#      nPairs <- dbGetQuery(rpairs@con, "select count(*) as c from Wdata")$c
#    } else
#    {
#
#      if (withProgressBar)
#      {
#        expPairs <- getExpectedSize(rpairs)
#        pgb <- txtProgressBar(max=expPairs)
#      }
#      nPairs <- 0
#      n <- 10000
#      links <- matrix(nrow=0, ncol=2)
#      possibleLinks <- matrix(nrow=0, ncol=2)
#      on.exit(clear(rpairs))
#      rpairs <- begin(rpairs)
#      while(nrow(slice <- nextPairs(rpairs, n)) > 0)
#      {
#        Wdata <- .fsWeightsBackend(slice[,-c(1,2,ncol(slice)),drop=FALSE],
#          m=m, u=u, cutoff=cutoff)$Wdata
#        if (any(is.na(Wdata)))
#          warning("Some weights have illegal values. Check error rate and frequencies!")
#        links <- rbind(links, as.matrix(slice[Wdata >= threshold.upper,1:2]))
#        possibleLinks <- rbind(possibleLinks,
#          as.matrix(slice[Wdata >= threshold.lower & Wdata < threshold.upper, 1:2]))
#        nPairs <- nPairs + nrow(slice)
#        if (withProgressBar)
#        {
#          setTxtProgressBar(pgb, nPairs)
#          flush.console()
#        }
#      }
#      if (withProgressBar) close(pgb)
#    }
#
#    new("RLResult", data = rpairs, links = as.matrix(links),
#      possibleLinks = as.matrix(possibleLinks),
#      nPairs = nPairs)
#  }
#) # end of setMethod
#
#

