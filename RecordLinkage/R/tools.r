# tools.r: various utility functions

# estimate the number of record pairs
setGeneric(
  name = "getExpectedSize",
  def = function(object, ...) standardGeneric("getExpectedSize")
)

setMethod(
  f = "getExpectedSize",
  signature = "data.frame",
  definition = function(object, blockfld=list())
  {
    if(!is.list(blockfld)) blockfld = list(blockfld)
    con <- dbConnect(dbDriver("SQLite"))
    dbWriteTable(con, "data", object)
    nData <- nrow(object)
    nAll <- nData * (nData - 1) / 2
    if (length(blockfld)==0) return(nAll)
    coln <- make.db.names(con, colnames(object))

    # ergibt Wahrscheinlichkeit, dass mit gegebenen Blockingfeldern
    # ein Paar nicht gezogen wird
    blockelemFun <- function(blockelem)
    {
      if(is.character(blockelem)) blockelem <- match(blockelem, colnames(object))
      freq <- dbGetQuery(con,
        sprintf("select count(*) as c from data group by %s having c > 1 and %s",
          paste("\"", coln[blockelem], "\"", sep="", collapse=", "),
          paste(
            sapply(coln[blockelem], sprintf, fmt = "\"%s\" is not null"),
            collapse = " and "
          )
        )
      )
      if (nrow(freq)==0) return(1)
      1 - (sum(sapply(freq,  function(x) x * (x-1) /2)) / nAll)
    }
    res <- nAll * (1-prod(sapply(blockfld, blockelemFun)))
    dbDisconnect(con)
    round(res)
  }
)

setMethod(
  f = "getExpectedSize",
  signature = "RLBigDataDedup",
  definition = function(object)
  {
    # With ff pairs, this function is actually not necessary. For compatibility,
    # the actual number of record pairs is returned
    nrow(object@pairs)
  }
)

setMethod(
  f = "getExpectedSize",
  signature = "RLBigDataLinkage",
  definition = function(object)
  {
    # With ff pairs, this function is actually not necessary. For compatibility,
    # the actual number of record pairs is returned
    nrow(object@pairs)
  }
)



# Subscript operator for Record Linkage objects
"[.RecLinkData" <- function(x,i)
{
  ret <- x
  ret$pairs <- x$pairs[i,]
  ret$Wdata <- x$Wdata[i]
  ret
}

"[.RecLinkResult" <- function(x,i)
{
  ret <- "[.RecLinkData"(x, i)
  ret$prediction <- x$prediction[i]
  ret
}

"[.RLBigData" <- function(x,i)
{
  ret <- x
  # this should be recoded more efficiently!
  ret@pairs <- as.ffdf(x@pairs[i,])
  ret@Wdata <- ff(x@Wdata[i])
  ret@WdataInd <- ff(x@WdataInd[i])
  ret
}

"[.RLResult" <- function(x,i)
{
  ret <- x
  ret@data <- x@data[i]
  ret@prediction <- ff(x@prediction[i])
  ret
}

# append two data sets into one

setGeneric(
  name = "%append%",
  def = function(x, y) standardGeneric("%append%")
)

setMethod(
  f = "%append%",
  signature = c("RecLinkData", "RecLinkData"),
  definition = function(x, y)
  {
    ret <- x
    ret$pairs <- rbind(x$pairs, y$pairs)
    # look if weights are stored in both, one or none of the objects
    withoutWeights <- is.null(x$Wdata) + is.null(y$Wdata)
    if(withoutWeights==1) warning("Weights in one of the objects only!")
    if(withoutWeights==0) ret$Wdata <- c(x$Wdata, y$Wdata)
    ret
  }
)


setMethod(
  f = "%append%",
  signature = c("RecLinkResult", "RecLinkResult"),
  definition = function(x, y)
  {
    ret <- selectMethod("%append%", c("RecLinkData", "RecLinkData"))(x, y)
    ret$prediction <- c(x$prediction, y$prediction)
    class(ret$prediction) <- "factor"
    levels(ret$prediction) <- c("N", "P", "L")
    ret
  }
)

# determine position with smallest x > or >= given threshold by binary search
# x must be sorted in ascending order or order be given in o
.searchThreshold <- function(x, threshold, inclusive=TRUE, o=NA)
{
  if (missing(o)) getX <- function(i) x[i] else getX <- function(i) x[o[i]]
  # choose appropriate comparison operator
  compFun <- if (inclusive) `>=` else `>`
  # return NA if all data are below threshold <=> the largest element of x does
  # not satisfy the condition
  if (!compFun(getX(length(x)), threshold)) return(NA)
  # if smallest value satisfies the condition, all of x do
  if (compFun(getX(1), threshold)) return(1)
  else
  {
    lower <- 1
    upper <- length(x)
    minInd <- round(length(x) / 2)
    # search until current value satisfies condition and its lower neighbour
    # does not
    while(!(compFun(getX(minInd), threshold) && !compFun(getX(minInd-1), threshold)))
    {
#      message(minInd)
      if (!compFun(getX(minInd), threshold)) # zu far left
      {
        lower <- minInd
        minInd <- ceiling((minInd + upper) / 2)
      } else
      {
        upper <- minInd
        minInd <- floor((lower + minInd) / 2)
      }
    }
  }
  minInd
}

# backend function for weight-based classification of big data sets
.ffWeightClassify <- function(rpairs, threshold.upper, threshold.lower)
{
    # get breakpoints between links, possible links and non-links
    # in sorted vector of weights
    minMatchInd <- .searchThreshold(rpairs@Wdata, threshold.upper, o=rpairs@WdataInd)
    minPossibleInd <- .searchThreshold(rpairs@Wdata, threshold.lower, o=rpairs@WdataInd)
    # if NA is returned (threshold too large), assign length + 1
    if (is.na(minMatchInd)) minMatchInd <- length(rpairs@Wdata) + 1
    if (is.na(minPossibleInd)) minPossibleInd <- length(rpairs@Wdata) + 1
    # determine number of links etc. by positions of breakpoints
    nMatches <- length(rpairs@Wdata) - minMatchInd + 1
    nPossibles <- length(rpairs@Wdata) - nMatches - minPossibleInd + 1
    nNonMatches <- length(rpairs@Wdata) - nMatches - nPossibles
    # index ranges for non-matches, possible matches, matches
    threshInd <- c(0, minPossibleInd-1, minMatchInd-1, length(rpairs@Wdata))
    # fill prediciton vector with the value that occurs most, fill in the others
    lev <- c("N", "P", "L")
    predOrder <- order(c(nNonMatches, nPossibles, nMatches))
    prediction <- ff(lev[predOrder[3]], levels=lev, length=length(rpairs@Wdata))
    for (i in 1:2)
    {
      # calculate Indices of lowest and highest weight of the current classification
      ind1 <- threshInd[predOrder[i]] + 1
      ind2 <- threshInd[predOrder[i]+1]
      if(ind2 >= ind1)
        prediction[rpairs@WdataInd[ind1:ind2]] <- lev[predOrder[i]]
    }
    new("RLResult", data = rpairs, prediction = prediction)
}
