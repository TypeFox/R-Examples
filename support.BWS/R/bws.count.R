bws.count <- function(data)
{
# Name : bws.count
# Title: Calculating count-based BW scores
# Arguments:
#  data   a data frame created using bws.dataset()


# set variables, vectors, and matrices

  numItems       <- attributes(data)$nitems
  freqencyItem   <- attributes(data)$fitem
  variableNames  <- colnames(data)[c(8:(7 + numItems))]
  id             <- c(subset(data, data$RES == TRUE, select = "ID"))
  uniqueId       <- unique(id[[1]])
  numRespondents <- length(uniqueId)


# create BEST matrix (B) and WORST matrix (W)
  
  B <- data.matrix(subset(data, data$RES == TRUE, select = variableNames))
  rownames(B) <- NULL
  W <- B

  B[which(B == -1, arr.ind = TRUE)] <- 0
  W[which(W ==  1, arr.ind = TRUE)] <- 0
  W[which(W == -1, arr.ind = TRUE)] <- 1

  colnames(B) <- variableNames
  colnames(W) <- variableNames

  B <- cbind(ID = id[[1]], B)
  W <- cbind(ID = id[[1]], W)

  
# calculate various BW scores

  # disaggregated scores
  disaggreB <- do.call(rbind, by(B[, 2:(1 + numItems)], B[, "ID"], colSums))
  disaggreW <- do.call(rbind, by(W[, 2:(1 + numItems)], W[, "ID"], colSums))
  diffBW     <- disaggreB - disaggreW
  std.diffBW <- diffBW / freqencyItem
  
  # aggregated scores
  aggreB  <- colSums(disaggreB)
  aggreW  <- colSums(disaggreW)
  aggreBW     <- aggreB - aggreW
  std.aggreBW <- aggreBW / (numRespondents * freqencyItem)
  sqrt.aggreBW     <- sqrt(aggreB / aggreW)
  std.sqrt.aggreBW <- sqrt.aggreBW / max(sqrt.aggreBW)


# format and return output

  rtn <- list(disaggregate = list(ID           = uniqueId,
                                  B            = disaggreB,
                                  W            = disaggreW,
                                  BW           = diffBW,
                                  stdBW        = std.diffBW),

              aggregate    = data.frame(
                                  B            = aggreB,
                                  W            = aggreW,
                                  BW           = aggreBW, 
                                  stdBW        = std.aggreBW,
                                  sqrtBW       = sqrt.aggreBW,
                                  std.sqrtBW   = std.sqrt.aggreBW),

              information  = list(nrespondents = numRespondents,
                                  nitems       = numItems,
                                  fitem        = freqencyItem,
                                  vnames       = variableNames))

  class(rtn) <- "bws.count"

  return(rtn)
}

