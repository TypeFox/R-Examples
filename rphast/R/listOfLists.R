
fill.in.array.lol <- function(lol, arr) {
  if (is.list(lol)) {
    for (i in 1:length(lol)) {
      i1 <- rep(FALSE, length(lol))
      i1[i] <- TRUE
      idx <- array(rep(i1, prod(dim(arr))/length(lol)), dim=dim(arr))
      arr[idx] <- fill.in.array.lol(lol[[i]], array(arr[idx], dim=dim(arr)[2:length(dim(arr))]))
    }
    return(arr)
  }
  return(lol)
}
  


rphast.simplify.list <- function(lol, pointer.only=FALSE) {
  if (!is.list(lol)) return(lol)
  if (!is.null(lol$externalPtr)) return(lol)
  if (length(lol) == 1) 
    return(rphast.simplify.list(lol[[1]]))
  currClass <- attr(lol, "class")
  attr(lol, "class") <- NULL
  # this is a little ugly but I can't think of a better way to deal with special
  # conversion issues.  No way to assign "NA" in C so if frames are undefined
  # they are -1, set to NA here.
  isFeat <- (!is.null(currClass) && currClass=="feat")
  isArray <- (!is.null(currClass) && currClass=="array")
  if (isFeat) {
    if (!is.null(lol$frame)) {
      lol$frame[lol$frame < 0] <- NA
    }
    isDataFrame <- TRUE
    isMatrix <- FALSE
  } else {
    isMatrix <- (!is.null(currClass) && currClass=="matrix")
    isDataFrame <- (!is.null(currClass) && currClass=="data.frame")
  }
  if (isArray) {
    arr <- array(dim=lol$dim, dimnames=lol$dimnames)
    lol$dim <- NULL
    lol$dimnames <- NULL
    lol <- drop(fill.in.array.lol(lol, arr))
    return(lol)
  }
  if (isMatrix || isDataFrame) {
    if (!is.null(lol$row.names)) {
      rowNames <- lol$row.names
      lol$row.names <- NULL
    } else rowNames <- NULL
    # NOTE this is where we take the transpose of the C matrix
    if (isMatrix) lol <- t(as.matrix(as.data.frame(lol), check.names=FALSE))
    if (isDataFrame) lol <- as.data.frame(lol, check.names=FALSE)
    if (!is.null(rowNames)) row.names(lol) <- rowNames
  } else {
    for (i in 1:length(lol))
      lol[[i]] <- rphast.simplify.list(lol[[i]])
    if (!is.null(currClass)) attr(lol, "class") <- currClass
  }
  lol
}
