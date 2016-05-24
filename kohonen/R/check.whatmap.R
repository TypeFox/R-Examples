check.whatmap <- function(x, whatmap)
{
  checkpar <- NULL
  whatmap <- unique(whatmap)
  
  if (!is.null(x$codes)) {
    checkpar <- x$codes
  } else {
    if (!is.null(x$data)) {
      checkpar <- x$data
    } else {
      if (is.list(x)) # not foolproof!
        checkpar <- x
    }
  }
  if (is.null(checkpar))
    stop("no possibility to check argument 'whatmap'!")

  if (!is.list(checkpar)) return(0) #OK, no selection
  if (is.null(whatmap)) {
    if (is.null(x$whatmap)) {
      return(1:length(checkpar))  #OK, no selection
    } else {
      return(x$whatmap)
    }
  }

  if ((is.numeric(whatmap) && all(whatmap %in% 1:length(checkpar))) |
      all(whatmap %in% names(checkpar))) {
    ## if necessary, convert the whatmap argument to indices
    if (!is.numeric(whatmap))
      whatmap <- which(names(checkpar) %in% whatmap)

    return(whatmap) # valid selection
  }

  stop("incorrect whatmap argument") # invalid selection
}
