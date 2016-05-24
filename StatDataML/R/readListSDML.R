# used for list tags as well as for listdata tags
readListSDML <- function(x)
{
  if (is.null(x)) return(NULL)

  if (xmlName(x) == "list") {
    ## dimension
    dimension <- readDimensionSDML(x[["dimension"]])

    ## properties (NULL if none)
    attrib <- readProperties(x[["properties"]])

    ## fetch sublists
    thislist <- xmlChildren(x[["listdata"]])
    thislist <- thislist[sapply(thislist[seq_along(thislist)],
                                xmlName) != "text"]
    thislist <- lapply(thislist, readListSDML)
    if (!length(thislist)) return(thislist)

    ## set names:
    ### no names for arrays
    if(any(names(thislist) == "array")) names(thislist) <- NULL

    ### list with dim attribute?
    if(length(dimension$dim) > 1) {
      dim(thislist) <- dimension$dim
      dimnames(thislist) <- dimension$names
    } else
    ### only one dimension?
      if(length(dimension$names[[1]]) > 0)
        names(thislist) <- dimension$names[[1]]

    ## set/append properties
    atL <- attributes(thislist)
    if (!is.null(atL)) attrib <- c(attrib, atL)
    if (!is.null(attrib)) attributes(thislist) <- attrib

    return(thislist)
  }

  if (xmlName(x) == "array")
    return(readArraySDML(x))

  if (xmlName(x) == "empty")
    return(NULL)
}
