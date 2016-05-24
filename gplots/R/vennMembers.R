# Extract intersections
vennMembers <- function(l, universe=NA, names, ...)
{
  venn_object <- getVennCounts(l, universe, intersections=TRUE, ...)
  map <- attr(venn_object, "intersections")
  if(missing(names))
    names <- colnames(venn_object)[-1]

  if(is.matrix(l) || is.data.frame(l))
  {
    ids <- rownames(l)
    retval <- list()
    for(i in names(map))
      retval[[i]] <- ids[map[[i]]]
  }
  else if(is.list(l))
    retval <- map


  flags <- do.call(rbind, strsplit(names(map), character(0), fixed=TRUE))
  rownames(flags) <- names(map)
  colnames(flags) <- names
  nameList <- list()
  for(i in 1:nrow(flags))     nameList[[i]] <- ifelse(flags[i,]=="1", colnames(flags), "")
  nameList <- do.call(data.frame,nameList)
  nameList <- apply(nameList, 2, paste, collapse=":")
  nameList <- gsub('::+', ':', nameList)
  nameList <- gsub('^:+', '',  nameList)
  nameList <- gsub(':+$', '',  nameList)

  names(retval) <- nameList

  sortTab <- cbind(sapply(nameList, nchar), nameList)
  ord     <- order(sortTab[,1], sortTab[,2])

  retval <- retval[ord]

  retval <- lapply(retval, as.character)

  retval
}
