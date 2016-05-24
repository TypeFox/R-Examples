##' Sample blockwise from clustered data
##' 
##' @title Block sampling
##' @param data Data frame
##' @param idvar Column defining the clusters
##' @param size Size of samples
##' @param replace Logical indicating wether to sample with replacement
##' @param \dots additional arguments to lower level functions
##' @return \code{data.frame}
##' @author Klaus K. Holst
##' @keywords models utilities
##' @export
##' @examples
##' 
##' d <- data.frame(x=rnorm(5), z=rnorm(5), id=c(4,10,10,5,5), v=rnorm(5))
##' (dd <- blocksample(d,size=20)) 
##' attributes(dd)$id
##' 
##' \dontrun{
##' blocksample(data.table::data.table(d),1e6)
##' }
blocksample <- function(data, size, idvar="id", replace=TRUE, ...) {
  if (length(idvar)==nrow(data)) {
      id0 <- idvar
  } else {
      if (inherits(data,"data.table")) {          
          id0 <- as.data.frame(data[,idvar,with=FALSE])[,1]
      } else id0 <- data[,idvar]
  }
  ii <-  cluster.index(id0)
  size <- ifelse(missing(size),ii$uniqueclust,size)
  ids <- sample(seq(ii$uniqueclust), size=size,replace=replace)  
  idx <- na.omit(as.vector(t(ii$idclustmat[ids,])))+1
  newid <- rep(seq(size), ii$cluster.size[ids])
  oldid <- id0[idx]
  res <- data[idx,]; res[,idvar] <- newid
  attributes(res)$id <- oldid
  return(res)
}
