##' Generic method for labeling elements of an object
##'
##' @title Label elements of object
##' @param x Object
##' @param \dots Additional arguments
##' @author Klaus K. Holst
##' @export
`baptize` <- function(x,...) UseMethod("baptize")

###{{{ baptize.lvm

##' @export
baptize.lvm <- function(x,labels,overwrite=FALSE,unique=FALSE,...) {
  p <- describecoef(x, mean=TRUE)
  sym <- lava.options()$sym
  MeanFix <- intfix(x)
  RegFix <- regfix(x)
  CovFix <- covfix(x)
  count <- 0
  curlab <- parlabels(x)
  coef(x)
  for (i in seq_along(p)) {
    p0 <- p[[i]]
    if (attributes(p0)$type=="reg") {
      curfix <- RegFix$values[p0[2],p0[1]]
      curlab <- RegFix$labels[p0[2],p0[1]]
      if (all(is.na(c(curfix,curlab))) | overwrite) {
        count <- count+1
##        st <- ifelse(missing(labels),paste0("p",count),labels[count])
        st <- ifelse(missing(labels),paste(p0[1],p0[2],sep=sym[1]),labels[count])
        regfix(x,from=p0[2],to=p0[1]) <- st
      }
    } else if (attributes(p0)$type=="cov") {
      curfix <- CovFix$values[p0[2],p0[1]]
      curlab <- CovFix$labels[p0[2],p0[1]]
      if (all(is.na(c(curfix,curlab))) | overwrite) {
        count <- count+1
##        st <- ifelse(missing(labels),paste0("p",count),labels[count])
##        st <- paste0("p",count)
        st <- ifelse(missing(labels),paste(p0[1],p0[2],sep=sym[2]),labels[count])
        covfix(x,p0[2],p0[1],exo=FALSE) <- st
      }
    } else { ## Mean parameter
      curfix <- MeanFix[[p0]]
      if (length(curfix)>0)
      if (is.na(curfix) | overwrite) {
        count <- count+1
        st <- ifelse(missing(labels),p0,labels[count])
##        st <- ifelse(missing(labels),paste0("m",count),labels[count])
        intfix(x,p0) <- st
      }
    }
  }
  if (index(x)$npar.ex>0) {
    x$exfix[is.na(x$exfix)] <- names(x$exfix)[is.na(x$exfix)]
    index(x) <- reindex(x)
  }
  return(x)
}

###}}}
