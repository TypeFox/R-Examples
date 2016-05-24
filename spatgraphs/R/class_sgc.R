#' Creator for sgc
#'
#' @param clusters list of clusters as point indices
#' @param type type
#' @param pars parameters
#' @param note notes
#'
#' @export
as.sgc<-function(clusters, type="?",pars=NULL,note=NULL)
{
  e <- as.sg(clusters, type, pars, note)
  e$parameters<-pars
  e$nclusters<-length(clusters)
  e$N<-max(unlist(lapply(clusters, max)))
  names(e)[1] <- "clusters"
  class(e) <- "sgc"
  if(!is.null(note)) e$note <- note
  e
}
###############################################################################
#' sgc print method
#'
#' @param x sgc object
#' @param ... ignored
#'
#' @export
print.sgc<-function(x,...)
{
  nam<-names(x$parameters)
  p<-"?"
  p<-paste(", par=(",paste(x$parameters,collapse=","),")",sep="")

  cat(paste("'Spatgraphs' cluster/component list-of-lists:",
            "\ngraph type '",x$type,"'",p,", ",x$nclusters," component",ifelse(x$N>1,"s","")," (",x$N," points).\n",sep=""))
  if(!is.null(x$note))cat(paste("Note: ", x$note,".\n",sep=""))
}

#####################################################################################
#' sgc summary
#'
#' @param object sgc object
#' @param ... ignored
#'
#' @export
summary.sgc<-function(object, ...)
{
  args<-list(...)
  print(object)
  cls<-sapply(object$clusters, length)
  cat("Isolated points:",sum(cls==1),"\n")
  cat("Cluster size stats:\n")
  print(summary(cls))
}

