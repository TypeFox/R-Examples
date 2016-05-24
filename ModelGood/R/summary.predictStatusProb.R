

#' @S3method summary predictStatusProb
summary.predictStatusProb <- function(object,cuts=seq(0,100,25),...){
    table("Predicted risk"=cut(round(100*object,1),
              cuts,
              include.lowest=TRUE,
              labels=paste(paste(cuts[-length(cuts)],cuts[-1],sep="-"),"%",sep="")))
}
