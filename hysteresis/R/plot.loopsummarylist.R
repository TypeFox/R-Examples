plot.loopsummarylist <- function(x,main=NULL,values=NULL,...){
  a <- x
  if (length(dim(a$models))<=1) mapply(plot.loopsummary,a$models,main=names(a$models),MoreArgs=list(values=values,...))
  else {
    length.values <- length(a$models[[1]]$values[,1])
    thenames <- do.call(paste,lapply(seq_len(ncol(a$Boot.Estimates)-length.values), function(i) a$Boot.Estimates[,i]))
    mapply(plot.loopsummary,a$models,main=thenames,MoreArgs=list(values=values,...))
  }
  if (!is.null(main)) mtext(paste(main),side=3,outer=T)
}

plot.loopsummarylist2r <- function(x,main=NULL,values=NULL,...){
  a <- x
  if (length(dim(a$models))<=1) mapply(plot.loop2rsummary,a$models,main=names(a$models),MoreArgs=list(values=values,...))
  else {
    length.values <- length(a$models[[1]]$values[,1])
    thenames <- do.call(paste,lapply(seq_len(ncol(a$Boot.Estimates)-length.values), function(i) a$Boot.Estimates[,i]))
    mapply(plot.loop2rsummary,a$models,main=thenames,MoreArgs=list(values=values,...))
  }
  if (!is.null(main)) mtext(paste(main),side=3,outer=T)
}
