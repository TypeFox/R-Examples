plot.fittedlooplist<-function(x,main=NULL,values=NULL,...){
  a <- x
   if (length(dim(a$models))<=1) mapply(plot.fittedloop,a$models,main=names(a$models),MoreArgs=list(values=values,...))
else {
  length.values <- length(a$models[[1]]$values)
  thenames <- do.call(paste,lapply(seq_len(ncol(a$Estimates)-length.values), function(i) a$Estimates[,i]))
  mapply(plot.fittedloop,a$models,main=thenames,MoreArgs=list(values=values,...))
}
  if (!is.null(main)) mtext(paste(main),side=3,outer=T)
  }

plot.fittedlooplist2r<-function(x,main=NULL,values=NULL,...){
  a <- x
   if (length(dim(a$models))<=1) mapply(plot.loop2r,a$models,main=names(a$models),MoreArgs=list(values=values,...))
else {
  length.values <- length(a$models[[1]]$values)
  thenames <- do.call(paste,lapply(seq_len(ncol(a$Estimates)-length.values), function(i) a$Estimates[,i]))
  mapply(plot.loop2r,a$models,main=thenames,MoreArgs=list(values=values,...))
}
  if (!is.null(main)) mtext(paste(main),side=3,outer=T)
  }
