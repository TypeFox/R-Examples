plot.ellipsefitlist<-function(x,main=NULL,values=NULL,...){
  a <- x
   if (length(dim(a$models))<=1) mapply(plot.ellipsefit,a$models,main=names(a$models),MoreArgs=list(values=values,...))
else {
  length.values <- length(a$models[[1]]$values)-any(names(a$models[[1]]$values)=="n")
  thenames <- do.call(paste,lapply(seq_len(ncol(a$Estimates)-length.values), function(i) a$Estimates[,i]))
  mapply(plot.ellipsefit,a$models,main=thenames,MoreArgs=list(values=values,...))
}
  if (!is.null(main)) mtext(paste(main),side=3,outer=T)
  }
