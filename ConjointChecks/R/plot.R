plot.checks<-function(x, items=NULL, item.labels=TRUE, ...) {
  #graphical args
  modify.args<-function(...,items) {
    list(...)->args
    if (is.null(args$ylim)) args$ylim<-c(0,1)
    if (is.null(args$xlab)) args$xlab<-""
    if (is.null(args$ylab)) args$ylab<-"Proportion Violations"
    if (is.null(args$type)) args$type<-"l"
    if (is.null(args$lty)) args$lty<-rep(1,length(items))
    if (is.null(args$col)) args$col<-"black"
    if (is.null(args$xaxt)) args$xaxt<-"n"
    args
  }
  #
  if (is.null(items)) items <- 1:ncol(x@tab)
  #
  modify.args(...,items=items)->args
  x@tab[,items]->dat
  c(list(y=dat),args)->list.of.args
  do.call("matplot",list.of.args)
  mtext(side=1,line=1,paste("Increasing Ability"))
  if (item.labels) {
    if (length(items)==1) mtext(side=1,line=2,paste("Item",items)) else {
      apply(x@tab,2,which.max)->maxes.x
      apply(x@tab,2,max,na.rm=TRUE)->maxes.y
      for (i in 1:length(items)) text(maxes.x[items],maxes.y[items],items,pos=3)
    }
  }
}
