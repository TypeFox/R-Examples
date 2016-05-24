print.RV.multcomp <- function (x,digits=4,...) {
  cat(paste("\n        Pairwise comparisons using ",x$method,"\n\n",sep=""))
  if (exists("model",where=x)) {
    cat(paste("model:  ",paste(x$model,collapse="\n\t"),"\n\n",sep=""))
  }
  if (exists("data.name",where=x)) {
    cat(paste("data:  ",x$data.name,"\n\n",sep=""))
  }
  rown <- if (is.null(rownames(x))) {FALSE} else {TRUE}
  print(x$p.value,digits=digits,na.print="-",row.names=rown)
  cat(paste("\nP value adjustment method: ",x$p.adjust.method,"\n",sep=""))
  invisible(x)
}
