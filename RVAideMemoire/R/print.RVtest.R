print.RVtest <- function (x,digits=4,quote=TRUE,prefix="",...) {
  cat("\n")
  cat(strwrap(x$method.test,prefix="\t"),sep="\n")
  cat("\n")
  cat("data: ",x$data.name,"\n")
  out <- character()
  if (!is.null(x$statistic)) {
    out <- c(out,paste(names(x$statistic),"=",format(round(x$statistic,4))))
  }
  if (!is.null(x$parameter)) {
    out <- c(out,paste(names(x$parameter),"=",format(round(x$parameter,3))))
  }
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value,digits=digits)
    out <- c(out,paste("p-value",if (substr(fp,1L,1L)=="<") {fp} else {paste("=",fp)}))
  }
  cat(strwrap(paste(out,collapse=", ")),sep="\n")
  if (!is.null(x$alternative)) {
    cat("alternative hypothesis: ")
    if (!is.null(x$null.value)) {
	if (length(x$null.value)==1L) {
	  alt.char <- switch(x$alternative,two.sided="not equal to", 
	    less="less than",greater="greater than")
	  cat("true",names(x$null.value),"is",alt.char,x$null.value,"\n")
	} else {
	  cat(x$alternative,"\nnull values:\n")
	  print(x$null.value,...)
      }
    }
  } else {
    cat(x$alternative,"\n")
  }
  if (!is.null(x$conf.int)) {
    cat(format(100*attr(x$conf.int,"conf.level")),"percent confidence interval:\n", 
	format(c(x$conf.int[1L],x$conf.int[2L])),"\n")
    }
  if (!is.null(x$estimate)) {
    cat("sample estimates:\n")
    print(x$estimate,...)
  }
  cat("\n")
  if (x$p.value<x$alpha) {
    cat(paste("        Pairwise comparisons ",ifelse(exists("method.multcomp",where=x),
	paste("using ",x$method.multcomp,sep=""),""),"\n\n",sep=""))
    print(x$p.value.multcomp,digits=digits,na.print="-",row.names=FALSE)
    cat(paste("\nP value adjustment method: ",x$p.adjust.method,"\n",sep=""))
  } else {
    if (exists("common",where=x)) {
	cat(x$common.name,"\n\n")
	print(x$common,digits=digits,row.names=FALSE)
    }
  }
  invisible(x)
}
