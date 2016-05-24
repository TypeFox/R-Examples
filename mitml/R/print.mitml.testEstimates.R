print.mitml.testEstimates <- function(x,...){
# print method for MI estimates

  cl <- x$call
  est <- x$estimates
  vc <- x$var.comp
  m <- x$m
  adj <- x$adj.df
  dfc <- x$df.com

  # header
  cat("\nCall:\n", paste(deparse(cl)), sep="\n")

  cat("\nFinal parameter estimates and inferences obtained from",m,"imputed data sets.\n")

  cat("\n")

  # print results
  if(!is.null(est)){
    val <- sprintf("%.3f",est)
    w <- max(sapply(c(colnames(est),val),nchar))
    out <- matrix("",nrow(est)+1,ncol(est)+1)
    out[,1] <- format(c("",rownames(est)))
    out[1,-1] <- format(colnames(est),justify="right",width=w)
    out[-1,-1] <- format(val,justify="right",width=w)
  
    for(i in 1:nrow(out)) cat(out[i,],"\n")
  }

  if(!is.null(vc)){
    if(!is.null(est)) cat("\n")

    val <- sprintf("%.3f",vc)
    w <- max(sapply(c("Estimate",val),nchar))
    out <- matrix("",nrow(vc)+1,2)
    out[,1] <- format(c("",rownames(vc)))
    out[1,-1] <- format("Estimate",justify="right",width=w)
    out[-1,-1] <- format(val,justify="right",width=w)
    
    for(i in 1:nrow(out)) cat(out[i,],"\n")
  }

  cat(if(adj){c("\nHypothesis test adjusted for small samples with",
              paste("df=[",paste(dfc,collapse=","),"]\ncomplete-data degrees of freedom.",sep=""))
      }else{"\nUnadjusted hypothesis test as appropriate in larger samples."},"\n")

  cat("\n")

  invisible()
}

summary.mitml.testEstimates <- function(object,...){
# summary method for objects of class mitml.testEstimates

  print.mitml.testEstimates(object,...)

}
