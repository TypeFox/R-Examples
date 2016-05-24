summary.ndlCrossvalidate <- function(object, ...)
{ n.fits = object$k
  fits = object$fits
  if(n.fits>=2)
    { names.fits <- names(fits[[1]])
      names.fits <- names.fits[-which(names.fits %in% c("df.null","df.model","crosstable","recall.predicted","precision.predicted"))]
      statistics.all <- data.frame(sapply(names.fits, function(statistic) sapply(fits, function(x) x[[statistic]])))
      rownames(statistics.all) <- 1:n.fits
      statistics.summary <- data.frame(matrix(c(apply(statistics.all,2,mean),apply(statistics.all,2,min),apply(statistics.all,2,max)),3,byrow=TRUE,dimnames=list(c("Mean","Minimum","Maximum"),names.fits)))
      statistics.summary[c("Minimum","Maximum"),"n.test"] <- NA
      crosstables <- array(,dim=c(dim(fits[[1]]$crosstable),n.fits), dimnames=c(dimnames(fits[[1]]$crosstable),NULL))
      recall.predicted <- NULL
      precision.predicted <- NULL
      for(i in 1:n.fits)
         { crosstables[,,i] <- fits[[i]]$crosstable
           recall.predicted <- rbind(recall.predicted, fits[[i]]$recall.predicted)
           precision.predicted <- rbind(precision.predicted, fits[[i]]$precision.predicted)
         }
      crosstable.summary <- apply(crosstables,c(1,2),mean)
      recall.predicted.summary <- apply(recall.predicted,2,mean)
      precision.predicted.summary <- apply(precision.predicted,2,mean)

      sumry = list(call=object$call, formula=object$formula, statistics.summary=statistics.summary, crosstable.summary=crosstable.summary, recall.predicted.summary = recall.predicted.summary, precision.predicted.summary = precision.predicted.summary, statistics.all=statistics.all, k=n.fits, n.total=object$n.total, n.train=object$n.train, n.test=object$n.test, ...)

    }
  else
    sumry=fits

  class(sumry) <- "summary.ndlCrossvalidate"  
  sumry;

}

print.summary.ndlCrossvalidate <- function(x, digits=max(3,getOption("digits")-3), ...)
{
  if(!is.null(x$digits))
    digits=x$digits

  cat("\nCross-validation summary statistics\n\n")

  cat("Call:\n")
  print(x$call)
  cat("\nFormula:\n")
  print(x$formula)
  cat(c("\nNumber of folds: ",x$k))
  cat(c("\nN(total): ",x$n.total))
  cat(c("\nN(train): ",x$n.train))
  cat(c("\nN(test):  ",x$n.test))
  cat("\n\n")

  sumry <- data.frame(format(apply(signif(x$statistics.summary[-which(names(x$statistics.summary) %in% c("n.test","d.lambda.prediction","d.tau.classification","p.lambda.prediction","p.tau.classification"))],digits),c(1,2),as.character),justify="right"))
  rownames(sumry) <- format(rownames(sumry),justify="left")
  colnames(sumry) <- format(colnames(sumry),justify="left")

  print(t(sumry),quote=FALSE)
  cat("\n")

  invisible(x)
}
