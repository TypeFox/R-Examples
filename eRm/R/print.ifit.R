`print.ifit` <-
function(x, visible=TRUE, ...)
# print method for itemfit
# x...object of class "ifit" from (itemfit)
{
  pvalues <- 1-pchisq(x$i.fit,x$i.df-1)  # df correction rh 10-01-20
  coef.table <- cbind(round(x$i.fit,3),x$i.df-1,round(pvalues,3),round(x$i.outfitMSQ,3),round(x$i.infitMSQ,3),round(x$i.outfitZ,2),round(x$i.infitZ,2))
  colnames(coef.table) <- c("Chisq","df","p-value","Outfit MSQ", "Infit MSQ", "Outfit t", "Infit t" )
  rownames(coef.table) <- names(x$i.fit)
  if (visible){       # added rh 10-01-20
    cat("\nItemfit Statistics: \n")
    print(coef.table)
    cat("\n")
  }
  invisible(coef.table)
}

