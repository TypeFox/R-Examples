`print.pfit` <-
function(x, visible=TRUE, ...)
# print method for personfit
# x...object of class "pfit" from (personfit)
{
  pvalues <- 1-pchisq(x$p.fit,x$p.df-1)  # df correction rh 10-01-20
  coef.table <- cbind(round(x$p.fit,3),x$p.df-1,round(pvalues,3),round(x$p.outfitMSQ,3),round(x$p.infitMSQ,3),round(x$p.outfitZ,2),round(x$p.infitZ,2))
  colnames(coef.table) <- c("Chisq","df","p-value","Outfit MSQ", "Infit MSQ", "Outfit t", "Infit t" )
  rownames(coef.table) <- names(x$p.fit)
  if (visible){       # added rh 10-01-20
     cat("\nPersonfit Statistics: \n")
     print(coef.table)
     cat("\n")
  }
  invisible(coef.table)
}

