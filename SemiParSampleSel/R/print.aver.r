print.aver <- function(x, ...){

es <- format(x$res, digits = 3, trim=TRUE)

cat("\nEstimated average with ",(1-x$sig.lev)*100,"% confidence interval:\n\n",sep="")

cat(es[2]," (",es[1],",",es[3],")\n\n",sep="")

invisible(x)

}



