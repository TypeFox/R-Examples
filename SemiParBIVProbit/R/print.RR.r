print.RR <- function(x, ...){


if(x$mar2 %in% x$bl){ #  && x$type != "naive"){

es <- format(x$res, digits = 3, trim=TRUE)

cat("\nRisk ratio with ",(1-x$prob.lev)*100,"% interval:\n\n",sep="")

cat(es[2]," (",es[1],",",es[3],")\n\n",sep="")

}


if( !(x$mar2 %in% x$bl) && x$type != "naive"){

es <- format(x$Ratios, digits = 3, trim=TRUE)

#cat("\nRisk ratios with ",(1-x$prob.lev)*100,"% intervals:\n\n",sep="")

cat("\n")
print(es)
cat("\n")

}

invisible(x)

}
