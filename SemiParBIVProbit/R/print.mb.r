print.mb <- function(x, ...){

if(x$Model == "B") nf <- nf1 <- "ATE"
if(x$Model == "BSS"){ nf <- "prevalence"; nf1 <- "Prevalence"}

if(is.null(x$IV)) mt <- "\nWorst-Case" else mt <- "\nIV"


cat(mt," Manski's bounds for ",nf," (%):\n","Lower Bound: ",x$LB,"\n","Upper Bound: ",x$UB,"\n",sep="")
cat("\nImbens&Manski's CI for ",nf," (%): [",x$CI[1],",",x$CI[2],"]","\n",sep="")

if( (x$UB - x$LB) < 0 ) cat("\nWARNING: The lower and upper bounds do not make sense!\n")

cat("\n",nf1," assuming random assignment (%): ",x$av.p,"\n\n",sep="")

invisible(x)

}



