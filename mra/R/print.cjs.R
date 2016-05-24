print.cjs <- function( x, alpha=c(0.05, 0.01), ... ){

nx <-	x$aux$nx 
ny <-	x$aux$ny

cap.coef <- round(x$capcoef, 5)
sur.coef <- round(x$surcoef, 5)
se.cap <- round( x$se.capcoef, 5 )
se.sur <- round( x$se.surcoef, 5 )

if(nx > ny){
	sur.coef<- c(sur.coef, rep("", nx-ny))
	se.sur  <- c(se.sur, rep("", nx-ny))
} else {
	cap.coef<- c(cap.coef, rep("", ny-nx))
	se.cap  <- c(se.cap, rep("", ny-nx))	
}

cat("Call:\n")
print(x$aux$call)
cat("\n")

cat(paste( format( c(" Capture var", names(cap.coef))), 
	format( c(" Est", cap.coef) ),  
	format( c(" SE", se.cap) ),
    "  ",
	format( c(" Survival var", names(sur.coef))), 
	format( c(" Est", sur.coef) ), 
	format( c(" SE", se.sur) ),
	"\n", sep="  "))

cat(paste("\nMessage =", x$message[2] ))
cat(paste("\nLink =", x$aux$link ))
cat(paste("\nModel df = ", x$df))
cat(paste("\nStd Errors and QAIC adjusted for C_hat = ", round(x$vif,6), "on", x$vif.df, "df"))
cat(paste("\nLog likelihood = ", x$loglike))
cat(paste("\nDeviance = ", x$dev))
cat(paste("\nAIC = ", x$aic))
cat(paste("\nAICc = ", x$aicc))
cat(paste("\nQAIC = ", x$qaic))
cat(paste("\nQAICc = ", x$qaicc))
#cat(paste("\nEBC = ", x$ebc, "\n"))
cat("\n")


if( inherits(x,"cjsgof") ){
	# This object has goodness of fit results. Print them.
	or.stars <- ifelse( x$or.pvalue <= alpha[2], "**", ifelse( x$or.pvalue <= alpha[1], "*", " "))
	HL.stars <- ifelse( x$HL.pvalue <= alpha[2], "**", ifelse( x$HL.pvalue <= alpha[1], "*", " "))
	mr.stars <- ifelse( x$gof.pvalue <= alpha[2], "**", ifelse( x$gof.pvalue <= alpha[1], "*", " "))

	cat(paste("\nOverall goodness of fit results:\n"))
	cat(paste( format( c("            ", "    Osius-Rojek:", "Hosmer-Lemeshow:", "  M&R ChiSquare:") ), 
		   format( c("Statistic", round(x$or.chi,  4), round(x$HL.chi,  4), round(x$gof.chi,  4) )),
		   format( c("df", x$or.df, x$HL.df, x$gof.df )),
		   format( c("p", round(x$or.pvalue,  4), round(x$HL.pvalue,  4), round(x$gof.pvalue,  4) )),
		   format( c("  ", or.stars, HL.stars, mr.stars )),
		   "\n", sep="\t")) 


	roc.desc <- ifelse( x$roc < 0.7, "marginal discrimination", 
			ifelse( (0.7 <= x$roc) & (x$roc < 0.8), "acceptable discrimination", 
			ifelse( (0.8 <= x$roc) & (x$roc < 0.9), "excellent discrimination", 
			"outstanding discrimination")))
	cat(paste("    ROC = ", round(x$roc, 4), ":", sep=""))
	cat(paste("\t", roc.desc, "\n", sep=""))

	cat(paste("\nTargeted lack of fit tests:\n"))
	cat(paste( format( c("         ", "  Occasion (Test 4):", "Individual (Test 5):") ), 
		   format( c("ChiSquare", round(c(x$t4.chi, x$t5.chi), 3) )),
		   format( c("Df", x$t4.df, x$t5.df) ),
		   format( c("p", round(c(x$t4.pvalue, x$t5.pvalue), 4) )),
		   "\n", sep="\t")) 

}

if(!is.null(x$n.hat)){
	cat("\nPopulation Size Estimates (se):\n")
	cat(paste( "N", 2:x$aux$ns , "=", round(x$n.hat[-1]), " (", round(x$se.n.hat[-1],2), "), ", sep=""))
	cat("\n")
}

invisible()

}


