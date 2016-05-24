print.kSamples <-
function (x, ...) 
{
######################################################
#
# This is a print function for objects of class kSamples,
# as they are produced by ad.test, kw.test,
# ad.combined.test, kw.combined.test, contingency2xt,
# contingency2xt.comb, Steel.test, SteelConfInt,
# and JT.test.
#
# Fritz Scholz, August 2015
#
#######################################################
	if(names(x)[2]=="k"){# checking whether the object x 
                         #came from ad.test or qn.test of JT.test
		if(x$test.name=="Steel"){
			cat("\nSteel Multiple Comparison Wilcoxon Test:\nk treatments against a common control (1st sample)\n\n")
		}else{
	    	 	cat(paste("\n\n",x$test.name,"k-sample test.\n"))
		}
	    	cat(paste("\nNumber of samples: ", x$k))
	    	cat("\nSample sizes: ",paste(x$ns,collapse=", "))
	    	cat(paste("\nNumber of ties:", x$n.ties))
	    	if(x$test.name == "Anderson-Darling"){
	    		cat(paste("\n\nMean of ",x$test.name," Criterion:", 
	        		x$k-1))
	    		cat(paste("\nStandard deviation of ",x$test.name," Criterion:", 
	        		x$sig))
	    		cat(paste("\n\nT.AD = (",x$test.name," Criterion - mean)/sigma"))
		}
		if(x$test.name != "Jonckheere-Terpstra" ){
			cat("\n\nNull Hypothesis: All samples come from a common population.\n\n")
		}else{
			cat("\n\nNull Hypothesis: All samples come from a common population.\n")
			cat("Alternative: Samples indicate a positive trend.\n\n")
		}
    		
    		if(x$method=="simulated") cat(paste("Based on Nsim =",x$Nsim,"simulations\n\n"))
    		if(x$test.name == "Anderson-Darling" ){print(signif(x$ad,4))}
		if(x$test.name == "van der Waerden scores" ) print(signif(x$qn,4))
		if(x$test.name == "Kruskal-Wallis" ) print(signif(x$qn,4))
		if(x$test.name == "normal scores" ) print(signif(x$qn,4))
 		if(x$test.name == "Jonckheere-Terpstra" ){
			print(signif(x$JT,4))
		}
		if (x$warning) {
        		cat("\n\nWarning: At least one sample size is less than 5,\n")
			cat("  asymptotic p-values may not be very accurate.\n")
    		}
    		invisible(x)
	}
  	if(names(x)[2]=="M"){# checking whether the object x came from ad.combined.test
   		cat(paste("Combination of",x$test.name,"K-Sample Tests.\n"))
    		cat(paste("\nNumber of data sets =", x$M,"\n"))
    		cat("\nSample sizes within each data set:\n")
    		ns <- NULL
    		k <- length(x$n.samples)
    		d.sets <- paste("Data set",1:k)
    		for(i in 1:k){
      			cat(d.sets[i],": ",x$n.samples[[i]])
       			cat("\n")
    		}
    		if(k>3) AD.name=paste("AD.1","...",paste("AD.",k,sep=""),sep="+")
    		if(k==2)AD.name=paste("AD.1+AD.2")
    		if(k==3)AD.name=paste("AD.1+AD.2+AD.3")
    		if(k>3) QN.name=paste("QN.1","...",paste("QN.",k,sep=""),sep="+")
    		if(k==2)QN.name=paste("QN.1+QN.2")
    		if(k==3)QN.name=paste("QN.1+QN.2+QN.3")
    		cat("Total sample size per data set: ")
    		cat(x$nt,"\n")
    		cat("Number of unique values per data set: ")
    		cat(x$nt-x$n.ties,"\n")
		if(x$test.name=="Anderson-Darling"){ 
    			cat(paste("\nAD.i =",x$test.name,"Criterion for i-th data set\n"))
    			cat("Means:",x$mu,"\n")
    			cat("Standard deviations:", x$sig,"\n")
			cat("\nT.i = (AD.i - mean.i)/sigma.i\n")
		}
    		cat("\nNull Hypothesis:\nAll samples within a data set come from a common distribution.\n")
    		cat("The common distribution may change between data sets.\n\n")
    		if(x$test.name=="Anderson-Darling"){ 
     			nx <- length(x$ad.list)
        		if(x$method=="simulated") cat(paste("Based on Nsim =",x$Nsim,"simulations\n\n"))
	    		for(i in 1:nx){
				cat(paste("for data set",i,"we get\n"))
    				print(signif(x$ad.list[[i]],4))
            			cat("\n")
    			}
			cat("Combined Anderson-Darling Criterion: AD.comb =",AD.name,"\n")
	    		cat("Mean =",x$mu.c,"   Standard deviation =",round(x$sig.c,5),"\n")
    			cat("\nT.comb = (AD.comb - mean)/sigma\n")
	    		cat("\n")
        		if(x$method=="simulated") cat(paste("Based on Nsim =",x$Nsim,"simulations\n\n"))
			ad.c <- x$ad.c
	    		print(signif(ad.c,4))
		}

    		if(x$test.name == "van der Waerden scores" | 
		   x$test.name == "normal scores" |
		   x$test.name == "Kruskal-Wallis" ){ 
    			nx <- length(x$qn.list)
        		if(x$method=="simulated") cat(paste("Based on Nsim =",x$Nsim,"simulations\n\n"))
	   	 	for(i in 1:nx){
				cat(paste("for data set",i,"we get\n"))
    				print(signif(x$qn.list[[i]],4))
            			cat("\n")
    			}
			cat("Combined Criterion: QN.combined =",QN.name,"\n")
			cat("\n")
        		if(x$method=="simulated") cat(paste("Based on Nsim =",x$Nsim,"simulations\n\n"))
	    		print(signif(x$qn.c,4))
		}
    		if (x$warning) {
        		cat("\n\nWarning: At least one sample size is less than 5,\n")
			cat("  asymptotic p-values may not be very accurate.\n")
    		}
		cat("\n")
		invisible(x)
  	}
  	if(x$test.name == "2 x t Contingency Table"){
   	 	cat(paste("\n      Kruskal-Wallis Test for 2 x",x$t,"Contingency Table\n\n"))
   		if(x$method=="simulated") cat(paste("      Based on Nsim =", x$Nsim,"simulations\n\n"))
		print(signif(x$KW.cont,4))
		cat("\n")
		invisible(x)
	}
  	if(x$test.name == "Combined 2 x t Contingency Tables"){
   	 	cat("\n   Combined Kruskal-Wallis Tests for 2 x t Contingency Tables\n\n")
   		if(x$method=="simulated") cat(paste("      Based on Nsim =", x$Nsim,"simulations\n\n"))
		nx <- length(x$kw.list)
		for( i in 1:nx){
				cat(paste("for data set",i,"we get\n"))
    				print(signif(x$kw.list[[i]],4))
            			cat("\n")
		}
		if(nx>3) KW.name=paste("KW.1","...",paste("KW.",k,sep=""),sep="+")
    		if(nx==2)KW.name=paste("KW.1+KW.2")
    		if(nx==3)KW.name=paste("KW.1+KW.2+KW.3")
		cat("Combined Criterion: KW.combined =",KW.name,"\n")
			cat("\n")
        		if(x$method=="simulated") cat(paste("Based on Nsim =",x$Nsim,"simulations\n\n"))
	    		print(signif(x$kw.c,4))
			cat("\n")
		invisible(x)
	}
  	if(x$test.name == "Steel"){
                print(signif(x$st,4))
		invisible(x)
	}
	if(x$test.name == "Steel.bounds"){
	cat("\nSteel Multiple Comparison Confidence Bounds for Shift Parameters\nBased on Wilcoxon Tests, k Treatments against a Common Control\n\n")
		cat("size of control sample: ",x$n0,"\n")
		if(length(x$ns) > 1){
			cat("sizes of treatment samples: ",paste(x$ns,collapse=", "),"\n")
		}else{
			cat("size of treatment sample: ",x$ns,"\n")
		}
		if(x$n.ties > 0){
			cat("number of ties: ",x$n.ties,"\n")
			cat("intervals should be widened on each end by the rounding epsilon\n")
			cat("to conservatively maintain the stated joint confidence level\n")
			    
		}
		

		cat("\nconservative bounds based on asymptotics\n\n")
		print(x$bounds[[1]])
		
		cat("\nbounds based on asymptotics,\n")
		cat("with level closest to nominal\n\n")
		print(x$bounds[[2]])
		
		if(x$method=="simulated"){
			cat("\nconservative bounds based on simulation\n\n")
			print(x$bounds[[3]])
			
			cat("\nbounds based on simulation,\n")
			cat("with level closest to nominal\n\n")
			print(x$bounds[[4]])
			
		}
		invisible(x)
	}
        
	
}

