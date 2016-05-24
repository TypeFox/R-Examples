################################################################################
# print method for objects of class "summary.din"                              #
################################################################################
print.summary.din <-
function(x, ...){

# Call: generic
# Input: object of class summary.din
# Print: prints the named list, of an object of class summary.din

################################################################################
# console output                                                               #
################################################################################
if(is.null(x$log.file)){
  d <- utils::packageDescription("CDM")
  base::packageStartupMessage(paste(d$Package," " , d$Version," (Built ",d$Date,")",sep=""))
#   cat("Call:\n",  x$CALL, "\n")
  cat("Call:\n",  x$call, "\n\n")
	  
	cat( "Date of Analysis:" , paste( x$end.analysis ) , "\n" )
	cat("Computation Time:" , print( 
			x$end.analysis - x$start.analysis),"\n\n")
	  
	  
#      "\nItem discrimination index:\n")
#      print(as.table(x$IDI))    
#      cat("\nSummary of skill pattern distribution:\n")
 #     print(x$SKILL.CLASS.PROB)
     cat("\nDeviance = ", x$deviance, " |   Log-Likelihood =" , 
			round( x$din.object$loglike,3)  , "\n") 	 
	
	#*** iterations
	cat( "\nNumber of iterations:" , x$din.object$iter , "\n")			
	if ( ! x$din.object$converged ){ cat("Maximum number of iterations was reached.\n") }
	#***		
			
	 cat( "\nNumber of item parameters:" , x$Npars[,1] , "\n")
	 cat( "Number of skill class parameters:" , x$Npars[,2] , "\n")	 
  cat("\nInformation criteria:",
      "\n  AIC = " , x$AIC,
      "\n  BIC = " , x$BIC, "\n")
  cat("\nMean of RMSEA item fit:" , 
	round( x$din.object$mean.rmsea ,3 ) , "\n")
  cat("\nItem parameters\n")
  obji <- x$item
  rownames(obji) <- NULL
	print( obji , digits=3 )
    cat("\nMarginal skill probabilities:\n")
    print(x$din.object$skill.patt , digits= 4)	  
	# tetrachoric skill correlations
	if( ncol(x$din.object$q.matrix ) > 1 ){ 
		obji <- skill.cor(x$din.object)$cor.skills	
		cat("\nTetrachoric correlations among skill dimensions\n")
		print( obji , digits=4 )
			}
	cat("\nSkill Pattern Probabilities \n\n")
	xt <- round( x$din.object$attribute.patt[,1] , digits=5 )
	names(xt) <- rownames( x$din.object$attribute.patt )
	print(xt)			
}else{

################################################################################
# logfile output                                                               #
################################################################################
    # tr <- try({sink(file = x$log.file)})
	tr <- try({sink(file = paste0( x$log.file , "__SUMMARY.Rout") )})

    if(is.null(tr)){
      gowidth <- getOption("width")
      options(width = 10000)
	  cat("#.......................................................\n")
	  d <-  utils::packageDescription("CDM")
#	  c1 <- citation("CDM")
	  cat(paste( "This is CDM package version " , d$Version , " (" , d$Date , ")\n",sep=""))
#	  print(c1)
	  cat("\n")	
      cat("#-------------------------\n")
      cat("# SUMMARY OF ANALYSIS\n")
	  cat( "Start:" , paste(x$start.analysis ))
	  cat( "\nEnd  :" , paste(x$end.analysis ))
      cat("\n#-------------------------\n\n")
      cat("Model rule:",x$display, "\n")
      cat("Number of observations:",nrow(x$data), "\n")
      cat("Number of items:",nrow(x$q.matrix), "\n")
      cat("Labels of items:", paste(rownames(x$q.matrix), collapse=", "), "\n")
      cat("Number of skills:",ncol(x$q.matrix), "\n")
      cat("Labels of skills:", paste(attributes(x$q.matrix)$skill.labels, collapse=", "), "\n")
      cat("Q-Matrix:\n\n")
      print(data.frame(x$q.matrix))
      
      cat("\n#-------------------------\n")
      cat("# SUMMARY OF MODEL FIT\n")
      cat("#-------------------------\n\n")
      cat("Loglikelihood:", x$loglike, "\n")
      cat("AIC:", x$AIC, "\n")
      cat("BIC:", x$BIC, "\n\n")
      cat("Item discrimination index:\n\n")
      print(as.table(x$IDI))
      cat("\nSummary of skill pattern distribution:\n\n")
      print(x$SKILL.CLASS.PROB)
      
      cat("\n#-------------------------\n")
      cat("SUMMARY OF MODEL RESULTS\n")
      cat("#-------------------------\n\n")
      cat("Item parameter estimates:\n\n")
      print(x$coef)
      cat("\nSkill probability:\n\n")
      print(x$skill.patt)
      cat("\nSkill pattern occurrence probability:\n\n")
      print(x$attribute.patt)
      cat("\nSkill class assignment and skill assignment probabilities for respective response pattern:\n\n")
      print(cbind("freq" = as.vector(table(x$pattern[,1])), unique(x$pattern)))
    
      sink()
      options(width = gowidth)
      cat("\nExtensive summary written to log file:\n", x$log.file,"\n")
    }else cat("\nError while trying to write summary to log file:\n", tr[1])          
  }
}
