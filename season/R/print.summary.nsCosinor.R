## print.summary.nsCosinor.R
print.summary.nsCosinor <- function(x, ...){

  if (class(x)!="summary.nsCosinor"){
    stop("Object must be of class 'summary.nsCosinor'")
  } 

  k<-length(x$cycles);
  cat("Statistics for non-stationary cosinor based on MCMC chains\n")
  cat("Number of MCMC samples = ",x$niters-x$burnin+1,"\n",sep="", ...)
  cat("\nStandard deviations\n")
  cat("Residual, mean=",x$stats$errorstats[1],", 95% CI [",x$stats$errorstats[2],", ",x$stats$errorstats[3],"]\n",sep="", ...)
  if(k==1){
    cat("Cycle=",x$cycles,"\n",sep="", ...)
    cat("Season, mean=",x$stats$wstats[1],", 95% CI [",x$stats$wstats[2],", ",x$stats$wstats[3],"]\n",sep="", ...)
    cat("\nPhase and amplitude\n")
    cat("Cycle=",x$cycles,"\n",sep="", ...)
    cat("Amplitude, mean=",x$stats$ampstats[1],", 95% CI [",x$stats$ampstats[2],", ",x$stats$ampstats[3],"]\n",sep="", ...)
    cat("Phase (radians), mean=",x$stats$phasestats[1],", 95% CI [",x$stats$phasestats[2],", ",x$stats$phasestats[3],"]\n",sep="", ...)
  } # end of k==1
  if(k>=2){
    for (index in 1:k){ 
      cat("Cycle=",x$cycles[index],"\n",sep="", ...)
      cat("Season, mean=",x$stats$wstats[index,1],", 95% CI [",x$stats$wstats[index,2],", ",x$stats$wstats[index,3],"]\n",sep="", ...)
    }
    cat("\nPhase and amplitude\n")
    for (index in 1:k){ 
      cat("Cycle=",x$cycles[index],"\n",sep="", ...)
      cat("Amplitude, mean=",x$stats$ampstats[index,1],", 95% CI [",x$stats$ampstats[index,2],", ",x$stats$ampstats[index,3],"]\n",sep="", ...)
      cat("Phase (radians), mean=",x$stats$phasestats[index,1],", 95% CI [",x$stats$phasestats[index,2],", ",x$stats$phasestats[index,3],"]\n",sep="", ...)
    }
  } # end of k>=2
}
