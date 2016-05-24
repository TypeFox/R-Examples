`print.runJags` <-
function(x, ... )
{

  if (class(x) != "runJags")
    stop("'x' must be of class 'runJags'")

##  cat("'jagsControl':", deparse(substitute(jags.control)), "\n")
  
  cat("CMD File:", x$cmd.file,"\nJAGS started at",x$start.time,"\n")
  if (x$exit==0) {
    cat("JAGS run completed successfully at", x$end.time,"\n")
  } else {
    cat("Error: JAGS run failed at", x$end.time,"\n")
  }

 if (length(x$elapsed.time)>0){
   cat("Elapsed times:\n")
   print(structure(x$elapsed.time,
         names=c("user","R.system","R.elapsed","Jags.user","Jags.system")))
 }
}

