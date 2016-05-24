print.Sprop <-
function(x,...) {
 cat('\nSprop object: Sample proportion estimate\n')
 if(x$call$N == Inf){ cat("Without ")}
 else{ cat("With ")}
 cat("finite population correction: N=",x$call$N,"\n\n")
 cat("Proportion estimate: ",round(x$p,4),"\n" )
 cat("Standard error: ",round(x$se,4),"\n\n")
 if(x$call$N < Inf){
  cat(100 * x$call$level,"% approximate hypergeometric confidence interval: \n",sep="")
  cat(" proportion: [",round(x$ci$approx[1],4),",",round(x$ci$approx[2],4),"]\n",sep="")
  cat(" number in population: [",x$nr$approx[1],",",x$nr$approx[2],"]\n",sep="")
  cat(100 * x$call$level,"% exact hypergeometric confidence interval: \n",sep="")
  cat(" proportion: [",round(x$ci$exact[1],4),",",round(x$ci$exact[2],4),"]\n",sep="")
  cat(" number in population: [",x$nr$exact[1],",",x$nr$exact[2],"]\n",sep="")
 }
 else{
 	cat(100 * x$call$level,"% asymptotic confidence interval:\n",sep="")
  cat(" proportion: [",round(x$ci$bin[1],4),",",round(x$ci$bin[2],4),"]\n",sep="")
	cat(100 * x$call$level,"% asymptotic confidence interval with correction by Wilson:\n",sep="")
  cat(" proportion: [",round(x$ci$ac[1],4),",",round(x$ci$ac[2],4),"]\n",sep="")
  cat(100 * x$call$level,"% exact confidence interval by Clopper-Pearson:\n",sep="")
  cat(" proportion: [",round(x$ci$cp[1],4),",",round(x$ci$cp[2],4),"]\n\n",sep="")
 }
 invisible(x)
}
