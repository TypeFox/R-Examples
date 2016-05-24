"print.binDesign" <-
function(x,...)

{

args<-list(...)
if(is.null(args$digits))
 {digits<-4}
else{digits<-args$digits}

 if(x$alternative=="less")
  {alt.hyp <- paste("true proportion is less than",x$p.hyp)
   ptrue <- paste("assumed true proportion","=", x$p.hyp - x$delta)}

 if(x$alternative=="greater")
  {alt.hyp <- paste("true proportion is greater than", x$p.hyp)
   ptrue <- paste("assumed true proportion","=", x$p.hyp + x$delta)}

 if(x$alternative=="two.sided")
  {alt.hyp <- paste("true proportion is not equal to",x$p.hyp )
   ptrue <- paste("assumed true proportion","=", x$p.hyp - x$delta,"or", x$p.hyp + x$delta)}


 if(x$power.reached==TRUE )
  {cat("Power was reached for n","=",x$nout,"\n")
   cat("with power","=",signif(x$powerout, digits),"\n")

   cat("alternative hypothesis:", alt.hyp ,"\n")
   cat( ptrue, "\n")
  }

 if(x$power.reached==FALSE && x$powerout!=0)
  {cat("Power of", x$power,"was not reached in the range of n","=",min(x$nit),",",max(x$nit),"\n")
   cat("Maximal power was reached for n","=",x$nout,"with power","=",signif(x$powerout, digits),"\n")

   cat("alternative hypothesis:", alt.hyp ,"\n")
   cat( ptrue, "\n")
  }


 if(x$power.reached==FALSE && x$powerout==0)
  {cat("Power of", x$power,"can not be reached in the range of n = ",min(x$nit),", ",max(x$nit),",\n")
   cat("null hypothesis can not be rejected for any sample size n in this range.","\n")

   cat("alternative hypothesis:", alt.hyp ,"\n")
   cat( ptrue, "\n")
  }

invisible(x)

}

