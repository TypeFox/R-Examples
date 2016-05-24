grchk <- function(xpar, ffn, ggr, trace=0, testtol=(.Machine$double.eps)^(1/3), ...) {
   # check gradient code in ggr for function in ffn
   # assume function ffn already passes checks in fnchk
   if (is.null(trace)) trace <- 0 else trace<-trace
   if (is.null(testtol)) testtol<-(.Machine$double.eps)^(1/3) else testtol<- testtol
   gname <- deparse(substitute(ggr)) # ?? May be embedded in ugr
   if (trace>1) cat("Gradient test with tolerance = ",testtol,"\n")
   if (trace>1) cat("Analytic gradient uses function ",gname,"\n")
   fval <- ffn(xpar, ...)
   if (trace>1) {
      cat("function at parameters = ", fval," with attributes:\n")
      print(attributes(fval))
   }
   if (trace>1) cat("Compute analytic gradient\n")
   ga<-ggr(xpar, ...)
   if (trace>1) print(ga)
   if (trace>1) cat("Compute numeric gradient\n")
   gn <- grad(func=ffn, x=xpar, ...) # numerically approximated gradient
   if (trace>1) print(gn)
   # Now test for equality (090612: ?? There may be better choices for the tolerances.
   if (trace>0) {
      cat("gradient test tolerance = ",testtol,"  fval=",fval,"\n")
      cat(" compare to max(abs(gn-ga))/(1+abs(fval)) = ",max(abs(gn-ga))/(1+abs(fval)),"\n")
   }
   if (max(abs(gn-ga))/(1 + abs(fval)) >= testtol) { 
      if(trace>0) cat("Gradient function might be wrong - check it! \n")
      gradchk<-FALSE
   } else gradchk<-TRUE
   attr(gradchk, "ga")<-ga
   attr(gradchk, "gn")<-gn
   gradchk
} # end gradchk

