fnchk <- function(xpar, ffn, trace=0, ... ) {
# fnchk <- function(xpar, ffn, cctrl=list(trace=0), ... )
#  A function to check the nonlinear optimization file that is "ffn", with gradient gr
#  The intention is to automatically test the gradient, hessian, Jacobian, Jacobian second derivatives,
#    as well as bounds
#
#  This function can take-in multiple starting values
#
# Input:
#  xpar = a vector of starting values
#  ffn = objective function (assumed to be sufficiently differentiable)
#  cctrl = a list of control information FOR THE CHECKING PROGRAM. See Details.
#          The name has been changed from control to avoid confusion with control list in optim/optimx
#  ...     = other arguments to the function identified by fname
#
#  NOTE: bounds do NOT appear here.
#
# Output:
#      fval
#      infeasible
#      excode
#      msg
#
#  Author:  John Nash
#  Date: Sept 18, 2011
#################################################################
  maxard10<-function(one, two) { # get maximum absolute relative difference scaled by 10. in denominator
# This internal function is used to make comparisons using a relative difference, but avoiding zero divide
    result<-max(abs((one-two)/(abs(one)+abs(two)+10.0)))
    return(result)
  }
#########
   if (trace>3) {
      cat("fnchk: ffn =\n")
      print(ffn)
      cat("xpar:")
      print(xpar)
      cat("dots:")
      print(list(...))
   }
   infeasible<-FALSE # set value OK, then alter if not feasible later
   excode <- 0 # ditto
   msg <- "fnchk OK" # ditto
   if (trace>1) {
      cat("about to call ffn(xpar, ...)\n")
      cat("ffn:")
      print(ffn)
      cat("xpar & dots:")
      print(xpar)
      print(list(...))
   }
   test<-try(fval<-ffn(xpar, ...)) # 
   if (trace>1) {
      cat("test in fnchk:")
      print(test)
   }
   # Note: This incurs one EXTRA function evaluation because optimx wraps other methods
   if (inherits(test, "try-error") ) {
      fval<-NA
      attr(fval, "inadmissible")<-TRUE
   }
   if (trace > 0) {
      cat("Function value at supplied parameters =")
      print(fval) # Use "print" rather than "cat" to allow extra structure to be displayed
      print(str(fval))
      print(is.vector(fval))
   }
   if (!is.null(attr(fval,"inadmissible")) && (attr(fval, "inadmissible"))) {
      infeasible <- TRUE
      excode <- -1
      msg <- "Function evaluation returns INADMISSIBLE"
      if (trace>0) cat(msg,"\n")
   }

   # Also check that it is returned as a scalar
   if (is.vector(fval)) {
      if (length(fval)>1) { # added 120411
        excode <- -4
        msg <- "Function evaluation returns a vector not a scalar"
        infeasible <- TRUE
        if (trace>0) cat(msg,"\n")
      }
   }

   if (is.list(fval)) {
      excode <- -4
      msg <- "Function evaluation returns a list not a scalar"
      infeasible <- TRUE
      if (trace>0) cat(msg,"\n")
   }

   if (is.matrix(fval)) {
      excode <- -4
      msg <- "Function evaluation returns a matrix list not a scalar"
      infeasible <- TRUE
      if (trace>0) cat(msg,"\n")
   }

   if (is.array(fval)) {
      excode <- -4
      msg <- "Function evaluation returns an array not a scalar"
      infeasible <- TRUE
      if (trace>0) cat(msg,"\n")
   }

   if ((length(fval)!=1) && !(is.vector(fval))) { #this may never get executed
      excode <- -4
      msg <- "Function returned not length 1, despite not vector, matrix or array"
      infeasible <- TRUE
      if (trace>0) cat(msg,"\n")
   }

   if ( ! (is.numeric(fval)) ) {
      excode <- -1 
      msg <- "Function evaluation returned non-numeric value"
      infeasible <- TRUE
      if (trace>0) cat(msg,"\n")
   }

   if (is.infinite(fval) || is.na(fval)) {
      excode <- -1 
      msg <- "Function evaluation returned Inf or NA (non-computable)"
      infeasible <- TRUE
      if (trace>0) cat(msg,"\n")
   }
   if (trace>0) cat("Function at given point=",fval,"\n")
   answer <- list(fval=fval, infeasible=infeasible, excode=excode, msg=msg)
}
### end of fnchk ***

