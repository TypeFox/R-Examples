grfwd <- function(par, userfn, fbase=NULL, env=optsp, ...) {
   # Forward different gradient approximation
   eps<-env$deps
   if (is.null(fbase)) fbase <- userfn(par, ...)  # ensure we have function value at par
   df <- rep(NA, length(par))
   teps <- eps * (abs(par) + eps)
   for (i in 1:length(par)) {
      dx <- par
      dx[i] <- dx[i] + teps[i]
      df[i] <- (userfn(dx, ...) - fbase)/teps[i]
   }
   df
}

