grcentral <- function(par, userfn, fbase=NULL, env=optsp, ...) {
   # Central difference gradient approximation
   if (is.null(fbase)) fbase <- userfn(par, ...)  # ensure we function value at par
   eps <- env$deps
   df <- rep(NA, length(par))
   teps <- 0.5 * eps * (abs(par) + eps)
   for (i in 1:length(par)) {
      dp <- par
      dp[i]<-dp[i]+teps[i]
      dm <- par
      dm[i]<-dm[i]-teps[i]
      df[i] <- 0.5*(userfn(dp, ...) - userfn(dm,...))/teps[i]
   }
   df
}
 
