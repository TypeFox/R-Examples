#' Parameter Estimation for Persistence and Search Proficiency
#' 
#' Finds Maximum Likelihood Estimates Weibull persistence parameters,
#' and for exponentially decreasing search proficiency.
#' 
#' @param rd list output from acme::read.data()
#' @param fname file name to which the output parameters are saved
#' @return \code{acme.est} returns a list with the following components:
#' \item{params}{5-element vector: alpha and rho parameters for the 
#' Weibull persistence distribution, a and b parameters for the 
#' exponentially decreasing search proficiency, and bt as the bleed-through
#' rate}
#' \item{info}{list of select system information inherited from 
#' acme::read.data() output}
#' 

#    Model Implementation:
#    acme.est()    Return and store in "acme.est" vector of param ests
##################################################################
# Five Parameter Estimates
#
acme.est <- function(rd, fname="acme.est") {
  if(missing(rd))
    print("Usage: acme.est(rd, <output.file>)");
  scav <- mle.wei(rd,  v=TRUE);
  srch <- mle.srch(rd, v=TRUE);
  est  <- list(params=c(scav$alp, scav$rho,
               a=srch$a.hat, b=srch$b.hat, bt=srch$bt.hat),
               info=rd$Info);
  save(est, file=fname);
  return(est);
}