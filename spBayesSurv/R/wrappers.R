
.DDPplots <- function ( xpred_, ygrid_, beta_, sigma2_, w_, probs_ ) {
  .Call("DDPplots", xpred_, ygrid_, beta_, sigma2_, w_, probs_, PACKAGE = "spBayesSurv")
}

.DistMat <- function ( si_, sj_ ) {
  .Call("DistMat", si_, sj_, PACKAGE = "spBayesSurv")
}

.CoxPHplots <- function ( xpred_, tgrid_, beta_, h_, d_, probs_ ) {
  .Call("CoxPHplots", xpred_, tgrid_, beta_, h_, d_, probs_, PACKAGE = "spBayesSurv")
}

baseline = function (...) 
{
  words <- as.character((match.call())[-1]);
  allf <- cbind(...);
  nterms <- ncol(allf);
  if (is.null(names(allf))){
    argname <- words[1:nterms]
  }else{
    argname <- ifelse(names(allf)=="", words[1:nterms], names(allf));
  }
  allf
}