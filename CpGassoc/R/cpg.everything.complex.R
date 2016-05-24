cpg.everything.complex <-
function(x,first=TRUE,logit.transform,Problems,beta.vales=NULL,logitperm=FALSE,...) {
  warn<-c("\nThere are beta values <= 0 or >=1\n- these have been set to missing by the logit transform.", 
                   "\nThere are beta values <0 or >1!",
                   "\n10% of beta values <0 or >1.\nAre you sure that the data have been entered correctly?")
  ranwarn<-c("The random effects model failed to converge at ",
                   " sites; test statistics for these\n",
                   "sites were set to NA.  This convergence failure can be avoided by\nusing a fixed effects model.")
  if(logitperm==FALSE) {
  if(first) {
    if(logit.transform & length(Problems)!=0) {
       warning(warn[1])
               }
    if(!logit.transform) {
      if(length(Problems!=0)) {
        warning(warn[2])
        }
      if(length(Problems)/length(as.matrix(beta.vales)) >.1) {
        warning(warn[3])
        }}}
  else {
    warning(ranwarn[1],Problems,ranwarn[2:3])
    }}}
