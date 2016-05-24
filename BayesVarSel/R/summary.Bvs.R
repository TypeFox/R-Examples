summary.Bvs <-
function(object,...){
  
  #we use object because it is requiered by S3 methods
  z <- object
  p <- z$p
  if (!inherits(object, "Bvs")) 
    warning("calling summary.Bvs(<fake-Bvs-x>) ...")
ans<-list()
#ans$coefficients <- z$betahat
#dimnames(ans$coefficients) <- list(names(z$lm$coefficients),"Estimate")

HPM <- z$HPMbin
MPM <- as.numeric(z$inclprob>=0.5)
astHPM <- matrix(" ",ncol=1,nrow=(p))
astMPM <- matrix(" ",ncol=1,nrow=(p))
astHPM[HPM==1] <- "*"
astMPM[MPM==1] <- "*"

incl.prob<-z$inclprob

summ.Bvs <- as.data.frame(cbind(round(incl.prob,digits=4),astHPM,astMPM))
dimnames(summ.Bvs)<- list(z$variables,c("Incl.prob.","HPM","MPM"))

ans$summary<-summ.Bvs  
ans$method <- z$method
ans$call<-z$call
class(ans) <- "summary.Bvs"
ans

}
