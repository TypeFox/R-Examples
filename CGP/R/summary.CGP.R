summary.CGP <-
function(object,...){
  
  val<-list(call=object$call,
            Lambda=object$lambda,
            Theta=object$theta,
            Alpha=object$alpha,
            Bandwidth=object$bandwidth,
            rmscv=object$rmscv,
            mu=object$mu,
            tau2=object$tau2,
            best.start=object$beststart,
            objval=object$objval)
  class(val)<-"summary.CGP"
  return(val)
  
}
