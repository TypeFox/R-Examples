#' Summarize a ds.mixture object
#'
#' Summarize a ds.mixture object. The function provides information on parameter
#' estimates, estimates of the abundance in the covered area and the average 
#' detectability and their respective standard errors and coefficients of 
#' variation. 
#'
#' @aliases summary.ds.mixture
#' @method summary ds.mixture
#' @S3method summary ds.mixture
#' @export
#' 
#' @param object A fitted mixture model detection function object.
#' @param ... Anything, but it will be ignored.
#' @return a summary of a \code{\link{ds.mixture}} object.
#'
#' @references
#' Miller, D.L. and L. Thomas (in prep.). Mixture model distance sampling detection functions.
#'
#' @author David L. Miller
#' @examples
#' library(mmds)
#' set.seed(0)
#' ## simulate some line transect data from a 2 point mixture
#' sim.dat<-sim.mix(c(-0.223,-1.897,inv.reparam.pi(0.3)),2,100,1)
#' ## fit the model
#' fit.sim.dat<-fitmix(sim.dat,1,2)
#' ## what happened?
#' summary(fit.sim.dat) 
#'
summary.ds.mixture<-function(object,...){
# See print.summary.ds for how this gets printed

   model<-object

   ans <- list()
   # Number of observations
   ans$n <- length(model$distance)

   # number of mixture components
   ans$mix.terms<-model$mix.terms

   ans$coeff<-list()

   # Parameter estimates and their standard errors
   ans$coeff$pars<-data.frame(estimate=model$pars,se=model$pars.se)

   # separate mixture proportions
   gp<-getpars(model$pars,model$mix.terms,model$zdim,model$z)
   ans$coeff$mix.prop<-gp$mix.prop
   names(ans$coeff$mix.prop)<-paste("pi_",1:model$mix.terms,sep="")

   # AIC
   ans$aic <- model$aic

   # Truncation distance
   ans$width <- model$width
  
   # point or line transects?
   if(model$pt){
      ans$ttype<-"point"
   }else{
      ans$ttype<-"line"
   }

   # average p
   ans$average.p<-model$pa
   ans$average.p.se<-model$pa.se
   ans$average.p.cv<-model$pa.se/model$pa
   # abundance
   ans$Nhat <- model$N
   ans$Nhat.se <- model$N.se
   ans$Nhat.cv <- model$N.se/model$N

   # set the class and return
   class(ans) <- "summary.ds.mixture"
   return(ans)
}
