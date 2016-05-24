#' Stepwise selection of mixture components
#'
#' Uses AIC to select the number of mixture components.
#'
#' @param ds.object \code{\link{ds.mixture}} object.
#' @param max.terms Maximum number of mixture components to attempt to fit. 
#'    Default 4.
#' @return a \code{\link{ds.mixture}} model object of the best (AIC-wise) model.
#'
#' @section Details:
#' This routine is most useful during model building. Setting up a basic 1-point
#' mixture model and then running this on the object returns a model with the 
#' lowest AIC.
#' 
#' Progress will be printed to the screen.
#'
#' @author David L. Miller
#'
#' @export
#'
#' @examples
#' library(mmds)
#' set.seed(0)
#' ## simulate some line transect data from a 2 point mixture
#' sim.dat<-sim.mix(c(-0.223,-1.897,inv.reparam.pi(0.3)),2,100,1)
#' ## fit the model
#' fit.sim.dat.1<-fitmix(sim.dat,1,1)
#' ## find best AIC model
#' step.ds.mixture(fit.sim.dat.1)
step.ds.mixture<-function(ds.object,max.terms=4){
   # equivalent of step() for ds.mixture objects
   # adds 1 to the mix.terms every time, starts from
   # wherever it was set.

   cat("\nStepwise mixture component selection.\n")

   mix.terms<-ds.object$mix.terms
   aics<-ds.object$aic
   i<-1   
   delaic<-1e6

   model.list<-list()
   model.list[[1]]<-ds.object

   cat("\nMixture terms   AIC\n")

   cat("   ",mix.terms,"        ",ds.object$aic,"\n")
   mix.terms<-mix.terms+1

   while( (mix.terms<max.terms)){
#data,width,mix.terms=1,pt=FALSE,model.formula="~1",initialvalues=NULL,showit=0,ctrl.options=c(maxit=10000),opt.method="BFGS+SANN",usegrad=TRUE,ftype="hn"){
      this.model<-try(fitmix(data=ds.object$data,
                             ftype=ds.object$ftype,
                             model.formula=ds.object$model.formula,
                             mix.terms=mix.terms,
                             initialvalues=NULL,
                             width=ds.object$width,
                             showit=ds.object$showit,
                             ctrl.options=ds.object$ctrl.options, 
                             opt.method=ds.object$opt.method,
                             usegrad=ds.object$usegrad,
                             pt=ds.object$pt))

      if(class(this.model)!="try-error"){
         cat("   ",mix.terms,"        ",this.model$aic,"\n")

         model.list[[i+1]]<-this.model
         aics<-c(aics,this.model$aic)
         mix.terms<-mix.terms+1
         i<-i+1
         if(abs(this.model$aic-aics[i-1]) < 3) break
      }else{
         break
      }
   }

   best<-which.min(aics)

   cat("\nBest model was a ",model.list[[best]]$mix.terms,"-point mixture.\n\n",sep="")

   invisible(model.list[[best]])

}
