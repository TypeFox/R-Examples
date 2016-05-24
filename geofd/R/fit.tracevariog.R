fit.tracevariog <-
function(emp.trace.vari, models, sigma2.0, phi.0, fix.nugget=FALSE, nugget=0, fix.kappa=TRUE, kappa=0.5, max.dist.variogram=NULL){

  # Argument validation
  if(missing(models) || is.null(models)){
    models <- c("spherical","exponential","gaussian","matern")
  }
  if(missing(max.dist.variogram) || is.null(max.dist.variogram)){
    max.dist.variogram <- max(emp.trace.vari$u)
  }
  if(!identical(class(emp.trace.vari),"variogram")) stop("the parameter emp.trace.vari must be of class variogram")
  if(!identical(dim(colors),dim(models))) stop("dimensions of parameters colors and models must be identical")
  if(!is.numeric(sigma2.0)) stop("the partial sill paramter must be a number")
  if(!is.numeric(phi.0)) stop("the range parameter must be a number")
  if(!is.logical(fix.nugget)) stop("the parameter fix.nugget must be logical")
  if(!is.numeric(nugget)) stop("the nugget parameter must be a number")
  if(!is.logical(fix.kappa)) stop("the parameter fix.kappa must be logical")
  if(!is.numeric(kappa)) stop("the kappa parameter must be a number")
  if(!is.numeric(max.dist.variogram)) stop("the max.dist.variogram parameter must be a number")

  # A comparison variable for choosing the best variogram model is loaded
  trace.vari <- c(NULL)
  trace.vari["value"] <- Inf
  # Is created an array where the variofit results will be pushed
  trace.vari.array <- c()
  # A loop is made in order to prove each model and choose the best using as criterion the variofit minimised sum of squares
  for(cont in c(1:length(models))){
    # The variofit function is called
    trace.vari.tmp <- variofit(emp.trace.vari, ini.cov.pars=c(sigma2.0,phi.0), max.dist=max.dist.variogram, fix.nugget=fix.nugget, nugget=nugget, fix.kappa=fix.kappa, kappa=kappa, cov.model=models[cont], messages=FALSE)
    # Each calculated trace variogram is pushed inside trace.vari.array
    trace.vari.array[[length(trace.vari.array)+1]] <- trace.vari.tmp
    # The last caculated variofit is compared against the last optimal variofit
    # If it is more optimal, then it is loaded in the variable wich will be returned
    if(as.numeric(trace.vari.tmp["value"])<as.numeric(trace.vari["value"])){
      trace.vari <- trace.vari.tmp
    }
  }

##################################################################
# Return:
##################################################################

  return(list(best=trace.vari, fitted=trace.vari.array))

}
