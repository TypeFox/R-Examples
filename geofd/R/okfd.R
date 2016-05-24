okfd <-
function(new.coords, coords, data,
    smooth.type=NULL, nbasis=max(50,dim(data)[1]), argvals=seq(0, 1, len = dim(data)[1]), lambda=0,
    cov.model=NULL, fix.nugget=FALSE, nugget=0, fix.kappa=TRUE, kappa=0.5, max.dist.variogram=NULL)
# argnames=c("argument", "sites", "values"),
{

  # Loading required libraries
  ## require(fda)
  ## require(geoR)

  # Argument validation
  smooth.type <- match.arg(smooth.type, c("bsplines","fourier"))
  # The argument cov.model is validated if it is NOT null
  if(!is.null(cov.model)){
    cov.model <- match.arg(cov.model, c("spherical","exponential","gaussian","matern"))
  }
  if(is.null(new.coords)) stop("new.coords is not an optional parameter")
  if(ncol(new.coords)!=2) stop("new.coords must be an n x 2 matrix")
  # nbasis, argvals and lambda are validated in runtime on their corresponding using function
  # max.dist, fix.nugget and nugget does not seem to be validated

  # Argument type conversion
  new.coords <- as.matrix(new.coords)
  coords <- as.matrix(coords)

  # Number of sites
  s <- dim(data)[2]

  fdmodel <- .simple.fdmodel(new.coords, coords, data, smooth.type, nbasis, argvals, lambda, cov.model, fix.nugget, nugget, fix.kappa, kappa, max.dist.variogram)

  ##################################################################
  # Doing prediction
  ##################################################################

  prediction <- .okfd.predict(argvals, fdmodel$fdobjects$datafd, coords, new.coords, fdmodel$trace.vari.objects$best, fdmodel$emp.trace.vari$Eu.d)

  ##################################################################
  # Return
  ##################################################################

  return.list <- list(
    coords=coords,
    data=data,
    argvals=argvals,
    nbasis=nbasis,
    lambda=lambda,
    new.coords=new.coords,
    emp.trace.vari=fdmodel$emp.trace.vari,
    trace.vari=fdmodel$trace.vari.objects$best,
#    In parameter $u there are the distances 
#    Eu.d=Eu.d, 
    new.Eu.d=prediction$new.Eu.d,
    functional.kriging.weights=prediction$functional.kriging.weights,
    krig.new.data=prediction$pred,
    pred.var=prediction$var,
    trace.vari.array=fdmodel$trace.vari.objects$fitted,
    datafd=fdmodel$fdobjects$datafd
  )
  #return.list$argnames <- argnames
  class(return.list) <- "geofd"

  return(return.list)

}
