okfd.cv <-
function(coords, data, argnames=c("argument", "sites", "values"), one.model=TRUE, smooth.type=NULL, array.nbasis=max(50,dim(data)[1]), argvals=seq(0,1,len=dim(data)[1]), array.lambda=0, cov.model=NULL, fix.nugget=FALSE, nugget=0, fix.kappa=TRUE, kappa=0.5, max.dist.variogram=NULL)
{
  # Argument validation
  smooth.type <- match.arg(smooth.type, c("bsplines","fourier"))
  if(!is.vector(array.nbasis)) stop("the argument \"array.nbasis\" must be a vector")
  if(!is.vector(array.lambda)) stop("the argument \"array.lambda\" must be a vector")
  if(!is.vector(argnames) || length(argnames)!=3) stop("the argument \"argnames\" must be a character vector of length 3")
  if(!is.logical(one.model)) stop("the argument \"one.model\" must be logical, i.e. TRUE or FALSE")

  # Init local variables 
  m <- dim(data)[1]
  diff.argvals <- diff(argvals)
  s <- dim(coords)[1]
  mse.cv <- matrix(0,nrow=length(array.nbasis),ncol=length(array.lambda))
  krig.cv <- array(0,dim=c(length(array.nbasis),length(array.lambda),s,m),
  dimnames <- c("nbasis","lambda",argnames[2],argnames[1]))

  k <- 0
  k.opt <- 1
  l.opt <- 1
  fdmodels <- c()

  # Loop over all number of basis functions parameters
  for (nbasis in array.nbasis){
    k <- k+1
    l <- 0
    # Loop over all smoothing penalization parameters
    for (lambda in array.lambda){
      l <- l+1
      # Option 1, one model is estimated using all sites
      if(one.model==TRUE){
        fdmodel <- .simple.fdmodel(new.coords, coords, data, smooth.type, nbasis, argvals, lambda, cov.model, fix.nugget, nugget, fix.kappa, kappa, max.dist.variogram)
        fdmodels[[length(fdmodels)+1]] <- fdmodel
      }
      # Loop over all sites number
      for (i in 1:s){
        # Option 2, one model is estimated using 's-i' sites
        if(one.model==FALSE){
          new.coords <- as.data.frame(coords)[i,]
          coords <- as.data.frame(coords)[-i,]
          data <- data[,-i]
          fdmodel <- .simple.fdmodel(new.coords, coords, data, smooth.type, nbasis, argvals, lambda, cov.model, fix.nugget, nugget, fix.kappa, kappa, max.dist.variogram)
          fdmodels[[length(fdmodels)+1]] <- fdmodel
        }
        # Prediction for the site 'i' using the other 's-i' sites,
        # the model and smoothed data specified by fdmodel
        res.okfd <- .okfd.predict(argvals=argvals, datafd=fdmodel$fdobjects$datafd, coords=coords, new.coords=as.data.frame(coords)[i,], trace.vari=fdmodel$trace.vari.objects$best, Eu.d=fdmodel$emp.trace.vari$Eu.d)
        # Predicted data is saved 
        krig.cv[k,l,i,] <- res.okfd$krig.new.data
        # An error measure is calculated
        aux <- (data[,i]-res.okfd$krig.new.data)^2
        aux <- diff.argvals * (aux[1:(m-1)]+aux[2:m])/2
        mse.cv[k,l] <- mse.cv[k,l] + sum( aux ) 
      }
      if (mse.cv[k,l] <= mse.cv[k.opt,l.opt]){
        k.opt <- k
        l.opt <- l
      }
    }
  }

  mse.cv.opt <- mse.cv[k.opt,l.opt]

  return(list(
    k.opt= k.opt,
    l.opt= l.opt,
    krig.cv= krig.cv,
    mse.cv= mse.cv,
    mse.cv.opt= mse.cv.opt,
    fdmodels=fdmodels
  ))

}
