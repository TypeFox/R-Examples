"gp" <-
function(gridsize=c(64,64),specdens="matern.specdens",specdens.param=c(1,4),variance.param=1,const.fixed=FALSE){
  # creates a gp object, based on the chosen grid, the spectral density of the desired correlation function and parameters for the correlation function
  if(min(gridsize)<2 || sum(log(gridsize,2)%%1)){
    stop("Must have a grid of power of 2 in each dimension")
  }

  object = new.env(parent=globalenv())
  
  object$gridsize=gridsize
  object$d=length(gridsize)
  object$specdens.param=specdens.param
  object$variance.param=variance.param
  object$const.fixed=const.fixed     # if fixed, coefficient for constant basis function set to zero
  if(is.character(specdens)){
    object$specdens=get(specdens,envir=parent.frame())
    object$specdens.name=specdens
  } else{  # specdens arg should be a function
    if(!is.function(specdens)){
      stop("Required specdens function argument is not a valid function.")
    }
    object$specdens=specdens
    object$specdens.name=substitute(specdens)
  }
  d=object$d
  if(d==1){
    object$gridsize=c(object$gridsize,1)
  }
  gridsize=object$gridsize
  object$coeff=matrix(0,nrow=gridsize[1],ncol=gridsize[2])
  if(d>2){
    stop("Only one- and two-dimensional processes are enabled;\nplease specify gridsize as a scalar or vector of length two")
  }
  if(d==2 && gridsize[1]!=gridsize[2]){
    warning("Unequal grid sizes imply differing resolution in the different dimensions.\n  Use with caution unless you know what you are doing.\nThis is not a good way to deal with non-square regions -\ninstead, just map your domain to a portion of the\nsquare domain needed for the spectral approach.\n")
  }
  # specifying the Fourier frequencies
  omega1=seq(0,gridsize[1]/2,len=gridsize[1]/2+1)
  omega1=c(omega1,-(rev(omega1))[2:(length(omega1)-1)])
  if(d==2){
    omega2=seq(0,gridsize[2]/2,len=gridsize[2]/2+1)
    omega2=c(omega2,-(rev(omega2))[2:(length(omega2)-1)])
    object$omega=expand.grid(omega1=omega1,omega2=omega2)
  } else{
    object$omega=matrix(omega1,ncol=1)
  }
  class(object)="gp"
  calc.variances(object) # calculates prior variances based on correlation function and its parameters
  zero.coeff(object)
  updateprocess(object)
  object$num.blocks=0
  cat("Note that the spatial range parameter is interpreted based\non the process living on a (0,1)^d grid\n")
  object
}
