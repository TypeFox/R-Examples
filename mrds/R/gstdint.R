#' Integral of pdf of distances
#'
#' Computes the integral of distpdf with scale=1 (stdint=TRUE) or specified
#' scale (stdint=FALSE)
#'
#' @param x lower,upper value for integration
#' @param ddfobj distance detection function specification
#' @param index specific data row index
#' @param select logical vector for selection of data values
#' @param width truncation width
#' @param standardize logical used to decide whether to divide through by the
#'   function evaluated at 0
#' @param point logical to determine if point count(TRUE) or line
#'   transect(FALSE)
#' @param stdint if TRUE, scale=1 otherwise specified scale used
#' @param doeachint if TRUE perform integration using \code{\link{integrate}}
#' @return vector of integral values of detection function
#' @note This is an internal function that is not intended to be invoked
#'   directly.
#' @author Jeff Laake and David L Miller
#' @keywords utility
#' @importFrom stats pnorm smooth.spline
gstdint <- function(x, ddfobj, index=NULL,select=NULL,width,
                   standardize=TRUE, point=FALSE, stdint=TRUE, doeachint=FALSE){

  if(!is.matrix(x)){
    x <- matrix(x, ncol=2)
  }
  ## NB this is not the integral of the PDF OR the detection function but
  ##    rather the integral of:
  ##       g(x)/w             for line transects
  ##       2*r*g(r)/width^2   for point transects

  # Set of observations for computation of detection function can
  # be specified with logical (select) and numeric (index) values.
  # Either or both can be specified although the latter is unlikely.
  if(is.null(select)){
    # use all
    if(is.null(index)){
      scale.dm <- ddfobj$scale$dm
    }else{
      # use only those with specific indices
      scale.dm <- ddfobj$scale$dm[index,,drop=FALSE]
    }
  }else{
    # Use those with select=TRUE
    if(is.null(index)){
      scale.dm <- ddfobj$scale$dm[select,,drop=FALSE]
    }else{
      # use the numeric index within those with select=TRUE
      scale.dm <- ddfobj$scale$dm[select,,drop=FALSE][index,,drop=FALSE]
    }
  }

  # when we have half-normal, key only use the exact analytic expression
  # for the integral using the error function/analytic expression
  if(ddfobj$type=="hn" & is.null(ddfobj$adjustment) & !doeachint){


    key.scale <- scalevalue(ddfobj$scale$parameters,scale.dm)
    if(point){
      # analytic expression for integral of 2*r*g(r)/width^2 when
      #  g(r) is half-normal
      int <- (2*(key.scale^2*exp(-x[,1]^2/(2*key.scale^2))-
                 key.scale^2*exp(-x[,2]^2/(2*key.scale^2))))/width^2
    }else{
      # define the error function in terms of pnorm
      erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
      # analytic expression for integral of g(x)/w when g(x) is half-normal

      # let's speed the above up a bit
      # when the left limit is zero (i.e. no left truncation) then just
      # need one call to erf
      if(all(x[,1]==0)){
        int <- (1/width)*(sqrt(pi/2)*key.scale*(erf(x[,2]/(key.scale*sqrt(2)))))

      }else{
        int <- (1/width)*sqrt(pi/2)*key.scale*(-erf(x[,1]/(key.scale*sqrt(2)))+
                                      erf(x[,2]/(key.scale*sqrt(2))))
      }
    }
    return(int)
  }else if(ddfobj$type=="hr" & is.null(ddfobj$adjustment) & !doeachint & !point){

    # do the spline shortcut
    # note that this will only work for a given shape parameter

    # Documentation by Jeff Laake from a previous incarnation of this code:
    # uses the cgftab which is a spline fitted to a table of standardized
    # integrals and the value is interpolated from the spline for each
    # observation.
    # This used to speed up integration of the detection function with scale
    # covariates. The detection function is integrated at a series of points
    # from 0 to W and then a spline is fitted to the computed values which are
    # cumulative (integral from 0 to x < integral from 0 to x+dx fr dx>0). The
    # spline is then used to predict values of the integral which depend on the 
    # scale which can depend on observation specific covariates.

    # set up the integration grid (these are values of x/sigma)
    xx <- exp(0.05*(1:100))-1
    xx <- cbind(c(0,xx[1:99]),xx)
    # fix scale and shape matrices
    ddfobj$scale$dm <- matrix(1, nrow=nrow(xx))
    if(!is.null(ddfobj$shape)) ddfobj$shape$dm <- matrix(1, nrow=nrow(xx))

    # Create cumulative sums of values of g(x) integrals from grid (xx)
    y <- cumsum(gstdint(xx, ddfobj=ddfobj, index=1:nrow(xx), width=width,
                        standardize=standardize, point=point, doeachint=TRUE))
    # Return smoothed spline of cumulative integral values
    spp <- smooth.spline(c(0,xx[,2]), c(0, y))

    # get the scale parameter
    xscale <- scalevalue(ddfobj$scale$parameters, scale.dm)

    # do the integration via predict.smooth.spline
    integrals <- predict(spp, as.vector(x[,2]/xscale))$y -
                 predict(spp, as.vector(x[,1]/xscale))$y

    # rescale for point or line
    if(!point){
      integrals <- xscale*integrals
    }else{
      integrals <- xscale^2*integrals
    }

    return(integrals)
  }else{
    # loop over the integration ranges, calculating integrals
    res <- rep(NA, nrow(x))
    for(i in 1:nrow(x)){
      res[i] <- integrate(distpdf, lower = x[i,1], upper = x[i,2], width=width,
                          ddfobj=ddfobj, select=select[i], index=index[i],
                          rel.tol = 1e-7,standardize=standardize,
                          stdint=stdint,point=point)$value
    }
    return(res)
  }
}
