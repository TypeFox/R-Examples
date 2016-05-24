#'  Detection functions
#'
#' Various functions used to specify key and adjustment functions for
#' detection functions.
#'
#' Multi-covariate detection functions (MCDS) are represented by a function
#' \eqn{g(x,w,\theta)} where x is distance, z is a set of covariates and
#' \eqn{\theta} is the parameter vector.  The functions are defined such that
#' \eqn{g(0,w,\theta)=1} and the covariates modify the scale \eqn{(x/\sigma)}
#' where a log link is used to relate \eqn{\sigma} to the covariates,
#' \eqn{\sigma=exp(\theta*w)}. A CDS function is obtained with a constant
#' \eqn{\sigma} which is equivalent to an intercept  design matrix, z.
#'
#' \code{detfct} will call either a gamma, half-normal, hazard-rate or uniform
#' function only returning the probability of detection at that distance. In
#' addition to the simple model above, we may specify adjustment terms to fit
#' the data better. These adjustments are either Cosine, Hermite and simple
#' polynomials.
#' These are specified as arguments to \code{detfct}, as detailed below.
#'
#' \code{detfct} function which calls the others and assembles the final result using either key(x)[1+series(x)] or (key(x)[1+series(x)])/(key(0)[1+series(0)]) (depending on the value of \code{standardize}).
#'
#' \code{keyfct.*} functions calculate key function values and \code{adjfct.*} calculate adjustment term values.
#'
#' \code{scalevalue} for either detection function it computes the scale with
#' the log link using the parameters and the covariate design matrix
#'
#' \code{fx}, \code{fr} non-normalized probability density for line transects and point
#' counts respectively
#'
#' @aliases detfct adjfct.cos adjfct.herm hermite.poly adjfct.poly keyfct.hn
#'  keyfct.hz keyfct.gamma scalevalue fx fr distpdf
#' @usage detfct(distance, ddfobj, select=NULL, index=NULL, width=NULL,
#'               standardize = TRUE, stdint=FALSE)
#'
#' adjfct.cos(distance, scaling = 1, adj.order, adj.parm = NULL, adj.exp=FALSE)
#'
#' adjfct.poly(distance, scaling = 1, adj.order, adj.parm = NULL, adj.exp=FALSE)
#'
#' adjfct.herm(distance, scaling = 1, adj.order, adj.parm = NULL, adj.exp=FALSE)
#'
#' scalevalue(key.scale, z)
#'
#' keyfct.hn(distance, key.scale)
#'
#' keyfct.hz(distance, key.scale, key.shape)
#'
#' keyfct.gamma(distance, key.scale, key.shape)
#'
#' fx(distance,ddfobj,select=NULL,index=NULL,width=NULL,
#'    standardize=TRUE,stdint=FALSE)
#'
#' fr(distance,ddfobj,select=NULL,index=NULL,width=NULL,
#'    standardize=TRUE,stdint=FALSE)
#'
#' distpdf(distance,ddfobj,select=NULL,index=NULL,width=NULL,standardize=TRUE,
#'            stdint=FALSE,point=FALSE)
#'
#' @param distance  vector of distances
#' @param ddfobj distance sampling object (see \code{\link{create.ddfobj}})
#' @param z design matrix for scale function
#' @param select logical vector for selection of data values
#' @param index specific data row index
#' @param key.scale vector of scale values
#' @param key.shape vector of shape values
#' @param adj.order vector of adjustment orders
#' @param adj.parm vector of adjustment parameters
#' @param width truncation width
#' @param standardize logical used to decide whether to divide through by the
#'  function evaluated at 0
#' @param scaling the scaling for the adjustment terms
#' @param stdint logical used to decide whether integral is standardized
#' @param point if TRUE, point counts; otherwise line transects
#' @param adj.exp if TRUE uses exp(adj) for adjustment to keep f(x)>0
#' @return
#' For \code{detfct}, the value is a vector of detection probabilities for the
#'   input set of x and z.
#' For \code{keyfct.*}, vector of detection probability for that
#'   key function at x.
#' For \code{adjfct.*}, vector of the value of the
#'   adjustment term at x.
#' For \code{scalevalue}, the value is a vector of the computed scales for the
#'   design matrix z.
#' @author Jeff Laake, David Miller
#' @seealso  \code{\link{mcds}},  \code{\link{cds}}
#' @references  Marques and Buckland 2004
#' Laake and Borchers 2004. in Buckland et al 2004.
#' Becker, E. F. and P. X. Quang. 2009. A gamma-shaped detection function for
#' line transect surveys with mark-recapture and covariate data. Journal of
#' Agricultural Biological and Environmental Statistics 14:207-223.
distpdf <- function(distance,ddfobj,select=NULL,index=NULL,width=NULL,
                    standardize=TRUE,stdint=FALSE,point=FALSE){

 if(!point){
   return(fx(distance=distance,ddfobj=ddfobj,select=select,index=index,
             width=width,standardize=standardize,stdint=stdint))
 }else
   return(fr(distance=distance,ddfobj=ddfobj,select=select,index=index,
             width=width,standardize=standardize,stdint=stdint))
}


fx <- function(distance,ddfobj,select=NULL,index=NULL,width=NULL,
               standardize=TRUE,stdint=FALSE){
  return(detfct(distance,ddfobj,select,index,width,standardize,stdint)/width)
}

fr <- function(distance,ddfobj,select=NULL,index=NULL,width=NULL,
               standardize=TRUE,stdint=FALSE){
  return(detfct(distance,ddfobj,select,index,width,standardize,stdint)*
         2*distance/width^2)
}


detfct <- function(distance,ddfobj,select=NULL,index=NULL,width=NULL,
    standardize=TRUE,stdint=FALSE){

  # Set of observations for computation of detection function can
  # be specified with logical (select) and numeric (index) values.
  # Either or both can be specified although the latter is unlikely.
  if(is.null(select)){
    # use all
    if(is.null(index)){
      scale.dm <- ddfobj$scale$dm
      shape.dm <- ddfobj$shape$dm
    }else{
      # use only those with specific indices
      scale.dm <- ddfobj$scale$dm[index,,drop=FALSE]
      shape.dm <- ddfobj$shape$dm[index,,drop=FALSE]
    }
  }else{
    # Use those with select=TRUE
    if(is.null(index)){
      scale.dm <- ddfobj$scale$dm[select,,drop=FALSE]
      shape.dm <- ddfobj$shape$dm[select,,drop=FALSE]
    }else{
      # use the numeric index within those with select=TRUE
      scale.dm <- ddfobj$scale$dm[select,,drop=FALSE][index,,drop=FALSE]
      shape.dm <- ddfobj$shape$dm[select,,drop=FALSE][index,,drop=FALSE]
    }
  }

  # Key function
  key <- ddfobj$type

  # calculate the key scale
  if(stdint){
    if(is.null(index)){
      key.scale <- rep(1,nrow(scale.dm))
    }else{
      key.scale <- 1
    }
  }else{
    if(!is.null(ddfobj$scale)){
      key.scale <- scalevalue(ddfobj$scale$parameters,scale.dm)
    }
  }

  # calculate the key shape
  if(!is.null(ddfobj$shape)){
    key.shape <- scalevalue(ddfobj$shape$parameters,shape.dm)
  }

  # for gamma shape parameter must be >1, see Becker and Quang (2009) p 213
  if(key=="gamma"){
    key.shape <- key.shape+1
    key.shape[key.shape==1] <- key.shape[key.shape==1]+0.000001
  }

  # 19-Jan-06 jll; added proper standardize code to get std integral.
  #  I left standardize code below in case it is needed for adjustment fcts
  #
  # Decide on keyfct.* and run.
  if(key == "hn"){
    key.vals <- keyfct.hn(distance, key.scale)
  }else if(key == "hr"){
    key.vals <- keyfct.hz(distance, key.scale, key.shape)
  }else if(key == "unif"){
    key.vals <- rep(1/width,length(distance))
  }else if(key == "gamma"){
    key.vals <- keyfct.gamma(distance, key.scale, key.shape)
  }else if(key == "th1"){
    key.vals <- keyfct.th1(distance, key.scale, key.shape)
  }else if(key == "th2"){
    key.vals <- keyfct.th2(distance, key.scale, key.shape)
  }

  # Adjustment functions
  # If we are using adjustment terms.
  if(!is.null(ddfobj$adjustment)){
    adj.series <- ddfobj$adjustment$series
    adj.scale <- ddfobj$adjustment$scale
    adj.order <- ddfobj$adjustment$order
    adj.parm <- ddfobj$adjustment$parameters
    adj.exp <- ddfobj$adjustment$exp

    # Find out if we are scaling by width or by key scale
    if(adj.scale == "width"){
      scaling <- width
    }else{
      scaling <- key.scale
    }

    ## Decide on adjustment term and run.
    # If we have simple polynomials
    if(adj.series == "poly"){
      adj.vals <- adjfct.poly(distance,scaling,adj.order,adj.parm,adj.exp)
    # Hermite polynomials
    }else if(adj.series == "herm"){
      adj.vals <- adjfct.herm(distance,scaling,adj.order,adj.parm,adj.exp)
    # Cosine series
    }else if(adj.series == "cos"){
      adj.vals <- adjfct.cos(distance,scaling,adj.order,adj.parm,adj.exp)
    }

    # If we have adjustment terms then we need to standardize the detection
    # function. So find the values for the key and adjustment terms at 0

    # dlm 25-Aug-05  This causes a division by zero error in the optimization
    #    so lets only do it when we need to, it cancels in the
    #    likelihood anyway.
    if(standardize == TRUE){
      if(key == "hn"){
        key.val.0 <- keyfct.hn(rep(0,length(distance)), key.scale)
      }else if(key == "hr"){
        key.val.0 <- keyfct.hz(rep(0,length(distance)), key.scale, key.shape)
      }else if(key == "gamma"){
        # for the gammma, use apex.gamma to find the apex first, then eval
        # need to update the scale to be +1 in this apex call
        ddfobj$shape$parameters <- log(exp(ddfobj$shape$parameters)+1)
        g.apex <- as.vector(apex.gamma(ddfobj))[1]
        key.val.0 <- keyfct.gamma(g.apex,key.scale, key.shape)
      }else if(key == "unif"){
        key.val.0 <- rep(1/width,length(distance))
      }else if(key == "th1"){
        key.val.0 <- keyfct.th1(rep(0,length(distance)), key.scale, key.shape)
      }else if(key == "th2"){
        key.val.0 <- keyfct.th2(rep(0,length(distance)),key.scale, key.shape)
      }
      if(adj.series == "poly"){
        adj.val.0 <- adjfct.poly(rep(0,length(distance)),scaling,
                                 adj.order,adj.parm,adj.exp)
      }else if(adj.series == "herm"){
        adj.val.0 <- adjfct.herm(rep(0,length(distance)),scaling,
                                 adj.order,adj.parm,adj.exp)
      }else if(adj.series == "cos"){
        adj.val.0 <- adjfct.cos(rep(0,length(distance)),scaling,
                                adj.order,adj.parm,adj.exp)
      }

      # Now return the standardized value of the detection function
      return((key.vals*(1+adj.vals))/(key.val.0*(1+adj.val.0)))

    }else{
      return(key.vals*(1+adj.vals))
    }
  }else{
    # If we have no adjustment terms then just return the key value.
    return(key.vals)
  }
}
