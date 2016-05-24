#' Variance propogation for DSM models
#'
#' Rather than use a bootstrap to calculate the variance in a \code{dsm} model,
#' use the clever variance propogation trick from Williams et al. (2011).
#'
#' The idea is to refit the spatial model but including the Hessian of the 
#' offset as an extra term. Variance estimates using this new model can then 
#' be used to calculate the variance of abundance estimates which incorporate 
#' detection function uncertainty. Further mathematical details are given in 
#' the paper in the references below.
#'
#' Many prediction grids can be supplied by supplying a list of 
#' \code{data.frame}s to the function.
#'
#' Note that this routine is only useful if a detection function has been used in the DSM.
#'
#' Based on (much more general) code from Mark Bravington and Sharon Hedley.
#'
#' @inheritParams dsm.var.gam
#' @return a list with elements
#'         \tabular{ll}{\code{model} \tab the fitted model object\cr
#'                      \code{pred.var} \tab variance of each region given
#'                      in \code{pred.data}.\cr
#'                      \code{bootstrap} \tab logical, always \code{FALSE}\cr
#'                      \code{pred.data} \tab as above\cr
#'                      \code{off.set} \tab as above\cr
#'                      \code{model}\tab the fitted model with the extra term\cr
#'                      \code{dsm.object} \tab the original model, as above\cr
#'                      \code{model.check} \tab simple check of subtracting the coefficients of the two models to see if there is a large difference\cr
#'                      \code{deriv} \tab numerically calculated Hessian of the offset\cr
#'                      }
#' @author Mark V. Bravington, Sharon L. Hedley. Bugs added by David L. Miller.
#' @references
#' Williams, R., Hedley, S.L., Branch, T.A., Bravington, M.V., Zerbini, A.N. and Findlay, K.P. (2011). Chilean Blue Whales as a Case Study to Illustrate Methods to Estimate Abundance and Evaluate Conservation Status of Rare Species. Conservation Biology 25(3), 526-535.
#' @export
#' @importFrom stats as.formula
#' @examples
#' \dontrun{
#'  library(Distance)
#'  library(dsm)
#'
#'  # load the Gulf of Mexico dolphin data (see ?mexdolphins)
#'  data(mexdolphins)
#'
#'  # fit a detection function and look at the summary
#'  hr.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance),
#'                 key = "hr", adjustment = NULL)
#'  summary(hr.model)
#'
#'  # fit a simple smooth of x and y
#'  mod1 <- dsm(N~s(x,y), hr.model, mexdolphins$segdata, mexdolphins$obsdata)
#'
#'  # Calculate the offset...
#'  off.set <- 444*1000*1000
#'
#'  # Calculate the variance
#'  mod1.var <- dsm.var.prop(mod1, mexdolphins$pred, off.set)
#'
#'  # summary over the whole area in mexdolphins$pred
#'
#'  # Plot a map of the CV
#'  #   need to format the prediction data with split
#'  mod1.var.map <- dsm.var.prop(mod1,
#'                  split(mexdolphins$pred,1:nrow(mexdolphins$pred)), off.set)
#'  plot(mod1.var.map)
#' }
#'
dsm.var.prop<-function(dsm.obj, pred.data,off.set,
    seglen.varname='Effort', type.pred="response") {

  is.gamm <- FALSE
  # if we have a gamm, then just pull out the gam object
  if(any(class(dsm.obj)=="gamm")){
    dsm.obj <- dsm.obj$gam
    is.gamm <- TRUE
  }

  # strip dsm class so we can use gam methods
  class(dsm.obj) <- class(dsm.obj)[class(dsm.obj)!="dsm"]

  # if all the offsets are the same then we can jsut supply 1 and rep it
  if(length(off.set)==1){
    if(is.null(nrow(pred.data))){
      off.set <- rep(list(off.set),length(pred.data))
    }else{
      off.set <- rep(off.set,nrow(pred.data))
    }
  }

  # make sure if one of pred.data and off.set is not a list we break
  # if we didn't have a list, then put them in a list so everything works
  if(is.data.frame(pred.data) & is.vector(off.set)){
    pred.data <- list(pred.data)
    off.set <- list(off.set)
  }else if(is.list(off.set)){
    if(length(pred.data)!=length(off.set)){
      stop("pred.data and off.set don't have the same number of elements")
    }
  }

  # pull out the ddf object
  ddf.obj <- dsm.obj$ddf

  # if there is no ddf object, then we should stop!
  # thanks to Adrian Schiavini for spotting this
  if(is.null(ddf.obj)){
    stop("No detection function in this analysis, use dsm.var.gam")
  }

  # this function changes the parameters in the ddf object
  tweakParams <- function(object, params) {
    if(missing(params)){
      return(object$par)
    }
    object$par <- params
    object$ds$aux$ddfobj <- mrds:::assign.par(object$ds$aux$ddfobj,params)
    return(object)
  }

  # function to find the derivatives of the offset
  funco <- function(p){
    # set the parameters to be p
    ipo <- tweakParams(ddf.obj, p)
    # calculate the offset
    ret <- log(2*as.vector(unique(predict(ipo,esw=TRUE, compute=TRUE)$fitted))*
            fo2data[[seglen.varname]])
    return(ret)
  }

  # pull out the data
  fo2data <- dsm.obj$data

  # find the derivatives
  p0 <- tweakParams(ddf.obj) # returns the parameters to numderiv
  firstD <- numderiv( funco, p0)

  # if the derivatives were zero, throw an error
  if(all(firstD==0)){
    stop('Doffset/Dpars==0... really??!!')
  }

  # now construct the extra term...
  formo <- dsm.obj$formula
  dmat.name <- '.D1'
  names.to.avoid <- unique( c( all.names( formo), names( fo2data)))
  while( dmat.name %in% names.to.avoid){
    dmat.name <- paste('.',dmat.name,sep="")
  }


  if(!is.gamm){
    # pull out the call
    formo[[3]] <- call( '+', formo[[3]], as.symbol(dmat.name))

    # put together the paraPen terms
    paraterm<-list(list(ddf.obj$hess))
    names(paraterm) <- dmat.name
    callo <- dsm.obj$call
    callo$paraPen <- c(callo$paraPen, paraterm)

    # insert the extra data into the frame
    fo2data[[ dmat.name]] <- firstD
  }else{
    # pull out the call
    callo <- eval(dsm.obj$gamm.call.list)

    # get the formula and make it a string
    formo <- paste(formo[[2]],formo[[1]],as.character(formo)[[3]],collapse="")

    # need to reparametrise
    S <- as.matrix(ddf.obj$hess)
    S.e <- eigen(S)
    if(is.matrix(firstD)){
      sqrt.D <- diag(sqrt(S.e$values))
      firstD <- firstD%*%sqrt.D

      rand.list <- list()

      # add the extra random effect term to the formula
      for(i in 1:ncol(firstD)){
        this.dmat.name <- paste0(dmat.name,i)
#        formo <- paste(formo," + s(",this.dmat.name,",bs=\"re\")",collapse="")
rand.list[[this.dmat.name]]<-pdMat(form=~1)
        fo2data[[ this.dmat.name]] <- firstD[,i]
      }
callo$random <- rand.list
    }else{
      sqrt.D <- sqrt(S.e$values)
      firstD <- firstD*sqrt.D

      # add the extra random effect term to the formula
#      formo <- paste(formo," + s(",dmat.name,",bs=\"re\")",collapse="")
rand.list <- list()
rand.list[[dmat.name]]<-pdMat(form=~1)
      fo2data[[ dmat.name]] <- firstD
callo$random <- rand.list
    }

    # make the formula a formula
    formo <- as.formula(formo)
  }

  # put the right objects into the call object
  callo$formula <- formo
  callo$family <- dsm.obj$family
  callo$data <- fo2data
  callo$knots <- dsm.obj$knots

  ## run the model
  if(!is.gamm){
    fit.with.pen <- with(dsm.obj,withCallingHandlers(eval(callo),
                                         warning=matrixnotposdef.handler))
  }else{
    fit.with.pen <- do.call("gamm",callo)
    fit.with.pen <- fit.with.pen$gam
  }
  # strip dsm class so we can use gam methods
  class(fit.with.pen) <- class(fit.with.pen)[class(fit.with.pen)!="dsm"]

  # Diagnostic from Mark
  # check that the fitted model isn't too different, used in summary()
  model.check <- summary(fitted(fit.with.pen) - fitted(dsm.obj))

  cft <- coef(fit.with.pen)
  preddo <- list(length(pred.data))
  dpred.db <- matrix(0, length(pred.data), length(cft))

  # depending on whether we have response or link scale predictions...
  if(type.pred=="response"){
    tmfn <- dsm.obj$family$linkinv
    dtmfn <- function(eta){vapply(eta, numderiv, numeric(1), f=tmfn)}
  }else if(type.pred=="link"){
    tmfn <- identity
    dtmfn <- function(eta){1}
  }

  # loop over the prediction grids
  for(ipg in seq_along(pred.data)){
    # if we have a single paramter model (e.g. half-normal) need to be careful
    if(is.matrix(firstD)){
      if(!is.gamm){
        pred.data[[ipg]][[dmat.name]] <- matrix(0,
                                                nrow(pred.data[[ipg]]),
                                                ncol(firstD))
      }else{
        for(i in 1:ncol(firstD)){
          this.dmat.name <- paste0(dmat.name,i)
          pred.data[[ipg]][[this.dmat.name]] <- matrix(0,
                                                       nrow(pred.data[[ipg]]),
                                                       ncol(firstD))
        }
      }
    }else{
      pred.data[[ipg]][[dmat.name]] <- rep(0, nrow(pred.data[[ipg]]))
    }

    ### fancy lp matrix stuff
    # set the offset to be zero here so we can use lp
    pred.data[[ipg]]$off.set<-rep(0,nrow(pred.data[[ipg]]))

    lpmat <- predict(fit.with.pen, newdata=pred.data[[ ipg]], type='lpmatrix')
    lppred <- lpmat %**% cft

    # if the offset is just one number then repeat it enough times
    if(length(off.set[[ipg]])==1){
      this.off.set <- rep(off.set[[ipg]],nrow(pred.data[[ipg]]))
    }else{
      this.off.set <- off.set[[ipg]]
    }

    preddo[[ipg]] <-  this.off.set %**% tmfn(lppred)
    dpred.db[ipg,] <- this.off.set %**% (dtmfn(lppred)*lpmat)
    # explanation of the above line and why we find this derivative
    # BTW in 'varpred', there is a decoy option 'vmethod' which at the
    # moment has to be '"delta"', for how to deal with nonlinearity in the
    # "link" of 'fitobj'. Could be done by simulation instead, and that would
    # be more accurate (if you did enough). However, in my limited experience:
    # once you've got a CV so big that the delta-method doesn't work, then
    # your estimate is officially Crap and there is not much point in
    # expending extra effort to work out exactly how Crap!
  }

  # "'vpred' is the covariance of all the summary-things." - MVB
  # so we want the diagonals if length(pred.data)>1
  # A B A^tr
  vpred <- dpred.db %**% tcrossprod(vcov(fit.with.pen), dpred.db)

  if(is.matrix(vpred)){
    vpred <- diag(vpred)
  }

  result <- list(pred.var = vpred,
                 bootstrap = FALSE,
                 var.prop = TRUE,
                 pred.data = pred.data,
                 pred = preddo,
                 off.set = off.set,
                 model = fit.with.pen,
                 dsm.object = dsm.obj,
                 model.check = model.check,
                 deriv = firstD,
                 seglen.varname=seglen.varname,
                 type.pred=type.pred
                )

  class(result) <- "dsm.var"

  return(result)
}

####### this is all utility stuff below here, taken from Mark's packages

# from Mark Bravington's handy2
numderiv<-function (f, x0, eps = 1e-04, TWICE. = TRUE, param.name = NULL,
    ..., SIMPLIFY = TRUE)
{
    if (is.null(param.name)) 
        ff <- function(x, ...) f(x, ...)
    else ff <- function(x, ...) {
        ll <- c(list(x), list(...))
        names(ll)[1] <- param.name
        do.call("f", ll)
    }
    f0 <- ff(x0, ...)
    n <- length(x0)
    m <- matrix(0, length(f0), n)
    for (i in 1:n) {
        this.eps <- eps * if (x0[i] == 0)
            1
        else x0[i]
        m[, i] <- (ff(x0 + this.eps * (1:n == i), ...) - f0)/this.eps
    }
    if (!is.null(dim(f0)))
        dim(m) <- c(dim(f0), n)
    if (TWICE.) {
        mc <- match.call()
        mc$eps <- -eps
        mc$TWICE. <- FALSE
        m <- 0.5 * (m + eval(mc, sys.frame(sys.parent())))
    }
    if (any(dim(m) == length(m)) && SIMPLIFY)
        m <- c(m)
    return(m)
}

# from mvbutils
"%**%"<-function(x, y){
    dimnames(x) <- NULL
    dimnames(y) <- NULL
    if (length(dim(x)) == 2 && length(dim(y)) == 2 && dim(x)[2] ==
        1 && dim(y)[1] == 1)
        return(c(x) %o% c(y))
    if ((!is.null(dim(x)) && any(dim(x) == 1)))
        dim(x) <- NULL
    if ((!is.null(dim(y)) && any(dim(y) == 1)))
        dim(y) <- NULL
    if (is.null(dim(x)) && is.null(dim(y))) {
        if (length(x) == length(y))
            x <- x %*% y
        else {
            if ((length(x) != 1) && (length(y) != 1))
                stop(paste("lengths of x (",length(x),") and y (",
                  length(y),") are incompatible",sep=""))
            else x <- x * y
        }
    }
    else x <- x %*% y
    if ((!is.null(dim(x)) && any(dim(x) == 1)))
        dim(x) <- NULL
    x
}
