############################################
## S3-METHOD THAT PREDICTS FROM A STmodel ##
############################################
##Functions in this file:
## predict.STmodel         EX:ok
## print.predictSTmode     EX:ok
## plot.predictSTmodel     EX:ok

##' Compute the conditional expectations (i.e. predictions) at the unobserved
##' space-time locations. Predictions are computed for the space-time locations in
##' \code{object} and/or \code{STdata}, conditional on the observations (and
##' temporal trends) in \code{object} and parameters given in \code{x}.
##' 
##' In addition to computing the conditional expectation at a number of
##' space-time locations the function also computes predictions based on only
##' the regression part of the model as well as the latent beta-fields.
##'
##' Prediction are computed as the conditional expectation of a latent field
##' given observations. This implies that \code{E(X_i| Y_i) != Y_i}, with the
##' difference being due to smoothing over the nugget. Further two possible
##' variance can be computed (see below), \code{V(X_i|Y_i)} and
##' \code{V(X_i|Y_i)+nugget_i}. Here the nugget for unobserved locations needs
##' to be specified as an additional argument \code{nugget.nobs}. The two
##' variances correspond, losely, to confidence and prediction intervals.
##' 
##' Variances are computed if \code{pred.var=TRUE} point-wise variances for the
##' predictions (and the latent beta-fields) are 
##' computed. If instead \code{pred.covar=TRUE} the full covariance matrices for
##' each predicted time series is computed; this implies that the covariances between
##' temporal predictions at the same location are calculated but \emph{not}, due
##' to memory restrictions, any covariances between locations.
##' \code{beta.covar=TRUE} gives the full covariance matrices for the latent
##' beta-fields.
##'
##' If \code{transform!="none"} the field is assumed to be log-Gaussian and
##' expectations are transformed, and if \code{pred.var=TRUE} the mean squared
##' prediction errors are given.
##' 
##' @title Computes Conditional Expectation at Unobserved Locations
##' 
##' @param object \code{STmodel} object for which to compute predictions.
##' @param x Model parameters for which to compute the conditional
##'   expectation. Either as a vector/matrix or an \code{estimateSTmodel} from
##'   \code{\link{estimate.STmodel}}.
##' @param STdata \code{STdata}/\code{STmodel} object with locations/times for
##'   which to predict. If not given predictions are computed for locations/times
##'   in \code{object}
##' @param Nmax Limits the size of matrices constructed when computing
##'   expectations. Use a smaller value if memory becomes a problem.
##' @param only.pars Compute only the regression parameters (using GLS) along
##'   with the related variance.
##' @param nugget.unobs Value of nugget at unonserved locations, either a scalar
##'   or a vector with one element per unobserved site. \strong{NOTE:} All sites in
##'   \code{STdata} are considered unobserved!
##' @param only.obs Compute predictions at only locations specified by
##'   observations in \code{STdata}. Used to limit computations when doing
##'   cross-validation.
##'   \code{only.obs=TRUE} \emph{implies} \code{pred.covar=FALSE} and
##'   \code{combine.data=FALSE}.
##'   Further \code{\link{createSTmodel}} will be called on any \code{STdata}
##'   input, possibly \emph{reordering the observations.}
##' @param pred.var,pred.covar Compute point-wise prediction variances; or
##'   compute covariance matrices for the predicted time series at each location.
##'   \code{pred.covar=TRUE} \emph{implies} \code{pred.var=TRUE} and sets
##'   \code{Nmax} equal to the number of timepoints.
##' @param beta.covar Compute the full covariance matrix for the latent
##'  beta-fields, otherwise only the diagonal elements of V(beta|obs) are
##'  computed. 
##' @param combine.data Combine \code{object} and \code{STdata} and predict for
##'  the joint set of points, see \code{\link{c.STmodel}}.
##' @param type A single character indicating the type of prediction to
##'  compute. Valid options are "f", "p", and "r", for \emph{full},
##'  \emph{profile} or \emph{restricted maximum likelihood} (REML). For profile
##'  and full the predictions are computed assuming that \emph{both} covariance
##'  parameters and regression parameters are known,
##'  e.g. \code{E(X|Y,cov_par,reg_par)}; for REML predictions are compute
##'  assuming \emph{only} covariance parameters known,
##'  e.g. \code{E(X|Y,cov_par)}. The main difference is that REML will have
##'  \emph{larger} variances due to the additional uncertainty in the
##'  regression parameters.
##' @param transform Regard field as log-Gaussian and apply exponential
##'  transformation to predictions. For the final expectations two options
##'  exist, either a unbiased prediction or the (biased) mean-squared error
##'  predictions.
##' @param LTA Compute long-term temporal averages. Either a logical value or a
##'  list; if \code{TRUE} then averages at each location (and variances if
##'  \code{pred.var=TRUE}) are computed; otherwise this should be a list with
##'  elements named after locations and each element containing a vector (or
##'  list of vectors) with dates over which to compute averages. If
##'  \code{only.obs=TRUE} averages are computed over only the observations.
##' @param ... Ignored additional arguments.
##' 
##' @return The function returns a list containing (objects not computed
##'   will be missing):
##'   \item{opts}{Copy of options used in the function call.}
##'   \item{pars}{A list with regression parameters and related variances.
##'               \code{pars} contain \code{gamma.E} and \code{alpha.E} with
##'               regression coefficients for the spatio-temporal model and
##'               land-use covaraiates; variances are found in \code{gamma.V}
##'               and \code{alpha.V}; cross-covariance between gamma and alpha in
##'               \code{gamma.alpha.C}.}
##'   \item{beta}{A list with estimates of the beta-fields, including the
##'               regression mean \code{mu}, conditional expectations \code{EX},
##'               possibly variances \code{VX}, and the full covariance matrix
##'               \code{VX.full}.} 
##'   \item{EX.mu}{predictions based on the regression parameters, geographic
##'                covariates, and temporal trends. I.e. only the deterministic
##'                part of the spatio-temporal model.}
##'   \item{EX.mu.beta}{Predictions based on the latent-beta fields, but excluding
##'                     the residual nu field.}
##'   \item{EX}{Full predictions at the space-time locations in
##'             \code{object} and/or \code{STdata}.}
##'   \item{EX.pred}{Only for \code{transform!="none"}, full predictions
##'                  including bias correction for prediction error.}
##'   \item{VX,VX.pred}{Pointwise variances and prediction variances (i.e. incl.
##'                     contribution from \code{nugget.unobs}) for all locations in \code{EX}.}
##'   \item{VX.full}{A list with (number of locations) elements, each element is a
##'                  (number of timepoints) - by - (number of timepoints) temporal
##'                  covariance matrix for the timeseries at each location.}
##'   \item{MSPE,MSPE.pred}{Pointwise mean-square prediction errors for the
##'                         log-Gaussian fields.}
##'   \item{log.EX,log.VX.pred,log.VX}{Pointwise predictions and variances for
##'           the un-transformed fields when \code{transform!="none"}}
##'   \item{LTA}{A data.frame with temporal averages for locations specified by
##'              \code{LTA}. }
##'   \item{I}{A vector with the locations of the observations in \code{object} or
##'            \code{STdata}. To extract predictions at the observations locations use
##'            \code{EX[I]}.}
##' 
##' @example Rd_examples/Ex_predict_STmodel.R
##' 
##' @author Johan Lindström
##' 
##' @family STmodel methods
##' @family predictSTmodel methods
##' @importFrom stats predict
##' @method predict STmodel
##' @export
predict.STmodel <- function(object, x, STdata=NULL, Nmax=1000, only.pars=FALSE,
                            nugget.unobs=0, only.obs=FALSE, pred.var=TRUE,
                            pred.covar=FALSE, beta.covar=FALSE,
                            combine.data=FALSE, type="p", LTA=FALSE, 
                            transform=c("none","unbiased","mspe"), ...){
##################################
### INITIAL SETUP AND CHECKING ###
  ##check class belongings
  stCheckClass(object, "STmodel", name="object")
  if( !is.null(STdata) ){
    stCheckClass(STdata, c("STdata","STmodel"), name="STdata")
  }
  ##first ensure that type is lower case
  type <- tolower(type)
  ##check if type is valid
  stCheckType(type)
  ##args for transform
  transform <- match.arg(transform)
  if( transform=="none" ){ transform <- NULL }

  if( inherits(x,"estimateSTmodel") ){
    x <- coef(x,"all")$par
  }##if( inherits(x,"estimateSTmodel") )
  
  ##figure out a bunch of dimensions
  ##dimensions$L!=0 is used to check for spatio-temporal covariates
  dimensions <- loglikeSTdim(object)
  
  ##check length of input parameters
  if( type=="f" && length(x)!=dimensions$nparam ){
    stop( paste("type=f, requires", dimensions$nparam,
                "parameters but length(x) =", length(x)) )
  }
  if( type!="f" && length(x)==dimensions$nparam ){
    ##drop the regression parameters
    x <- x[(dimensions$nparam-dimensions$nparam.cov+1):dimensions$nparam]
  }
  if( type!="f" && length(x)!=dimensions$nparam.cov ){
    stop( paste("type!=f, requires", dimensions$nparam.cov,
                "parameters but length(x) =", length(x)) )
  }
  ##only.pars is strange when type=="f"
  if( only.pars && type=="f"){
    warning("only.pars=TRUE and type=(f)ull only returns KNOWN parameters.",
            immediate. = TRUE)
  }
  ##only.pars, and we can ignore a bunch of things
  if( only.pars ){
    only.obs <- pred.covar <- beta.covar <- combine.data <- LTA <-FALSE
    transform <- NULL
  }else{
    ##check if only.obs is valid
    if( only.obs && is.null(STdata) ){
      stop("only.obs=TRUE requires STdata.")
    }
    ##do we have an object to combine with
    if( combine.data && is.null(STdata) ){
      warning("No data to combine with; predicting for 'object'", immediate. = TRUE)
      combine.data <- FALSE
    }
    ##can't combine data if we're predicting at only obs.
    if(only.obs && combine.data){
      warning("only.obs=TRUE implies combine.data=FALSE.", immediate. = TRUE)
      combine.data <- FALSE
    }
    ##computing prediction covariates require !only.obs and pred.var
    if(pred.covar && !pred.var){
      warning("pred.covar=TRUE implies pred.var=TRUE.", immediate. = TRUE)
      pred.var <- TRUE
    }
    if(pred.covar && only.obs){
      warning("only.obs=TRUE implies pred.covar=FALSE.", immediate. = TRUE)
      pred.covar <- FALSE
    }
    ##computing beta covariances requires pred.var
    if(beta.covar && !pred.var){
      warning("beta.covar=TRUE implies pred.var=TRUE.", immediate. = TRUE)
      pred.var <- TRUE
    }
    if( !is.null(transform) ){
      if( type!="r" && transform!="unbiased" ){
        warning("transform 'unbiased' and 'mspe' are equivalent when type!='r'",
                immediate. = TRUE)
        transform <- "unbiased"
      }
      if(pred.covar){
        warning("transform implies pred.covar=FALSE.", immediate. = TRUE)
        pred.covar <- FALSE
      }
    }##if( !is.null(transform) )
  }##if( only.pars ){...}else{...}
  
  ##create the STdata used for predictions
  if( is.null(STdata) ){
    ##pure copy, predict in the data set
    STdata <- object
  }else if(combine.data){
    ##combine the two datasets, and use for predictions (expands object with
    ##ONLY locations from STdata)
    STdata <- c(object, STdata)
  }else{
    if( is.null(object$trend.fnc) && dim(object$trend)[2]!=1 ){
      ##trend.fnc missing and trend is more than a constant
      ##(backwards compability)
      object$trend.fnc <- internalCreateTrendFnc(object$trend)
    }
    
    if( !inherits(STdata,"STmodel") ){
      ##predict only at STdata, and STdata not of class STmodel: need to cast.
      ##First drop ST covariate if no covariate in object
      if(dimensions$L==0){
        STdata$SpatioTemporal <- NULL
      }
      ##and fix the trend (right no of trends, names and dates)
      suppressMessages(STdata <- updateTrend(STdata, fnc=object$trend.fnc,
                                             extra.dates=STdata$trend$date))
      
      ##since we're not using covariances for the prediction locations
      ##(specified seperately in nugget.unobs), we just pick a simple covariance
      ##structure, avoiding problems with missing covariates/levels.
      cov.nu <- object$cov.nu
      cov.nu$nugget <- TRUE
      ##Create an STmodel from STdata
      STdata <- createSTmodel(STdata, LUR=object$LUR.list, ST=object$ST.list,
                              cov.beta=object$cov.beta, cov.nu=cov.nu,
                              locations=object$locations.list,
                              scale=!is.null(object$scale.covars),
                              scale.covars=object$scale.covars)
    }else{
      ##STdata is an STmodel object ->
      ##test for consistent covariates and scaling (only equal scaling allowed).
      areSTmodelsConsistent(object, STdata, "STdata")
    }##if( !inherits(STdata,"STmodel") ){...}{else{...}

    ##STdata is now of type STmodel, lets fix the trend
    suppressMessages(STdata <- updateTrend.STmodel(STdata, fnc=object$trend.fnc))
  }##if( is.null(STdata) ){...}else if(...){...}else{...}

  ##Should we compute LTAs?
  if( !is.list(LTA) && LTA ){
    ##yes, create list
    if( only.obs ){
      ##averages over observtions only,
      ##i.e. split observation dates by their locations.
      LTA.list <- split(STdata$obs$date, STdata$obs$ID)
    }else{
      ##averages over all predictions
      LTA.list <- lapply(STdata$locations$ID,
                         function(x){ return(STdata$trend$date) })
      names(LTA.list) <- STdata$locations$ID
    }
  }else if( is.list(LTA) ){
    ##yes, list given
    LTA.list <- LTA
    LTA <- TRUE
  }else{
    ##no, empty list and set LTA=FALSE
    LTA.list <- NULL
    LTA <- FALSE
  }##if( !is.list(LTA) && LTA ){...}else if(...){...}else{...}
  
  ##check the validity of LTA.list
  if( LTA ){
    ##check validity of sites
    LTA.missing <- !(names(LTA.list) %in% STdata$locations$ID)
    if( any(LTA.missing) ){
      warning("Removed ", sum(LTA.missing),
              " elements from LTA not found in STdata$locations$ID",
              immediate. = TRUE)
      LTA.list <- LTA.list[ !LTA.missing ]
    }

    if( length(LTA.list)!=0 ){
      ##convert components to lists
      LTA.list <- lapply(LTA.list,
                         function(x){ switch(is.list(x)+1, list(x), x) })
      ##check validity of dates
      LTA.tot <- 0
      for(i in 1:length(LTA.list)){
        for(j in length(LTA.list[[i]])){
          LTA.missing <- !(LTA.list[[i]][[j]] %in% STdata$trend$date)
          LTA.list[[i]][[j]] <- LTA.list[[i]][[j]][ !LTA.missing ]
          LTA.tot <- LTA.tot + sum(LTA.missing)
        }
        ##drop empty points
        LTA.list[[i]] <- LTA.list[[i]][ sapply(LTA.list[[i]], length)!=0 ]
      }##for(i in 1:length(LTA.list))
      if( LTA.tot!=0 ){
        warning("Removed ", LTA.tot,
                " dates not found in STdata$trend$date from all LTA elements",
                immediate. = TRUE)
      }
      ##drop empty points
      LTA.list <- LTA.list[ sapply(LTA.list, length)!=0 ]
    }##if( length(LTA.list)!=0 )

    ##nothing left -> don't compute LTA's
    if( length(LTA.list)==0 ){
      LTA <- FALSE
      LTA.list <- NULL
    }
  }##if( LTA )

#########################################
### COMPUTE COVARIANCE MATRICES FOR Y ###
  ##extract parameters from x
  tmp <- loglikeSTgetPars(x, object)
  if(type=="f"){
    gamma <- tmp$gamma
    alpha <- tmp$alpha
  }
  cov.pars.beta <- tmp$cov.beta
  cov.pars.nu <- tmp$cov.nu

  ##the following only needed for type!="f", or predictions
  if(type!="f" || !only.pars){
    ##nugget for unobserved sites
    nugget.unobs <- internalFixNuggetUnobs(nugget.unobs, STdata,
                                           cov.pars.nu$nugget)
    ##extract sparse matrices
    Fobs <- expandF(object$F, object$obs$idx, n.loc=dimensions$n.obs)
  
    ##Create the Xtilde = [M FX] matrix
    Xtilde <- as.matrix( Fobs %*% Matrix::bdiag(object$LUR) )
    ##Add the spatio-temporal covariate (if it exists)
    if( dimensions$L!=0 ){
      Xtilde <- cbind(object$ST, Xtilde)
    }
    ##create covariance matrices, beta-field
    i.sigma.B <- makeSigmaB(cov.pars.beta$pars, dist = object$D.beta,
                            type = object$cov.beta$covf,
                            nugget = cov.pars.beta$nugget)
    ##and nu-field
    i.sigma.nu <- makeSigmaNu(cov.pars.nu$pars, dist = object$D.nu,
                              type = object$cov.nu$covf,
                              nugget = cov.pars.nu$nugget,
                              random.effect = cov.pars.nu$random.effect,
                              blocks1 = object$nt, ind1 = object$obs$idx,
                              sparse=TRUE)
    ##calculate block cholesky factor of the matrices, in-place to conserve memory
    i.sigma.B <- makeCholBlock(i.sigma.B, n.blocks=dimensions$m)
    i.sigma.nu <- chol(i.sigma.nu)
    ##invert the matrices, in-place to conserve memory
    i.sigma.B <- invCholBlock(i.sigma.B, n.blocks=dimensions$m)
    i.sigma.nu <- chol2inv(i.sigma.nu)

    ##F'*inv(sigma.nu)
    tF.iS <- t(Fobs) %*% i.sigma.nu
    ##F'*inv(sigma.nu)*F
    tF.iS.F <- as.matrix(tF.iS %*% Fobs)
    ##inv(sigma.B|Y) = F'*inv(sigma.nu)*F + inv(sigma.B)
    R.i.sigma.B.Y <- tF.iS.F + i.sigma.B
    ##calculate cholesky factor of inv(sigma.B|Y)
    R.i.sigma.B.Y <- chol(R.i.sigma.B.Y)
    
    ##compute iSigma.nu*Xtilde and iSigma.nu*Xtilde
    iS.X <- i.sigma.nu %*% Xtilde
    iS.Y <- i.sigma.nu %*% object$obs$obs
  
    ##compute F'*iSigma.nu*Xtilde and F'*iSigma.nu*Y
    tF.iS.X <- as.matrix(t(Fobs) %*% iS.X)
    tF.iS.Y <- as.matrix(t(Fobs) %*% iS.Y)
    
    ##compute inv(sigma.B|Y)*F'*iSigma.nu*Xtilde
    iSBY.tF.iS.X <- solveTriBlock(R.i.sigma.B.Y, tF.iS.X, transpose=TRUE)
    iSBY.tF.iS.X <- solveTriBlock(R.i.sigma.B.Y, iSBY.tF.iS.X, transpose=FALSE)
  }#if(type!="f" || !only.pars)
  
###############################################
### REGRESISON PARAMETERS (alpha and gamma) ###
  if(type=="f"){
    ##regression parameters given
    gamma.E <- c(tmp$gamma)
    alpha.E <- unlist(tmp$alpha)
    ##known parameters, set variances to zero
    gamma.V <- matrix(0,length(gamma.E),length(gamma.E))
    alpha.V <- matrix(0,length(alpha.E),length(alpha.E))
    gamma.alpha.C <- matrix(0,length(gamma.E),length(alpha.E))
    ##create a combined vector
    gamma.alpha <- as.matrix(c(gamma.E,alpha.E))
  }else{
    ##compute inv(Xtilde'*Sigma^-1*Xtilde) (and force symmetric)
    i.XSX <- as.matrix( t(Xtilde) %*% iS.X )
    i.XSX <- symmpart(i.XSX - t(tF.iS.X) %*% iSBY.tF.iS.X)
    
    ##compute t(Xtilde'*Sigma^-1*Y)
    tmp2 <- as.matrix( object$obs$obs %*% iS.X )
    tmp2 <- tmp2 - t(tF.iS.Y) %*% iSBY.tF.iS.X
    
    ##compute alpha and gamma
    gamma.alpha <- solve(i.XSX, t(tmp2))

    ##extract computed parameters
    if(dimensions$L!=0){
      gamma.E <- c(gamma.alpha[1:dimensions$L])
      alpha.E <- c(gamma.alpha[-(1:dimensions$L)])
    }else{
      gamma.E <- double(0)
      alpha.E <- c(gamma.alpha)
    }
    ##compute the variances
    i.XSX <- solve(i.XSX)
    if(dimensions$L!=0){
      gamma.V <- i.XSX[1:dimensions$L,1:dimensions$L,drop=FALSE]
      alpha.V <- i.XSX[-c(1:dimensions$L),-c(1:dimensions$L),drop=FALSE]
      gamma.alpha.C <- i.XSX[1:dimensions$L,-c(1:dimensions$L),drop=FALSE]
    }else{
      gamma.V <- matrix(0,0,0)
      alpha.V <- i.XSX
      gamma.alpha.C <- matrix(0,length(gamma.E),length(alpha.E))
    }
  }##if(type=="f"){...}else{...}
  
  ##force to matrices
  gamma.E <- as.matrix(gamma.E)
  alpha.E <- as.matrix(alpha.E)
  ##add names
  names.tmp <- loglikeSTnames(object, TRUE)
  names.tmp <- names.tmp[1:(dimensions$nparam-dimensions$nparam.cov)]
  if(dimensions$L!=0){
    ##first for gamma
    rownames(gamma.E) <- names.tmp[1:dimensions$L]
    colnames(gamma.V) <- rownames(gamma.V) <- rownames(gamma.E)
    rownames(gamma.alpha.C) <- rownames(gamma.E)
    rownames(alpha.E) <- names.tmp[-(1:dimensions$L)]
  }else{
    rownames(alpha.E) <- names.tmp
  }
  ##and then for alpha
  colnames(alpha.V) <- rownames(alpha.V) <- rownames(alpha.E)
  colnames(gamma.alpha.C) <- rownames(alpha.E)
  
###############################
### CONSTRUCT RETURN OBJECT ###
  out <- list()
  ##set class
  class(out) <- "predictSTmodel"
  ##save options used for predict.
  out$opts <- list(only.pars=only.pars, nugget.unobs=nugget.unobs,
                   only.obs=only.obs, pred.var=pred.var, pred.covar=pred.covar,
                   beta.covar=beta.covar, combine.data=combine.data, type=type,
                   transform=transform, LTA=LTA, LTA.list=LTA.list)
  ##collect the regression results
  out$pars <- list(gamma.E=gamma.E, alpha.E=alpha.E, 
                   gamma.V=gamma.V, alpha.V=alpha.V,
                   gamma.alpha.C=gamma.alpha.C)

  if( out$opts$only.pars ){
    ##return the parameters
    return( out )
  }
  ##remove variables not needed, regression parameters...
  rm(gamma.E, alpha.E, gamma.V, alpha.V, gamma.alpha.C, names.tmp)
  ##and options stored elsewhere
  rm(only.pars, nugget.unobs, only.obs, pred.var, pred.covar,
     beta.covar, combine.data, type, LTA, LTA.list, transform)
    
####################################
### COMPUTE inv(Sigma_oo)*(Y-mu) ###
  
  ##compute F'*iSigma.nu*(C-mu) and iSigma.nu*(C-mu)
  tF.iS.C <- tF.iS.Y - tF.iS.X %*% gamma.alpha
  iS.C <- iS.Y - iS.X %*% gamma.alpha
  
  ##compute inv(sigma.B|Y)*F'*iSigma.nu*(C-mu)
  iSBY.tF.iS.C <- solveTriBlock(R.i.sigma.B.Y, tF.iS.C, transpose=TRUE)
  iSBY.tF.iS.C <- solveTriBlock(R.i.sigma.B.Y, iSBY.tF.iS.C, transpose=FALSE)
  
  ##compute iSigma.tilde*(C-mu)
  ##this is the same for all the conditionals...
  iSoo.C <- as.matrix( iS.C - t( tF.iS ) %*% iSBY.tF.iS.C )

  ##for REML we also need iSigma.tilde*Xtilde
  if( out$opts$type=="r" ){
    iSoo.Xtilde <- as.matrix( iS.X - t(tF.iS) %*% iSBY.tF.iS.X )
  }

  ##remove variables not needed
  rm(iS.X, iS.Y, tF.iS.X, tF.iS.Y, iSBY.tF.iS.X, tF.iS.C, iS.C, iSBY.tF.iS.C)

###################################
### DEFINE PREDICTION LOCAITONS ###
  ##create matrices with site index and temporal trends for all
  ##OBS date has to be INCREASING, i.e. same sorting as for STdata$obs...
  if( out$opts$only.obs ){
    ##unobserved points
    idx.unobs <- STdata$obs$idx
    ##time of observation for each site
    T1 <- STdata$obs$date
    ##trend functions for all points matrix
    Funobs <- STdata$F
    if( dimensions$L!=0 ){
      ##extract the spatio-temporal trends
      ST.unobs <- STdata$ST
    }
    ##size of the unobserved locations
    N.unobs <- length( unique(idx.unobs) )
    T.unobs <- length( unique(T1) )
  }else{
    ##size of the unobserved locations
    N.unobs <- dim(STdata$locations)[1]
    T.unobs <- dim(STdata$trend)[1]
  
    ##unobserved points
    idx.unobs <- rep(1:N.unobs, each=T.unobs)
    ##time of observation for each site
    T1 <- rep(STdata$trend$date, N.unobs)
    ##trend functions for all points matrix
    Funobs <- matrix(1, (T.unobs*N.unobs), dimensions$m)
    for(i in (1:dimensions$m)){
      ##if not constant, find relevant row in $trend
      if(colnames(STdata$F)[i] != "const"){
        Funobs[,i] <- rep(STdata$trend[,colnames(STdata$F)[i]], N.unobs)
      }
    }##for(i in (1:dimensions$m))
    if( dimensions$L!=0 ){
      ##extract the spatio-temporal trends
      ST.unobs <- matrix(STdata$ST.all, (T.unobs*N.unobs), dimensions$L)
    }##if( dimensions$L!=0 )
  }##if( out$opts$only.obs ) ... else ...
  
  ##compute number of observations for each time-point,
  ##taking into account that there might be missing dates.
  date.all <- sort(unique( c(object$trend$date, STdata$trend$date) ))
  nt.unobs <- nt.obs <- double(length(date.all))
  for(i in c(1:length(date.all))){
    nt.obs[i] <- sum( object$obs$date==date.all[i] )
  }

  ##extract sparse matrices
  Funobs <- expandF(Funobs, idx.unobs, n.loc=N.unobs)
  ##compute LUR times the temporal trends [M F*X] for unobserved
  Xtilde.unobs <- as.matrix( Funobs %*% Matrix::bdiag(STdata$LUR.all) )
  if( dimensions$L!=0 ){
    Xtilde.unobs <- cbind(ST.unobs, Xtilde.unobs)
  }
  
###########################
### COMPUTE BETA FIELDS ###
  ##compute LUR*alpha (mean values for the beta fields)
  out$beta <- list()
  out$beta$mu <- matrix(Matrix::bdiag(STdata$LUR.all) %*% out$pars$alpha.E,
                        ncol=length(STdata$LUR.all))
  colnames(out$beta$mu) <- names(STdata$LUR.all)
  rownames(out$beta$mu) <- rownames( STdata$LUR.all[[1]] )
  
  ##create distance matrix between unobserved and observed elements
  ##unobserved locations
  loc.unobs.nu <- STdata$locations[, c("x.nu","y.nu"), drop=FALSE]
  loc.unobs.beta <- STdata$locations[, c("x.beta","y.beta"), drop=FALSE]

  ##observed locations
  I.obs <- match( colnames(object$D.nu), object$locations$ID)
  loc.obs.nu <- object$locations[I.obs, c("x.nu","y.nu"), drop=FALSE]
  loc.obs.beta <- object$locations[I.obs, c("x.beta","y.beta"), drop=FALSE]
  
  ##We also have that in the stripped data the indecies may not match so we need a
  ##second vector giving which idx in striped matches original idx.
  Ind.2.1 <- match(object$locations$ID[I.obs], STdata$locations$ID, nomatch=0)
    
  ##precompute the cross-covariance for the beta-fields
  sigma.B.C <- makeSigmaB(cov.pars.beta$pars,
                          dist = crossDist(loc.unobs.beta, loc.obs.beta),
                          type = object$cov.beta$covf,
                          nugget = cov.pars.beta$nugget,
                          ind2.to.1=Ind.2.1, sparse=TRUE)

  ##The right most part is F'*iSigma.tilde*(C-mu)
  out$beta$EX <- out$beta$mu + matrix(sigma.B.C %*% (t(Fobs) %*% iSoo.C),
                                      ncol=dim(out$beta$mu)[2])
  
  ##and variances
  if( out$opts$pred.var ){
    Sby.iSb.Sou <- as.matrix( i.sigma.B %*% t(sigma.B.C) )
    Sby.iSb.Sou <- solveTriBlock(R.i.sigma.B.Y, Sby.iSb.Sou, transpose=TRUE)
    Sby.iSb.Sou <- solveTriBlock(R.i.sigma.B.Y, Sby.iSb.Sou, transpose=FALSE)

    ##contribution from REML estimate
    if( out$opts$type=="r" ){
      var.beta.REML <- (cBind(rep(0,dimensions$L), Matrix::bdiag(STdata$LUR.all)) -
                        sigma.B.C %*% (t(Fobs)%*%iSoo.Xtilde))
      var.beta.REML <- as.matrix(var.beta.REML)
      if( out$opts$beta.covar ){
        ##full matrix
        V.REML <- (var.beta.REML %*% i.XSX) %*% t(var.beta.REML)
      }else{
        ##only diagonal elements
        V.REML <- rowSums( (var.beta.REML %*% i.XSX) * var.beta.REML )
      }
    }else{
      V.REML <- 0
    }##if( out$opts$type=="r" ){...}else{...}
    
    if( out$opts$beta.covar ){
      ##sigma.B for unobserved locations
      sigma.B.uu <- makeSigmaB(cov.pars.beta$pars,
                               dist = crossDist(loc.unobs.beta),
                               type = object$cov.beta$covf,
                               nugget = cov.pars.beta$nugget)
      ##compute variance matrix
      tmp <- sigma.B.uu - sigma.B.C %*% tF.iS.F %*% Sby.iSb.Sou + V.REML
      
      out$beta$VX.full <- list()
      for(i in 1:dim(out$beta$EX)[2]){
        Ind <- (1:dim(out$beta$EX)[1]) + (i-1)*dim(out$beta$EX)[1]
        out$beta$VX.full[[i]] <- as.matrix(tmp[Ind, Ind, drop=FALSE])
        rownames(out$beta$VX.full[[i]]) <- rownames(out$beta$EX)
        colnames(out$beta$VX.full[[i]]) <- rownames(out$beta$EX)
      }
      names(out$beta$VX.full) <- colnames(out$beta$EX)
      tmp <- diag(tmp)
    }else{
      ##compute only diagonal elements of the STATIONARY sigma.beta matrix
      sigma.B.uu <- makeSigmaB(cov.pars.beta$pars, dist = matrix(0,1,1),
                               type = object$cov.beta$covf,
                               nugget = cov.pars.beta$nugget)
      ##expand to full matrix
      sigma.B.uu <- matrix(diag(sigma.B.uu), ncol=dim(sigma.B.uu)[1],
                           nrow=dim(loc.unobs.beta)[1], byrow=TRUE)
      ##and compute the relevant covariance
      tmp <- (c(sigma.B.uu) - rowSums(sigma.B.C * t(tF.iS.F %*% Sby.iSb.Sou)) +
              V.REML)
    }##if( out$opts$beta.covar ){...}else{...}
    
    out$beta$VX <- matrix(tmp, ncol=dim(out$beta$EX)[2])
    dimnames(out$beta$VX) <- dimnames(out$beta$EX)
    
    ##remove variables not needed
    rm(tmp, Sby.iSb.Sou, sigma.B.uu, V.REML)
  }##  if( out$opts$pred.var )
  ##remove variables not needed
  rm(i.sigma.B, tF.iS.F)

######################################
### MATRICES HOLDING RETURN VALUES ###
  ##mean value of the field
  out$EX.mu <- as.matrix( Xtilde.unobs %*% gamma.alpha )
  ##compute predicitons based on the beta fields
  out$EX.mu.beta <- as.matrix( Funobs %*% c(out$beta$EX) )
  ##Add the spatio-temporal covariates, if any
  if( dimensions$L!=0 ){
    out$EX.mu.beta <- out$EX.mu.beta + (ST.unobs %*% out$pars$gamma.E)
  }
  
  if( !out$opts$only.obs ){
    ##reshape out$EX.mu and out$EX.mu.beta to match T-by-N
    dim(out$EX.mu.beta) <- dim(out$EX.mu) <- c(T.unobs, N.unobs)
    ##and add names
    colnames(out$EX.mu) <- STdata$locations$ID
    rownames(out$EX.mu) <- as.character(STdata$trend$date)
    dimnames(out$EX.mu.beta) <- dimnames(out$EX.mu)
  }else{
    dimnames(out$EX.mu.beta) <- dimnames(out$EX.mu) <- NULL
  }
  ##lets save the contribution from the mean.
  out$EX <- out$EX.mu
  ##additional fields for log-transformation
  if( !is.null(out$opts$transform) ){
    ##temporary variables for transformed data
    EX.trans <- out$EX.mu
    EX.trans.pred <- out$EX.mu
    ##reserve places in the out structure
    out$log.EX <- out$EX.pred <- NA
  }
  
  ##Create matrices containing pointwise variances
  ##needed if a) variances requested OR b) log-transform
  if( out$opts$pred.var || !is.null(out$opts$transform) ){
    out$VX <- matrix(NA, dim(out$EX)[1], dim(out$EX)[2])
    out$VX.pred <- matrix(NA, dim(out$EX)[1], dim(out$EX)[2])
    ##add names
    if( !out$opts$only.obs ){
      dimnames(out$VX) <- dimnames(out$VX.pred) <- dimnames(out$EX)
    }
    ##and a list of covariances.
    if( out$opts$pred.covar ){
      out$VX.full <- list()
    }
    ##additional fields for log-transformation
    ##(only if variances requested AND log-transform)
    if( out$opts$pred.var && !is.null(out$opts$transform) ){
      out$MSPE <- out$MSPE.pred <- out$VX
    }
  }##if( out$opts$pred.var )

  ##list holding LTA elements
  if( out$opts$LTA ){
    out$LTA <- vector("list", length(out$opts$LTA.list))
    names(out$LTA) <- names(out$opts$LTA.list)
  }
  
################################
### COMPUTE PREDICITONS OF X ###
  ##cross distance for the nu coordinates
  cross.D.nu <- crossDist(loc.unobs.nu, loc.obs.nu)

  ## sigma.B.C * F'
  sigma.B.C.tF <- sigma.B.C %*% t(Fobs)
  
  ##now we need to split the conditional expectation into parts to
  ##reduce the memory footprint
  if( out$opts$pred.covar || out$opts$LTA ){
    Ind.list <- split(1:length(out$EX), idx.unobs)
  }else{
    Ind.list <- split(1:length(out$EX),
                      rep(1:ceiling(length(out$EX)/Nmax),
                          each=Nmax)[1:length(out$EX)])
  }
  
  ##loop over the parts
  for( i in 1:length(Ind.list) ){
    ##index of the points we want to get conditional expectations for
    Ind <- Ind.list[[i]]
    ##compute number of unobserved per block for this subset
    T1.Ind <- T1[Ind]
    for(j in c(1:length(date.all))){
      nt.unobs[j] <- sum( T1.Ind==date.all[j] )
    }
    ##order T1.Ind for computation of cross.covariance
    T1.order <- order(T1.Ind)
    
    ##create full cross-covariance matrix for all the observations
    sigma.B.full.C <- as.matrix(Funobs[Ind,,drop=FALSE] %*% sigma.B.C.tF)
    ##parameter order is sill, nugget, range (always predict with nugget=0)
    sigma.nu.C <- makeSigmaNu(cov.pars.nu$pars, dist = cross.D.nu,
                              type = object$cov.nu$covf, nugget = 0,
                              random.effect = cov.pars.nu$random.effect,
                              ind1 = (idx.unobs[Ind])[T1.order],
                              ind2 = object$obs$idx,
                              blocks1 = nt.unobs, blocks2 = nt.obs,
                              ind2.to.1=Ind.2.1)
    sigma.nu.C[T1.order,] <- sigma.nu.C
    sigma.nu.C <- sigma.nu.C + sigma.B.full.C
    ##calculate conditional expectation
    out$EX[Ind] <- out$EX.mu[Ind] + sigma.nu.C %*% iSoo.C

    if( out$opts$LTA ){
      ##which site are we at and does it have a matching LTA request?
      ##asked here to avoid unnescesary VX.full computations.
      ID.tmp <- STdata$locations$ID[unique(idx.unobs[Ind])]
      ##sanity check
      if( length(ID.tmp)!=1 ){
        stop("Something wrong with LTA-computations.")
      }
      LTA.tmp <- out$opts$LTA.list[[ID.tmp]]
    }else{
      LTA.tmp <- NULL
    }
    
    ##calculate pointwise variance
    if( out$opts$pred.var || !is.null(out$opts$transform) ){
      ##pick out the observed locations in this iteration
      I.loc <- sort(unique( idx.unobs[Ind] ))
      I.loc2 <- (rep(I.loc,dimensions$m) +
                 rep((0:(dimensions$m-1)) * N.unobs, each=length(I.loc)))
      ##compute relevant part of the covariance matrix
      unobs.D.beta <- crossDist(loc.unobs.beta[I.loc,,drop=FALSE])
      sigma.B.uu <- makeSigmaB(cov.pars.beta$pars, dist = unobs.D.beta,
                               type = object$cov.beta$covf,
                               nugget = cov.pars.beta$nugget, sparse=TRUE)
      ##first the unobserved covariance matrix
      V.uu <- (Funobs[Ind,I.loc2,drop=FALSE] %*%
               Matrix(sigma.B.uu %*% t(Funobs[Ind,I.loc2,drop=FALSE])))
      ##distance matrices for the unobserved nu locations
      unobs.D.nu <- crossDist(loc.unobs.nu[I.loc,,drop=FALSE])
      V.uu <- V.uu + makeSigmaNu(cov.pars.nu$pars, dist = unobs.D.nu,
                                 type = object$cov.nu$covf, nugget = 0,
                                 random.effect = cov.pars.nu$random.effect,
                                 ind1 = idx.unobs[Ind]-min(I.loc)+1,
                                 blocks1 = nt.unobs)
      
      ##compute iSigma.nu*Sigma.ou
      iS.Sou <- i.sigma.nu %*% t(sigma.nu.C)
      ##compute F'*iSigma.nu*Sigma.ou
      tF.iS.Sou <- as.matrix(tF.iS %*% t(sigma.nu.C))
      ##compute chol(inv(sigma.B.Y))' * (F'*iSigma.nu*Sigma.ou)
      Sby.tF.iS.Sou <- solveTriBlock(R.i.sigma.B.Y, tF.iS.Sou, transpose=TRUE)

      ##contribution from REML estimate
      if( out$opts$type=="r" ){
        tmp <- Xtilde.unobs[Ind,,drop=FALSE] - sigma.nu.C %*% iSoo.Xtilde
        if( out$opts$pred.covar || !is.null(LTA.tmp) ){
          ##full matrix
          V.REML <- (tmp %*% i.XSX) %*% t(tmp)
          if( !is.null(out$opts$transform) ){
            ##also compute LAMBDA correction for transform
            lambda.full <- (tmp %*% i.XSX) %*% t(Xtilde.unobs[Ind,,drop=FALSE])
            lambda <- diag(lambda.full)
          }
        }else{
          ##only diagonal elements
          V.REML <- rowSums( (tmp %*% i.XSX) * tmp)
          if( !is.null(out$opts$transform) ){
            ##also compute LAMBDA correction for transform
            lambda <- rowSums((tmp %*% i.XSX) * (Xtilde.unobs[Ind,,drop=FALSE]))
          }
        }
      }else{
        V.REML <- 0
        lambda <- 0
      }

      ##combine parts, either to covariance matrix or VX for each estimate.
      if( out$opts$pred.covar || !is.null(LTA.tmp) ){
        ##full matrix
        tmp <- - sigma.nu.C %*% iS.Sou + t(Sby.tF.iS.Sou) %*% Sby.tF.iS.Sou
        V.cond <- V.cond.0 <- as.matrix(V.uu + tmp + V.REML)
        ##add nugget to the diagonal
        diag(V.cond) <- diag(V.cond) + out$opts$nugget.unobs[idx.unobs[Ind]]
        
        ##extract diagonal elements and store in output structure
        out$VX[Ind] <- diag(V.cond.0)
        out$VX.pred[Ind] <- diag(V.cond)
        if(out$opts$pred.covar){
          out$VX.full[[i]] <- V.cond.0
        }
        
      }else{
        ##only diagonal elements
        tmp <- (- colSums(t(sigma.nu.C) * iS.Sou) +
                colSums(Sby.tF.iS.Sou * Sby.tF.iS.Sou))
        out$VX[Ind] <- diag(V.uu) + tmp + V.REML
        out$VX.pred[Ind] <- out$VX[Ind] + out$opts$nugget.unobs[idx.unobs[Ind]]
      }

      ##ensure positive variances
      out$VX[Ind] <- pmax(out$VX[Ind], 0)
      out$VX.pred[Ind] <- pmax(out$VX.pred[Ind], 0)
    }else{
      V.cond <- V.cond.0 <- NULL
    }##if( out$opts$pred.var || !is.null(out$opts$transform) ){...}else{...}
    
    ##compute log-transform, store in EX.trans, EX.trans.pred
    if( !is.null(out$opts$transform) ){
      ##use -lambda or -2*lambda, lambda=0 if type!="r"
      if( out$opts$transform=="unbiased" ){ c <- 1 }else{ c <- 2 }

      ##expectations
      EX.trans[Ind] <- exp( out$EX[Ind] + out$VX[Ind]/2 - c*lambda )
      EX.trans.pred[Ind] <- exp( out$EX[Ind] + out$VX.pred[Ind]/2 - c*lambda )

      ##compute MSPE
      if( out$opts$pred.var ){
        ##first compute/extract the unbiased estimate for each location
        if( out$opts$transform=="mspe" ){
          EX.ub <- exp( out$EX[Ind] + out$VX[Ind]/2 - lambda )
          EX.ub.pred <- exp( out$EX[Ind] + out$VX.pred[Ind]/2 - lambda )
        }else{
          EX.ub <- EX.trans[Ind]
          EX.ub.pred <- EX.trans.pred[Ind]
        }
        
        ##prediction variance for unobserved locations (add nugget)
        V.uu.pred <- V.uu + diag(out$opts$nugget.unobs[idx.unobs[Ind]])
        
        if( !is.null(LTA.tmp) ){
          ##LTA to be computed, we'll need the full matrix

          ##then find the adjustment factor for the second term, only relevant for
          ##"unbiased" and type="r" (set to 1 o.w.)
          if(out$opts$transform=="unbiased" && out$opts$type=="r"){
            Z <- exp(lambda.full)
            Z <- Z+t(Z) - exp(lambda.full + t(lambda.full))
          }else{
            Z <- 1
          }
          MSPE.part <- exp(V.uu) * (1 - exp(-V.cond.0)*Z)
          MSPE.part.pred <- exp(V.uu.pred) * (1 - exp(-V.cond)*Z)
          ##Compute MSPE for EX and EX.pred
          out$MSPE[Ind] <- (EX.ub*EX.ub * diag(MSPE.part))
          out$MSPE.pred[Ind] <- (EX.ub.pred*EX.ub.pred * diag(MSPE.part.pred))
        }else{
          ##then find the adjustment factor for the second term, only relevant for
          ##"unbiased" and type="r" (set to 1 o.w.)
          if(out$opts$transform=="unbiased" && out$opts$type=="r"){
            Z <- 2*exp(lambda)-exp(2*lambda)
          }else{
            Z <- 1
          }
          ##Compute MSPE for EX and EX.pred
          out$MSPE[Ind] <- (EX.ub*EX.ub * exp( diag(V.uu) ) *
                            ( 1 - exp(-out$VX[Ind])*Z ))
          out$MSPE.pred[Ind] <- (EX.ub.pred*EX.ub.pred * exp( diag(V.uu.pred) ) *
                                 ( 1 - exp(-out$VX.pred[Ind])*Z ))
        }##if( !is.null(LTA.tmp) ){...}else{...}
        
        ##ensure positive MSPE
        out$MSPE[Ind] <- pmax(out$MSPE[Ind], 0)
        out$MSPE.pred[Ind] <- pmax(out$MSPE.pred[Ind], 0)
      }else{
        MSPE.part <- MSPE.part.pred <- NULL
      }##if( out$opts$pred.var ){...}else{...}

      if( !is.null(LTA.tmp) ){
        ##observatiosn for this location
        EX.tmp <- cbind(exp(out$EX.mu[Ind]), exp(out$EX.mu.beta[Ind]),
                        EX.trans[Ind], EX.trans.pred[Ind])
        colnames(EX.tmp) <- c("EX.mu", "EX.mu.beta", "EX", "EX.pred")
        ##compute mean value
        out$LTA[[ID.tmp]] <- internalComputeLTA(LTA.tmp, EX.tmp, T1[Ind],
                                                V=MSPE.part, 
                                                V.pred=MSPE.part.pred,
                                                E.vec=EX.ub, 
                                                E.vec.pred=EX.ub.pred)
      }##if( !is.null(LTA.tmp) )
    }else{
      ##if not transform, seperate code for LTA (differences are too big...)
      if( !is.null(LTA.tmp) ){
        ##observatiosn for this location
        EX.tmp <- cbind(out$EX.mu[Ind], out$EX.mu.beta[Ind], out$EX[Ind])
        colnames(EX.tmp) <- c("EX.mu", "EX.mu.beta", "EX")
        ##compute mean value
        out$LTA[[ID.tmp]] <- internalComputeLTA(LTA.tmp, EX.tmp, T1[Ind],
                                                V=V.cond.0, V.pred=V.cond)
      }##if( !is.null(LTA.tmp) )
    }##if( !is.null(out$opts$transform) ){...} else {...}
    
  }##for( i in 1:length(Ind.list) )

  ##combine LTA results
  if(out$opts$LTA){
    ##bind results
    out$LTA <- do.call(cbind, out$LTA)
    ##and add names
    tmp <- sapply(out$opts$LTA.list, length)
    if( all(tmp==1) ){
      colnames(out$LTA) <- names(tmp)
    }else{
      colnames(out$LTA) <- paste(rep(names(tmp),times=tmp),
                                 unlist(lapply(tmp,seq)), sep=".")
    }
    ##convert to data.frame
    out$LTA <- as.data.frame(t(out$LTA))
  }##if(out$opts$LTA)

  if( !is.null(out$opts$transform) ){
    ##move EX of log-field
    out$log.EX <- out$EX
    ##store EX of transformed field
    out$EX <- EX.trans
    out$EX.pred <- EX.trans.pred
    ##transform partial fields
    out$EX.mu <- exp(out$EX.mu)
    out$EX.mu.beta <- exp(out$EX.mu.beta)
    ##rename VX/VX.pred to log-fields
    I <- grep("^VX",names(out))
    names(out)[I] <- paste("log.",names(out)[I],sep="")
  }

  ##add names to full prediction variances
  if( out$opts$pred.covar ){
    names(out$VX.full) <- colnames(out$EX)
    for(i in 1:length(out$VX.full))
      colnames(out$VX.full[[i]]) <- rownames(out$VX.full[[i]]) <- rownames(out$EX)
  }##if(out$opts$pred.covar)
  
  ##index for the unobserved values into the predicted values
  if( length(STdata$obs$obs)!=0 ){
    if( out$opts$only.obs ){
      I <- 1:length(out$EX)
    }else{
      I <- (match(STdata$obs$ID,colnames(out$EX))-1)*dim(out$EX)[1] +
        match(STdata$obs$date, STdata$trend$date)
    }
    out$I <- data.frame(I=I, date=STdata$obs$date, ID=STdata$obs$ID,
                        stringsAsFactors=FALSE)
  }##if( length(STdata$obs$obs)!=0 )

  return( out )
}##function predict.STmodel


###################################
## S3 methods for predictSTmodel ##
###################################
##' \code{\link[base:print]{print}} method for class \code{predictSTmodel}.
##'
##' @title Print details for \code{predictSTmodel} object
##' @param x \code{predictSTmodel} object to print information for.
##' @param ... Ignored additional arguments.
##' @return Nothing
##'
##' @author Johan Lindström
##'
##' @examples
##'   ##load data
##'   data(pred.mesa.model)
##'   print(pred.mesa.model)
##'   
##' 
##' @family predictSTmodel methods
##' @method print predictSTmodel
##' @export
print.predictSTmodel <- function(x, ...){
  ##check class belonging
  stCheckClass(x, "predictSTmodel", name="x")

  cat("Prediction for STmodel.\n\n")
  if( x$opts$only.pars ){
    cat("Only computed regression parameters:\n")
  }else{
    cat("Regression parameters:\n")
  }
  if( length(x$pars$gamma.E)==0 ){
    cat("\t","No spatio-temporal covariate.\n")
  }else{
    cat("\t", length(x$pars$gamma), "Spatio-temporal covariate(s).\n")
  }
  cat("\t", length(x$pars$alpha.E),
      "beta-fields regression parameters in x$pars.\n\n")
  if( x$opts$only.pars ){
    return(invisible())
  }

  if( x$opts$type=="r" ){
    cat("Regression parameters are assumed to be unknown and\n")
    cat("\tprediction variances include uncertainties\n")
    cat("\tin regression parameters.\n\n")
  }else{
    cat("Regression parameters are assumed to be known and\n")
    cat("\tprediction variances do NOT include\n")
    cat("\tuncertainties in regression parameters.\n\n")
  }
  
  cat("Prediction of beta-fields, (x$beta):\n")
  str(x$beta,1)
  cat("\n")
  
  if( !is.null(x$opts$transform) ){
    cat("Predictions for log-Gaussian field of type:",
        x$opts$transform, "\n\n")
  }

  if( x$opts$only.obs ){
    cat("Predictions only for", length(x$EX), "observations.\n")
  }else{
    cat("Predictions for", dim(x$EX)[1], "times at",
        dim(x$EX)[2], "locations.\n")
  }
  str(x[grep("EX",names(x))],1)
  cat("\n")
  
  if( length(grep("VX",names(x)))!=0 ){
    cat("Variances have been computed.\n")
    str(x[grep("VX",names(x))],1)
  }else{
    cat("Variances have NOT been computed.\n")
  }
  cat("\n")

  if( !is.null(x$opts$transform) ){
    if( length(grep("MSPE",names(x)))!=0 ){
      cat("Mean squared prediciton errors have been computed.\n")
      str(x[grep("MSPE",names(x))],1)
    }else{
      cat("Mean squared prediciton errors NOT been computed.\n")
    }
    cat("\n")
  }
  
  if( isTRUE(x$opts$LTA) ){
    cat(dim(x$LTA)[1], "temporal averages have been compute.\n")
    str(x["LTA"],1)
  }
  cat("\n")
  
  return(invisible())
}##function print.predictSTmodel


##' \code{\link[graphics:plot]{plot}} method for classes \code{predictSTmodel}
##' and \code{predCVSTmodel}. Provides several different plots of the
##' data.  
##'
##' @title Plots for \code{predictSTmodel} and \code{predCVSTmodel} Objects
##' 
##' @param x \code{predictSTmodel} or \code{predCVSTmodel} object to plot.
##' @param y Plot predictions as a function of either \code{"time"} or
##'   \code{"obs"}ervations. 
##' @param STmodel \code{STdata}/\code{STmodel} object containing observations
##'   with which to compare the predictions (not used for
##'   \code{plot.predCVSTmodel}). 
##' @param ID The location for which we want to plot predictions. A
##'   string matching names in \code{colnames(x$EX)} (or \code{x$I$ID},
##'   number(s) which are used as \code{ID = colnames(x$EX)[ID]}, or
##'   \code{"all"} in which case all predictions are used.
##'   If several locations are given (or \code{"all"}) then
##'   \code{y} must be \code{"obs"}.
##' @param col A vector of three colours: The first is the colour of the
##'   predictions, second for the observations and third for the polygon
##'   illustrating the confidence bands. \cr For \code{y="obs"} the colours are
##'   1) colour of the points, 2) colour of the 1-1 line, and 3) colour of the
##'   polygon. If \code{ID="all"}, picking \code{col[1]="ID"} will colour code
##'   the observations-prediction points by site ID.
##' @param pch,cex,lty,lwd Vectors with two elements giving the point type,
##'   size, line type and line width to use when plotting the predictions and
##'   observations respectively. Setting a value to \code{NA} will give no
##'   points/lines for the predictions/observations. \cr
##'   When plotting predictions
##'   as a function of observations \code{lty[2]} is used for the addition of
##'   \code{\link[graphics:abline]{abline}(0,1, lty=lty[2], col=col[2],
##'   lwd=lwd[2])}; \code{pch[2]} and \code{cex[2]} are ignored. 
##' @param p Width of the plotted confidence bands (as coverage percentage,
##'   used to find appropriate two-sided normal quantiles).
##' @param pred.type Which type of prediction to plot, one of
##'   \code{"EX"}, \code{"EX.mu"}, \code{"EX.mu.beta"}, or \code{"EX.pred"};
##'   see the output from \code{\link{predict.STmodel}}
##' @param pred.var Should we plot confidence bands based on prediction (TRUE)
##'   or confidence intrevalls (FALSE), see \code{\link{predict.STmodel}}.
##'   Only relevant if \code{pred.type="EX"} or \code{pred.type="EX.pred"}. \cr
##'   \strong{NOTE:} \emph{The default differs for \code{plot.predictSTmodel}
##'   and \code{plot.predCVSTmodel}!}
##' @param add Add to existing plot?
##' @param ... Additional parameters passed to
##'   \code{\link[graphics:plot]{plot}}.
##' 
##' @return Nothing
##'
##' @example Rd_examples/Ex_plot_predictSTmodel.R
##'
##' @author Johan Lindström
##' 
##' @family predictSTmodel methods
##' @method plot predictSTmodel
##' @export
plot.predictSTmodel <- function(x, y="time", STmodel=NULL, ID=x$I$ID[1],
                                col=c("black","red","grey"), pch=c(NA,NA),
                                cex=c(1,1), lty=c(1,1), lwd=c(1,1), p=0.95,
                                pred.type="EX", pred.var=FALSE,
                                add=FALSE, ...){
  ##check class belonging
  stCheckClass(x, "predictSTmodel", name="x")
  if( !is.null(STmodel) ){
    stCheckClass(STmodel, c("STmodel","STdata"), name="STmodel")
  }
  ##check for observations
  if( x$opts$only.pars ){
    message("Prediction structure contains only parameters.")
    ##return
    return(invisible())
  }
  ##we have to use y, cast to resonable name
  tmp <- internalPlotPredictChecks(y, pred.type, pred.var, x$opts$transform)
  plot.type <- tmp$plot.type
  pred.type <- tmp$pred.type
  pred.var <- tmp$pred.var

  ##check ID
  ID <- internalFindIDplot(ID, colnames(x$EX))
  ID.all <- internalCheckIDplot(ID, y)

  ##pick out data for plotting
  if( !ID.all ){
    ##pick out predictions, two cases
    if( x$opts$only.obs ){
      I1 <- x$I$ID==ID
      I2 <- 1
      date <- x$I$date[I1]
    }else{
      I1 <- 1:dim(x[[ pred.type[1] ]])[1]
      I2 <- ID
      date <- rownames(x[[ pred.type[1] ]])
    }
    sd <- x[[pred.var]][I1,I2]
    if( is.null(sd) ){
      sd <- NA
    }else{
      sd <- sqrt(sd)
    }
    pred <- data.frame(x=x[[ pred.type[1] ]][I1, I2], sd=sd, date=date,
                       x.log=NA, stringsAsFactors=FALSE)
    if( length(pred.type)==2 ){
      pred$x.log <- x[[ pred.type[2] ]][I1, I2]
    }
    ##pick out observations, if any
    obs <- STmodel$obs[STmodel$obs$ID==ID,,drop=FALSE]
    if( !is.null(obs) ){
      tmp <- obs
      obs <- matrix(NA, dim(pred)[1], 1)
      rownames(obs) <- as.character(pred$date)
      obs[as.character(tmp$date),] <- tmp$obs
    }
  }else{
    ##all only relevant when we have plot by observations
    ##pick out predictions, two cases
    if(length(ID)==1 && ID=="all"){
      ID <- unique(x$I$ID)
    }
    I.ID <- x$I$ID %in% ID
    I <- x$I$I[I.ID]
    sd <- x[[pred.var]][I]
    if( is.null(sd) ){
      sd <- NA
    }else{
      sd <- sqrt(sd)
    }
    pred <- data.frame(x=x[[ pred.type[1] ]][I], sd=sd,
                       date=x$I$date[I.ID], ID=x$I$ID[I.ID],
                       x.log=NA, stringsAsFactors=FALSE)
    if( length(pred.type)==2 ){
      pred$x.log <- x[[ pred.type[2] ]][I]
    }
    
    ##pick out observations, if any
    obs <- createDataMatrix(STmodel)
    I <- (match(x$I$ID[I.ID],colnames(obs))-1)*dim(obs)[1] +
        match(as.character(x$I$date[I.ID]), rownames(obs))
    obs <- obs[I]
    if( dim(pred)[1]!=length(obs) ){
      stop("Number of observed locations do not match no. predicted.")
    }
  }##if( !ID.all ){...}else{...}

  if( !is.null(x$opts$transform) && !is.null(obs) ){
    ##log-field, first fix observations
    obs <- exp(obs)
  }
  
  internalPlotPredictions(plot.type, ID, pred, obs, col, pch, cex, lty, lwd, p,
                          add, x$opts$transform, ...)
  
  ##return
  return(invisible())
}##function plot.predictSTmodel
