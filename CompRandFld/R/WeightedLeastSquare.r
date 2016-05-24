####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Emails: simone.padoan@unibocconi.it,
### moreno.bevilacqua@uv.cl
### Institutions: Department of Decision Sciences,
### University Bocconi of Milan and
### Departamento de Estadistica
### Universidad de Valparaiso
### File name: WeightedLeastSquare.r
### Description:
### This file contains a set of procedures in order
### to estimate the parameters of some covariance
### function models for a given dataset.
### Last change: 28/03/2013.
####################################################


### Procedures are in alphabetical order.

print.WLS <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    if(x$model=='Gaussian'){process <- x$model
                            model <- x$model}
    if(x$model=='ExtGauss'){process <- 'Max-Stable'
                            model <- 'Extremal Gaussian'}
    if(x$model=='BrowResn'){process <- 'Max-Stable'
                            model <- 'Brown-Resnick'}
    if(x$model=='ExtT'){process <- 'Max-Stable'
                        model <- 'Extremal T'}
    if(x$weighted)
      method <- 'Weighted Least Squares'
    else
      method <- method <- 'Least Squares'

    cat('\n##############################################################')
    cat('\nResults:', method,'Fitting of', process, 'Random Fields.\n')
    cat('\nModel used from the', method, ':', model, '\n')
    cat('\nCovariance model:', x$corrmodel, '\n')
    cat('Number of spatial coordinates:', x$numcoord, '\n')
    cat('Number of dependent temporal realisations:', x$numtime, '\n')
    cat('Number of replicates of the random field:', x$numrep, '\n')
    cat('Number of estimated parameters:', length(x$param), '\n')
    cat('The value of the', method, 'at the minimum:',-x$wls,'\n')
    cat('Number of spatial bins', length(x$bins),'\n')
    cat('Number of temporal bins', length(x$bint),'\n')
    cat('Min and max spatial distances:', x$srange,'\n')
    if(length(x$coordt)>1) cat('Min and max temporal interval:', x$trange,'\n')

    cat('\nEstimated parameters:\n')
    print.default(x$param, digits = digits, print.gap = 2,
                  quote = FALSE)

    cat('\n##############################################################\n')
    invisible(x)
  }


WlsInit <- function(coordx, coordy, coordt, corrmodel, data, distance, fcall, fixed, grid,
                    likelihood, margins, maxdist, maxtime, model, numblock, param, parscale,
                    paramrange, replicates, start, taper, tapsep, threshold, type, varest, vartype,
                    weighted, winconst, winstp)
  {
    ### Initialization parameters:
    initparam <- InitParam(coordx, coordy, coordt, corrmodel, data, distance, fcall, fixed,
                           grid, likelihood, margins, maxdist, maxtime,model, numblock,
                           param, parscale, paramrange, replicates, start, taper, tapsep,
                           threshold, "WLeastSquare", type, varest, vartype,
                           weighted, winconst, winstp)


    if(!is.null(initparam$error))
      stop(initparam$error)
    ### Set the initial type of likelihood objects:
    initparam$type <- CheckType(type)
    if(substr(model,1,6)=='Binary') return(initparam)
    if(is.null(start)) start <- NA else start <- unlist(start)
    if(is.null(fixed)) fixed <- NA else fixed <- unlist(fixed)
    ### Checks if all the starting values have been passed by the user:
    if(initparam$numstart==initparam$numparam)
      {### If full or pairwise then checks the mean parameter:
        if(model=='Gaussian' & (type %in% c('Standard','Pairwise','Tapering'))){
          if(is.na(fixed["mean"])){
              if(is.na(start["mean"])) initparam$param <- c(initparam$fixed["mean"], initparam$param)
              else initparam$param <- c(start["mean"], initparam$param)
            initparam$namesparam <- sort(names(initparam$param))
            initparam$param <- initparam$param[initparam$namesparam]
            initparam$numparam <- initparam$numparam+1
            initparam$flagnuis['mean'] <- 1
            initparam$numfixed <- initparam$numfixed-1
            if(initparam$numfixed > 0) initparam$fixed <- fixed
            else initparam$fixed <- NULL}
          else initparam$fixed['mean'] <- fixed["mean"]}
        return(initparam)}
    ### Updates if there are some starting values passed by the user:
    if(initparam$numstart > 0){
      initparam$param <- initparam$param[seq(1,initparam$numparam)[-pmatch(initparam$namesstart,initparam$namesparam)]]
      initparam$fixed <- c(initparam$fixed, initparam$start)}
    ### Define an internal function for the model fitting by least squares method:
    Lsquare <- function(bins, bint, corrmodel, fixed, fun, lenbins, moments,
                        namescorr, namesnuis, numbins, numbint, param)
      {
        param <- c(param, fixed)#set the parameters set:
        paramcorr <- param[namescorr]#set the correlation parameters:
        nuisance <- param[namesnuis]#set the nuisance parameters:
        #computes the weighted least squares:
        result <- .C(fun, as.double(bins), as.double(bint), as.integer(corrmodel),
                     as.double(lenbins), as.double(moments), as.integer(numbins),
                    as.integer(numbint), as.double(nuisance), as.double(paramcorr),
                     res=double(1), PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)$res
        return(result)
      }
    fname <- NULL ### the function name (variogram and least square)
    ###### ----------- START Estimation of the empirical variogram ---------- #####
    numbins <- as.integer(13) # number of spatial bins
    bins <- double(numbins) # vector of spatial bins
    numvario <- numbins-1
    moments <- double(numvario) # vector of spatial moments
    lenbins <- integer(numvario) # vector of spatial bin sizes
    ### Checks the type of variogram
    if(model=='ExtGauss' || model=='BrowResn' || model=='ExtT'){
      fname <- 'Binned_Madogram'
      data <- Dist2Dist(data, to='Uniform')}
    if(initparam$spacetime){### Computes the spatial-temporal variogram:
      numbint <- initparam$numtime-1 # number of temporal bins
      bint <- double(numbint)        # vector temporal bins
      momentt <- double(numbint)     # vector of temporal moments
      lenbint <- integer(numbint)    # vector of temporal bin sizes
      numbinst <- numvario*numbint    # number of spatial-temporal bins
      binst <- double(numbinst)      # spatial-temporal bins
      momentst <- double(numbinst)   # vector of spatial-temporal moments
      lenbinst <- integer(numbinst)  # vector of spatial-temporal bin sizes
      if(type=="Tapering") fname <- 'Binned_Variogram_st2'
      else                 fname <- 'Binned_Variogram_st'
      # Compute the spatial-temporal moments:
      EV=.C(fname, bins=bins, bint=bint, as.double(initparam$data), lenbins=lenbins,
      lenbinst=lenbinst, lenbint=lenbint,moments=moments, momentst=momentst,
      momentt=momentt, as.integer(numbins), as.integer(numbint),
         PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
         bins=EV$bins;bint=EV$bint;
         lenbins<-EV$lenbins;lenbint<-EV$lenbint;lenbinst<-EV$lenbinst;
         moments<-EV$moments;momentst<-EV$momentst;momentt<-EV$momentt
      indbin <- lenbins>0
      bins <- bins[indbin]
      moments <- moments[indbin]
      lenbins <- lenbins[indbin]
      numbins <- sum(indbin)
      indbint <- lenbint>0
      bint <- bint[indbint]
      momentt <- momentt[indbint]
      lenbint <- lenbint[indbint]
      numbint <- sum(indbint)
      indbinst <- lenbinst>0
      momentst <- momentst[indbinst]
      lenbinst <- lenbinst[indbinst]
      numbinst <- sum(indbinst)
      # Set the moment vectors and their sizes:
      moment <- matrix(momentst,nrow=numbins,ncol=numbint,byrow=TRUE)
      lenbin <- matrix(lenbinst,nrow=numbins,ncol=numbint,byrow=TRUE)
      moment <- rbind(momentt, moment)
      moment <- cbind(c(0,moments),moment)
      lenbin <- rbind(lenbint, lenbin)
      lenbin <- cbind(c(1,lenbins),lenbin)
      moments <- moment
      lenbins <- lenbin
      bins <- c(-bins[1],bins)
      bint <- c(0,bint)
      numbins <- numbins+1
      numbint <- numbint+1
    }
    else{### Computes the spatial variogram:
      if(type=="Tapering") fname <- 'Binned_Variogram2'
      else                 fname <- 'Binned_Variogram'
      numbint <- 1 # number of temporal bins
      bint <- double(numbint) # vector temporal bins
      momentt <- double(1) # vector of temporal moments
      momentst <- double(1)   # vector of spatial-temporal moments
      lenbint <- integer(1) # vector of temporal bin sizes
      lenbinst <- integer(1)  # vector of spatial-temporal bin sizes
      EV=.C(fname, bins=bins, as.double(data), lenbins=lenbins, moments=moments, as.integer(numbins),
         PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
      bins<-EV$bins;lenbins<-EV$lenbins;moments<-EV$moments
      centers <- bins[1:numvario]+diff(bins)/2
      indbin <- lenbins>0
      bins <- bins[indbin]
      centers <- centers[indbin]
      moments <- moments[indbin]
      lenbins <- lenbins[indbin]
      numbins <- sum(indbin)
      variogram <- moments/lenbins
      if(!is.na(initparam$param['scale']))
         initparam$param['scale'] <- centers[max(variogram) == variogram]}
    ###### ----------- END Estimation of the empirical variogram ---------- #####
    ###### ---------- START model fitting by least squares method ----------######
    if(model=='Gaussian') # Gaussian spatial or spatial-temporal case:
      fname <- 'LeastSquare_G'
    if(model=='ExtGauss') # max-stable case (extremal Gaussian):
      fname <- 'LeastSquare_MEG'
    if(model=='BrowResn') # max-stable case (Brown-Resnick):
      fname <- 'LeastSquare_MBR'
    if(model=='ExtT') # max-stable case (extremal-t):
      fname <- 'LeastSquare_MET'
    ### Fit the model:
    fitted <- optim(initparam$param, Lsquare, bins=bins, bint=bint, corrmodel=initparam$corrmodel,
                    fixed=initparam$fixed, fun=fname, lenbins=lenbins, moments=moments,
                    namescorr=initparam$namescorr, namesnuis=initparam$namesnuis,
                    numbins=numbins, numbint=numbint, control=list(fnscale=-1, reltol=1e-14,
                    maxit=1e8), hessian=FALSE)
    ###### ---------- END model fitting by least squares method ----------######
    ### Updates the parameter estimates:
    initparam$param <- fitted$par
    ### Updates with the starting values passed by the user:
    if(initparam$numstart > 0){
      initparam$param <- c(initparam$param,initparam$start)
      initparam$fixed <- initparam$fixed[seq(1,initparam$numstart+initparam$numfixed)[-pmatch(initparam$namesstart,names(initparam$fixed))]]
      initparam$namesparam <- sort(names(initparam$param))
      initparam$param <- initparam$param[initparam$namesparam]}
    #    for(i in 1 : initparam$numstart){
    #        initparam$param <- c(initparam$param, initparam$start[initparam$namesstart[i]])
    #        initparam$fixed <- initparam$fixed[!initparam$namesfixed==initparam$namesstart[i]]}
    #  initparam$param <- initparam$param[initparam$namesparam]}
    ### If full or pairwise then checks the mean parameter:
    if(model=='Gaussian' & (type %in% c('Standard','Pairwise','Tapering'))){
        if(is.na(fixed["mean"])){
          if(is.na(start["mean"])) initparam$param <- c(initparam$fixed["mean"], initparam$param)
          else initparam$param <- c(start["mean"], initparam$param)
          initparam$namesparam <- sort(names(initparam$param))
          initparam$param <- initparam$param[initparam$namesparam]
          initparam$numparam <- initparam$numparam + 1
          initparam$flagnuis["mean"] <- 1
          initparam$numfixed <- initparam$numfixed - 1
          if(initparam$numfixed > 0) initparam$fixed <- fixed
          else initparam$fixed <- NULL}
        else initparam$fixed["mean"] <- fixed["mean"]}
   return(initparam)
  }


WLeastSquare <- function(data, coordx, coordy=NULL, coordt=NULL, corrmodel, distance="Eucl",
                         fixed=NULL,grid=FALSE, maxdist=NULL, maxtime=NULL, model='Gaussian',
                         optimizer='Nelder-Mead', numbins=NULL, replicates=1, start=NULL,
                         weighted=FALSE)
  {
    ### Check first if the model is not binary:
    if(substr(model,1,6)=='Binary') stop("The weighted least squares method can not be used with binary data")

    call <- match.call()
    ### Check the parameters given in input:
    checkinput <- CheckInput(coordx, coordy, coordt, corrmodel, data, distance,"Fitting", fixed, grid, 'None',
                             "Frechet", maxdist, maxtime, model, NULL, optimizer, NULL, replicates, start, NULL,
                             NULL, NULL, 'WLeastSquare', FALSE, 'SubSamp', weighted)

    if(!is.null(checkinput$error))
      stop(checkinput$error)
    # check the number of bins:
    if(!is.null(numbins) & !is.integer(numbins))
      if(numbins < 0)
        stop('insert a positive integer value for the number of bins\n')
    # set the default number of spatial bins:
    if(is.null(numbins))
      numbins <- 13
    ### Define the object function for the weighted least squares method:
    WLsquare <- function(bins, bint, corrmodel, fixed, fun, lenbins, moments,
                         namescorr, namesnuis, numbins, numbint, param)
      {
        param <- c(param, fixed)#set the parameters set:
        paramcorr <- param[namescorr]#set the correlation parameters:
        nuisance <- param[namesnuis]#set the nuisance parameters:
        #computes the weighted least squares:
        result <- .C(fun, as.double(bins), as.double(bint), as.integer(corrmodel),
                     as.double(lenbins), as.double(moments), as.integer(numbins),
                     as.integer(numbint), as.double(nuisance), as.double(paramcorr),
                     res=double(1), PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)$res
       return(result)
      }
    ### Initializes global variables:
    WLeastSquare <- NULL
    fname <- NULL
    variogramt <- NULL
    variogramst <- NULL
    ### Initializes the parameter values:
    parscale <- NULL
    initparam <- InitParam(coordx, coordy, coordt, corrmodel, data, distance, "Fitting", fixed, grid,
                           'None', "Frechet", maxdist, maxtime, model, NULL, NULL,
                           parscale, optimizer=='L-BFGS-B', replicates, start,NULL, NULL, NULL,
                           'WLeastSquare', 'WLeastSquare', FALSE, 'SubSamp', FALSE, 1, 1)
    if(!is.null(initparam$error))
      stop(initparam$error)
    ###### ----------- START Estimation of the empirical variogram ---------- #####
    numvario <- numbins-1
    bins <- double(numbins) # vector of spatial bins
    moments <- double(numvario) # vector of spatial moments
    lenbins <- integer(numvario) # vector of spatial bin sizes
    ### Checks the type of variogram:
    fname <- 'Binned_Variogram'
    if(model=='ExtGauss' || model=='BrowResn' || model=='ExtT'){
      initparam$data <- Dist2Dist(initparam$data, to='Uniform')
      fname <- 'Binned_Madogram'}
    if(initparam$spacetime){### Computes the spatial-temporal variogram:
      numbint <- initparam$numtime-1 # number of temporal bins
      bint <- double(numbint)        # vector temporal bins
      momentt <- double(numbint)     # vector of temporal moments
      lenbint <- integer(numbint)    # vector of temporal bin sizes
      numbinst <- numvario*numbint    # number of spatial-temporal bins
      binst <- double(numbinst)      # spatial-temporal bins
      momentst <- double(numbinst)   # vector of spatial-temporal moments
      lenbinst <- integer(numbinst)  # vector of spatial-temporal bin sizes
      fname <- 'Binned_Variogram_st'
      # Compute the spatial-temporal moments:
      EV=.C(fname, bins=bins, bint=bint, as.double(initparam$data), lenbins=lenbins, lenbinst=lenbinst,
      lenbint=lenbint,moments=moments,momentst=momentst,momentt=momentt,
      as.integer(numbins), as.integer(numbint),PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
      bins <- EV$bins
      binst <- EV$binst
      bint <- EV$bint;
      lenbins<-EV$lenbins;lenbinst<-EV$lenbinst;lenbint<-EV$lenbint;
      moments<-EV$moments;momentt<-EV$momentt;momentst<-EV$momentst;
      indbin <- lenbins>0
      centers <- bins[1:numvario]+diff(bins)/2
      bins <- bins[indbin]
      centers <- centers[indbin]
      moments <- moments[indbin]
      lenbins <- lenbins[indbin]
      # Computes the spatial marginal variogram:
      variograms <- moments/lenbins
      numbins <- sum(indbin)
      indbint <- lenbint>0
      bint <- bint[indbint]
      momentt <- momentt[indbint]
      lenbint <- lenbint[indbint]
      numbint <- sum(indbint)
      # Computes the temporal marginal variogram:
      variogramt <- momentt/lenbint
      indbinst <- lenbinst>0
      momentst <- momentst[indbinst]
      lenbinst <- lenbinst[indbinst]
      numbinst <- sum(indbinst)
      # Computes the spatial-temporal variogram:
      variogramst <- momentst/lenbinst
      # Set the moment vectors and their sizes:
      moment <- matrix(momentst,nrow=numbins,ncol=numbint,byrow=TRUE)
      lenbin <- matrix(lenbinst,nrow=numbins,ncol=numbint,byrow=TRUE)
      moment <- rbind(momentt, moment)
      moment <- cbind(c(0,moments),moment)
      lenbin <- rbind(lenbint, lenbin)
      lenbin <- cbind(c(1,lenbins),lenbin)
      moments <- moment
      lenbins <- lenbin
      bins <- c(-bins[1],bins)
      bint <- c(0,bint)
      numbins <- numbins+1
      numbint <- numbint+1}
      # Set an initial value for the scale parameter:
    #if(!is.null(initparam$param['scale_s']))
          #initparam$param['scale_s'] <- bins[max(variograms)==variograms]}
    else{### Computes the spatial variogram:
      numbint <- 1 # number of temporal bins
      bint <- double(numbint) # vector temporal bins
      momentt <- double(1) # vector of temporal moments
      momentst <- double(1)   # vector of spatial-temporal moments
      lenbint <- integer(1) # vector of temporal bin sizes
      lenbinst <- integer(1)  # vector of spatial-temporal bin sizes
      EV=.C(fname, bins=bins, as.double(initparam$data), lenbins=lenbins, moments=moments, as.integer(numbins),
         PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
      bins<-EV$bins;lenbins<-EV$lenbins;moments<-EV$moments
      indbin <- lenbins>0
      centers <- bins[1:numvario]+diff(bins)/2
      bins <- bins[indbin]
      centers <- centers[indbin]
      moments <- moments[indbin]
      lenbins <- lenbins[indbin]
      numbins <- sum(indbin)
      variograms <- moments/lenbins
      # Set an initial value for the scale parameter:
      if(!is.null(initparam$param['scale']))
        initparam$param['scale'] <- bins[max(variograms)==variograms]}
    ###### ----------- END Estimation of the empirical variogram ---------- #####

    ###### ---------- START model fitting by weighted least squares method ----------######
    if(model=='Gaussian') # Gaussian random field:
      if(weighted) fname <- 'WLeastSquare_G'
      else fname <- 'LeastSquare_G'
    # Max-stable random field (extremal Gaussian):
    if(model=='ExtGauss')
      if(weighted) fname <- 'WLeastSquare_MEG'
      else fname <- 'LeastSquare_MEG'
    if(model=='BrowResn')
      if(weighted) fname <- 'WLeastSquare_MBR'
      else fname <- 'LeastSquare_MBR'
    if(model=='ExtT')
      if(weighted) fname <- 'WLeastSquare_MET'
      else fname <- 'LeastSquare_MET'
    ### Computes estimates by the weighted least squares method:
    if(optimizer=='L-BFGS-B')
      fitted <- optim(initparam$param, WLsquare, bins=bins, bint=bint, corrmodel=initparam$corrmodel,
                      fixed=initparam$fixed, fun=fname, lenbins=lenbins, method=optimizer,
                      moments=moments, namescorr=initparam$namescorr, namesnuis=initparam$namesnuis,
                      numbins=numbins, numbint=numbint, control=list(fnscale=-1, factr=1, pgtol=1e-14,
                      maxit = 1e8), lower=initparam$lower, upper=initparam$upper, hessian=FALSE)
    else
      fitted <- optim(initparam$param, WLsquare, bins=bins, bint=bint, corrmodel=initparam$corrmodel,
                      fixed=initparam$fixed, fun=fname, lenbins=lenbins, method=optimizer,
                      moments=moments, namescorr=initparam$namescorr, namesnuis=initparam$namesnuis,
                      numbins=numbins, numbint=numbint, control=list(fnscale=-1, reltol=1e-14, maxit=1e8),
                      hessian=FALSE)
    ###### ---------- END model fitting by weighted least squares method ----------######
    ### Removes the global variobales:
    .C('DeleteGlobalVar', PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
    ### Set the output:
    WLeastSquare <- list(bins=bins,
                         bint=bint,
                         centers=centers,
                         coordx = initparam$coordx,
                         coordy = initparam$coordy,
                         coordt = initparam$coordt,
                         convergence = fitted$convergence,
                         corrmodel = corrmodel,
                         data = initparam$data,
                         fixed = initparam$fixed,
                         grid = grid,
                         iterations = fitted$counts,
                         message = fitted$message,
                         model=model,
                         numcoord=initparam$numcoord,
                         numrep=initparam$numrep,
                         numtime=initparam$numtime,
                         param = fitted$par,
                         srange=initparam$srange,
                         trange=initparam$trange,
                         variograms = variograms,
                         variogramt = variogramt,
                         variogramst = variogramst,
                         weighted = weighted,
                         wls = fitted$value)
    structure(c(WLeastSquare, call = call), class = c("WLS"))
  }
