####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Emails: simone.padoan@unibocconi.it,
### moreno.bevilacqua@uv.cl
### Institutions: Department of Decision Sciences,
### University Bocconi of Milan and
### Departamento de Estadistica
### Universidad de Valparaiso
### File name: EmpiricalEstimators.r
### Description:
### This file contains a set of procedures in order
### to estimate the empircal covariance or the extreme
### dependence structures for a given dataset.
### Last change: 28/03/2013.
####################################################

### Procedures are in alphabetical order.

EVariogram <- function(data, coordx, coordy=NULL, coordt=NULL, cloud=FALSE, distance="Eucl",
                       grid=FALSE, gev=c(0,1,0), maxdist=NULL, maxtime=NULL, numbins=NULL,
                       replicates=1, type='variogram')
  {
    call <- match.call()
    corrmodel <- 'gauss'
    ### Check the parameters given in input:
    if(is.null(type))
      type <- 'variogram'
    # Checks if its a variogram or a madogram
    if(!is.null(type) & all(type!="lorelogram", type!='variogram', type!='madogram', type!='Fmadogram'))
      stop('the admitted types are: variogram, madogram (Fmadogram) or lorelogram\n')
    if(type=="lorelogram" & cloud)
        stop("there is not cloud version of lorelogram\n")
    # Set the type of model:
    if(type=='variogram'){
        model <- 'Gaussian'
        fname <- 'Binned_Variogram'}
    if(type=="lorelogram"){
        model <- "BinaryGauss"
        fname <- "Binned_Lorelogram"}
    if(type=="madogram" || type=="Fmadogram") model <- 'ExtGauss'
    # Checks if its a spatial or spatial-temporal random field:
    if(!is.null(coordt))
      if(is.numeric(coordt))
        if(length(coordt)>1) corrmodel <- 'gneiting'
    # Checks the input:
    checkinput <- CheckInput(coordx, coordy, coordt, corrmodel, data, distance, "Fitting", NULL, grid,
                             'None', "Frechet", maxdist, maxtime, model, NULL, 'Nelder-Mead', NULL,
                             replicates, NULL, NULL,NULL, 0, 'WLeastSquare', FALSE, 'SubSamp', FALSE)
    # Checks if there are errors in the input:
    if(!is.null(checkinput$error))
      stop(checkinput$error)
    ### START -- Specific checks of the Empirical Variogram:
    if(!is.null(cloud) & !is.logical(cloud))
      stop('insert a logical value (TRUE/FALSE) for the cloud parameter\n')
    if(type=='madogram' & (!is.numeric(gev) || is.null(gev) || !length(gev)==3))
      stop('insert a numeric vector with the three GEV parameters\n')
    if(!is.null(numbins) & !is.integer(numbins))
      if(numbins < 0)
        stop('insert a positive integer value for the number of bins\n')
    if(is.null(numbins))
      numbins <- 13
    ### END -- Specific checks of the Empirical Variogram

    ### Initialization parameters:
    initparam <- InitParam(coordx, coordy, coordt, corrmodel, data,distance, "Fitting",
                           NULL, grid, 'None', "Frechet", maxdist,
                           maxtime, model, NULL, NULL, FALSE, FALSE, replicates,
                           NULL, NULL, NULL, 0, 'WLeastSquare', 'WLeastSquare', FALSE,
                           'SubSamp', FALSE, 1, 1)
    # Checks if there are inconsistences:
    if(!is.null(initparam$error))
      stop(initparam$error)
    numvario <- numbins-1
    if(cloud){
        numbins <- numvario <- initparam$numpairs
        fname <- 'Cloud_Variogram'}
    ### Estimation of the empirical spatial or spatial-temporal variogram:
    bins <- double(numbins) # spatial bins
    moments <- double(numvario) # vector of spatial moments
    lenbins <- integer(numvario) # vector of spatial bin sizes
    bint <- NULL
    lenbinst <- NULL
    lenbint <- NULL
    variogramst <- NULL
    variogramt <- NULL
    if(type=='madogram'){ # Trasform to GEV margins:
      if(cloud) fname <- 'Cloud_Madogram' else fname <- 'Binned_Madogram'
        if(gev[3]>=1)
          stop('the shape parameter can not be greater or equal to 1')
        initparam$data <- Dist2Dist(initparam$data, to='Gev',
                                    loc=rep(gev[1], initparam$numcoord),
                                    scale=rep(gev[2], initparam$numcoord),
                                    shape=rep(gev[3], initparam$numcoord))}
    if(type=='Fmadogram'){ # Transform to Uniform margins:
      if(cloud) fname <- 'Cloud_Madogram' else fname <- 'Binned_Madogram'
      initparam$data <- Dist2Dist(initparam$data, to='Uniform')}
    if(initparam$spacetime){
      numbint <- initparam$numtime-1 # number of temporal bins
      bint <- double(numbint)        # temporal bins
      momentt <- double(numbint)     # vector of temporal moments
      lenbint <- integer(numbint)    # vector of temporal bin sizes
      numbinst <- numvario*numbint   # number of spatial-temporal bins
      binst <- double(numbinst)      # spatial-temporal bins
      momentst <- double(numbinst)   # vector of spatial-temporal moments
      lenbinst <- integer(numbinst)  # vector of spatial-temporal bin sizes
      if(cloud) fname <- 'Cloud_Variogram_st' else fname <- 'Binned_Variogram_st'
      if(type=="lorelogram") fname <- "Binned_Lorelogram_st"
      # Compute the spatial-temporal moments:
      EV=.C(fname, bins=bins, bint=bint, as.double(initparam$data),
            lenbins=lenbins, lenbinst=lenbinst, lenbint=lenbint,moments=moments, momentst=momentst,momentt=momentt,
            as.integer(numbins), as.integer(numbint),PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
      bins<-EV$bins;bint<-EV$bint;
      moments<-EV$moments;  momentst<-EV$momentst;  momentt<-EV$momentt;
      lenbins<-EV$lenbins;lenbint<-EV$lenbint;lenbinst<-EV$lenbinst;
      centers <- bins[1:numvario]+diff(bins)/2
      if(type=="lorelogram"){
      elorel <- matrix(momentst,nrow=numvario,ncol=numbint,byrow=TRUE)
      elorel <- rbind(c(0,momentt),cbind(moments,elorel))
      a <- rowSums(elorel)==0
      b <- colSums(elorel)==0
      d <- rowSums(elorel)!=0
      f <- colSums(elorel)!=0
      if(sum(a)){
          elorel <- elorel[-which(a),]
          centers <- c(0,centers)[which(d)]
          centers <- centers[-1]}
      if(sum(b)){
          elorel <- elorel[,-which(b)]
          bint <- c(0,bint)[which(f)]
          bint <- bint[-1]}
      elorel[elorel==0] <- NA
      moments <- as.vector(elorel[,1][-1])
      momentt <- as.vector(elorel[1,][-1])
      momentst <- c(t(elorel[-1,-1]))
      lenbins <- rep(1,length(moments))
      lenbint <- rep(1,length(momentt))
      lenbinst <- rep(1,length(momentst))}
      indbin <- lenbins>0
      indbint <- lenbint>0
      indbinst <- lenbinst>0
      bins <- bins[indbin]
      bint <- bint[indbint]
      centers <- centers[indbin]
      moments <- moments[indbin]
      lenbins <- lenbins[indbin]
      momentt <- momentt[indbint]
      lenbint <- lenbint[indbint]
      momentst <- momentst[indbinst]
      lenbinst <- lenbinst[indbinst]
      # Computes the spatial marginal variogram:
      variograms <- moments/lenbins
      # Computes the temporal marginal variogram:
      variogramt <- momentt/lenbint
      # Computes the spatial-temporal variogram:
      variogramst <- momentst/lenbinst}
    else{# Computes the spatial moments
      EV=.C(fname, bins=bins, as.double(initparam$data), lenbins=lenbins, moments=moments, as.integer(numbins),
         PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
       bins<-EV$bins
       moments<-EV$moments
       lenbins<-EV$lenbins
      # Computes the spatial variogram:
      indbin <- lenbins>0
      bins <- bins[indbin]
      numbins <- length(bins)
      # check if cloud or binned variogram:
      if(cloud) centers <- bins else centers <- bins[1:(numbins-1)]+diff(bins)/2
      moments <- moments[indbin]
      lenbins <- lenbins[indbin]
      variograms <- moments/lenbins}
    # Start --- compute the extremal coefficient
    extcoeff <- NULL
    if(type == 'madogram'){ # Check if its a madogram:
      if(gev[3] == 0)
        extcoeff <- exp(variograms / gev[2])
      else
        extcoeff <- Dist2Dist(gev[1] + variograms / gamma(1 - gev[3]))
      extcoeff[extcoeff>2] <- NA}
    if(type == 'Fmadogram'){ # Check if its a Fmadogram:
      extcoeff <- (1 + 2 * variograms) / (1 - 2 * variograms)
      extcoeff[extcoeff>2] <- NA}

    .C('DeleteGlobalVar', PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)

    EVariogram <- list(bins=bins,
                       bint=bint,
                       cloud=cloud,
                       centers=centers,
                       extcoeff=extcoeff,
                       lenbins=lenbins,
                       lenbinst=lenbinst,
                       lenbint=lenbint,
                       srange=initparam$srange,
                       variograms=variograms,
                       variogramst=variogramst,
                       variogramt=variogramt,
                       trange=initparam$trange,
                       type=type)

    structure(c(EVariogram, call = call), class = c("Variogram"))

  }

