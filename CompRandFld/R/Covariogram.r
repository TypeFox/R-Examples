####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Emails: simone.padoan@unibocconi.it,
### moreno.bevilacqua@uv.cl
### Institutions: Department of Decision Sciences,
### University Bocconi of Milan and
### Departamento de Estadistica
### Universidad de Valparaiso
### File name: Covariogram.r
### Description:
### This file contains a set of procedures
### to compute and plot the estimated covariance
### function and the variogram after fitting a
### random field by composite-likelihood.
### Last change: 28/03/2013.
####################################################


### Procedures are in alphabetical order.

### Compute and plot the (estimated) covariance function and the variogram
### from a fitted model obtain from the FitComposite or the WLeastSquare procedure
Covariogram <- function(fitted, lags=NULL, lagt=NULL, answer.cov=FALSE, answer.vario=FALSE,
                        answer.extc=FALSE, answer.range=FALSE, fix.lags=NULL, fix.lagt=NULL,
                        show.cov=FALSE, show.vario=FALSE, show.extc=FALSE, show.range=FALSE,
                        add.cov=FALSE, add.vario=FALSE, add.extc=FALSE, pract.range=95,
                        vario=NULL, ...)
  {
    result <- NULL
    # define the bivariate Gaussian distribution:
    vpbnorm <- function(corrmodel,lags,lagt,nuisance,numlags,numlagt,param,threshold)
    {
        rho=double(numlags*numlagt)
        p=.C("vpbnorm",  as.integer(corrmodel), as.double(lags), as.double(lagt),
           as.integer(numlags), as.integer(numlagt), as.double(nuisance),
           as.double(param), rho=double(numlags*numlagt), as.double(threshold), PACKAGE='CompRandFld',
           DUP=TRUE, NAOK=TRUE)
        return(p$rho)
    }
    # define the correlation function:
    CorrelationFct <- function(corrmodel, lags, lagt, numlags, numlagt, param)
    {
        corr=double(numlags*numlagt)
        p=.C('VectCorrelation', corr=double(numlags*numlagt), as.integer(corrmodel), as.double(lags),
           as.integer(numlags), as.integer(numlagt), as.double(param),
           as.double(lagt), PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
        return(p$corr)
    }
    # define the extremal coefficient function:
    ExtremalCoeff <- function(corrmodel, lags, model, nuisance, numlags, param)
    {
        

        p=.C('ExtCoeff', as.integer(corrmodel), extc=double(numlags), as.double(lags),
           as.integer(model), as.integer(numlags), as.double(nuisance),
           as.double(param), PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
        return(p$extc)
     }
    # Pratical range in the Gaussian case:
    PracRangeNorm <- function(corrmodel, lags, lagt, nuisance, numlags, numlagt, param, pract.range)
        return(nuisance["nugget"]+nuisance["sill"]*CorrelationFct(corrmodel, lags, lagt, numlags, numlagt, param)-
               (nuisance["nugget"]+nuisance["sill"]*(1-pract.range)))
    # Pratical range in the max-stable case:
    PracRangeExt <- function(corrmodel, lags, model, nuisance, numlags, param, pract.range)
                return(ExtremalCoeff(corrmodel, lags, model, nuisance, numlags, param)-(2*pract.range))

    isextc <- FALSE
    isvario <- !is.null(vario) # is empirical variogram is passed?
    if(isvario) isextc <- !is.null(vario$extcoeff)
    ispatim <- fitted$numtime>1# is a space-time random field?
    # START ---- check input --- #
    if(!class(fitted)=='FitComposite' & !class(fitted)=='WLS')
        stop("Enter an object obtained FitComposite or WLeastSquare")

    if(isvario & !class(vario)=='Variogram')
        stop("Enter an object obtained from the function EVariogram")

    if(!is.numeric(pract.range) & answer.range)
        stop("Enter a number for the parameter % of sill")
    else{
        if(pract.range < 0 || pract.range > 100)
            stop("Entered an incorrect value for the % of sill")
        else
          pract.range <- pract.range / 100}
    # Set the type of process:
    model <- CheckModel(fitted$model)
    gaussian <- model==1
    binary <- model==2
    maxstab <- model>=3
    # determine the spatial(temporal) distances:
    if(binary) {slow <- 1e-3; zero <- NA} else {slow <- 0; zero <- 0}
    if(is.null(lags)) lags <- seq(slow,fitted$srange[2],length=100)
    # determine the spatio-temporal distances:
    if(is.null(lagt)) lagt <- seq(0,fitted$trange[2],by=1)
    if(isvario & ispatim){
        lags <- c(0,vario$centers)
        lagt <- c(0,vario$bint)}
    numlags <- length(lags)
    numlagt <- length(lagt)
    # check the fixed lags:
    if(!is.null(fix.lags)){
        if(!is.numeric(fix.lags) || fix.lags<0 || fix.lags>numlags){
            cat("invalid spatial fixed lag\n")
            return(result)}}
    if(!is.null(fix.lagt))
        if(!is.numeric(fix.lagt) || fix.lagt<0 || fix.lagt>numlagt){
            cat("invalid temporal fixed lag\n")
            return(result)}
    # END ---- check input --- #
    # set range (for practical range)
    if(fitted$corrmodel=="matern"){
        lower <- 1e-10
        upper <- 1e20}
    else{
        lower <- 0
        upper <- 1e100}
    param <- c(fitted$fixed,fitted$param)[CorrelationParam(fitted$corrmodel)]
    nuisance <- c(fitted$fixed,fitted$param)[NuisanceParam(fitted$model)]
    corrmodel <- CheckCorrModel(fitted$corrmodel)
    correlation <- CorrelationFct(corrmodel, lags, lagt, numlags, numlagt, param)
    # Gaussian random field:
    if(gaussian){
        # covariance function
        covariance <- nuisance["nugget"]+nuisance["sill"]*correlation
        # variogram function
        variogram <- nuisance["nugget"]+nuisance["sill"]*(1-correlation)
        vario.main <- "Spatial variogram"
        vario.ylab <- "Variogram"
        if(ispatim){
            dim(covariance) <- c(numlags, numlagt)
            dim(variogram) <- c(numlags, numlagt)
            vario.main <- "Space-time variogram"
            vario.zlab <- "Variogram"
        }
        # compute the practical range:
        if(show.range || answer.range){
            Range <- uniroot(PracRangeNorm, c(lower, upper), corrmodel=corrmodel,
                             nuisance=nuisance, numlags=1, numlagt=1, lagt=lagt,
                             param=param, pract.range=pract.range)$root}}
    # max-stable random field:
    if(maxstab){
        # covariance function
        covariance <- 1-nuisance["sill"]+nuisance["sill"]*correlation
        # extremal coefficient function
        extrcoeff <- ExtremalCoeff(corrmodel, lags, model, nuisance, numlags, param)
        # compute the practical range:
        Range <- uniroot(PracRangeExt,c(lower, upper),corrmodel=corrmodel,model=model,
                         nuisance=nuisance,numlags=1,param=param,pract.range=pract.range)$root}
    # binary random field:
    if(binary){
        covariance <- nuisance["nugget"]+nuisance["sill"]*correlation
        p <- pnorm(nuisance["mean"])
        q <- vpbnorm(corrmodel, lags, lagt, nuisance,
                     numlags, numlagt, param, fitted$threshold)
        variogram <- log(q*(1-2*p+q)/(p-q)^2)
        vario.main <- "Spatial lorelogram"
        vario.ylab <- "Lorelogram"
        if(ispatim){
            dim(covariance) <- c(numlags, numlagt)
            dim(variogram) <- c(numlags, numlagt)
            vario.main <- "Space-time lorelogram"
            vario.zlab <- "Lorelogram"}}
    # display the covariance function
    if(show.cov){
        if(ispatim){# spatio-temporal case:
            # build the covariance matrix:
            plagt <- !is.null(fix.lags)
            plags <- !is.null(fix.lagt)
            numplot <- 1+plags+plagt
            par(mfrow=c(1,numplot))
            # temporal section
            par(mai=c(.2,.2,.2,.2))
            persp(lags, lagt, covariance, xlab="Distance", ylab="Time",
                  zlab="Covariance", ltheta=90,
                  shade=0.75, ticktype="detailed", phi=30,
                  theta=30,main="Space-time covariance",
                  , cex.axis=.8, cex.lab=.8,zlim=c(0,max(covariance))) #
            if(plagt){
                par(mai=c(.5,.5,.3,.3),mgp=c(1.6,.6,0))
                plot(lagt, covariance[fix.lags,], xlab="Time",
                     ylab="Covariance", type="l",cex.axis=.8,cex.lab=.8,
                     main="Space-time cov: temporal profile",...)}
            if(plags){
                par(mai=c(.5,.5,.3,.3),mgp=c(1.6,.6,0))
                plot(lags, covariance[,fix.lagt], xlab="Distance",
                     ylab="Covariance", type="l",cex.axis=.8,cex.lab=.8,
                     main="Space-time cov: spatial profile",...)}
        }
        else{# spatial case:
            if(add.cov & dev.cur()!=1){
                lines(lags, covariance,...)
                if(show.range) abline(v=Range)}
            else{
                plot(lags, covariance, type='l', ylim=c(min(covariance),
                     max(covariance)), main="Spatial covariance",
                     xlab="Distance", ylab="Covariance",...)
                if(show.range) abline(v=Range)}}}
    # display the variogram function
    if(show.vario){
        if(ispatim){# spatio-temporal case:
            plagt <- !is.null(fix.lags)
            plags <- !is.null(fix.lagt)
            nrowp <- 1
            ncolp <- 1
            if(isvario){
                ncolp <- ncolp+1
                if(plags || plagt) nrowp <- nrowp+1}
            else ncolp <- ncolp+plags+plagt
            par(mfrow=c(nrowp,ncolp))
            sup <- 0
            tup <- 0
            if(isvario){
                nbins <- length(vario$centers)
                nbint <- length(vario$bint)
                evario <- matrix(vario$variogramst,nrow=nbins,
                                 ncol=nbint,byrow=TRUE)
                evario <- rbind(c(zero,vario$variogramt),cbind(vario$variograms,evario))
                evario.grid <- expand.grid(c(0,vario$centers),c(0,vario$bint))
                scatterplot3d::scatterplot3d(evario.grid[,1],evario.grid[,2], c(evario),
                              type="h",highlight.3d=TRUE,cex.axis=.7,cex.lab=.7,
                              main=paste("Empirical",vario.main),xlab="Distance",
                              ylab="Time",zlab=vario.zlab,mar=c(2,2,2,2),mgp=c(0,0,0))
                if(plagt) tup <- max(evario[fix.lags,],na.rm=TRUE)
                if(plags) sup <- max(evario[,fix.lagt],na.rm=TRUE)
            }
            par(mai=c(.2,.2,.2,.2))
            persp(lags, lagt, variogram, xlab="Distance",
                  ylab="Time", zlab=vario.zlab, ltheta=90,
                  shade=0.75, ticktype="detailed", phi=30,
                  theta=30,main=vario.main, cex.axis=.8,
                   cex.lab=.8)  #zlim=c(0,max(variogram))
            if(plagt){
                par(mai=c(.5,.5,.3,.3),mgp=c(1.6,.6,0))
                plot(lagt, variogram[fix.lags,], xlab="Time",cex.axis=.8,cex.lab=.8,
                     ylab=vario.ylab, type="l", ylim=c(min(0,min(variogram)),max(nuisance["nugget"]+
                     nuisance["sill"],tup)), main=paste(vario.ylab,": temporal profile",
                     sep=""),...)
                if(isvario) points(lagt, evario[fix.lags,],...)}
            if(plags){
                par(mai=c(.5,.5,.3,.3),mgp=c(1.6,.6,0))
                plot(lags, variogram[,fix.lagt], xlab="Distance",cex.axis=.8,cex.lab=.8,
                     ylab=vario.ylab, type="l", ylim=c(min(0,min(variogram)),max(nuisance["nugget"]+
                     nuisance["sill"],sup)), main=paste(vario.ylab,": spatial profile",
                     sep=""),...)
                if(isvario) points(lags, evario[,fix.lagt],...)}}
        else{# spatial case:
            if(add.vario & dev.cur()!=1){
                points(vario$centers, vario$variograms,...)
                lines(lags, variogram,...)
                if(show.range) abline(v=Range)}
            else{
                bnds <- range(variogram)
                bnds[1] <- min(bnds[1], vario$variograms)
                bnds[2] <- max(bnds[2], vario$variograms)
                plot(lags, variogram, type='l', ylim=c(bnds[1], bnds[2]),
                     main=vario.main,xlab="Distance",
                     ylab=vario.ylab,...)
                points(vario$centers, vario$variograms,...)
                if(show.range) abline(v=Range)}}}
    # display the extremal function
    if(show.extc){
        if(add.extc & dev.cur()!=1){
            if(isvario & isextc) points(vario$centers, vario$extcoeff,...)
            lines(lags, extrcoeff,...)
            if(show.range) abline(v=Range)}
        else{
            bnds <- range(extrcoeff)
            if(isvario & isextc){
                bnds[1] <- min(bnds[1], vario$extcoeff)
                bnds[2] <- max(bnds[2], vario$extcoeff)}
            plot(lags, extrcoeff, type='l',ylim=c(bnds[1], bnds[2]),
                 main="Spatial extremal coefficient",ylab="Extremal coefficient",
                 xlab="Distance",...)
            if(isvario & isextc)
                points(vario$centers, vario$extcoeff,...)
            if(show.range) abline(v=Range)}}
    # return the estimated covariance function
    if(answer.cov) result <- list(lags=lags, covariance=covariance)
    # return the estimated variogram function
    if(answer.vario)
        if(!is.list(result)) result <- list(lags=lags, variogram=variogram)
        else result$variogram <- variogram
    # return the range
    if(answer.range)
        if(!is.list(result)) result <- list(range=Range)
        else result$range <- Range
    # return the estimated extremal coefficient function
    if(answer.extc)
        if(!is.list(result)) result <- list(lags=lags, extrcoeff=extrcoeff)
        else result$extrcoeff <- extrcoeff
    if(!is.null(result))
    return(result)
  }
