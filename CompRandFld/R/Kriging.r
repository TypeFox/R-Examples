####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Emails: simone.padoan@unibocconi.it,
### moreno.bevilacqua@uv.cl
### Institutions: Department of Decision Sciences,
### University Bocconi of Milan and
### Departamento de Estadistica
### Universidad de Valparaiso
### File name: Kriging.r
### Description:
### This file contains a set of procedures
### for computing simple (tapered) and ordinary kriging
### predictor  at an unknown space (time) locations.
### Last change: 28/03/2013.
####################################################

Kri<- function(data, coordx, coordy=NULL, coordt=NULL, corrmodel, distance="Eucl", grid=FALSE, loc, maxdist=NULL,
               maxtime=NULL, param, taper=NULL, tapsep=NULL, time=NULL, type="Standard", type_krig="Simple")

{
    CheckDistance<- function(distance)
    {
        CheckDistance <- NULL
        CheckDistance <- switch(distance,
                                eucl=0,
                                Eucl=0,
                                chor=1,
                                Chor=1,
                                geod=2,
                                Geod=2,
                                proj=3,
                                Proj=3)
        return(CheckDistance)
    }
    ### computing covariance matrix
    covmatrix <- Covmatrix(coordx, coordy, coordt, corrmodel, distance, grid, TRUE, maxdist, maxtime, param, taper, tapsep, type)
     ### control on input kriging
    checkinput <- CheckInput(loc, NULL, time, covmatrix, data, NULL, "Kriging",NULL,
                             NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL,
                             NULL, NULL,NULL, NULL, type_krig, NULL, NULL, NULL)
    if(!is.null(checkinput$error)) stop(checkinput$error)
    corrmodel <- CheckCorrModel(covmatrix$corrmodel)
    distance <- CheckDistance(covmatrix$distance)
    if(!is.null(covmatrix$tapmod)) tapmod <-CheckCorrModel(covmatrix$tapmod)
    else    tapmod <-NULL
    corrparam <- unlist(covmatrix$param[covmatrix$namescorr])# selecting the correlation parametrs
    loc <- matrix(loc,ncol=2)
    numloc <- nrow(loc)
    tloc <- length(time);if(!tloc) tloc <- 1
    locx <- loc[,1];locy <- loc[,2]
    fix <- covmatrix$numcoord*covmatrix$numtime*numloc
    pred <- NULL
    varpred <- NULL
    k <- 0
    if(type=="Standard") {
    ## Computing correlation between the locations to predict and the locations observed
    cc=.C('Corr_c',corri=double(covmatrix$numcoord*covmatrix$numtime*numloc*tloc), as.double(covmatrix$coordx),as.double(covmatrix$coordy),
    as.double(covmatrix$coordt),as.integer(corrmodel),as.integer(covmatrix$grid),as.double(locx),as.double(locy),as.integer(covmatrix$numcoord),
    as.integer(numloc),as.integer(tloc),as.integer(covmatrix$numtime),as.double(corrparam),as.integer(covmatrix$spacetime),
    as.double(time),as.integer(distance),PACKAGE='CompRandFld',DUP=TRUE,NAOK=TRUE)
    corri<-cc$corri
    ## cholesky decomposition
    cholvarcov <- try(chol(covmatrix$covmatrix),silent=TRUE)
    if(!is.matrix(cholvarcov)) return("Covariance matrix is not positive definite")
    invcov<-chol2inv(cholvarcov)
    for(i in 1:tloc)  {
    ck <- window(corri,start=1+k,end=fix+k)
    cc <- t(matrix(ck*covmatrix$param["sill"],nrow=covmatrix$numcoord*covmatrix$numtime,ncol=numloc))
    krig_weights <- cc%*%invcov
    totvar <- as.numeric(covmatrix$param["sill"]+covmatrix$param["nugget"])

  if(type_krig=='Simple'||type_krig=='simple')  {
     pp <- param$mean + krig_weights %*% (c(data)-param$mean)   ## simple kriging predictor
     vv <- diag(as.matrix(totvar - krig_weights%*%t(cc)))       ## simple kriging predictor variance
     }
  if(type_krig=='Ordinary'||type_krig=='ordinary')  {
     suminv <- sum(invcov)
     mm <- sum(invcov%*%c(data))/suminv           # (1^t%*% C^-1 %*%data)/(1^t%*% C^-1 %*%1)   mean estimation
     pp <- mm +krig_weights %*% (c(data)-mm)      ## ordinary  kriging predictor
     vv <-  diag(as.matrix( totvar - krig_weights %*% t(cc) + (1-rowSums(krig_weights))^2/suminv))  ## ordinary kriging predictor variance
     }
     pred=rbind(pred,c(pp))
     varpred=rbind(varpred,c(vv))
     k=k+fix
    }}

    if(type=="Tapering")  {
    
    ct=.C('Corr_c_tap',corri=double(covmatrix$numcoord*covmatrix$numtime*numloc*tloc), corri_tap=double(covmatrix$numcoord*covmatrix$numtime*numloc*tloc), as.double(covmatrix$coordx),as.double(covmatrix$coordy),as.double(covmatrix$coordt),
    as.integer(corrmodel),as.integer(tapmod),as.integer(covmatrix$grid),as.double(locx),as.double(locy),
    as.integer(covmatrix$numcoord),as.integer(numloc),as.integer(tloc),as.integer(covmatrix$numtime),as.double(corrparam),
    as.integer(covmatrix$spacetime),as.double(time),as.integer(distance),PACKAGE='CompRandFld',DUP=TRUE,NAOK=TRUE)
    corri_tap<-ct$corri_tap
    corri<-ct$corri
    cholvarcov <- try(spam::chol.spam(covmatrix$covmatrix),silent=TRUE)
    if(class(cholvarcov)=="try-error") return("Covariance matrix is not positive definite")
    invcov<- spam::solve.spam(cholvarcov)
    for(i in 1:tloc){
    ck <- window(corri,start=1+k,end=fix+k)
    ck_tap <- window(corri_tap,start=1+k,end=fix+k)
    cc <- t(matrix(ck*covmatrix$param["sill"],nrow=covmatrix$numcoord*covmatrix$numtime,ncol=numloc))
    cc_tap <- t(matrix(ck_tap,nrow=covmatrix$numcoord*covmatrix$numtime,ncol=numloc))
    krig_weights_tap <- cc_tap%*%as.matrix(invcov)
    totvar <- as.numeric(covmatrix$param["sill"]+covmatrix$param["nugget"])
    if(type_krig=='Simple'||type_krig=='simple')  {
     pp <- param$mean + krig_weights_tap %*% (c(data)-param$mean)                                                               ## simple tapering kriging predictor
     vv <- diag(as.matrix(totvar - 2*krig_weights_tap%*%t(cc))+krig_weights_tap%*%covmatrix$covmatrix%*%t(krig_weights_tap) )   ## simple tapering kriging predictor variance
     }
     # if(type_krig=='Ordinary'||type_krig=='ordinary')  {
     #suminv <- sum(invcov)
     #mm  <- sum(invcov%*%c(data))/suminv                #   (1^t%*% C^-1 %*%data)/(1^t%*% C^-1) %*%1 mean estimation
     #pp <- mm +krig_weights_tap %*% (c(data)-mm)                                     ## ordinary  kriging predictor
     #vv <-  diag(as.matrix( totvar - krig_weights %*% t(cc) + (1-rowSums(krig_weights))^2/suminv))  ## ordinary kriging predictor variance
     #}
     pred <- rbind(pred,c(pp))
     varpred <- rbind(varpred,c(vv))
     k <- k+fix
      }
    }
     # Delete the global variables:
    .C('DeleteGlobalVar', PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)

   if(tloc==1)  {c(pred);c(varpred)}
    # Return the objects list:
    Kg <- list(    coordx = covmatrix$coordx,
                   coordt = covmatrix$coordt,
                   covmatrix=covmatrix$covmatrix,
                   corrmodel = corrmodel,
                   data=data,
                   distance = distance,
                   grid=covmatrix$grid,
                   loc=loc,
                   nozero=covmatrix$nozero,
                   numcoord = covmatrix$numcoord,
                   numloc= numloc,
                   numtime = covmatrix$numtime,
                   numt = tloc,
                   maxdist=maxdist,
                   maxtime=maxtime,
                   param = covmatrix$param,
                   pred=pred,
                   spacetime = covmatrix$spacetime,
                   tapmod=tapmod,
                   time=time,
                   type=type,
                   type_krig=type_krig,
                   varpred=varpred)
    structure(c(Kg, call = call), class = c("Kg"))
}

