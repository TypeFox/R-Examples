    ####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Emails: simone.padoan@unibocconi.it,
### moreno.bevilacqua@uv.cl
### Institutions: Department of Decision Sciences,
### University Bocconi of Milan and
### Departamento de Estadistica
### Universidad de Valparaiso
### File name: Covmatrix.r
### Description:
### This file contains a set of procedures
### for computing a covariance (tapered) matrix for a given
### space(time) covariance model.
### Last change: 28/03/2013.
####################################################


### Procedures are in alphabetical order.
Covmatrix <- function(coordx, coordy=NULL, coordt=NULL, corrmodel, distance="Eucl", grid=FALSE,
                      iskrig=FALSE, maxdist=NULL, maxtime=NULL, param, taper=NULL, tapsep=NULL, type="Standard")

{
    Cmatrix <- function(corrmodel, dime, nuisance, numpairs, numpairstot, paramcorr, setup, spacetime, type)
    {
       if(type=="Tapering")  {
       
        if(!spacetime) fname <- "CorrelationMat_tap" else fname <- "CorrelationMat_st_tap"
        p=.C(fname, corr=double(numpairs), as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),
           PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
        corr<-p$corr
        vcov <- corr*nuisance['sill']
        vcov[vcov==(nuisance['sill'])] <- nuisance['sill']+nuisance['nugget']
        varcov <- new("spam",entries=vcov*setup$taps,colindices=setup$ja,
                         rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
         }
        else{
        corr <- double(numpairstot)
        if(!spacetime) fname <- "CorrelationMat" else fname <- "CorrelationMat_st"
        p=.C(fname, corr=double(numpairs), as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),
           PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
        corr<-p$corr
        # Builds the covariance matrix:
        varcov <- (nuisance['nugget'] + nuisance['sill']) * diag(dime)
        corr <- corr * nuisance['sill']
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr
        }
      return(varcov)
    }
    ### END
    # Check the user input
    if(is.null(param$mean)) param$mean<-0
    if(is.null(param$nugget)) param$nugget<-0


    checkinput <- CheckInput(coordx, coordy, coordt, corrmodel, NULL, distance, "Simulation",
                             NULL, grid, NULL, NULL, maxdist, maxtime,  model='Gaussian', NULL, NULL,
                             param, 1, NULL, taper, tapsep, NULL, type, NULL, NULL, NULL)

    if(!is.null(checkinput$error)) stop(checkinput$error)
    # Initialising the parameters:
    initparam <- InitParam(coordx, coordy, coordt, corrmodel, NULL, distance, "Simulation",
                           NULL, grid, NULL, NULL, maxdist, maxtime, 'Gaussian', NULL,
                           param, NULL, NULL, 1, NULL, taper, tapsep, NULL, type, type,
                           NULL, NULL, FALSE, NULL, NULL)
    dime=initparam$numcoord*initparam$numtime
    numpairstot=dime*(dime-1)*0.5
    if(!is.null(initparam$error)) stop(initparam$error)
    setup<-initparam$setup
    if(initparam$type=="Tapering"){
    if(initparam$spacetime) fname= "CorrelationMat_st_tap"
    if(!initparam$spacetime) fname= "CorrelationMat_tap"
    corr <- double(initparam$numpairs)
    tapmod <- setup$tapmodel
    tp=.C(fname, tapcorr=double(initparam$numpairs),as.integer(tapmod),as.double(c(0,0,1)),
       as.double(1),PACKAGE='CompRandFld',DUP=TRUE,NAOK=TRUE)
    setup$taps<-tp$tapcorr
    }
    covmatrix<- Cmatrix(initparam$corrmodel,dime,initparam$param[initparam$namesnuis],initparam$numpairs,numpairstot,
                        initparam$param[initparam$namescorr],setup,initparam$spacetime,initparam$type)
    initparam$param=initparam$param[names(initparam$param)!='mean']

    # Delete the global variables:
    if(!iskrig)  .C('DeleteGlobalVar', PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)

    # Return the objects list:
    CovMat <- list(coordx = initparam$coordx,
                   coordy = initparam$coordy,
                   coordt = initparam$coordt,
                   covmatrix=covmatrix,
                   corrmodel = corrmodel,
                   distance = distance,
                   grid=   grid,
                   nozero=initparam$setup$nozero,
                   maxdist = maxdist,
                   maxtime = maxtime,
                   namescorr = initparam$namescorr,
                   numcoord = initparam$numcoord,
                   numtime = initparam$numtime,
                   param = initparam$param,
                   spacetime = initparam$spacetime,
                   tapmod=taper)
    structure(c(CovMat, call = call), class = c("CovMat"))
}

