####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Emails: simone.padoan@unibocconi.it,
### moreno.bevilacqua@uv.cl
### Institutions: Department of Decision Sciences,
### University Bocconi of Milan and
### Departamento de Estad?stica
### Universidad de Valparaiso
### File name: Simulation.r
### Description:
### This file contains a set of procedures
### for the simulation of Gaussian random fields and
### related functions.
### Last change: 28/03/2013.
####################################################


### Procedures are in alphabetical order.

# Simulate spatial and spatio-temporal random felds:
RFsim <- function(coordx, coordy=NULL, coordt=NULL, corrmodel, distance="Eucl", grid=FALSE,
                  model='Gaussian', numblock=NULL, param, replicates=1, threshold=NULL)
{
    ### START -- exact direct simulation based on the Cholesky decomposition
    CholRFsim<- function(corrmodel, grid,  model, nuisance, numcoord, numcoordx,
                         numcoordy, numtime, paramcorr, replicates, spacetime,
                         threshold)
    {
        # Compute the correlation function:
        if(!spacetime) fname <- "CorrelationMat" else fname <- "CorrelationMat_st"
        CV=.C(fname, cr=double(.5*(numcoord*numtime)*(numcoord*numtime-1)), as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),
           PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
        # Builds the covariance matrix:
        corr<-CV$cr
        varcov <- (nuisance['nugget'] + nuisance['sill']) * diag(numcoord*numtime)
        corr <- corr * nuisance['sill']
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr
        # Simulation based on the Cholesky decomposition:
        if(grid){
            sim <- array(double(numcoord*numtime*replicates), c(numcoordx, numcoordy, numtime, replicates))
            for(i in 1:replicates) {
                 cholvarcov <- try(chol(varcov),silent=TRUE)
                 if(!is.matrix(cholvarcov)) return("Covariance matrix is not positive definite")
                 da=nuisance['mean']+t(cholvarcov)%*%rnorm(numcoord*numtime)
                 l=0
                 for(k in 1:(numtime)) {
                 sim[,,l+1,i]=da[seq(1+l,numcoord*numtime*replicates+l,numtime)]
                 l=l+1 }}

            if(replicates==1)
                if(!spacetime) sim <- array(sim, c(numcoordx, numcoordy))
                else sim <- array(sim, c(numcoordx, numcoordy, numtime))
            else if(!spacetime) sim <- array(sim, c(numcoordx, numcoordy, replicates))
            }
        else{
            sim <- array(double(numcoord*numtime*replicates), c(numtime, numcoord, replicates))
            for(i in 1:replicates){
                cholvarcov <- try(chol(varcov),silent=TRUE)
                if(!is.matrix(cholvarcov)) return("Covariance matrix is not positive definite")
                sim[,,i] <- nuisance['mean']+t(chol(varcov))%*%rnorm(numcoord*numtime)}
            if(replicates==1)
                if(!spacetime) sim <- c(sim)
                else sim <- matrix(sim, nrow=numtime, ncol=numcoord)
            else if(!spacetime) sim <- matrix(sim, nrow=replicates, ncol=numcoord, byrow=TRUE)}
        # If a Binary-Gaussian field is required, the Gaussian field is transformed into a Binary version:
        if(model==2){
            simdim <- dim(sim)
            sim <- as.numeric(sim>threshold)
            dim(sim) <- simdim}
      return(sim)
    }
    ### END -- exact direct simulation based on the Cholesky decomposition
    # Check the user input
    checkinput <- CheckInput(coordx, coordy, coordt, corrmodel, NULL, distance, "Simulation",
                             NULL, grid, NULL, NULL, NULL, NULL, model, numblock, NULL, param,
                             replicates, NULL, NULL, NULL, threshold, "Standard", NULL, NULL, NULL)

    if(!is.null(checkinput$error)) stop(checkinput$error)
    # Initialising the simulation parameters:

    initparam <- InitParam(coordx, coordy, coordt, corrmodel, NULL, distance, "Simulation",
                           NULL, grid, NULL, NULL, NULL, NULL, model, numblock,
                           param, NULL, NULL, replicates, NULL, NULL, NULL, threshold, NULL, "Standard",
                           NULL, NULL, FALSE, NULL, NULL)
    if(!is.null(initparam$error)) stop(initparam$error)
    ### Simulation of Gaussian or Binary-Gaussian random fields:
    if(model %in% c("Gaussian","BinaryGauss")){
        # Direct method based on Cholesky decomposition:
        sim <- CholRFsim(initparam$corrmodel,grid,initparam$model,initparam$param[initparam$namesnuis],
                         initparam$numcoord,initparam$numcoordx,initparam$numcoordy,initparam$numtime,
                         initparam$param[initparam$namescorr],initparam$numrep,initparam$spacetime,
                         initparam$threshold)}

    if(model=="ExtGauss"){
        # Simulate directly Extremal Gaussian random fields from Random Felds:
        sim <- RandomFields::MaxStableRF(x=initparam$coordx, y=initparam$coordy, model=corrmodel, grid=grid,
                                         maxstable="extr",param=initparam$param[initparam$namessim],
                                         n=replicates,pch='')
        # Formatting of output:
        if(grid){
            if(replicates==1) sim <- array(sim, c(initparam$numcoordx, initparam$numcoordy))
            else sim <- array(sim, c(initparam$numcoordx, initparam$numcoordy, replicates))}
        else{
            if(replicates==1) sim <- c(sim)
            else sim <- matrix(sim, nrow=replicates,ncol=initparam$numcoord, byrow=TRUE)}}
    if(model=="ExtT" || model=="BrowResn"){
        # Set simulation variables:
        sim <- NULL
        RandomFields::DeleteAllRegisters()
        RandomFields::RFparameters(Storing=TRUE, PrintLevel=1)
        for(i in 1:replicates){
            maxima <- double(initparam$numcoord)
            # Simulate underlying Gaussian random fields (to get then Student-t fields):
            onesim <- RandomFields::GaussRF(x=initparam$coordx, y=initparam$coordy,model=corrmodel,
                                            param=initparam$param[initparam$namessim],grid=grid,
                                            n=initparam$numblock,pch='')
            # Compute compotentwise maxima:
              mx=.C("ComputeMaxima",as.double(initparam$param["df"]),mx=maxima,
                 as.integer(initparam$model),as.integer(initparam$numblock),
                 as.integer(initparam$numcoord),onesim,
                 PACKAGE='CompRandFld',DUP=TRUE,NAOK=TRUE)
            # Update:
            sim <- c(sim,mx$mx)}
        # Formatting of output:
        if(grid){
            if(replicates==1) sim <- array(sim, c(initparam$numcoordx, initparam$numcoordy))
            else sim <- array(sim, c(initparam$numcoordx, initparam$numcoordy, replicates))}
        else{
            if(replicates==1) sim <- c(sim)
            else sim <- matrix(sim, nrow=replicates,ncol=initparam$numcoord, byrow=TRUE)}}
    # Delete the global variables:
     .C('DeleteGlobalVar', PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
    # Return the objects list:
    RFsim <- list(coordx = initparam$coordx,
                  coordy = initparam$coordy,
                  coordt = initparam$coordt,
                  corrmodel = corrmodel,
                  data = sim,
                  distance = distance,
                  grid = grid,
                  model = model,
                  numcoord = initparam$numcoord,
                  numtime = initparam$numtime,
                  param = initparam$param,
                  randseed=.Random.seed,
                  replicates = initparam$replicates,
                  spacetime = initparam$spacetime,
                  threshold = initparam$threshold)

    structure(c(RFsim, call = call), class = c("RFsim"))
}

