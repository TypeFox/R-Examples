##' survspat function
##'
##' A function to run a Bayesian analysis on censored spatial survial data assuming a proportional hazards model using an adaptive Metropolis-adjusted
##' Langevin algorithm.
##'
##' @param formula the model formula in a format compatible with the function flexsurvreg from the flexsurv package 
##' @param data a SpatialPointsDataFrame object containing the survival data as one of the columns
##' @param dist choice of distribution function for baseline hazard. Current options are: exponentialHaz, weibullHaz, gompertzHaz, makehamHaz, tpowHaz
##' @param cov.model an object of class covmodel, see ?covmodel ?ExponentialCovFct or ?SpikedExponentialCovFct
##' @param mcmc.control mcmc control parameters, see ?mcmcpars
##' @param priors an object of class Priors, see ?mcmcPriors
##' @param shape when data is a data.frame, this can be a SpatialPolygonsDataFrame, or a SpatialPointsDataFrame, used to model spatial variation at the small region level. The regions are the polygons, or they represent the (possibly weighted) centroids of the polygons.
##' @param ids named list entry shpid character string giving name of variable in shape to be matched to variable dataid in data. dataid is the second entry of the named list.
##' @param control additional control parameters, see ?inference.control
##' @return an object inheriting class 'mcmcspatsurv' for which there exist methods for printing, summarising and making inference from.
##' @seealso \link{tpowHaz}, \link{exponentialHaz}, \link{gompertzHaz}, \link{makehamHaz}, \link{weibullHaz}, 
##' \link{covmodel}, link{ExponentialCovFct}, \code{SpikedExponentialCovFct},
##' \link{mcmcpars}, \link{mcmcPriors}, \link{inference.control} 
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor. Auxiliary Variable Markov Chain Monte Carlo for Spatial Survival and Geostatistical Models. Benjamin M. Taylor. Submitted. \url{http://arxiv.org/abs/1501.01665}
##' }
##' @export


survspat <- function(   formula,
                        data,
                        dist,
                        cov.model,
                        mcmc.control,
                        priors,
                        shape=NULL,
                        ids=list(shpid=NULL,dataid=NULL),
                        control=inference.control(gridded=FALSE)){
                        
    formula <- as.formula(formula)                    
    
    # initial checks
    if(!(inherits(data,"SpatialPointsDataFrame")|inherits(data,"data.frame"))){
        stop("'data' must be of class 'SpatialPointsDataFrame' or 'data.frame'.")
    }   
    if(inherits(data,"data.frame")&is.null(shape)){
        stop("If inherits(data,'data.frame')==TRUE, shape cannot be NULL")
    }   

    latentmode <- "points"
    if(inherits(data,"data.frame")){
        latentmode <- "polygons"
    }
    if(inherits(cov.model,"SPDEmodel")){
        latentmode <- "SPDE"
    }

    if(inherits(data,"SpatialPointsDataFrame")&!is.null(shape)){
        if(latentmode!="SPDE"){
            warning("Non NULL shape, ignoring shape.",immediate.=TRUE)
            Sys.sleep(2)
        }
    }

    if(latentmode=="polygons" & control$gridded){
        stop("Cannot have control$gridded==TRUE and !is.null(shape)==TRUE at the same time")
    }

    if(latentmode=="SPDE" & control$gridded){
        stop("Cannot use SPDE mode and have !is.null(shape)==TRUE at the same time. Please provide a polygon object on which to produce spatial predictions.")
    }

    responsename <- as.character(formula[[2]])
    if(latentmode=="points"){
        survivaldata <- data@data[[responsename]]
    }
    else{
        survivaldata <- data[[responsename]]
    }
    checkSurvivalData(survivaldata) 

    if(latentmode=="SPDE"){
        if(proj4string(data)!=proj4string(shape)){
            stop("'shape' and 'data' must have the same proj4string.")
        } 
        if(is.null(control$cellwidth)){
            stop("Must specify 'cellwidth' in inference.control in SPDE mode.")
        }   
    }
    
    
    # if(!inherits(data,"SpatialPointsDataFrame")){
    #     stop("data must be an object of class SpatialPointsDataFrame")
    # }                     

    # okay, start the MCMC!

    # start timing, maybe
    if(!control$timeonlyMCMC){
        start <- Sys.time()
    }

    if(latentmode=="points" | latentmode=="SPDE"){
        coords <- coordinates(data)
    }
    else{
        coords <- coordinates(shape)
    }

    control$dist <- dist
    
    funtxt <- ""    
    if(control$gridded){
        funtxt <- "_gridded"
    }
    if(latentmode=="polygons"){
        funtxt <- "_polygonal"
    }
    if(latentmode=="SPDE"){
        funtxt <- "_SPDE"
    }
    
    gridobj <- NULL
   
    if(control$gridded){
        gridobj <- FFTgrid(spatialdata=data,cellwidth=control$cellwidth,ext=control$ext)
    	del1 <- gridobj$del1
    	del2 <- gridobj$del2
    	Mext <- gridobj$Mext
    	Next <- gridobj$Next
    	mcens <- gridobj$mcens
    	ncens <- gridobj$ncens    	
    	## COMPUTE GRID DISTANCES ##
    	x <- gridobj$mcens
        y <- gridobj$ncens    
        xidx <- rep(1:Mext,Next)
        yidx <- rep(1:Next,each=Mext)
        dxidx <- pmin(abs(xidx-xidx[1]),Mext-abs(xidx-xidx[1]))
        dyidx <- pmin(abs(yidx-yidx[1]),Next-abs(yidx-yidx[1]))
        u <- sqrt(((x[2]-x[1])*dxidx)^2+((y[2]-y[1])*dyidx)^2)
        
        spix <- grid2spix(xgrid=mcens,ygrid=ncens,proj4string=CRS(proj4string(data)))
        
        control$fftgrid <- gridobj
        control$idx <- over(data,geometry(spix))
        control$Mext <- Mext
        control$Next <- Next
        control$uqidx <- unique(control$idx)
        
        cat("Output grid size: ",Mext/control$ext," x ",Next/control$ext,"\n")
    }
    else{
        if(latentmode=="polygons"){
            u <- as.vector(as.matrix(dist(coords)))
            control$idx <- match(data[,ids$dataid],shape@data[,ids$shpid])
            control$n <- nrow(shape)
            control$uqidx <- unique(control$idx)
        }
        else if(latentmode=="SPDE"){
            
            matobj <- setupPrecMatStruct(shape=shape,cellwidth=control$cellwidth,no=cov.model$order)
            control$precmat <- matobj$f
            control$grid <- matobj$grid

            # require(INLA)
            # now reorder the indices for faster computation
            # tempprec <- control$precmat(SPDEprec(0.1,cov.model$order))            
            # reord <- inla.qreordering(tempprec)
            # env <- environment(control$precmat)
            # image.plot(as.matrix(tempprec))
            # image.plot(as.matrix(tempprec)[reord$reordering,reord$reordering])
            # rord <- function(oldidx,newidx){
            #     ans <- rep(NA,length(oldidx))
            #     for(i in 1:length(newidx)){
            #         ans[oldidx==i] <- newidx[i]
            #     }
            #     return(ans)
            # }
            # env$index[,1] <- rord(env$index[,1],reord$reordering)
            # env$index[,2] <- rord(env$index[,2],reord$reordering)
            # shuff <- function(x){
            #     if(x[2]>x[1]){
            #         return(c(x[2:1],x[3]))
            #     }
            #     else{
            #         return(x)
            #     }    
            # }
            # env$index <- t(apply(env$index,1,shuff))
            # env$grid <- env$grid[reord$reordering]
            # tempprec1 <- control$precmat(SPDEprec(0.1,cov.model$order)) 
            #browser()

            control$idx <- over(data,geometry(control$grid))
            control$n <- length(control$grid)
            control$uqidx <- unique(control$idx)
        }
        else{
            u <- as.vector(as.matrix(dist(coords)))
        }
    }

    DATA <- data    
    
    if(latentmode=="points" | latentmode=="SPDE"){
        data <- data@data  
    }
                        
    ##########
    # This chunk of code borrowed from flexsurvreg    
    ##########  
                      
    call <- match.call()
    indx <- match(c("formula", "data"), names(call), nomatch = 0)
    if (indx[1] == 0){ 
        stop("A \"formula\" argument is required")
    }
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    m <- eval(temp, parent.frame())

    Terms <- attr(m, "terms")
    X <- model.matrix(Terms, m)
    
    ##########
    # End of borrowed code    
    ##########
 
    X <- X[, -1, drop = FALSE]                                   
                        
    mcmcloop <- mcmcLoop(N=mcmc.control$nits,burnin=mcmc.control$burn,thin=mcmc.control$thin,progressor=mcmcProgressTextBar)                            
              
    info <- distinfo(dist)()
    
    control$omegatrans <- info$trans
    control$omegaitrans <- info$itrans
    control$omegajacobian <- info$jacobian # used in computing the derivative of the log posterior with respect to the transformed omega (since it is easier to compute with respect to omega) 
    control$omegahessian <- info$hessian


    #######

    control$censoringtype <- attr(survivaldata,"type")

    if(control$censoringtype=="left" | control$censoringtype=="right"){
        control$censored <- survivaldata[,"status"]==0
        control$notcensored <- !control$censored

        control$Ctest <- any(control$censored)
        control$Utest <- any(control$notcensored)

        control$idxi <- list()
        control$idxicensored <- list()
        control$idxinotcensored <- list()
        lapply(control$uqidx,function(i){control$idxi[[i]] <<- which(control$idx==i)})
        lapply(control$uqidx,function(i){control$idxicensored[[i]] <<- control$censored & control$idx==i})
        lapply(control$uqidx,function(i){control$idxinotcensored[[i]] <<- control$notcensored & control$idx==i})
        
        if(control$Ctest){
            control$idxicensored <- lapply(control$idxicensored,function(x){try(which(x),silent=TRUE)})
        }
        if(control$Utest){
            control$sumidxinotcensored <- lapply(control$idxinotcensored,function(x){try(sum(x),silent=TRUE)})
            control$idxinotcensored <- lapply(control$idxinotcensored,function(x){try(which(x),silent=TRUE)})
        }
    }
    else{
        control$rightcensored <- survivaldata[,"status"] == 0
        control$notcensored <- survivaldata[,"status"] == 1
        control$leftcensored <- survivaldata[,"status"] == 2
        control$intervalcensored <- survivaldata[,"status"] == 3

        control$Rtest <- any(control$rightcensored)        
        control$Utest <- any(control$notcensored) 
        control$Ltest <- any(control$leftcensored)
        control$Itest <- any(control$intervalcensored)

        control$idxirightcensored <- list()
        control$idxinotcensored <- list()
        control$idxileftcensored <- list()
        control$idxiintervalcensored <- list()

        lapply(control$uqidx,function(i){control$idxirightcensored[[i]] <<- control$rightcensored & control$idx==i})
        lapply(control$uqidx,function(i){control$idxinotcensored[[i]] <<- control$notcensored & control$idx==i})
        lapply(control$uqidx,function(i){control$idxileftcensored[[i]] <<- control$leftcensored & control$idx==i})
        lapply(control$uqidx,function(i){control$idxiintervalcensored[[i]] <<- control$intervalcensored & control$idx==i})

        if(control$Rtest){
            control$idxirightcensored <- lapply(control$idxirightcensored,function(x){try(which(x),silent=TRUE)})
        }
        if(control$Utest){
            control$idxinotcensored <- lapply(control$idxinotcensored,function(x){try(which(x),silent=TRUE)})
        }
        if(control$Ltest){
            control$idxileftcensored <- lapply(control$idxileftcensored,function(x){try(which(x),silent=TRUE)})
        }
        if(control$Itest){
            control$idxiintervalcensored <- lapply(control$idxiintervalcensored,function(x){try(which(x),silent=TRUE)})
        }
    }

    #######

    
    cat("\n","Getting initial estimates of model parameters using maximum likelihood on non-spatial version of the model","\n")
    mlmod <- maxlikparamPHsurv(surv=survivaldata,X=X,control=control)
    estim <- mlmod$par
    print(mlmod)
    cat("Done.\n")
    #browser() 
    
    cat("Calibrating MCMC algorithm and finding initial values ...\n")
    
    betahat <- estim[1:ncol(X)]
    omegahat <- estim[(ncol(X)+1):length(estim)]
    
    if(latentmode!="SPDE"){  
        control$sigmaidx <- match("sigma",cov.model$parnames)
        if(is.na(control$sigmaidx)){
            stop("At least one of the parameters must be the variance of Y, it should be named sigma")
        }
    }
    
    Yhat <- estimateY(  X=X,
                        betahat=betahat,
                        omegahat=omegahat,
                        surv=survivaldata,
                        control=control)
                        
    calibrate <- get(paste("proposalVariance",funtxt,sep=""))    
       
    other <- calibrate( X=X,
                        surv=survivaldata,
                        betahat=betahat,
                        omegahat=omegahat,
                        Yhat=Yhat,
                        priors=priors,
                        cov.model=cov.model,
                        u=u,
                        control=control) 
    
    #gammahat <- other$gammahat
    etahat <- other$etahat                                                                        
    SIGMA <- other$sigma 

    beta <- betahat
    omega <- omegahat
    eta <- etahat 
 
    gamma <- rep(0,nrow(X))
    if(control$gridded){
        gamma <- matrix(0,control$Mext,control$Next)
    }
    if(latentmode=="polygons" | latentmode=="SPDE"){
        gamma <- rep(0,control$n)
    }
        
    lenbeta <- length(beta)
    lenomega <- length(omega)
    leneta <- length(eta)
    lengamma <- length(gamma)
    
    
    
    npars <- lenbeta + lenomega + leneta + lengamma
    
    SIGMA[1:(lenbeta+lenomega),1:(lenbeta+lenomega)] <- (1.65^2/((lenbeta+lenomega)^(1/3)))*SIGMA[1:(lenbeta+lenomega),1:(lenbeta+lenomega)]
    SIGMA[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)] <- 0.4*(2.38^2/leneta)* SIGMA[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta),(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)]
    SIGMA[(lenbeta+lenomega+leneta+1):(lenbeta+lenomega+leneta+lengamma),(lenbeta+lenomega+leneta+1):(lenbeta+lenomega+leneta+lengamma)] <- (1.65^2/(lengamma^(1/3)))*SIGMA[(lenbeta+lenomega+leneta+1):(lenbeta+lenomega+leneta+lengamma),(lenbeta+lenomega+leneta+1):(lenbeta+lenomega+leneta+lengamma)]   
    
    if(control$gridded){
        matidx <- matrix(0,control$Mext,control$Next)
        matidx[1:(control$Mext/control$ext),1:(control$Next/control$ext)] <- 1
        matidx <- as.logical(matidx) # used to select which Y's to save 
    }
  
    
    diagidx <- 1:npars
    diagidx <- matrix(diagidx,nrow=npars,ncol=2)
    SIGMApars <- as.matrix(SIGMA[1:(lenbeta+lenomega+leneta),1:(lenbeta+lenomega+leneta)])
    
    SIGMAparsINV <- solve(SIGMApars)
    cholSIGMApars <- t(chol(SIGMApars))  
  
    SIGMAgamma <- SIGMA[diagidx][(lenbeta+lenomega+leneta+1):npars]
    SIGMAgammaINV <- 1/SIGMAgamma
    cholSIGMAgamma <- sqrt(SIGMAgamma)
    
     
    
    cat("Running MCMC ...\n")
    
    h <- 1
    
    
    
    LOGPOST <- get(paste("logPosterior",funtxt,sep=""))
    
    oldlogpost <- LOGPOST(  surv=survivaldata,
                            X=X,
                            beta=beta,
                            omega=omega,
                            eta=eta,
                            gamma=gamma,
                            priors=priors,
                            cov.model=cov.model,
                            u=u,
                            control=control,
                            gradient=TRUE)
                            
                          
                                                        
    
    betasamp <- c()
    omegasamp <- c()
    etasamp <- c()
    Ysamp <- c()
    
    gamma <- c(gamma) # turn gamma into a vector 
    
    tarrec <- oldlogpost$logpost
    
    print(SIGMA[1:8,1:8])
    
    loglik <- c()
    
    gammamean <- 0
    count <- 1

    # start timing, maybe
    if(control$timeonlyMCMC){
        start <- Sys.time()
    }

    bad <-  c()
    
    
    while(nextStep(mcmcloop)){

        stuffpars <- c(beta,omega,eta)
        propmeanpars <- stuffpars + (h/2)*SIGMApars%*%oldlogpost$grad[1:(lenbeta+lenomega+leneta)]
        newstuffpars <- propmeanpars + sqrt(h)*cholSIGMApars%*%rnorm(lenbeta+lenomega+leneta)
        
 
        propmeangamma <- gamma + (h/2)*SIGMAgamma*oldlogpost$grad[(lenbeta+lenomega+leneta+1):npars]
        newstuffgamma <- propmeangamma + sqrt(h)*cholSIGMAgamma*rnorm(lengamma)
        ngam <- newstuffgamma
        if(control$gridded){
            ngam <- matrix(ngam,control$Mext,control$Next)
        }
                                       
        newlogpost <- LOGPOST(  surv=survivaldata,
                                X=X,
                                beta=newstuffpars[1:lenbeta],
                                omega=newstuffpars[(lenbeta+1):(lenbeta+lenomega)],
                                eta=newstuffpars[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)],
                                gamma=ngam,
                                priors=priors,
                                cov.model=cov.model,
                                u=u,
                                control=control,
                                gradient=TRUE)

        revmeanpars <- newstuffpars + (h/2)*SIGMApars%*%newlogpost$grad[1:(lenbeta+lenomega+leneta)]
        revmeangamma <- newstuffgamma + (h/2)*SIGMAgamma*newlogpost$grad[(lenbeta+lenomega+leneta+1):npars]       

        revdiffpars <- as.matrix(stuffpars-revmeanpars)
        forwdiffpars <- as.matrix(newstuffpars-propmeanpars)
        revdiffgamma <- as.matrix(gamma-revmeangamma)
        forwdiffgamma <- as.matrix(newstuffgamma-propmeangamma)


        logfrac <- newlogpost$logpost - oldlogpost$logpost - 
                            (0.5/h)*t(revdiffpars)%*%SIGMAparsINV%*%revdiffpars + 
                            (0.5/h)*t(forwdiffpars)%*%SIGMAparsINV%*%forwdiffpars -
                            (0.5/h)*sum(revdiffgamma*SIGMAgammaINV*revdiffgamma) + 
                            (0.5/h)*sum(forwdiffgamma*SIGMAgammaINV*forwdiffgamma)
        
        ac <- min(1,exp(as.numeric(logfrac)))
        if(is.na(ac) | is.nan(ac)){
            ac <- 0
            bad <- c(bad,iteration(mcmcloop))
            warning("An acceptance probability could not be calculated for this iteration, this is likely because the spatial decay parameter was too big for this choice of 'ext'. Either increase ext, or tighten prior on spatial decay parameter. At the end of the run check $bad to see which iterations this affected. Stop the run if this problem persists.",immediate.=TRUE)
        }
        
        if(ac>runif(1)){
            beta <- newstuffpars[1:lenbeta]
            omega <- newstuffpars[(lenbeta+1):(lenbeta+lenomega)]
            eta <- newstuffpars[(lenbeta+lenomega+1):(lenbeta+lenomega+leneta)]
            gamma <- newstuffgamma

            oldlogpost <- newlogpost
        }
        
        h <- exp(log(h) + (1/(iteration(mcmcloop)^0.5))*(ac-0.574))
        if(iteration(mcmcloop)%%100==0){
            cat("\n","h =",h,"\n")
        }
        
        
        if(is.retain(mcmcloop)){
            betasamp <- rbind(betasamp,as.vector(beta))
            omegasamp <- rbind(omegasamp,as.vector(omega))
            etasamp <- rbind(etasamp,as.vector(eta))
            if(control$gridded){
                Ysamp <- rbind(Ysamp,as.vector(oldlogpost$Y[matidx]))
            }
            else{
                Ysamp <- rbind(Ysamp,as.vector(oldlogpost$Y))
            }
            tarrec <- c(tarrec,oldlogpost$logpost)
            loglik <- c(loglik,oldlogpost$loglik)
            gammamean <- ((count-1)/count)*gammamean + (1/count)*gamma
            count <- count + 1
        }
    }
    
    colnames(Ysamp) <- paste("Y",1:ncol(Ysamp),sep="")
    
    # Compute DIC
    Dhat <- -2*LOGPOST(  surv=survivaldata,
                                X=X,
                                beta=colMeans(betasamp),
                                omega=colMeans(omegasamp),
                                eta=colMeans(etasamp),
                                gamma=gammamean,
                                priors=priors,
                                cov.model=cov.model,
                                u=u,
                                control=control,
                                gradient=TRUE)$loglik
    pD <- -2*mean(loglik) - Dhat
    DIC <- Dhat + 2*pD

    retlist <- list()
    retlist$formula <- formula
    retlist$data <- DATA
    retlist$dist <- dist
    retlist$cov.model <- cov.model
    retlist$mcmc.control <- mcmc.control
    retlist$priors <- priors
    retlist$control <- control
    
    retlist$terms <- Terms
    retlist$mlmod <- mlmod
    
    ####
    #   Back transform for output
    ####
    
    if(length(omega)>1){
        omegasamp <- t(apply(omegasamp,1,control$omegaitrans))
    }
    else{
        omegasamp <- t(t(apply(omegasamp,1,control$omegaitrans)))
    }
    colnames(omegasamp) <- info$parnames
    
    if(length(eta)>1){
        etasamp <- t(apply(etasamp,1,cov.model$itrans))
    }
    else{
        etasamp <- t(t(apply(etasamp,1,cov.model$itrans)))
    }
    colnames(etasamp) <- cov.model$parnames     
    
    ####    
    
    colnames(betasamp) <- colnames(model.matrix(formula,data))[-1] #attr(Terms,"term.labels")
    retlist$betasamp <- betasamp
    retlist$omegasamp <- omegasamp
    retlist$etasamp <- etasamp
    retlist$Ysamp <- Ysamp
    
    retlist$loglik <- loglik
    retlist$Dhat <- Dhat
    retlist$pD <- pD
    retlist$DIC <- DIC
    
    retlist$X <- X
    retlist$survivaldata <- survivaldata
    
    retlist$gridded <- control$gridded
    if(latentmode=="SPDE"){
        retlist$precmat <- control$precmat
        retlist$grid <- control$grid
        retlist$idx <- control$idx
        retlist$uqidx <- control$uqidx
    }
    else if(control$gridded){
        retlist$M <- Mext/control$ext
        retlist$N <- Next/control$ext
        retlist$xvals <- mcens[1:retlist$M]
        retlist$yvals <- ncens[1:retlist$N]
        lookup <- as.vector(matrix(1:(Mext*Next),Mext,Next)[1:retlist$M,1:retlist$N])
        ref <- as.vector(matrix(1:(retlist$M*retlist$N),retlist$M,retlist$N))
        retlist$cellidx <- sapply(control$idx,function(i){ref[lookup==i]})
    }
    
    retlist$tarrec <- tarrec
    retlist$lasth <- h
    retlist$bad <- bad
    
    retlist$omegatrans <- control$omegatrans
    retlist$omegaitrans <- control$omegaitrans
    
    retlist$control <- control
    retlist$censoringtype <- attr(survivaldata,"type")

    retlist$shape <- shape
    retlist$ids <- ids    
    retlist$latentmode <- latentmode
    
    retlist$time.taken <- Sys.time() - start
    
    cat("Time taken:",retlist$time.taken,"\n")
    
    class(retlist) <- c("list","mcmcspatsurv")

    return(retlist)
}


