##' getleneta function
##'
##' A function to compute the length of eta
##'
##' @param cov.model a covariance model 
##' @return the length of eta
##' @export

getleneta <- function(cov.model){
    if(inherits(cov.model,"fromRandomFieldsCovarianceFct")){
        leneta <- 2
    }
    else if(inherits(cov.model,"fromUserFunction") | inherits(cov.model,"SPDEmodel")){
        leneta <- cov.model$npar
    }
    else{
        stop("Unknkown covariance type")
    }
    return(leneta)
} 


##' getparranges function
##'
##' A function to extract parameter ranges for creating a grid on which to evaluate the log-posterior, used in calibrating the MCMC. This function
##' is not intended for general use.
##'
##' @param priors an object of class mcmcPriors
##' @param leneta the length of eta passed to the function 
##' @param mult defaults to 1.96 so the grid formed will be mean plus/minus 1.96 times the standard deviation
##' @return an appropriate range used to calibrate the MCMC: the mean of the prior for eta plus/minus 1.96 times the standard deviation
##' @export

getparranges <- function(priors,leneta,mult=1.96){
    rgs <- list()
    if(length(priors$etaprior$mean)==1 & length(priors$etaprior$sd)==1){
        for(i in 1:leneta){
            rgs[[i]] <- c(priors$etaprior$mean-mult*priors$etaprior$sd,priors$etaprior$mean+mult*priors$etaprior$sd)
        }
    }
    else if(length(priors$etaprior$mean)==1 & length(priors$etaprior$sd)>1){
        for(i in 1:leneta){
            rgs[[i]] <- c(priors$etaprior$mean-mult*priors$etaprior$sd[i],priors$etaprior$mean+mult*priors$etaprior$sd[i])
        }
    }
    else if(length(priors$etaprior$mean)>1 & length(priors$etaprior$sd)==1){
        for(i in 1:leneta){
            rgs[[i]] <- c(priors$etaprior$mean[i]-mult*priors$etaprior$sd,priors$etaprior$mean[i]+mult*priors$etaprior$sd)
        }
    }
    else{
        for(i in 1:leneta){
            rgs[[i]] <- c(priors$etaprior$mean[i]-mult*priors$etaprior$sd[i],priors$etaprior$mean[i]+mult*priors$etaprior$sd[i])
        }
    } 
    return(rgs)
}    

##' gencens function
##'
##' A function to generate observed times given a vector of true survival times and a vector of censoring times. Used in the simulation of
##' survival data. 
##'
##' @param survtimes a vector of survival times 
##' @param censtimes a vector of censoring times for left or right censored data, 2-column matrix of censoring times for interval censoring (number of rows equal to the number of observations). 
##' @param type the type of censoring to generate can be 'right' (default), 'left' or 'interval' 
##' @return an object of class 'Surv', the censoring indicator is equal to 1 if the
##' event is uncensored and 0 otherwise for right/left censored data, or for interval censored data, the indicator is 0 uncensored, 1 right censored, 
##' 2 left censored, or 3 interval censored.
##' @export


gencens <- function(survtimes,censtimes,type="right"){

    if(!(type=="right" | type=="left" | type=="interval")){
        stop("type must be one of 'right', 'left' or 'interval'")    
    }
    
    n <- length(survtimes)
    if(type=="right" | type=="left"){
        if(length(survtimes)!=length(censtimes)){
            stop("survtimes and censtimes should have the same length")
        }
    }
    else{        
        if(ncol(censtimes)!=2 | nrow(censtimes)!=n){
            stop("censtimes should be a 2-column matrix with number of rows equal to length(survtimes)")
        }
    } 
    
    if(type=="right"){
        obstimes <- survtimes
        cens <- rep(1,n)   
        for(i in 1:n){
            if(censtimes[i]<survtimes[i]){
                obstimes[i] <- censtimes[i]
                cens[i] <- 0
            }
        }
        return(Surv(time=obstimes,event=cens,type="right"))
    }
    if(type=="left"){
        obstimes <- survtimes
        cens <- rep(1,n)   
        for(i in 1:n){
            if(censtimes[i]>survtimes[i]){
                obstimes[i] <- censtimes[i]
                cens[i] <- 0
            }
        }
        return(Surv(time=obstimes,event=cens,type="left"))
    }
    if(type=="interval"){
        obstimes <- cbind(survtimes,NA)
        cens <- rep(1,n)
        uncensidx <- sample(1:n,floor(n/2))
        alter <- rep(TRUE,n)
        alter[uncensidx] <- FALSE #
        censtimes <- t(apply(censtimes,1,sort))   
        for(i in 1:n){
            if(censtimes[i,1]<survtimes[i] & censtimes[i,2]>survtimes[i] & alter[i]){
                obstimes[i,1] <- censtimes[i,1]
                obstimes[i,2] <- censtimes[i,2]
                cens[i] <- 3 # interval censored
            }
            else if(censtimes[i,2]<survtimes[i] & alter[i]){
                obstimes[i,1] <- censtimes[i,2]
                cens[i] <- 0 # right censored
            }
            else if(censtimes[i,1]>survtimes[i] & alter[i]){
                obstimes[i,1] <- censtimes[i,1]
                cens[i] <- 2 # left censored
            }
        }
        return(Surv(time=obstimes[,1],time2=obstimes[,2],event=cens,type="interval"))
    }
    
  
}






##' plotsurv function
##'
##' A function to produce a 2-D plot of right censored spatial survival data.
##'
##' @param spp A spatial points data frame
##' @param ss A Surv object (with right-censoring) 
##' @param maxcex maximum size of dots default is equavalent to setting cex equal to 1
##' @param transform optional transformation to apply to the data, a function, for example 'sqrt'
##' @param background a background object to plot default is null, which gives a blamk background note that if non-null, the parameters xlim and ylim will be derived from this object.
##' @param eventpt The type of point to illustrate events, default is 19 (see ?pch) 
##' @param eventcol the colour of events, default is black
##' @param censpt The type of point to illustrate events, default is "+" (see ?pch)
##' @param censcol the colour of censored observations, default is red
##' @param xlim optional x-limits of plot, default is to choose this automatically 
##' @param ylim optional y-limits of plot, default is to choose this automatically
##' @param xlab label for x-axis
##' @param ylab label for y-axis
##' @param add logical, whether to add the survival plot on top of an existing plot, default is FALSE, which produces a plot in a new device
##' @param ... other arguments to pass to plot
##' @return Plots the survival data non-censored observations appear as dots and censored observations as crosses. The size of the dot is proportional to the observed time.
##' @export

plotsurv <- function(spp,ss,maxcex=1,transform=identity,background=NULL,eventpt=19,eventcol="red",censpt="+",censcol="black",xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,add=FALSE,...){
    crds <- coordinates(spp)
    if(is.null(xlim)){
        if(is.null(background)){
            xlim <- range(crds[,1])
        }
    }
    if(is.null(ylim)){
        if(is.null(background)){
            ylim <- range(crds[,2])
        }
    }
    
    if(is.null(xlab)){
        xlab <- "x-coordinates"
    }
    if(is.null(ylab)){
        ylab <- "y-coordinates"
    }
    
    stimes <- ss[,"time"]
    stimes <- transform(stimes)

    event <- ss[,"status"] == 1 # event indicator
    cexx <- maxcex* stimes / max(stimes)    
    
    if(!add){
        plot(background,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
    }
    points(crds[event,],pch=eventpt,col=eventcol,cex=cexx[event])
    points(crds[!event,],pch=censpt,col=censcol,cex=cexx[!event])
    
}



##' inference.control function
##'
##' A function to control inferential settings. This function is used to set parameters for more advanced use of spatsurv. 
##'
##' @param gridded logical. Whether to perform compuation on a grid. Default is FALSE.
##' @param cellwidth the width of computational cells to use 
##' @param ext integer the number of times to extend the computational grid by in order to perform compuitation. The default is 2.
##' @param optimcontrol a list of optional arguments to be passed to optim for non-spatial models
##' @param hessian whether to return a numerical hessian. Set this to TRUE for non-spatial models.
##' equal to the number of parameters of the baseline hazard
##' @param plotcal logical, whether to produce plots of the MCMC calibration process, this is a technical option and should onyl be set 
##' to TRUE if poor mixing is evident (the printed h is low), then it is also useful to use a graphics device with multiple plotting windows. 
##' @param timeonlyMCMC logical, whether to only time the MCMC part of the algorithm, or whether to include in the reported running time the time taken to calibrate the method (default)
##' @return returns parameters to be used in the function survspat
##' @seealso \link{survspat}
##' @export

inference.control <- function(gridded=FALSE,cellwidth=NULL,ext=2,optimcontrol=NULL,hessian=FALSE,plotcal=FALSE,timeonlyMCMC=FALSE){
    ans <- list()
    ans$gridded <- gridded
    ans$cellwidth <- cellwidth 
    ans$ext <- ext 
    ans$optimcontrol <- optimcontrol
    ans$hessian <- hessian
    ans$plotcal <- plotcal
    ans$timeonlyMCMC <- timeonlyMCMC
    class(ans) <- c("inference.control","list")
    return(ans)
}





##' getsurvdata function
##'
##' A function to return the survival data from an object of class mcmcspatsurv. This function is not intended for general use.
##'
##' @param x an object of class mcmcspatsurv 
##' @return the survival data from an object of class mcmcspatsurv
##' @export

getsurvdata <- function(x){
    responsename <- as.character(x$formula[[2]])
    return(x$data[[responsename]])
}







##' checkSurvivalData function
##'
##' A function to check whether the survival data to be passed to survspat is in the correct format
##'
##' @param s an object of class Surv, from the survival package 
##' @return if there are any issues with data format, these are returned with the data an error message explaining any issues with the data
##' @export

checkSurvivalData <- function(s){
    if(class(s)!="Surv"){
        stop("Survival data must be of class 'Surv', see ?Surv")
    }
    
    if(attr(s,"type")=="right" | attr(s,"type")=="left" | attr(s,"type")=="interval"){
        if(any(as.matrix(s)<0,na.rm=TRUE)){
            stop("Survival data must not contain negative times, please change the offset of your data so that all times are non-negative")
        } 
        
        if(attr(s,"type")=="left" | attr(s,"type")=="interval"){
            cat("\n #####################################################\n # WARNING: CODE FOR LEFT AND INTERVAL CENSORED DATA #\n #          IS UNDER DEVELOPMENT                     #\n #####################################################\n\n")
            warning("*** CODE UNDER DEVELOPMENT ***",immediate.=TRUE)
        }
           
    }
    else{
        stop("Survival data must be of type 'left', 'right', or 'interval', see ?Surv")
    }
}




##' setupHazard function
##'
##' A function to set up the baseline hazard, cumulative hazard and derivative functions for use in evaluating the log posterior. 
##' This fucntion is not intended for general use.
##'
##' @param dist an object of class 'basehazardspec'
##' @param pars parameters with which to create the functions necessary to evaluate the log posterior
##' @param grad logical, whetether to create gradient functions for the baseline hazard and cumulative hazard
##' @param hess logical, whetether to create hessian functions for the baseline hazard and cumulative hazard
##' @return a list of functions used in evaluating the log posterior
##' @export

setupHazard <- function(dist,pars,grad=FALSE,hess=FALSE){
    funlist <- list()
    
    funlist$h <- basehazard(dist)(pars)
    if(grad){
        funlist$gradh <- gradbasehazard(dist)(pars)
    }
    if(hess){
        funlist$hessh <- hessbasehazard(dist)(pars)
    }
    
    funlist$H <- cumbasehazard(dist)(pars)
    if(grad){
        funlist$gradH <- gradcumbasehazard(dist)(pars)
    }
    if(hess){
        funlist$hessH <- hesscumbasehazard(dist)(pars)
    }
    
    return(funlist) 
}


##' invtransformweibull function
##'
##' A function to transform estimates of the (alpha, lambda) parameters of the weibull baseline hazard function, so they are commensurate 
##' with R's inbuilt density functions, (shape, scale).
##'
##' @param x a vector of paramters
##' @return the transformed parameters. For the weibull model, this transforms 'shape' 'scale' (see ?dweibull) to 'alpha' and 'lambda' for the MCMC
##' @export

invtransformweibull <- function(x){
    a <- x[1] # shape
    b <- x[2] # scale
    alpha <- a
    lambda <- (1/b)^a

    ans <- c(alpha=alpha,lambda=lambda) # note this is logged later for use in the MCMC
    
    return(ans)
}



##' transformweibull function
##'
##' A function to back-transform estimates of the parameters of the weibull baseline hazard function, so they are commensurate with R's inbuilt density functions. 
##' Transforms from (shape, scale) to (alpha, lambda)
##'
##' @param x a vector of paramters
##' @return the transformed parameters. For the weibull model, this is the back-transform from 'alpha' and 'lambda' to 'shape' 'scale' (see ?dweibull).
##' @export

transformweibull <- function(x){

    alpha <- x[1]
    lambda <- x[2]
    
    shape <- alpha
    scale <- exp((-1/alpha)*log(lambda))

    ans <- c(shape=shape,scale=scale)    
    
    return(ans)
}



##' spatsurvVignette function
##'
##' Display the introductory vignette for the spatsurv package. 
##'
##' @return displays the vignette by calling browseURL
##' @export

spatsurvVignette <- function(){
    browseURL("www.lancaster.ac.uk/staff/taylorb1/preprints/spatsurv.pdf") 
}


##' allocate function
##'
##' A function to allocate coordinates to an observation whose spatial location is known to the regional level
##'
##' @param poly a SpatialPolygonsDataFrame, on which the survival data exist in aggregate form
##' @param popden a sub-polygon raster image of population density
##' @param survdat data.frame containing the survival data
##' @param pid name of the variable in the survival data that gives the region identifier in poly
##' @param sid the name of the variable in poly to match the region identifier in survdat to
##' @param n the number of different allocations to make. e.g. if n is 2 (the default) two candidate sets of locations are available.
##' @param wid The default is 2000, interpreted in metres ie 2Km. size of buffer to add to window for raster cropping purposes: this ensures that for each polygon, the cropped raster covers it completely. 
##' @return matrices x and y, both of size (number of observations in survdat x n) giving n potential candidate locations of points in the columns of x and y.
##' @export

allocate <- function(poly,popden,survdat,pid,sid,n=2,wid=2000){
    nr <- length(poly)
    X <- matrix(NA,nrow(survdat),n)
    Y <- matrix(NA,nrow(survdat),n)
    for(i in 1:nr){
        progressreport(i,nr)
        spol <- gBuffer(poly[i,],width=wid)
        win <- as(poly[i,],"owin")
        den <- asImRaster(crop(popden,spol))
        idx <- survdat[,sid]==poly@data[i,pid]
        ns <- sum(idx)
        if(ns==0){
            next
        }
        else{
            pts <- rpoint(ns*n,f=den,win=win)
            for(j in 1:n){              
                X[which(idx),] <- matrix(pts$x,ns,n)
                Y[which(idx),] <- matrix(pts$y,ns,n)
            }
        }       
    }
    return(list(x=X,y=Y))
}



