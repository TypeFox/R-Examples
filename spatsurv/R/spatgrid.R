##' FFTgrid function
##'
##' A function to generate an FFT grid and associated quantities including cell dimensions,
##' size of extended grid, centroids,
##'
##' @param spatialdata a SpatialPixelsDataFrame object
##' @param cellwidth width of computational cells
##' @param ext multiplying constant: the size of the extended grid: ext*M by ext*N
##' @return a list
##' @export

FFTgrid <- function(spatialdata,cellwidth,ext){

    bb <- bbox(spatialdata)
    xrange <- bb[1,]
    yrange <- bb[2,]
    
    M <- 2^ceiling(log(diff(xrange)/cellwidth,base=2)) # minimum gridsize in x-direction
    N <- 2^ceiling(log(diff(yrange)/cellwidth,base=2)) # minimum gridsize in x-direction
    
    xdim <- M * cellwidth
    ydim <- N * cellwidth
    xadd <- (xdim - diff(xrange))/2
    yadd <- (ydim - diff(yrange))/2    

    xrange <- xrange + c(-xadd,xadd)
    yrange <- yrange + c(-yadd,yadd)  
    

    del1 <- (xrange[2]-xrange[1])/M
    del2 <- (yrange[2]-yrange[1])/N 
    
    Mext <- ext*M	
    Next <- ext*N
    
    mcens <- xrange[1]+.5*del1+(0:(Mext-1))*del1
    ncens <- yrange[1]+.5*del2+(0:(Next-1))*del2	

    obj <- list(del1=del1,del2=del2,Mext=Mext,Next=Next,mcens=mcens,ncens=ncens)    
    class(obj) <- "FFTgrid"
    
    return(obj)   
}

##' grid2spix function
##'
##' A function to convert a regular (x,y) grid of centroids into a SpatialPixels object
##'
##' @param xgrid vector of x centroids (equally spaced)
##' @param ygrid vector of x centroids (equally spaced)
##' @param proj4string an optional proj4string, projection string for the grid, set using the function CRS
##' @return a SpatialPixels object
##' @export

grid2spix <- function(xgrid,ygrid,proj4string=CRS(as.character(NA))){
    return(SpatialPixels(SpatialPoints(expand.grid(xgrid,ygrid),proj4string=proj4string)))
} 

##' grid2spts function
##'
##' A function to convert a regular (x,y) grid of centroids into a SpatialPoints object
##'
##' @param xgrid vector of x centroids (equally spaced)
##' @param ygrid vector of x centroids (equally spaced)
##' @param proj4string an optional proj4string, projection string for the grid, set using the function CRS
##' @return a SpatialPoints object
##' @export

grid2spts <- function(xgrid,ygrid,proj4string=CRS(as.character(NA))){
    return(SpatialPoints(expand.grid(xgrid,ygrid),proj4string=proj4string))
} 


##' grid2spdf function
##'
##' A function to convert a regular (x,y) grid of centroids into a SpatialPoints object 
##'
##' @param xgrid vector of x centroids (equally spaced)
##' @param ygrid vector of x centroids (equally spaced)
##' @param proj4string an optional proj4string, projection string for the grid, set using the function CRS
##' @return a SpatialPolygonsDataFrame
##' @export

grid2spdf <- function(xgrid,ygrid,proj4string=CRS(as.character(NA))){
  m <- length(xgrid)
  n <- length(ygrid)
  spts <- SpatialPixels(SpatialPoints(expand.grid(xgrid,ygrid),proj4string=proj4string))
  sps <- as(spts,"SpatialPolygons")
  spdf <- SpatialPolygonsDataFrame(sps,data=data.frame(grid=1:(m*n)),match.ID=FALSE)
  return(spdf)
}



##' gridY function
##'
##' A function to put estimated individual Y's onto a grid
##'
##' @param Y estimate of Y
##' @param control control parameters
##' @return ...
##' @export

gridY <- function(Y,control){
    newy <- matrix(-var(Y)/2,control$fftgrid$Mext,control$fftgrid$Next)

    Ytemp <- sapply(control$uqidx,function(i){mean(Y[control$idx==i])})    
    
    newy[control$uqidx] <- Ytemp
    return(newy)
}


##' gridY_polygonal function
##'
##' A function to put estimated individual Y's onto a grid
##'
##' @param Y estimate of Y
##' @param control control parameters
##' @return ...
##' @export

gridY_polygonal <- function(Y,control){
    newy <- rep(-var(Y)/2,control$n) #matrix(-var(Y)/2,control$fftgrid$Mext,control$fftgrid$Next)

    Ytemp <- sapply(control$uqidx,function(i){mean(Y[control$idx==i])})    
    
    newy[control$uqidx] <- Ytemp
    return(newy)
}




##' GammafromY function
##'
##' A function to change Ys (spatially correlated noise) into Gammas (white noise). Used in the MALA algorithm.
##'
##' @param Y Y matrix
##' @param rootQeigs square root of the eigenvectors of the precision matrix
##' @param mu parameter of the latent Gaussian field
##' @return Gamma
##' @export

GammafromY <- function(Y,rootQeigs,mu){
    nc <- dim(rootQeigs)[2]
    nb <- length(Y)
    return((1/nb)*Re(fft(fft(Y-mu)*rootQeigs,inverse=TRUE)))    
}




##' YfromGamma function
##'
##' A function to change Gammas (white noise) into Ys (spatially correlated noise). Used in the MALA algorithm. 
##'
##' @param Gamma Gamma matrix
##' @param invrootQeigs inverse square root of the eigenvectors of the precision matrix
##' @param mu parameter of the latent Gaussian field
##' @return Y
##' @export

YfromGamma <- function(Gamma,invrootQeigs,mu){
    nc <- dim(invrootQeigs)[2]
    nb <- length(Gamma)
    return(mu + (1/nb)*Re(fft(invrootQeigs*fft(Gamma,inverse=TRUE))))
}



##' getOptCellwidth function
##'
##' A function to compute an optimal cellwidth close to an initial suggestion. This maximises the efficiency of the 
##' MCMC algorithm when in the control argument of the function survspat, the option gridded is set to TRUE
##'
##' @param dat any spatial data object whose bounding box can be computed using the function bbox.
##' @param cellwidth an initial suggested cellwidth
##' @param ext the extension parameter for the FFT transform, set to 2 by default 
##' @param plot whether to plot the grid and data to illustrate the optimal grid
##' @return the optimum cell width
##' @export


getOptCellwidth <- function(dat,cellwidth,ext=2,plot=TRUE){
    dataArea <- prod(apply(bbox(dat),1,diff))
    fun <- function(x){
        gridobj <- FFTgrid(spatialdata=dat,cellwidth=x,ext=ext)
        del1 <- gridobj$del1
        del2 <- gridobj$del2
        mcens <- gridobj$mcens
        ncens <- gridobj$ncens 
        xdiff <- diff(range(mcens)) + del1
        ydiff <- diff(range(ncens)) + del2
        gridArea <- xdiff * ydiff
        return(gridArea/dataArea)
    }
    opt <- optimise(fun,interval=c(cellwidth/2,cellwidth*2))

    gridobj <- FFTgrid(spatialdata=dat,cellwidth=opt$minimum,ext=ext)
    del1 <- gridobj$del1
    del2 <- gridobj$del2
    Mext <- gridobj$Mext
    Next <- gridobj$Next
    mcens <- gridobj$mcens
    ncens <- gridobj$ncens 

    cat("Suggest using cellwidth:",opt$minimum,"which gives an output grid size of",Mext/ext," by ",Next/ext," this should be rounded upwards, if at all.\n")

    if(plot){
        plot(expand.grid(mcens[1:(Mext/ext)],ncens[1:(Next/ext)]),pch="+",col="red",cex=0.5)
        sp::plot(dat,cex=0.5,add=TRUE,pch=19)
    }

    return(opt$minimum)
}


##' showGrid function
##'
##' A function to show the grid that will be used for a given cellwidth
##'
##' @param dat any spatial data object whose bounding box can be computed using the function bbox.
##' @param cellwidth an initial suggested cellwidth
##' @param ext the extension parameter for the FFT transform, set to 2 by default 
##' @return a plot showing the grid and the data. Ideally the data should only just fit inside the grid.
##' @export

showGrid <- function(dat,cellwidth,ext=2){
    gridobj <- FFTgrid(spatialdata=dat,cellwidth=cellwidth,ext=ext)
    del1 <- gridobj$del1
    del2 <- gridobj$del2
    Mext <- gridobj$Mext
    Next <- gridobj$Next
    mcens <- gridobj$mcens
    ncens <- gridobj$ncens 

    plot(expand.grid(mcens[1:(Mext/ext)],ncens[1:(Next/ext)]),pch="+",col="red",cex=0.5)
    sp::plot(dat,cex=0.5,add=TRUE,pch=19)
}




## GP function
##
## A function to store a realisation of a spatial gaussian process for use in MCMC algorithms that include Bayesian parameter estimation.
## Stores not only the realisation, but also computational quantities.
##
## @param gamma the transformed (white noise) realisation of the process
## @param fftgrid an object of class FFTgrid, see ?genFFTgrid
## @param covFunction an object of class function returning the spatial covariance
## @param covParameters an object of class CovParamaters, see ?CovParamaters
## @param d matrix of grid distances
## @return a realisation of a spatial Gaussian process on a regular grid
## @export#
#
#GP <- function(gamma,fftgrid,covmodel,covparameters,u){
#    if(!inherits(gamma,"matrix")){
#        stop("argument 'gamma' must be a matrix")
#    }
#    if(class(fftgrid)!="FFTgrid"){
#        stop("argument 'fftgrid' must be an object of class FFTgrid, see ?genFFTgrid")
#    }
#    if(!any(class(covFunction)=="CovFunction")){
#        stop("argument 'covFunction' must be an object of class CovFunction")
#    }
#    if(class(covParameters)!="CovParamaters"){
#        stop("argument 'covParamaters' must be an object of class CovParamaters, see ?CovParamaters")
#    }
#    
#    Mext <- length(fftgrid$mcens)
#    Next <- length(fftgrid$ncens)
#    covbase <- covFunction(d=d,CovParameters=covParameters)
#    covbase <- lapply(covbase,matrix,nrow=Mext,ncol=Next)
#    ################
#
#    gp <- list()
#    gp$gamma <- gamma
#    gp$CovFunction <- covFunction 
#    gp$CovParameters <- covParameters
#    gp$covbase <- covbase 
#    #gp$fftcovbase <- lapply(covbase,function(x){Re(fft(x))})
#    gp$rootQeigs <- sqrt(1/Re(fft(covbase[[1]])))#sqrt(1/gp$fftcovbase$eval)
#    gp$invrootQeigs <- 1/gp$rootQeigs
#    gp$Y <- YfromGamma(gamma,invrootQeigs=gp$invrootQeigs,mu=covParameters$mu)
#    gp$expY <- exp(gp$Y) # for computational purposes
#    gp$mcens <- fftgrid$mcens
#    gp$ncens <- fftgrid$ncens
#    
#    class(gp) <- "GPrealisation"    
#    return(gp)    
#    
#}