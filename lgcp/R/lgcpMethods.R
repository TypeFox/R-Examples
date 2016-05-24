###
# generic methods
###

##' expectation function
##'
##' Generic function used in the computation of Monte Carlo expectations.
##'
##' @param obj an object    
##' @param ... additional arguments  
##' @return method expectation
##' @export

expectation <- function(obj,...){
    UseMethod("expectation")
}



##' extract function
##'
##' Generic function for extracting information dumped to disk. See \link{extract.lgcpPredict} for further information.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method extract
##' @seealso \link{extract.lgcpPredict}
##' @export

extract <- function(obj,...){
    UseMethod("extract")
}



##' lgcpgrid function
##'
##' Generic function for the hadling of list objects where each element of the list is
##' a matrix. Each matrix is assumed to have the same dimension. Such objects arise from the 
##' various routines in the package lgcp.
##'
##' lgcpgrid objects are list objects with names len, nrow, ncol, grid, xvals, yvals, zvals. The first three elements of the list
##' store the dimension of the object, the fourth element, grid, is itself a list object consisting of matrices in which the 
##' data is stored. The last three arguments can be used to give what is effectively a 3 dimensional array a physical reference.
##'
##' For example, the mean of Y from a call to lgcpPredict, obj$y.mean for example, is stored in an lgcpgrid object. If several time points have been
##' stored in the call to lgcpPredict, then the grid element of the lgcpgrid object contains the output for each of the time points in succession. So 
##' the first element, obj$y.mean$grid[[1]],contains the output from the first time point and so on.
##'
##'
##' @param grid a list object with each member of the list being a numeric matrix, each matrix having the same dimension
##' @param ... other arguments     
##' @return method lgcpgrid
##' @seealso \link{lgcpgrid.list}, \link{lgcpgrid.array}, \link{lgcpgrid.matrix}
##' @export

lgcpgrid <- function(grid,...){
    UseMethod("lgcpgrid")
}



###
# Functions associated with output grids
###



##' lgcpgrid.list function
##'
##' Creates an lgcpgrid object from a list object plus some optional coordinates. Note that each element of the list should be a matrix,
##' and that each matrix should have the same dimension.
##'
##' @method lgcpgrid list
##' @param grid a list object with each member of the list being a numeric matrix, each matrix having the same dimension
##' @param xvals optional vector of x-coordinates associated to grid. By default, this is the cell index in the x direction.
##' @param yvals optional vector of y-coordinates associated to grid. By default, this is the cell index in the y direction.
##' @param zvals optional vector of z-coordinates (time) associated to grid. By default, this is the cell index in the z direction.
##' @param ... other arguments      
##' @return an object of class lgcpgrid
##' @seealso \link{lgcpgrid.array}, \link{as.list.lgcpgrid}, \link{print.lgcpgrid}, 
##' \link{summary.lgcpgrid}, \link{quantile.lgcpgrid}, \link{image.lgcpgrid}, \link{plot.lgcpgrid} 
##' @export

lgcpgrid.list <- function(grid,xvals=1:dim(grid[[1]])[1],yvals=1:dim(grid[[1]])[2],zvals=1:length(grid),...){
    obj <- list()
    obj$grid <- grid
    obj$len <- length(grid)
    for(i in 1:obj$len){
        verifyclass(grid[[i]],"matrix")
    }
    dims <- sapply(grid,dim)
    if (!all(dims[1,]==dims[1,1]) | !all(dims[2,]==dims[2,1])){
        stop("Matrices in list must have the same dimension")
    }
    obj$nrow <- nrow(grid[[1]])
    obj$ncol <- ncol(grid[[1]])
    obj$xvals <- xvals
    obj$yvals <- yvals
    obj$zvals <- zvals
    class(obj) <- c("lgcpgrid","lgcpobject")
    return(obj)
}



##' lgcpgrid.array function
##'
##' Creates an lgcp grid object from an 3-dimensional array.
##'
##' @method lgcpgrid array
##' @param grid a three dimensional array object
##' @param xvals optional vector of x-coordinates associated to grid. By default, this is the cell index in the x direction.
##' @param yvals optional vector of y-coordinates associated to grid. By default, this is the cell index in the y direction.
##' @param zvals optional vector of z-coordinates (time) associated to grid. By default, this is the cell index in the z direction.
##' @param ... other arguments      
##' @return an object of class lgcpgrid
##' @seealso \link{lgcpgrid.list}, \link{as.list.lgcpgrid}, \link{print.lgcpgrid}, 
##' \link{summary.lgcpgrid}, \link{quantile.lgcpgrid}, \link{image.lgcpgrid}, \link{plot.lgcpgrid} 
##' @export

lgcpgrid.array <- function(grid,xvals=1:dim(grid)[1],yvals=1:dim(grid)[2],zvals=1:dim(grid)[3],...){
    if(length(dim(grid))!=3){
        stop("Array object must have 3 dimensions.")
    }
    obj <- list()
    d <- dim(grid)
    obj$len <- d[3]  
    obj$nrow <- d[1]
    obj$ncol <- d[2]
    obj$grid <- unlist(apply(grid,3,function(x){return(list(x))}),recursive=FALSE)
    obj$xvals <- xvals
    obj$yvals <- yvals
    obj$zvals <- zvals
    class(obj) <- c("lgcpgrid","lgcpobject")
    return(obj)
}

##' lgcpgrid.matrix function
##'
##' Creates an lgcp grid object from an 2-dimensional matrix.
##'
##' @method lgcpgrid matrix
##' @param grid a three dimensional array object
##' @param xvals optional vector of x-coordinates associated to grid. By default, this is the cell index in the x direction.
##' @param yvals optional vector of y-coordinates associated to grid. By default, this is the cell index in the y direction.
##' @param ... other arguments      
##' @return an object of class lgcpgrid
##' @seealso \link{lgcpgrid.list}, \link{as.list.lgcpgrid}, \link{print.lgcpgrid}, 
##' \link{summary.lgcpgrid}, \link{quantile.lgcpgrid}, \link{image.lgcpgrid}, \link{plot.lgcpgrid} 
##' @export

lgcpgrid.matrix <- function(grid,xvals=1:nrow(grid),yvals=1:ncol(grid),...){
    arr <- array(dim=c(dim(grid),1))
    arr[,,1] <- grid
    return(lgcpgrid.array(arr,xvals=xvals,yvals=yvals,zvals=1))
}



##' as.list.lgcpgrid function
##'
##' Method to convert an lgcpgrid object into a list of matrices.
##'
##' @method as.list lgcpgrid
##' @param x an object of class lgcpgrid  
##' @param ... other arguments   
##' @return conversion from lgcpgrid to list
##' @seealso \link{lgcpgrid.list}, \link{lgcpgrid.array}, \link{print.lgcpgrid}, 
##' \link{summary.lgcpgrid}, \link{quantile.lgcpgrid}, \link{image.lgcpgrid}, \link{plot.lgcpgrid} 
##' @export

as.list.lgcpgrid <- function(x,...){
    return(x$grid)
}


##' as.array.lgcpgrid function
##'
##' Method to convert an lgcpgrid object into an array.
##'
##' @method as.array lgcpgrid
##' @param x an object of class lgcpgrid  
##' @param ... other arguments
##' @return conversion from lgcpgrid to array
##' @export
   
as.array.lgcpgrid <- function(x,...){
    arr <- array(dim=c(x$nrow,x$ncol,x$len))
    for (i in 1:x$len){
        arr[,,i] <- x$grid[[i]]
    }
    return(arr)
} 

##' raster.lgcpgrid function
##'
##' A function to convert lgcpgrid objects into either a raster object, or a RasterBrick object.
##'
##' @method raster lgcpgrid
##' @param x an lgcpgrid object
##' @param crs PROJ4 type description of a map projection (optional). See ?raster
##' @param transpose Logical. Transpose the data? See ?brick method for array 
##' @param ... additional arguments
##' @return ...
##' @export

raster.lgcpgrid <- function(x,crs=NA,transpose=FALSE,...){
    if (is.null(x$xvals)|is.null(x$yvals)){ # for backwards compatibility 
        x$xvals <- 1:nrow(x$grid[[1]])
        x$yvals <- 1:ncol(x$grid[[1]])
        x$zvals <- 1:length(x$grid)
    }
    M <- length(x$xvals)
    N <- length(x$yvals)
    dx <- diff(x$xvals[1:2])
    dy <- diff(x$yvals[1:2])
    if (length(x$zvals)==1){    
        return(raster(t(x$grid[[1]][,N:1]), xmn=x$xvals[1]-dx/2, xmx=x$xvals[M]+dx/2, ymn=x$yvals[1]-dy/2, ymx=x$yvals[N]+dy/2, crs=crs))
    }
    else{
        x$grid <- lapply(x$grid,function(xx){t(xx[,N:1])})
        ar <- as.array(x)        
        br <- brick(ar, xmn=x$xvals[1]-dx/2, xmx=x$xvals[M]+dx/2, ymn=x$yvals[1]-dy/2, ymx=x$yvals[N]+dy/2, crs=crs,transpose=transpose)
        names(br) <- paste("Time",x$zvals,sep="")
        return(br)
    }
}

suppressWarnings(setMethod("raster","lgcpgrid",raster.lgcpgrid))

##' as.SpatialPixelsDataFrame function
##'
##' Generic function for conversion to SpatialPixels objects. 
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method as.SpatialPixels
##' @seealso \link{as.SpatialPixelsDataFrame.lgcpgrid}
##' @export

as.SpatialPixelsDataFrame <- function(obj,...){
    UseMethod("as.SpatialPixelsDataFrame")
}



##' as.SpatialPixelsDataFrame.lgcpgrid function
##'
##' Method to convert lgcpgrid objects to SpatialPixelsDataFrame objects.
##' 
##'
##' @method as.SpatialPixelsDataFrame lgcpgrid
##' @param obj an lgcpgrid object
##' @param ... additional arguments to be passed to SpatialPoints, eg a proj4string
##' @return Either a SpatialPixelsDataFrame, or a list consisting of SpatialPixelsDataFrame objects.
##' @export

as.SpatialPixelsDataFrame.lgcpgrid <- function(obj,...){
    if (is.null(obj$xvals)|is.null(obj$yvals)){ # for backwards compatibility 
        obj$xvals <- 1:nrow(obj$grid[[1]])
        obj$yvals <- 1:ncol(obj$grid[[1]])
        obj$zvals <- 1:length(obj$grid)
    }
    spts <- SpatialPoints(expand.grid(obj$xvals,obj$yvals),...)
    dat <- list()
    for(i in 1:length(obj$grid)){
        dat[[i]] <- as.vector(obj$grid[[i]])
    }
    
    if(length(obj$grid)==1){
        names(dat) <- "out"
    }
    else{    
        names(dat) <- paste("outTime",obj$zvals,sep="")
    }
    dat <- as.data.frame(dat)
        
    return(SpatialPixelsDataFrame(spts,data=dat))
}





##' print.lgcpgrid function
##'
##' Print method for lgcp grid objects.
##'
##' @method print lgcpgrid
##' @param x an object of class lgcpgrid  
##' @param ... other arguments   
##' @return just prints out details to the console
##' @seealso \link{lgcpgrid.list}, \link{lgcpgrid.array}, \link{as.list.lgcpgrid},
##' \link{summary.lgcpgrid} \link{quantile.lgcpgrid} \link{image.lgcpgrid} \link{plot.lgcpgrid} 
##' @export

print.lgcpgrid <- function(x,...){
    cat("lgcpgrid object.\n")
    cat(paste("   Length: ",x$len,"\n",sep=""))
    cat(paste("  Matsize: [ ",x$nrow," , ",x$ncol," ]\n",sep=""))
}



##' summary.lgcpgrid function
##'
##' Summary method for lgcp objects. This just applies the summary function to each
##' of the elements of object$grid.
##'
##' @method summary lgcpgrid
##' @param object an object of class lgcpgrid  
##' @param ... other arguments   
##' @return Summary per grid, see ?summary for further options
##' @seealso \link{lgcpgrid.list}, \link{lgcpgrid.array}, \link{as.list.lgcpgrid}, \link{print.lgcpgrid}, 
##' \link{quantile.lgcpgrid}, \link{image.lgcpgrid}, \link{plot.lgcpgrid} 
##' @export

summary.lgcpgrid <- function(object,...){
    for(i in 1:object$len){
        cat(paste("Summary for grid ",i,":\n",sep=""))
        print(summary(as.vector(object$grid[[i]],...)))
        cat("\n")
    }
}



##' quantile.lgcpgrid function
##'
##' Quantile method for lgcp objects. This just applies the quantile function to each of
##' the elements of x$grid
##'
##' @method quantile lgcpgrid
##' @param x an object of class lgcpgrid  
##' @param ... other arguments   
##' @return Quantiles per grid, see ?quantile for further options
##' @seealso \link{lgcpgrid.list}, \link{lgcpgrid.array}, \link{as.list.lgcpgrid}, \link{print.lgcpgrid}, 
##' \link{summary.lgcpgrid}, \link{image.lgcpgrid}, \link{plot.lgcpgrid} 
##' @export

quantile.lgcpgrid <- function(x,...){
    for(i in 1:x$len){
        cat(paste("Quantiles for grid ",i,":\n",sep=""))
        print(quantile(x$grid[[i]],...))
        cat("\n")
    }
}



##' image.lgcpgrid function
##'
##' Produce an image plot of an lgcpgrid object.
##'
##' @method image lgcpgrid
##' @param x an object of class lgcpgrid
##' @param sel vector of integers between 1 and grid$len: which grids to plot. Default NULL, in which case all grids are plotted.
##' @param ask logical; if TRUE the user is asked before each plot  
##' @param ... other arguments   
##' @return grid plotting
##' @seealso \link{lgcpgrid.list}, \link{lgcpgrid.array}, \link{as.list.lgcpgrid}, \link{print.lgcpgrid}, 
##' \link{summary.lgcpgrid}, \link{quantile.lgcpgrid}, \link{plot.lgcpgrid} 
##' @export

image.lgcpgrid <- function(x,sel=1:x$len,ask=TRUE,...){
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    for(i in sel){
        if(is.null(x$xvals)|is.null(x$yvals)){
            image.plot(1:x$nrow,1:x$ncol,x$grid[[i]],...)
        }
        else{
            image.plot(x$xvals,x$yvals,x$grid[[i]],...)
        }
    }    
} 



##' plot.lgcpgrid function
##'
##' This is a wrapper function for image.lgcpgrid
##'
##' @method plot lgcpgrid
##' @param x an object of class lgcpgrid
##' @param sel vector of integers between 1 and grid$len: which grids to plot. Default NULL, in which case all grids are plotted.
##' @param ask logical; if TRUE the user is asked before each plot  
##' @param ... other arguments   
##' @return an image-type plot
##' @seealso \link{lgcpgrid.list}, \link{lgcpgrid.array}, \link{as.list.lgcpgrid}, \link{print.lgcpgrid}, 
##' \link{summary.lgcpgrid},\link{quantile.lgcpgrid}, \link{image.lgcpgrid}
##' @export

plot.lgcpgrid <- function(x,sel=1:x$len,ask=TRUE,...){
    image.lgcpgrid(x,sel=1:x$len,ask=TRUE,...)  
} 



###
# Functions associated with objects of class lgcpPredict
###



##' print.lgcpPredict function
##'
##' Print method for lgcpPredict objects.
##'
##' @method print lgcpPredict
##' @param x an object of class lgcpPredict    
##' @param ... additional arguments 
##' @return just prints information to the screen
##' @seealso \link{lgcpPredict}
##' @export

print.lgcpPredict <- function(x,...){
    cat("lgcpPredict object.\n")
    cat("\n")
    cat("  General Information\n")
    cat("  -------------------\n")
    #flag <- TRUE

    #if(is.null(x$spatialonly)){
    #    flag <- FALSE
    #}
    #else{
    #    flag <- x$spatialonly
    #}
    
    #if(!flag){ # is.null(x$spatialonly) for backwards compatibility
    #    cat(paste("      FFT Gridsize: [ ",x$ext*x$M," , ",x$ext*x$N," ]\n",sep=""))
    #    cat("\n")
    #
    #    cat(paste("              Data:\n",sep=""))
    #    neattable(matrix(c("Time |",x$aggtimes,"Counts |",sapply(x$nis,sum)),nrow=2,byrow=TRUE),indent=4)
    #    cat("\n")
    #    cat(paste("        Parameters: sigma=",round(x$sigma,3),", phi=",round(x$phi,3),", theta=",round(x$theta,3),"\n",sep=""))
    #}
    #else{ # spatial only
        cat(paste("      FFT Gridsize: [ ",x$ext*x$M," , ",x$ext*x$N," ]\n",sep=""))
        cat("\n")
    
    #    cat(paste("        Parameters: sigma=",round(x$sigma,3),", phi=",round(x$phi,3),"\n",sep=""))
    #}
    if (!is.null(x$gridfunction)){
        cat(paste("    Dump Directory: ",x$gridfunction$dirname,"\n",sep=""))
    }
    if (!is.null(x$gridaverage)){
        cat("\n")
        cat(paste("     Grid Averages:\n",sep=""))
        fname <- unlist(x$gridaverage[[1]])
        fclass <- sapply(x$gridaverage[[2]],function(xx){class(xx)})
        neattable(matrix(c("Function",fname,"Output Class",fclass),ncol=2),indent=4)
        cat("\n")
    }
    
    cat(paste("        Time taken: ",round(x$timetaken,3)," ",units(x$timetaken),"\n",sep=""))
    cat("\n")
    
    cat("  MCMC Information\n")
    cat("  ----------------\n")
    cat(paste("    Number Iterations: ",x$mcmcinfo$N,"\n",sep=""))
    cat(paste("              Burn-in: ",x$mcmcinfo$burnin,"\n",sep=""))
    cat(paste("             Thinning: ",x$mcmcinfo$thin,"\n",sep=""))
    cat(paste("      Mean Acceptance: ",round(mean(x$mcmcacc),3),"\n",sep=""))  
    cat(paste("      Adaptive Scheme: ",class(x$mcmcpars$adaptivescheme)[1],"\n",sep=""))
    cat(paste("               Last h: ",x$lasth,"\n",sep=""))
}  



##' neattable function
##'
##' Function to print right-aligned tables to the console.
##'
##' @param mat a numeric or character matrix object
##' @param indent indent
##' @return prints to screen with specified indent
##' @examples
##' mat <- rbind(c("one","two","three"),matrix(round(runif(9),3),3,3))
##' neattable(mat)
##' @export

neattable <- function(mat,indent=0){
    mat <- t(apply(mat,1,as.character))
    wd <- max(nchar(mat)) + 1
    sp <- wd - nchar(mat)    
    for (i in 1:nrow(mat)){
        cat(rep(" ",indent),sep="")
        for (j in 1:ncol(mat)){
            cat(rep(" ",sp[i,j]),sep="")
            cat(mat[i,j],sep="")
        }
        cat("\n")
    }
}   


###
# Accessor functions
### 

##' xvals.lgcpPredict function
##'
##' Gets the x-coordinates of the centroids of the prediction grid.
##'
##' @method xvals lgcpPredict
##' @param obj an object of class lgcpPredict
##' @param ... additional arguments
##' @return the x coordinates of the centroids of the grid
##' @seealso \link{lgcpPredict}
##' @export

xvals.lgcpPredict <- function(obj,...){
    return(obj$mcens)
}



##' yvals.lgcpPredict function
##'
##' Gets the y-coordinates of the centroids of the prediction grid.
##'
##' @method yvals lgcpPredict
##' @param obj an object of class lgcpPredict
##' @param ... additional arguments
##' @return the y coordinates of the centroids of the grid
##' @seealso \link{lgcpPredict}
##' @export

yvals.lgcpPredict <- function(obj,...){
    return(obj$ncens)
}



##' plot.lgcpPredict function
##'
##' Simple plotting function for objects of class \code{lgcpPredict}.
##'
##' @method plot lgcpPredict
##' @param x an object of class lgcpPredict
##' @param type Character string: what type of plot to produce. Choices are "relrisk" (=exp(Y)); "serr" (standard error of relative risk); or "intensity" (=lambda*mu*exp(Y)).
##' @param sel vector of integers between 1 and grid$len: which grids to plot. Default NULL, in which case all grids are plotted.
##' @param plotdata whether or not to overlay the data
##' @param ask logical; if TRUE the user is asked before each plot
##' @param clipWindow whether to plot grid cells outside the observation window  
##' @param ... additional arguments passed to image.plot
##' @return plots the Monte Carlo mean of quantities obtained via simulation. By default the mean relative risk is plotted.
##' @seealso \link{lgcpPredict}
##' @export

plot.lgcpPredict <- function(x,type="relrisk",sel=1:x$EY.mean$len,plotdata=TRUE,ask=TRUE,clipWindow=TRUE,...){#
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    for(i in sel){
        if (!clipWindow){
            if (type=="relrisk"){
                image.plot(x$mcens,x$ncens,x$EY.mean$grid[[i]],sub=paste("Relative Risk, time",x$aggtimes[i]),...)
            }
            else if (type=="serr"){
                serr <- serr(x)$grid
                image.plot(x$mcens,x$ncens,serr[[i]],sub=paste("S.E. Relative Risk, time",x$aggtimes[i]),...)
            }
            else if (type=="intensity"){
                image.plot(x$mcens,x$ncens,x$temporal[i]*x$grid[[i]][1:x$M,1:x$N]*x$EY.mean$grid[[i]],sub=paste("Poisson Intensity, time",x$aggtimes[i]),...)
            }
            else{
                stop("type must be 'relrisk', 'serr' or 'intensity'")            
            }
        }
        else{
            incl <- "centroid"
         
            if(!is.null(x$cellInside)){
                dm <- dim(x$cellInside)
                grinw <- array(as.logical(x$cellInside),dim=dm)
            }
            else{
                grinw <- gridInWindow(xvals,yvals,window,inclusion=incl)
            }
            grinw[grinw==0] <- NA
            if (type=="relrisk"){
                image.plot(x$mcens,x$ncens,grinw*x$EY.mean$grid[[i]],sub=paste("Relative Risk, time",x$aggtimes[i]),...)
            }
            else if (type=="serr"){
                serr <- serr(x)$grid
                image.plot(x$mcens,x$ncens,grinw*serr[[i]],sub=paste("S.E. Relative Risk, time",x$aggtimes[i]),...)
            }
            else if (type=="intensity"){
                image.plot(x$mcens,x$ncens,grinw*x$temporal[i]*x$grid[[i]][1:x$M,1:x$N]*x$EY.mean$grid[[i]],sub=paste("Poisson Intensity, time",x$aggtimes[i]),...)
            }
            else{
                stop("type must be 'relrisk', 'serr' or 'intensity'") 
            }
        }        
        
        if(plotdata){
            plot(x$xyt[as.integer(x$xyt$t)==x$aggtimes[[i]]],add=TRUE)
        }
        plot(x$xyt$window,lwd=2,add=TRUE)
    }    
}



##' meanfield function
##'
##' Generic function to extract the mean of the latent field Y.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method meanfield
##' @export

meanfield <- function(obj,...){
    UseMethod("meanfield")
}



##' meanfield.lgcpPredict function
##'
##' This is an accessor function for objects of class \code{lgcpPredict} and returns the mean of the 
##' field Y as an lgcpgrid object.
##'
##' @method meanfield lgcpPredict
##' @param obj an object of class lgcpPredict
##' @param ... additional arguments
##' @return returns the cell-wise mean of Y computed via Monte Carlo.
##' @seealso \link{lgcpPredict}, \link{lgcpgrid}
##' @export

meanfield.lgcpPredict <- function(obj,...){
    return(obj$y.mean)
}  



##' varfield function
##'
##' Generic function to extract the variance of the latent field Y.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method meanfield
##' @seealso \link{lgcpPredict}
##' @export

varfield <- function(obj,...){
    UseMethod("varfield")
}



##' varfield.lgcpPredict function
##'
##' This is an accessor function for objects of class \code{lgcpPredict} and returns the variance of the 
##' field Y as an lgcpgrid object.
##'
##' @method varfield lgcpPredict
##' @param obj an object of class lgcpPredict
##' @param ... additional arguments
##' @return returns the cell-wise variance of Y computed via Monte Carlo.
##' @seealso \link{lgcpPredict}
##' @export

varfield.lgcpPredict <- function(obj,...){
    return(obj$y.var)
}



##' rr function
##'
##' Generic function to return relative risk.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method rr
##' @seealso \link{lgcpPredict}, \link{rr.lgcpPredict}
##' @export

rr <- function(obj,...){
    UseMethod("rr")
}



##' rr.lgcpPredict function
##'
##' Accessor function returning the relative risk = exp(Y) as an lgcpgrid object. 
##'
##' @method rr lgcpPredict
##' @param obj an lgcpPredict object
##' @param ... additional arguments
##' @return the relative risk as computed my MCMC
##' @seealso \link{lgcpPredict}
##' @export

rr.lgcpPredict <- function(obj,...){
    return(obj$EY.mean)
}

##' serr function
##'
##' Generic function to return standard error of relative risk.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method serr
##' @seealso \link{lgcpPredict}, \link{serr.lgcpPredict}
##' @export

serr <- function(obj,...){
    UseMethod("serr")
}



##' serr.lgcpPredict function
##'
##' Accessor function returning the standard error of relative risk as an lgcpgrid object. 
##'
##' @method serr lgcpPredict
##' @param obj an lgcpPredict object
##' @param ... additional arguments
##' @return Standard error of the relative risk as computed by MCMC.
##' @seealso \link{lgcpPredict}
##' @export

serr.lgcpPredict <- function(obj,...){
    se <- lapply(obj$EY.var$grid,sqrt)
    return(lgcpgrid(se,xvals=obj$mcens,yvals=obj$ncens,zvals=obj$aggtimes))
}


##' intens function
##'
##' Generic function to return the Poisson Intensity.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method intens
##' @seealso \link{lgcpPredict}, \link{intens.lgcpPredict}
##' @export

intens <- function(obj,...){
    UseMethod("intens")
}



##' intens.lgcpPredict function
##'
##' Accessor function returning the Poisson intensity as an lgcpgrid object.
##'
##' @method intens lgcpPredict
##' @param obj an lgcpPredict object
##' @param ... additional arguments
##' @return the cell-wise mean Poisson intensity, as computed by MCMC.
##' @seealso \link{lgcpPredict}
##' @export

intens.lgcpPredict <- function(obj,...){
    intens <- list()
    MN <- dim(obj$EY.mean$grid[[1]])
    print(MN)
    cellarea <- diff(obj$mcens[1:2])*diff(obj$ncens[1:2]) 
    for(i in 1:obj$EY.mean$len){
        intens[[i]] <- cellarea*obj$temporal[i]*obj$grid[[i]][1:MN[1],1:MN[2]]*obj$EY.mean$grid[[i]]
    }
    return(lgcpgrid(intens,xvals=obj$mcens,yvals=obj$ncens,zvals=obj$aggtimes))
}



##' seintens function
##'
##' Generic function to return the standard error of the Poisson Intensity.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method seintens
##' @seealso \link{lgcpPredict}, \link{seintens.lgcpPredict}
##' @export

seintens <- function(obj,...){
    UseMethod("seintens")
}



##' seintens.lgcpPredict function
##'
##' Accessor function returning the  standard error of the Poisson intensity as an lgcpgrid object.
##'
##' @method seintens lgcpPredict
##' @param obj an lgcpPredict object
##' @param ... additional arguments
##' @return the cell-wise standard error of the Poisson intensity, as computed by MCMC.
##' @seealso \link{lgcpPredict}
##' @export

seintens.lgcpPredict <- function(obj,...){
    seintens <- list()
    MN <- dim(obj$EY.var$grid[[1]])
    cellarea <- diff(obj$mcens[1:2])*diff(obj$ncens[1:2]) 
    for(i in 1:obj$EY.var$len){
        seintens[[i]] <- cellarea*obj$temporal[i]*obj$grid[[i]][1:MN[1],1:MN[2]] * sqrt(obj$EY.var$grid[[i]])
    }
    return(lgcpgrid(seintens,xvals=obj$mcens,yvals=obj$ncens,zvals=obj$aggtimes))
}


##' discreteWindow function
##'
##' Generic function for extracting the FFT discrete window.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method discreteWindow
##' @seealso \link{discreteWindow.lgcpPredict}
##' @export

discreteWindow <- function(obj,...){
    UseMethod("discreteWindow")
}



##' discreteWindow.lgcpPredict function
##'
##' A function for extracting the FFT discrete window from an lgcpPredict object.
##'
##' @method discreteWindow lgcpPredict
##' @param obj an lgcpPredict object
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
##' @param ... additional arguments
##' @return ...
##' @export

discreteWindow.lgcpPredict <- function(obj,inclusion="touching",...){
    xv <- xvals(obj)
    yv <- yvals(obj)
    gr <- expand.grid(xv,yv)
        
    if(inclusion=="centroid"){
        return(matrix(inside.owin(gr[,1],gr[,2],obj$xyt$window),obj$M,obj$N))
    }
    else if(inclusion=="touching"){
        return(matrix(touchingowin(xv,yv,obj$xyt$window),obj$M,obj$N))
    }
    else{
        stop("Invlaid choice for argument 'inclusion'.")
    }
}




##' gridfun function
##'
##' A generic function for returning \code{gridfunction} objects. 
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method gridfun
##' @seealso \link{setoutput}, \link{lgcpgrid}
##' @export

gridfun <- function(obj,...){
    UseMethod("gridfun")
}



##' gridfun.lgcpPredict function
##' 
##' Accessor function for \code{lgcpPredict objects}: returns the \code{gridfunction} argument
##' set in the \code{output.control} argument of the function \code{lgcpPredict}.
##'
##' @method gridfun lgcpPredict
##' @param obj an object of class lgcpPredict
##' @param ... additional arguments
##' @return returns the output from the gridfunction option of the setoutput argument of lgcpPredict
##' @seealso \link{setoutput}, \link{lgcpgrid}
##' @export

gridfun.lgcpPredict <- function(obj,...){
    if (is.null(obj$gridfunction)){
        stop("No gridfunction was specified")
    }
    return(obj$gridfunction)
}



##' gridav function
##'
##' A generic function for returning \code{gridmeans} objects.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method gridav
##' @seealso \link{setoutput}, \link{lgcpgrid}
##' @export

gridav <- function(obj,...){
    UseMethod("gridav")
}



##' gridav.lgcpPredict function
##'
##' Accessor function for \code{lgcpPredict objects}: returns the \code{gridmeans} argument
##' set in the \code{output.control} argument of the function \code{lgcpPredict}.
##'
##' @method gridav lgcpPredict
##' @param obj an object of class lgcpPredict
##' @param fun an optional character vector of length 1 giving the name of a function to return Monte Carlo average of
##' @param ... additional arguments
##' @return returns the output from the gridmeans option of the setoutput argument of lgcpPredict
##' @seealso \link{setoutput}, \link{lgcpgrid}
##' @export

gridav.lgcpPredict <- function(obj,fun=NULL,...){
    if (is.null(obj$gridaverage)){
        stop("No Monte Carlo averages were requested.")
    }
    if(is.null(fun)){
        x <- list()
        x$names <- obj$gridaverage[[1]]
        x$output <- obj$gridaverage[[2]]
        class(x) <- c("gridaverage","lgcpobject")
        return(x)
    }
    else{
        if(length(fun)>1){
            stop("fun must be a character vector of length 1.")
        }
        if(sum(obj$gridaverage[[1]]==fun)==0){
            stop("The Monte Carlo average of that function was not requested.")
        }
        return(obj$gridaverage[[2]][[which(obj$gridaverage[[1]]==fun)]])
    }
}



##' print.gridaverage function
##'
##' Print method for \code{gridaverage} objects
##'
##' @method print gridaverage
##' @param x an object of class gridaverage  
##' @param ... other arguments   
##' @return just prints out details
##' @export

print.gridaverage <- function(x,...){
    cat("gridaverage object.\n")
    cat("  Functions:",paste(unlist(x$names),collapse=" "),"\n")
    cat("  $names returns function names, $output returns function output.\n")
}



##' hvals function
##'
##' Generic function to return the values of the proposal scaling \eqn{h}{h} in the MCMC algorithm.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method hvals
##' @export

hvals <- function(obj,...){
    UseMethod("hvals")
}



##' hvals.lgcpPredict function
##'
##' Accessor function returning the value of \eqn{h}{h}, the MALA proposal scaling constant over the iterations of the algorithm for 
##' objects of class \code{lgcpPredict}
##'
##' @method hvals lgcpPredict
##' @param obj an object of class lgcpPredict
##' @param ... additional arguments
##' @return returns the values of h taken during the progress of the algorithm
##' @seealso \link{lgcpPredict}
##' @export

hvals.lgcpPredict <- function(obj,...){
    return(obj$hrec)
}



##' window.lgcpPredict function
##'
##' Accessor function returning the observation window from objects of class \code{lgcpPredict}. Note that for
##' computational purposes, the window of an \code{stppp} object will be extended to accommodate the requirement that 
##' the dimensions must be powers of 2. The function \code{window.lgcpPredict} returns the extended window.
##'
##' @method window lgcpPredict
##' @param x an object of class lgcpPredict
##' @param ... additional arguments
##' @return returns the observation window used durign computation
##' @seealso \link{lgcpPredict}
##' @export

window.lgcpPredict <- function(x,...){
    return(x$xyt$window)
}



##' mcmctrace function
##'
##' Generic function to extract the information required to produce MCMC trace plots.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method mcmctrace
##' @export

mcmctrace <- function(obj,...){
    UseMethod("mcmctrace")
}



##' mcmctrace.lgcpPredict function
##'
##' If \code{MCMCdiag} was positive when \code{lgcpPredict} was called, then this retrieves information from the chains stored.
##'
##' @method mcmctrace lgcpPredict
##' @param obj an object of class lgcpPredict
##' @param ... additional arguments
##' @return returns the saved MCMC chains in an object of class \code{mcmcdiag}.
##' @seealso \link{lgcpPredict}, \link{plot.mcmcdiag}
##' @export

mcmctrace.lgcpPredict <- function(obj,...){
    if (is.null(obj$mcmcsamp)){
        stop("No information is available on the chain unless it was dumped to disk.")
    }
    tr <- list()
    tr$cellidx <- obj$ijcells
    tr$trace <- obj$mcmcsamp
    class(tr) <- "mcmcdiag"
    return(tr)
}

##' plot.mcmcdiag function
##'
##' The command \code{plot(trace(lg))}, where \code{lg} is an object of class \code{lgcpPredict} will plot the 
##' mcmc traces of a subset of the cells, provided they have been stored, see \code{mcmpars}.
##'
##' @method plot mcmcdiag
##' @param x an object of class mcmcdiag
##' @param idx vector of chain indices to plot, default plots all chains
##' @param ... additional arguments passed to plot
##' @return plots the saved MCMC chains
##' @seealso \link{mcmctrace.lgcpPredict}, \link{mcmcpars}, 
##' @export

plot.mcmcdiag <- function(x,idx=1:dim(x$trace)[2],...){
    plot(NULL,xlim=c(1,dim(x$trace)[1]),ylim=range(x$trace),xlab="Iteration",ylab="Y",...)
    n <- length(idx)
    cols <- rainbow(n)
    lab <- apply(x$cellidx,1,paste,collapse=",")[idx]
    for (i in 1:n){
        lines(x$trace[,idx[i]],col=cols[i])
    }
    legend("topleft",lty=rep("solid",n),col=cols,legend=lab,title="Cell Index")
}



##' plotExceed function
##'
##' A generic function for plotting exceedance probabilities.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return generic function returning method plotExceed
##' @seealso \link{plotExceed.lgcpPredict}, \link{plotExceed.array}
##' @export

plotExceed <- function(obj,...){
    UseMethod("plotExceed")
}



##' plotExceed.array function
##'
##' Function for plotting exceedance probabilities stored in array objects. Used in \code{plotExceed.lgcpPredict}.
##'
##' @method plotExceed array
##' @param obj an object
##' @param fun the name of the function used to compute exceedances (character vector of length 1). Note that the named function must be in memory.
##' @param lgcppredict an object of class lgcpPredict that can be used to supply an observation window and x and y coordinates
##' @param xvals optional vector giving x coords of centroids of cells
##' @param yvals optional vector giving y coords of centroids of cells
##' @param window optional obervation window
##' @param cases optional xy (n x 2) matrix of locations of cases to plot 
##' @param nlevel number of colour levels to use in plot, default is 64
##' @param ask whether or not to ask for a new plot between plotting exceedances at different thresholds.
##' @param mapunderlay optional underlay to plot underneath maps of exceedance probabilities. Use in conjunction with rainbow parameter 'alpha' (eg alpha=0.3) to set transparency of exceedance layer.
##' @param alpha graphical parameter takign values in [0,1] controlling transparency of exceedance layer. Default is 1.
##' @param sub optional subtitle for plot
##' @param ... additional arguments passed to image.plot
##' @return generic function returning method plotExceed
##' @seealso \link{plotExceed.lgcpPredict}
##' @export

plotExceed.array <- function(obj,fun,lgcppredict=NULL,xvals=NULL,yvals=NULL,window=NULL,cases=NULL,nlevel=64,ask=TRUE,mapunderlay=NULL,alpha=1,sub=NULL,...){
    if (!is.null(lgcppredict)){
        xvals <- lgcppredict$mcens
        yvals <- lgcppredict$ncens
        window <- lgcppredict$xyt$window
    }
    if (is.null(xvals)){
        xvals <- 1:dim(obj)[1]
    }
    if (is.null(yvals)){
        yvals <- 1:dim(obj)[2]
    }
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    
    relrisks <- attr(get(fun),"threshold")
    
    if (is.null(window)){
        grinw <- 1
    }
    else{
        incl <- "touching"
        #if(!inherits(try(obj$inclusion),"try-error")){
        #    if(!is.null(obj$inclusion)){
        #        incl <- obj$inclusion
        #    }
        #}
        
        if(!is.null(lgcppredict)){
            dm <- dim(lgcppredict$cellInside)
            grinw <- array(as.logical(lgcppredict$cellInside),dim=dm)
        }
        else{
            grinw <- gridInWindow(xvals,yvals,window,inclusion=incl)
        }
        grinw[as.numeric(grinw)==0] <- NA
    }
    
    for (i in 1:length(relrisks)){
        subt <- sub
        if(is.null(subt)){
            subt <- paste("Threshold:",relrisks[i])
        }
        if(is.null(mapunderlay)){
            image.plot(xvals,yvals,grinw*obj[,,i],col=rev(heat.colors(nlevel)),sub=subt,...)
        }
        else{
            plot(mapunderlay)
            image.plot(xvals,yvals,grinw*obj[,,i],col=rev(heat.colors(nlevel)),sub=subt,add=TRUE,...)        
        }        
        
        if (!is.null(window)){
            plot(window,lwd=2,add=TRUE)
        }
        if(!is.null(cases)){
            points(cases,pch="+",cex=0.5)
        }
    }
}



##' is.pow2 function
##'
##' Tests whether a number id
##'
##' @param num a numeric
##' @return logical: is num a power of 2?
##' @examples
##' is.pow2(128)  # TRUE
##' is.pow2(64.9) # FALSE
##' @export

is.pow2 <- function(num){
    return(isTRUE(all.equal(log(num,base=2),as.integer(log(num,base=2)))))
}



##' plotExceed.lgcpPredict function
##'
##' Function for plotting exceedance probabilities stored in \code{lgcpPredict} ojects.
##'
##' @method plotExceed lgcpPredict
##' @param obj an object
##' @param fun the name of the function used to compute exceedances (character vector of length 1). Note that the named function must be in memory.
##' @param nlevel number of colour levels to use in plot, default is 64
##' @param ask whether or not to ask for a new plot between plotting exceedances at different thresholds.
##' @param plotcases whether or not to plot the cases on the map
##' @param mapunderlay optional underlay to plot underneath maps of exceedance probabilities. Use in conjunction with rainbow parameter 'alpha' (eg alpha=0.3) to set transparency of exceedance layer.
##' @param alpha graphical parameter takign values in [0,1] controlling transparency of exceedance layer. Default is 1.
##' @param ... additional arguments passed to image.plot
##' @return plot of exceedances
##' @seealso \link{lgcpPredict}, \link{MonteCarloAverage}, \link{setoutput}
##' @examples
##' \dontrun{exceedfun <- exceedProbs(c(1.5,2,4))}
##' \dontrun{
##'     plot(lg,"exceedfun") # lg is an object of class lgcpPredict
##'                          # in which the Monte Carlo mean of 
##'                          # "exceedfun" was computed  
##'                          # see ?MonteCarloAverage and ?setoutput
##' }
##' @export

plotExceed.lgcpPredict <- function(obj,fun,nlevel=64,ask=TRUE,plotcases=FALSE,mapunderlay=NULL,alpha=1,...){
    if(!inherits(get(fun),"function")){
        stop("Argument fun should be the name of a function in memory")
    }
    if (is.null(attr(get(fun),"threshold"))){
        stop("Function fun has no threshold attribute: is fun the correct function?")
    }
    funidx <- which(obj$gridaverage[[1]]==fun)
    len <- 1
    if (inherits(obj$gridaverage[[2]][[funidx]],"list")){
        len <- length(obj$gridaverage[[2]][[funidx]])
        for (i in 1:len){
            if(plotcases){
                plotExceed.array(obj$gridaverage[[2]][[funidx]][[i]],fun=fun,xvals=obj$mcens,yvals=obj$ncens,window=obj$xyt$window,cases=obj$xyt[obj$xyt$t==obj$aggtimes[i]],nlevel=nlevel,ask=ask,main=paste("Time:",obj$aggtimes[i]),mapunderlay=mapunderlay,alpha=alpha,...)
            }
            else{
                plotExceed.array(obj$gridaverage[[2]][[funidx]][[i]],fun=fun,xvals=obj$mcens,yvals=obj$ncens,window=obj$xyt$window,nlevel=nlevel,ask=ask,main=paste("Time:",obj$aggtimes[i]),mapunderlay=mapunderlay,alpha=alpha,...)
            }
        }
    }
    else if (inherits(obj$gridaverage[[2]][[funidx]],"array")){
        warning("Labels on plot assume that this is the last time point",immediate.=TRUE)
        if(plotcases){
            plotExceed.array(obj$gridaverage[[2]][[funidx]],fun=fun,xvals=obj$mcens,yvals=obj$ncens,window=obj$xyt$window,cases=obj$xyt[obj$xyt$t==rev(obj$aggtimes)[1]],nlevel=nlevel,ask=ask,main=paste("Time:",rev(obj$aggtimes)[1]),mapunderlay=mapunderlay,alpha=alpha,...)
        }
        else{
            plotExceed.array(obj$gridaverage[[2]][[funidx]],fun=fun,xvals=obj$mcens,yvals=obj$ncens,window=obj$xyt$window,nlevel=nlevel,ask=ask,main=paste("Time:",rev(obj$aggtimes)[1]),mapunderlay=mapunderlay,alpha=alpha,...)
        }
    }
    else{
        stop("obj$gridaverage[[2]][[funidx]] should either be an array or a list object.")    
    }
}

##' gridInWindow function
##'
##' For the grid defined by x-coordinates, xvals, and y-coordinates, yvals, and an owin object W, this function just returns
##' a logical matrix M, whose [i,j] entry is TRUE if the point(xvals[i], yvals[j]) is inside the observation window. 
##'
##' @param xvals x coordinates
##' @param yvals y coordinates
##' @param win owin object
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
##' @return matrix of TRUE/FALSE, which elements of the grid are inside the observation window win
##' @export

gridInWindow <- function(xvals,yvals,win,inclusion="touching"){
    nx <- length(xvals)
    ny <- length(yvals)  
    if(inclusion=="centroid"){
        return(matrix(inside.owin(rep(xvals,ny),rep(yvals,each=nx),win),nx,ny))
    }    
    else if(inclusion=="touching"){
        return(matrix(touchingowin(xvals,yvals,win),nx,ny))
    }
    else{
        stop("Invlaid choice for argument 'inclusion'.")
    }
}



##' quantile.lgcpPredict function
##'
##' \bold{This function requires data to have been dumped to disk}: see \code{?dump2dir} and \code{?setoutput}. The routine \code{quantile.lgcpPredict} 
##' computes quantiles of functions of Y. For example, to get cell-wise quantiles of exceedance probabilities, set \code{fun=exp}. 
##' Since computign the quantiles is an expensive operation, the option to output the quantiles on a subregion of interest is also provided (by
##' setting the argument \code{inWindow}, which has a sensible default).
##'
##' @method quantile lgcpPredict
##' @param x an object of class lgcpPredict
##' @param qt a vector of the required quantiles
##' @param tidx the index number of the the time interval of interest, default is the last time point.
##' @param fun a 1-1 function (default the identity function) to be applied cell-wise to the grid. Must be able to evaluate sapply(vec,fun) for vectors vec.
##' @param inWindow an observation owin window on which to compute the quantiles, can speed up calculation. Default is x$xyt$window.
##' @param crop2parentwindow logical: whether to only compute the quantiles for cells inside x$xyt$window (the 'parent window')
##' @param startidx optional starting sample index for computing quantiles. Default is 1.   
##' @param sampcount number of samples to include in computation of quantiles after startidx. Default is all
##' @param ... additional arguments
##' @return an array, the [,,i]th slice being the grid of cell-wise quantiles, qt[i], of fun(Y), where Y is the MCMC output dumped to disk.
##' @seealso \link{lgcpPredict}, \link{dump2dir}, \link{setoutput}, \link{plot.lgcpQuantiles}
##' @export

quantile.lgcpPredict <- function(x,qt,tidx=NULL,fun=NULL,inWindow=x$xyt$window,crop2parentwindow=TRUE,startidx=1,sampcount=NULL,...){
    if (is.null(x$gridfunction$dirname)){
        stop("dump2dir not specified, MCMC output must have be dumped to disk to use this function. See ?dump2dir.")
    }
    if (length(tidx)>1){
        stop("tidx should either be NULL, or a vector of length 1")
    }
    if (is.null(fun)){
        fun <- function(a){return(a)}
    }
    fn <- paste(x$gridfunction$dirname,"simout.nc",sep="")
    ncdata <- nc_open(fn)
    datadim <- ncdata$var$simrun$varsize
    if (is.null(tidx)){
        tidx <- datadim[3]
    }
    if(!is.null(inWindow)){
        incl <- "centroid"
        if(!inherits(try(x$inclusion),"try-error")){
            if(!is.null(x$inclusion)){
                incl <- x$inclusion
            }
        }
        if (crop2parentwindow){
            grinw <- matrix(as.logical(gridInWindow(x$mcens,x$ncens,inWindow,inclusion=incl) * gridInWindow(x$mcens,x$ncens,x$xyt$window,inclusion=incl)),length(xvals(x)),length(yvals(x)))
        }
        else{
            grinw <- gridInWindow(x$mcens,x$ncens,inWindow,inclusion=incl)
        }
    }
    result <- array(dim=c(datadim[1],datadim[2],length(qt)))
    if(is.null(sampcount)){
        sampcount <- datadim[4]
    }
    pb <- txtProgressBar(min=1,max=datadim[1]*datadim[2],style=3)
    for (i in 1:datadim[1]){
        for(j in 1:datadim[2]){
            setTxtProgressBar(pb,j+(i-1)*datadim[2])
            if(!is.null(inWindow)){
                if (isTRUE(grinw[i,j])){
                    tr <- ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(i,j,tidx,startidx), count=c(1,1,1,sampcount))
                    tr <- sapply(tr,fun)
                    result[i,j,] <- quantile(tr,qt)
                }
                else{
                    result[i,j,] <- NA
                }
            }
            else{
                tr <- ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(i,j,tidx,1), count=c(1,1,1,-1))
                tr <- sapply(tr,fun)
                result[i,j,] <- quantile(tr,qt)
            }
        }
    }
    close(pb)
    nc_close(ncdata)
    attr(result,"quantiles") <- qt
    attr(result,"xcoords") <- x$mcens
    attr(result,"ycoords") <- x$ncens
    attr(result,"window") <- NULL
    if(!is.null(inWindow)){
        attr(result,"window") <- inWindow
    }
    attr(result,"ParentWindow") <- x$xyt$window
    class(result) <- c("lgcpQuantiles","array")
    return(result)   
}



##' plot.lgcpQuantiles function
##'
##' Plots \code{lgcpQuantiles} objects: output from \code{quantiles.lgcpPredict}
##'
##' @method plot lgcpQuantiles
##' @param x an object of class lgcpQuantiles
##' @param sel vector of integers between 1 and grid$len: which grids to plot. Default NULL, in which case all grids are plotted.
##' @param ask logical; if TRUE the user is asked before each plot  
##' @param crop whether or not to crop to bounding box of observation window
##' @param plotwin logical whether to plot the window attr(x,"window"), default is FALSE
##' @param ... other arguments  passed to image.plot  
##' @return grid plotting
##' This is a wrapper function for image.lgcpgrid
##' @seealso \link{quantile.lgcpPredict}
##' @examples
##' \dontrun{qtiles <- quantile(lg,qt=c(0.5,0.75,0.9),fun=exp)} 
##'                           # assumed that lg has class lgcpPredict
##' \dontrun{plot(qtiles)}
##' @export

plot.lgcpQuantiles <- function(x,sel=1:dim(x)[3],ask=TRUE,crop=TRUE,plotwin=FALSE,...){
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    for (i in sel){
        if(crop & !is.null(attr(x,"window"))){
            image.plot(attr(x,"xcoords"),attr(x,"ycoords"),x[,,i],xlim=attr(x,"window")$xrange,ylim=attr(x,"window")$yrange,sub=paste("quantile:",attr(x,"quantiles")[i]),...)
        }
        else{
            image.plot(attr(x,"xcoords"),attr(x,"ycoords"),x[,,i],xlim=attr(x,"ParentWindow")$xrange,ylim=attr(x,"ParentWindow")$yrange,sub=paste("quantile:",attr(x,"quantiles")[i]),...)
        }
        
        if (plotwin){
            if(!is.null(attr(x,"window"))){
                plot(attr(x,"window"),lwd=1,add=TRUE)
            }
        }
        plot(attr(x,"ParentWindow"),lwd=2,add=TRUE)    
    }
} 



##' identify.lgcpPredict function
##'
##' Identifies the indices of grid cells on plots of \code{lgcpPredict} objects. Can be used to identify
##' a small number of cells for further information eg trace or autocorrelation plots (provided data has been dumped to disk). On calling
##' \code{identify(lg)} for example (see code below), the user can click multiply with the left mouse button on the graphics device; once
##' the user has selected all points of interest, the right button is pressed, which returns them.
##'
##' @method identify lgcpPredict
##' @param x an object of class lgcpPredict
##' @param ... additional arguments 
##' @return a 2 x n matrix containing the grid indices of the points of interest, where n is the number of points selected via the mouse.
##' @seealso \link{lgcpPredict}, \link{loc2poly}
##' @examples
##' \dontrun{plot(lg)} # lg an lgcpPredict object
##' \dontrun{pt_indices <- identify(lg)}
##' @export

identify.lgcpPredict <- function(x,...){
    points(rep(x$mcens,x$N),rep(x$ncens,each=x$M),col=NA)
    id <- identify(x=rep(x$mcens,x$N),y=rep(x$ncens,each=x$M),plot=FALSE)
    retidx <- function(idx){
        y <- floor(idx/x$N) + 1
        x <- idx - y * x$N + x$M
        return(c(x,y))
    }
    return(t(sapply(id,retidx)))
}



##' identifygrid function
##'
##' Identifies the indices of grid cells on plots of objects.
##'
##' @param x the x grid centroids
##' @param y the y grid centroids 
##' @return a 2 x n matrix containing the grid indices of the points of interest, where n is the number of points selected via the mouse.
##' @seealso \link{lgcpPredict}, \link{loc2poly}, \link{identify.lgcpPredict}
##' @export

identifygrid <- function(x,y){
    M <- length(x)
    N <- length(y)
    points(rep(x,N),rep(y,each=M),col=NA)
    id <- identify(rep(x,N),y=rep(y,each=M),plot=FALSE)
    retidx <- function(idx){
        y <- floor(idx/N) + 1
        x <- idx - y * N + M
        return(c(x,y))
    }
    return(t(sapply(id,retidx)))
}



##' loc2poly function
##'
##' Converts a polygon selected via the mouse in a graphics window into an polygonal owin object. (Make sure the x and y scales are correct!)
##' Points must be selected traversing the required window in one direction (ie either clockwise, or anticlockwise), points must not be overlapping.
##' Select the sequence of edges via left mouse button clicks and store the polygon with a right click. 
##'
##' @param n the maximum number of points to locate
##' @param type same as argument type in function locator. see ?locator. Default draws lines
##' @param col colour of lines/points
##' @param ... other arguments to pass to locate
##' @return a polygonal owin object
##' @seealso \link{lgcpPredict}, \link{identify.lgcpPredict}
##' @examples
##' \dontrun{plot(lg)} # lg an lgcpPredict object
##' \dontrun{subwin <- loc2poly())}
##' @export

loc2poly <- function(n=512,type="l",col="black",...){
    l <- locator(n=n,type=type,col=col,...)
    win <- try(owin(poly=list(x=l$x,y=l$y)),silent=TRUE)
    if(class(win)=="try-error"){
        win <- try(owin(poly=list(x=rev(l$x),y=rev(l$y))),silent=TRUE) # try the points ordering in the other direction
        if(class(win)=="try-error"){
            owin(poly=list(x=l$x,y=l$y)) # if there is still a problem with the chosen set of points, then this should print out the underlying issue
        }
        else{
            return(win)
        }
    }
    else{
        return(win)    
    }
}



##' expectation.lgcpPredict function
##'
##' \bold{This function requires data to have been dumped to disk}: see \code{?dump2dir} and \code{?setoutput}. This function computes the 
##' Monte Carlo Average of a function where data from a run of \code{lgcpPredict} has been dumped to disk.
##'
##' A Monte Carlo Average is computed as:
##' \deqn{E_{\pi(Y_{t_1:t_2}|X_{t_1:t_2})}[g(Y_{t_1:t_2})] \approx \frac1n\sum_{i=1}^n g(Y_{t_1:t_2}^{(i)})}{E_{\pi(Y_{t_1:t_2}|X_{t_1:t_2})}[g(Y_{t_1:t_2})] \approx \frac1n\sum_{i=1}^n g(Y_{t_1:t_2}^{(i)})}
##' where \eqn{g}{g} is a function of interest, \eqn{Y_{t_1:t_2}^{(i)}}{Y_{t_1:t_2}^{(i)}} is the \eqn{i}{i}th retained sample from the target  
##' and \eqn{n}{n} is the total number of retained iterations. For example, to compute the mean of \eqn{Y_{t_1:t_2}}{Y_{t_1:t_2}} set,
##' \deqn{g(Y_{t_1:t_2}) = Y_{t_1:t_2},}{g(Y_{t_1:t_2}) = Y_{t_1:t_2},}
##' the output from such a Monte Carlo average would be a set of \eqn{t_2-t_1}{t_2-t_1} grids, each cell of which 
##' being equal to the mean over all retained iterations of the algorithm (NOTE: this is just an example computation, in
##' practice, there is no need to compute the mean on line explicitly, as this is already done by default in \code{lgcpPredict}).
##'
##' @method expectation lgcpPredict
##' @param obj an object of class lgcpPredict
##' @param fun a function accepting a single argument that returns a numeric vector, matrix or array object
##' @param maxit Not used in ordinary circumstances. Defines subset of samples over which to compute expectation. Expectation is computed using information from iterations 1:maxit, where 1 is the first non-burn in iteration dumped to disk. 
##' @param ... additional arguments  
##' @return the expectated value of that function
##' @seealso \link{lgcpPredict}, \link{dump2dir}, \link{setoutput}
##' @export

expectation.lgcpPredict <- function(obj,fun,maxit=NULL,...){
    if (is.null(obj$gridfunction$dirname)){
        stop("dump2dir not specified, MCMC output must have be dumped to disk to use this function.  See ?dump2dir.")
    }
    fn <- paste(obj$gridfunction$dirname,"simout.nc",sep="")
    ncdata <- nc_open(fn)
    if (!is.null(maxit)){
        if (maxit<2 | maxit>ncdata$dim$iter$len){
            stop(paste("maxit must be between 2 and",ncdata$dim$iter$len))
        }
        endloop <- maxit
    }
    else{
        endloop <- ncdata$dim$iter$len
    }
    Y <- lgcpgrid(ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(1,1,1,1), count=c(-1,-1,-1,1))) # first "simulation"
    result <- tryCatch(lapply(Y$grid,fun),finally=nc_close(ncdata))
    pb <- txtProgressBar(min=1,max=endloop,style=3)
    setTxtProgressBar(pb,1)
    for (i in 2:endloop){
        ncdata <- nc_open(fn)
        Y <- lgcpgrid(ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(1,1,1,i), count=c(-1,-1,-1,1)))
        result <- tryCatch(add.list(result,lapply(Y$grid,fun)),finally=nc_close(ncdata))
        setTxtProgressBar(pb,i)
    }
    close(pb)
    for (i in 1:length(result)){
        result[[i]] <- result[[i]] / endloop
    }
    return(result)
}


##' expectation.lgcpPredictSpatialOnlyPlusParameters function
##'
##' \bold{This function requires data to have been dumped to disk}: see \code{?dump2dir} and \code{?setoutput}. This function computes the 
##' Monte Carlo Average of a function where data from a run of \code{lgcpPredict} has been dumped to disk.
##'
##' @method expectation lgcpPredictSpatialOnlyPlusParameters
##' @param obj an object of class lgcpPredictSpatialOnlyPlusParameters
##' @param fun a function with arguments 'Y', 'beta', 'eta', 'Z' and 'otherargs'. See vignette("Bayesian_lgcp") for an example
##' @param maxit Not used in ordinary circumstances. Defines subset of samples over which to compute expectation. Expectation is computed using information from iterations 1:maxit, where 1 is the first non-burn in iteration dumped to disk.
##' @param ... additional arguments
##' @usage "expectation(obj,fun,maxit=NULL,...)"
##' @return the expectated value of that function
##' @export

expectation.lgcpPredictSpatialOnlyPlusParameters <- function(obj,fun,maxit=NULL,...){
    if (is.null(obj$gridfunction$dirname)){
        stop("dump2dir not specified, MCMC output must have be dumped to disk to use this function.  See ?dump2dir.")
    }
    fn <- paste(obj$gridfunction$dirname,"simout.nc",sep="")
    ncdata <- nc_open(fn)
    if (!is.null(maxit)){
        if (maxit<2 | maxit>ncdata$dim$iter$len){
            stop(paste("maxit must be between 2 and",ncdata$dim$iter$len))
        }
        endloop <- maxit
    }
    else{
        endloop <- ncdata$dim$iter$len
    }
    Y <- lgcpgrid(ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(1,1,1,1), count=c(-1,-1,-1,1))) # first "simulation"
    result <- tryCatch(lapply(Y$grid,fun,eta=obj$etarec[1,,drop=FALSE],beta=obj$betarec[1,,drop=FALSE],Z=obj$Z,otherargs=obj),finally=nc_close(ncdata))
    pb <- txtProgressBar(min=1,max=endloop,style=3)
    setTxtProgressBar(pb,1)
    for (i in 2:endloop){
        ncdata <- nc_open(fn)
        Y <- lgcpgrid(ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(1,1,1,i), count=c(-1,-1,-1,1)))
        result <- tryCatch(add.list(result,lapply(Y$grid,fun,eta=obj$etarec[i,,drop=FALSE],beta=obj$betarec[i,,drop=FALSE],Z=obj$Z,otherargs=obj)),finally=nc_close(ncdata))
        setTxtProgressBar(pb,i)
    }
    close(pb)
    for (i in 1:length(result)){
        result[[i]] <- result[[i]] / endloop
    }
    return(result)
}



##' extract.lgcpPredict function
##'
##' \bold{This function requires data to have been dumped to disk}: see \code{?dump2dir} and \code{?setoutput}. \code{extract.lgcpPredict}  
##' extracts chunks of data that have been dumped to disk. The subset of data can either be specified using an (x,y,t,s) box or (window,t,s) region
##' where window is a polygonal subregion of interest.
##'
##' @method extract lgcpPredict
##' @param obj an object of class lgcpPredict
##' @param x range of x-indices: vector (eg c(2,4)) corresponding to desired subset of x coordinates. If equal to -1, then all cells in this dimension are extracted
##' @param y range of y-indices as above
##' @param t range of t-indices: time indices of interest
##' @param s range of s-indices ie the simulation indices of interest
##' @param inWindow an observation owin window over which to extract the data (alternative to specifying x and y).
##' @param crop2parentwindow logical: whether to only extract cells inside obj$xyt$window (the 'parent window') 
##' @param ... additional arguments
##' @return extracted array
##' @seealso \link{lgcpPredict}, \link{loc2poly}, \link{dump2dir}, \link{setoutput}
##' @export

extract.lgcpPredict <- function(obj,x=NULL,y=NULL,t,s=-1,inWindow=NULL,crop2parentwindow=TRUE,...){
    if (is.null(obj$gridfunction$dirname)){
        stop("dump2dir not specified, MCMC output must have be dumped to disk to use this function.  See ?dump2dir.")
    }
    if(all(is.null(c(x,y,inWindow))) | (!is.null(inWindow)&any(!is.null(c(x,y))))){
        stop("either x and y OR inWindow must be given")
    }
    fn <- paste(obj$gridfunction$dirname,"simout.nc",sep="")
    ncdata <- nc_open(fn)
    
    if(is.null(inWindow)){
        if(((length(x)>1)&any(x==-1)) | ((length(y)>1)&any(y==-1)) | ((length(t)>1)&any(t==-1)) | ((length(s)>1)&any(s==-1))){
            stop("Error in one or more of x, y, t, or s: cannot mix -1's with positive indices")
        }
    
        if(x[1]!=-1 & any(x<1 | x>ncdata$var$simrun$varsize[1])){
            stop(paste("Supplied value(s) of x must be between 1 and",ncdata$var$simrun$varsize[1],"OR equal to -1, see ?extract.lgcpPredict"))
        }
        if(y[1]!=-1 & any(y<1 | y>ncdata$var$simrun$varsize[2])){
            stop(paste("Supplied value(s) of y must be between 1 and",ncdata$var$simrun$varsize[2],"OR equal to -1, see ?extract.lgcpPredict"))
        }
        if(t[1]!=-1 & any(t<1 | t>ncdata$var$simrun$varsize[3])){
            stop(paste("Supplied value(s) of t must be between 1 and",ncdata$var$simrun$varsize[3],"OR equal to -1, see ?extract.lgcpPredict"))
        }
        if(s[1]!=-1 & any(s<1 | s>ncdata$var$simrun$varsize[4])){
            stop(paste("Supplied value(s) of s must be between 1 and",ncdata$var$simrun$varsize[4],"OR equal to -1, see ?extract.lgcpPredict"))
        }    
    }    
    
    if (is.null(inWindow)){
        startidx <- c(min(x),min(y),min(t),min(s))
        mo <- (startidx==-1)
        startidx[mo] <- 1
        ct <- rep(-1,4)
        ct[which(!mo)] <- c(max(x),max(y),max(t),max(s))[which(!mo)] - c(min(x),min(y),min(t),min(s))[which(!mo)] + 1
        data <- ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=startidx, count=ct)
        attr(data,"mode") <- "block"
        attr(data,"xcoords") <- obj$mcens[min(x):max(x)]
        attr(data,"ycoords") <- obj$ncens[min(y):max(y)]
        class(data) <- c("lgcpExtract","array")
    }
    else{
        incl <- "centroid"
        if(!inherits(try(obj$inclusion),"try-error")){
            if(!is.null(obj$inclusion)){
                incl <- obj$inclusion
            }
        }
        if (crop2parentwindow){
            grinw <- matrix(as.logical(gridInWindow(obj$mcens,obj$ncens,inWindow,inclusion=incl) * gridInWindow(obj$mcens,obj$ncens,obj$xyt$window,inclusion=incl)),length(xvals(obj)),length(yvals(obj)))
        }
        else{
            grinw <- gridInWindow(obj$mcens,obj$ncens,inWindow,inclusion=incl)
        }
        n <- sum(grinw)
        idx <- which(grinw,arr.ind=TRUE)
        startidx <- c(min(t),min(s))
        mo <- (startidx==-1)
        startidx[mo] <- 1
        ct <- rep(-1,2)
        ct[which(!mo)] <- c(max(t),max(s))[which(!mo)] - c(min(t),min(s))[which(!mo)] + 1
        vs <- ncdata$var[[1]]$varsize[3:4]
        d <- ct
        d[d==-1] <- vs[d==-1]
        data <- array(dim=c(n,d))
        for (i in 1:n){
            data[i,,] <- ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(idx[i,],startidx), count=c(1,1,ct))
        }
        attr(data,"xcoords") <- obj$mcens
        attr(data,"ycoords") <- obj$ncens
        attr(data,"mode") <- "inWindow"
        attr(data,"window") <- inWindow
        attr(data,"mask") <- grinw # can be used to translate returned vector back into 2 dimensional object (idea is to save storage space)
        class(data) <- c("lgcpExtractInWindow","lgcpExtract","array")
    }
    nc_close(ncdata)   
    return(data)
}



##' add.list function
##'
##' This function adds the elements of two list objects together and returns the result in another list object.
##'
##' @param list1 a list of objects that could be summed using "+"
##' @param list2 a list of objects that could be summed using "+"
##' @return a list with ith entry the sum of list1[[i]] and list2[[i]]
##' @export

add.list <- function(list1,list2){
    l1 <- length(list1)
    if(length(list2)!=l1){
        stop("Lists must be of the same length")
    }
    lst <- list()
    if(l1>0){
        for (i in 1:l1){
            lst[[i]] <- list1[[i]] + list2[[i]]
        }
    }
    return(lst)    
}



##' multiply.list function
##'
##' This function multiplies the elements of two list objects together and returns the result in another list object.
##'
##' @param list1 a list of objects that could be summed using "+"
##' @param list2 a list of objects that could be summed using "+"
##' @return a list with ith entry the sum of list1[[i]] and list2[[i]]
##' @export

multiply.list <- function(list1,list2){
    l1 <- length(list1)
    if(length(list2)!=l1){
        stop("Lists must be of the same length")
    }
    lst <- list()
    if(l1>0){
        for (i in 1:l1){
            lst[[i]] <- list1[[i]] * list2[[i]]
        }
    }
    return(lst)    
}



##' smultiply.list function
##'
##' This function multiplies each element of a list by a scalar constant.
##'
##' @param list a list of objects that could be summed using "+"
##' @param const a numeric constant
##' @return a list with ith entry the scalar multiple of const * list[[i]]
##' @export

smultiply.list <- function(list,const){
    lst <- list()
    for (i in 1:length(list)){
        lst[[i]] <- list[[i]] * const
    }
    return(lst)    
}



##' showGrid function
##'
##' Generic method for displaying the FFT grid used in computation.
##'
##' @param x an object
##' @param ... additional arguments
##' @return generic function returning method showGrid
##' @seealso \link{showGrid.default}, \link{showGrid.lgcpPredict}, \link{showGrid.stppp}
##' @export

showGrid <- function(x,...){
    UseMethod("showGrid")
}



##' showGrid.default function
##'
##' Default method for printing a grid to a screen. Arguments are vectors giving the x any y coordinates of the
##' centroids.
##'
##' @method showGrid default
##' @param x an vector of grid values for the x coordinates
##' @param y an vector of grid values for the y coordinates
##' @param ... additional arguments passed to points
##' @return plots grid centroids on the current graphics device
##' @seealso \link{showGrid.lgcpPredict}, \link{showGrid.stppp}
##' @export

showGrid.default <- function(x,y,...){
     points(cbind(rep(x,length(y)),rep(y,each=length(x))),...)
}



##' showGrid.lgcpPredict function
##'
##' This function displays the FFT grid used on a plot of an \code{lgcpPredict} object.
##' First plot the object using for example \code{plot(lg)}, where \code{lg} is an object
##' of class \code{lgcpPredict}, then for any of the plots produced, a call to
##' \code{showGrid(lg,pch=="+",cex=0.5)} will display the centroids of the FFT grid.
##'
##' @method showGrid lgcpPredict
##' @param x an object of class lgcpPredict
##' @param ... additional arguments  passed to points
##' @return plots grid centroids on the current graphics device
##' @seealso \link{lgcpPredict}, \link{showGrid.default}, \link{showGrid.stppp}
##' @export

showGrid.lgcpPredict <- function(x,...){
    points(cbind(rep(xvals(x),x$N),rep(yvals(x),each=x$M)),...)
}



##' showGrid.stppp function
##'
##' If an stppp object has been created via simulation, ie using the function \code{lgcpSim}, then
##' this function will display the grid centroids that were used in the simulation
##'
##' @method showGrid stppp
##' @param x an object of class stppp. Note this function oly applies to SIMULATED data.
##' @param ... additional arguments  passed to points
##' @return plots grid centroids on the current graphics device. FOR SIMULATED DATA ONLY.
##' @seealso \link{lgcpSim}, \link{showGrid.default}, \link{showGrid.lgcpPredict}
##' @examples
##' \dontrun{xyt <- lgcpSim()}
##' \dontrun{plot(xyt)}
##' \dontrun{showGrid(xyt,pch="+",cex=0.5)}
##' @export

showGrid.stppp <- function(x,...){
    if (is.null(attr(x,"xvals"))){
        stop("No grid available")
    }    
    xv <- attr(x,"xvals")
    yv <- attr(x,"yvals")
    points(cbind(rep(xv,length(yv)),rep(yv,each=length(xv))),...)
}


##' autocorr function
##'
##' \bold{This function requires data to have been dumped to disk}: see \code{?dump2dir} and \code{?setoutput}. The routine \code{autocorr.lgcpPredict} 
##' computes cellwise selected autocorrelations of Y.
##' Since computing the quantiles is an expensive operation, the option to output the quantiles on a subregion of interest is also provided (by
##' setting the argument \code{inWindow}, which has a sensible default).
##'
##' @param x an object of class lgcpPredict
##' @param lags a vector of the required lags
##' @param tidx the index number of the the time interval of interest, default is the last time point.
##' @param inWindow an observation owin window on which to compute the autocorrelations, can speed up calculation. Default is x$xyt$window, set to NULL for full grid.
##' @param crop2parentwindow logical: whether to only compute autocorrelations for cells inside x$xyt$window (the 'parent window')
##' @param ... additional arguments
##' @return an array, the [,,i]th slice being the grid of cell-wise autocorrelations.
##' @seealso \link{lgcpPredict}, \link{dump2dir}, \link{setoutput}, \link{plot.lgcpAutocorr}, \link{ltar}, \link{parautocorr}, \link{traceplots}, \link{parsummary}, \link{textsummary},
##' \link{priorpost}, \link{postcov}, \link{exceedProbs}, \link{betavals}, \link{etavals}
##' @export

autocorr <- function(x,lags,tidx=NULL,inWindow=x$xyt$window,crop2parentwindow=TRUE,...){
    if (is.null(x$gridfunction$dirname)){
        stop("dump2dir not specified, MCMC output must have be dumped to disk to use this function. See ?dump2dir.")
    }
    if (length(tidx)>1){
        stop("tidx should either be NULL, or a vector of length 1")
    }
    fn <- paste(x$gridfunction$dirname,"simout.nc",sep="")
    ncdata <- nc_open(fn)
    datadim <- ncdata$var$simrun$varsize
    if (is.null(tidx)){
        tidx <- datadim[3]
    }
    if(!is.null(inWindow)){
        incl <- "centroid"
        if(!inherits(try(x$inclusion),"try-error")){
            if(!is.null(x$inclusion)){
                incl <- x$inclusion
            }
        }
        if (crop2parentwindow){
            if(identical(x$xyt$window,inWindow)){
                dm <- dim(x$cellInside)
                grinw <- array(as.logical(x$cellInside),dim=dm)
            }
            else{
                grinw <- matrix(as.logical(gridInWindow(x$mcens,x$ncens,inWindow,inclusion=incl) * gridInWindow(x$mcens,x$ncens,x$xyt$window,inclusion=incl)),length(xvals(x)),length(yvals(x)))
            }
        }
        else{
            if(identical(x$xyt$window,inWindow)){
                dm <- dim(x$cellInside)
                grinw <- array(as.logical(x$cellInside),dim=dm)
            }
            else{
                grinw <- gridInWindow(x$mcens,x$ncens,inWindow,inclusion=incl)
            }
        }
    }
    result <- array(dim=c(datadim[1],datadim[2],length(lags)))
    sampcount <- datadim[4]
    pb <- txtProgressBar(min=1,max=datadim[1]*datadim[2],style=3)
    trigger <- TRUE
    for (i in 1:datadim[1]){
        for(j in 1:datadim[2]){
            setTxtProgressBar(pb,j+(i-1)*datadim[2])
            if(!is.null(inWindow)){
                if (isTRUE(grinw[i,j])){
                    tr <- ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(i,j,tidx,1), count=c(1,1,1,-1))
                    result[i,j,] <- acf(tr,plot=FALSE)$acf[lags+1]
                }
                else{
                    result[i,j,] <- NA
                }
            }
            else{
                tr <- ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(i,j,tidx,1), count=c(1,1,1,-1))
                result[i,j,] <- acf(tr,plot=FALSE)$acf[lags+1]
            }
            
            if(trigger){
                ltst <- length(acf(ncvar_get(nc=ncdata, varid=ncdata$var[[1]], start=c(1,1,tidx,1), count=c(1,1,1,-1)),plot=FALSE)$acf)
                tst <- (lags+1)>ltst
                if(any(tst)){
                    stop(paste("Cannot return lag",ltst,"or above."))
                }
                trigger <- FALSE
            }
        }
    }
    close(pb)
    nc_close(ncdata)
    attr(result,"lags") <- lags
    attr(result,"xcoords") <- x$mcens
    attr(result,"ycoords") <- x$ncens
    attr(result,"window") <- NULL
    if(!is.null(inWindow)){
        attr(result,"window") <- inWindow
    }
    attr(result,"ParentWindow") <- x$xyt$window
    class(result) <- c("lgcpAutocorr","array")
    return(result)   
}

##' plot.lgcpAutocorr function
##'
##' Plots \code{lgcpAutocorr} objects: output from \code{autocorr}
##'
##' @method plot lgcpAutocorr
##' @param x an object of class lgcpAutocorr
##' @param sel vector of integers between 1 and grid$len: which grids to plot. Default NULL, in which case all grids are plotted.
##' @param ask logical; if TRUE the user is asked before each plot  
##' @param crop whether or not to crop to bounding box of observation window
##' @param plotwin logical whether to plot the window attr(x,"window"), default is FALSE
##' @param ... other arguments  passed to image.plot  
##' @return a plot
##' @seealso \link{autocorr}
##' @examples
##' \dontrun{ac <- autocorr(lg,qt=c(1,2,3))} 
##'                           # assumes that lg has class lgcpPredict
##' \dontrun{plot(ac)}
##' @export

plot.lgcpAutocorr <- function(x,sel=1:dim(x)[3],ask=TRUE,crop=TRUE,plotwin=FALSE,...){
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    for (i in sel){
        if(crop & !is.null(attr(x,"window"))){
            image.plot(attr(x,"xcoords"),attr(x,"ycoords"),x[,,i],xlim=attr(x,"window")$xrange,ylim=attr(x,"window")$yrange,sub=paste("lag:",attr(x,"lags")[i]),...)
        }
        else{
            image.plot(attr(x,"xcoords"),attr(x,"ycoords"),x[,,i],xlim=attr(x,"ParentWindow")$xrange,ylim=attr(x,"ParentWindow")$yrange,sub=paste("lag:",attr(x,"lags")[i]),...)
        }
        
        if (plotwin){
            if(!is.null(attr(x,"window"))){
                plot(attr(x,"window"),add=TRUE)
            }
        }
        plot(attr(x,"ParentWindow"),add=TRUE)    
    }
} 


##' touchingowin function
##'
##' A function to compute which cells are touching an owin or spatial polygons object
##'
##' @param x grid centroids in x-direction note this will be expanded into a GRID of (x,y) values in the function
##' @param y grid centroids in y-direction note this will be expanded into a GRID of (x,y) values in the function
##' @param w an owin or SpatialPolygons object 
##' @return vector of TRUE or FALSE according to whether the cell
##' @export

touchingowin <- function(x, y, w){
    if (inherits(w,"owin")){
        w <- as(w,"SpatialPolygons")
    }
    gri <- grid2spoly(x,y)
    int <- !gDisjoint(w,gri,byid=TRUE) & !gTouches(w,gri,byid=TRUE) # cells that have some internal points in common
    return(apply(int,1,any))
} 


##' grid2spoly function
##'
##' A function to convert a regular (x,y) grid of centroids into a SpatialPolygons object
##'
##' @param xgrid vector of x centroids (equally spaced)
##' @param ygrid vector of x centroids (equally spaced)
##' @param proj4string proj 4 string: specify in the usual way
##' @return a SpatialPolygons object
##' @export

grid2spoly <- function(xgrid,ygrid,proj4string=CRS(as.character(NA))){
    m <- length(xgrid)
    n <- length(ygrid)
    spts <- SpatialPixels(SpatialPoints(cbind(rep(xgrid,n),rep(ygrid,each=m))),proj4string=proj4string)
    return(as(spts,"SpatialPolygons"))
}


##' plot.lgcpZmat function
##'
##' A function to plot lgcpZmat objects
##'
##' @method plot lgcpZmat
##' @param x an lgcpZmat object, see ?getZmat
##' @param ask graphical parameter ask, see ?par
##' @param pow power parameter, raises the image values to this power (helps with visualisation, default is 1.)
##' @param main title for plot, default is null which gives an automatic title to the plot (the name of the covariate)
##' @param misscol colour to identify imputed grid cells, default is yellow
##' @param obswin optional observation window to add to plot using plot(obswin).
##' @param ... other paramters 
##' @return a sequence of plots of the interpolated covariate values
##' @export

plot.lgcpZmat <- function(x,ask=TRUE,pow=1,main=NULL,misscol="black",obswin=NULL,...){
    MAIN <- main
    mcens <- attr(x,"mcens")
    ncens <- attr(x,"ncens")
    M <- attr(x,"M")
    N <- attr(x,"N")
    cellInside <- attr(x,"cellInside")
    cellInside[cellInside==0] <- NA    
    
    npar <- ncol(x)
    parn <- colnames(x)
    for(i in 1:npar){
        if(is.null(MAIN)){
            main <- parn[i]
        }
        
        blah <- try(image.plot(mcens,ncens,(matrix(x[,i],M,N)*cellInside)^pow,main=main,...))
        
        if(!inherits(blah,"try-error")){        
            if(any(attr(x,"missingind")==1)|is.null(attr(x,"missingind"))){ # is.null(...) for backward compatibility
                egr <- expand.grid(mcens,ncens)
                idx <- which(attr(x,"missingind")==1)
                points(egr[idx,],col=misscol,pch="+")
            }
        }
        
        if(!is.null(obswin)){
            plot(obswin,add=TRUE)
        }
        
        if(ask){
            cat ("Press [enter] to continue")
            line <- readline()
        }
    }
    
    
}



##' lgcpvignette function
##'
##' Display the introductory vignette for the lgcp package. 
##'
##' @return displays the vignette by calling browseURL
##' @export

lgcpvignette <- function(){
    browseURL("http://www.jstatsoft.org/v52/i04/paper") 
}

##' lgcpbayes function
##'
##' Display the introductory vignette for the lgcp package. 
##'
##' @return displays the vignette by calling browseURL
##' @export

lgcpbayes <- function(){
    cat("This vignette is currently under review. Please contact the package author (b.taylor1-at-lancaster.ac.uk) to obtain a preprint.")
    #browseURL("... blah ...") 
}
