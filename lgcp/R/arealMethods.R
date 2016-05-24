##' stapp function
##'
##' Generic function for space-time aggregated point-process data 
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method stapp
##' @export

stapp <- function(obj,...){
    UseMethod("stapp")
}



##' stapp.SpatialPolygonsDataFrame function
##'
##' Construct a space-time aggregated point-process (stapp) object from a SpatialPolygonsDataFrame (along with some other info)
##'
##' @method stapp SpatialPolygonsDataFrame
##' @param obj an SpatialPolygonsDataFrame object
##' @param counts a (length(t) by N) matrix containing aggregated case counts for each of the geographical regions defined by the SpatialPolygonsDataFrame, where N is the number of regions
##' @param t vector of times, for each element of t there should correspond a column in the matrix 'counts'
##' @param tlim vector giving the upper and lower bounds of the temporal observation window
##' @param window the observation window, of class owin, see ?owin   
##' @param ... additional arguments
##' @return an object of class stapp
##' @export

stapp.SpatialPolygonsDataFrame <- function(obj,counts,t,tlim,window,...){
    if (is.null(obj$population)){
        stop("SpatialPolygonsDataFrame must have non-null obj$population: giving regional population.")
    }
    if (length(t)!=dim(counts)[2]){
        stop(paste("Vector t must have length",dim(counts)[2]))
    }
    if (length(tlim)!=2){
        stop("tlim must be a vector of length exactly 2 and sorted in ascending order")
    }
    if (tlim[1] >= tlim[2]){
        stop("Vector tlim must be sorted in ascending order")
    }
    if (any(as.integer(t)>as.integer(tlim[2])) | any(as.integer(t)<as.integer(tlim[1]))){
        stop("There is a t such that as.integer(t) is outside range of as.integer(tlim). Please amend tlim so all observations fall in observation time window.")
    }
    
    obj$counts <- counts # first tack on the count data   
    
    ret <- list()
    ret$spdf <- obj
    ret$t <- t
    ret$tlim <- tlim 
    ret$window <- window 
    class(ret) <- "stapp"
    return(ret) 
}



##' stapp.list function
##'
##' A wrapper function for stapp.SpatialPolygonsDataFrame
##'
##' Construct a space-time aggregated point-process (stapp) object from a list object. The first element of the list 
##' should be a SpatialPolygonsDataFrame, the second element of the list a counts matrix, the third element of the list a
##' vector of times, the fourth element a vector giving the bounds of the temporal observation window and the fifth element
##' a spatstat owin object giving the spatial observation window.
##'
##' @method stapp list
##' @param obj an list object as described above, see ?stapp.SpatialPolygonsDataFrame for further details on the requirements of the list
##' @param ... additional arguments
##' @return an object of class stapp
##' @export

stapp.list <- function(obj,...){
    stapp.SpatialPolygonsDataFrame(obj[[1]],obj[[2]],obj[[3]],obj[[4]],obj[[5]])
}



##' SpatialPolygonsDataFrame.stapp function
##'
##' A function to return the SpatialPolygonsDataFrame part of an stapp object
##'
##' @param from stapp object
##' @return an object of class SpatialPolygonsDataFrame
##' @export

SpatialPolygonsDataFrame.stapp <- function(from){
    return(from$spdf)
}



##' print.stapp function
##'
##' Print method for stapp objects
##'
##' @method print stapp
##' @param x an object of class stapp
##' @param printhead whether or not to print the head of the counts matrix
##' @param ... additional arguments
##' @return prints the stapp object x
##' @export 

print.stapp <- function(x,printhead=TRUE,...){
    cat("Space-time aggregated point pattern\n")
    cat("  Aggregate times : ",x$t,"\n")
    cat("      Time Window : [",x$tlim[1],",",x$tlim[2],"]\n")
    cat(" Aggregate counts :\n")
    if(printhead){
        print(head(slot(x$spdf,"data")$counts))
    }
    else{
        cat("    A matrix of size [",dim(slot(x$spdf,"data")$counts)[1],",",dim(slot(x$spdf,"data")$counts)[2],"]\n")
    }
    print(x$window)
}


##' as.owin.stapp function
##'
##' A function to extract the SpatialPolygons part of W and return it as an owin object.
##'

##' @method as.owin stapp
##' @param W see ?as.owin
##' @param ... see ?as.owin
##' @param fatal see ?as.owin
##' @return an owin object
##' @export

as.owin.stapp <- function (W,...,fatal=TRUE){ 
    return(as.owin(W$spdf,...,fatal=fatal))
}



##' as.stppp function
##'
##' Generic function for converting to stppp objects
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method as.stppp
##' @export

as.stppp <- function(obj,...){
    UseMethod("as.stppp")
}



##' as.stppp.stapp function
##'
##' A function to convert stapp objects to stppp objects for use in lgcpPredict. The regional counts in the stapp object are
##' assigned a random location within each areal region proportional to a population density (if that is available) else
##' the counts are distributed uniformly across the observation windows. 
##'
##' @method as.stppp stapp
##' @param obj an object of class stapp
##' @param popden a 'spatialAtRisk' of sub-class 'fromXYZ' object representing the population density, or for better results, lambda(s) can also
##' be used here. Cases are distributed across the spatial region according to popden. NULL by default, which has the effect of assigning counts uniformly.
##' @param n if popden is NULL, then this parameter controls the resolution of the uniform. Otherwise if popden is of class 'fromFunction', it controls the size of the imputation grid used for sampling. Default is 100.
##' @param dmin If any reginal counts are missing, then a set of polygonal 'holes' in the observation window will be computed for each. dmin is the parameter used to control the simplification of these holes (see ?simplify.owin). default is zero.
##' @param check logical. If any reginal counts are missing, then roughly speaking, check specifies whether to check the 'holes'. 
##' @param ... additional arguments
##' @return ...
##' @export


as.stppp.stapp <- function(obj,popden=NULL,n=100,dmin=0,check=TRUE,...){
    verifyclass(popden,"fromXYZ")
    #verifyclass(owinlist,"owinlist") no longer needed

    if (is.null(popden)){
        bb <- obj$spdf@bbox
        popden <- spatialAtRisk.fromXYZ(X=seq(bb[1,1],bb[1,2],length.out=n),Y=seq(bb[2,1],bb[2,2],length.out=n),Zm=matrix(1,n,n))
    }

    popden$Zm[is.na(popden$Zm)] <- 0     
    
    xgrid <- xvals(popden)
    ygrid <- yvals(popden)
    nx <- length(xgrid)
    ny <- length(ygrid)    
    xls <- rep(xgrid,ny)
    yls <- rep(ygrid,each=nx)     
    xdiv <- xgrid[2] - xgrid[1]
    ydiv <- ygrid[2] - ygrid[1]   
    
    #cellInside <- lapply(owinlist,function(w){inside.owin(x=xls,y=yls,w=w)}) no longer needed
    #EJP: olay <- overlay(SpatialPoints(cbind(xls,yls)),obj$spdf)
    olay <- over(SpatialPoints(cbind(xls,yls)), geometry(obj$spdf))
    olay[is.na(olay)] <- 0
    
    rcounts <- obj$spdf@data$counts
    popvec <- popden$Zm
    
    x <- c()
    y <- c()
    t <- c()
    for (i in 1:dim(rcounts)[1]){
        rcts <- rcounts[i,] # total counts in region
        ncts <- sum(rcts,na.rm=TRUE)
        if (ncts==0){
                next
            } 
        #idx <- sample(1:(nx*ny),ncts,replace=TRUE,prob=cellInside[[i]]*popvec) no longer needed
        idx <- sample(1:(nx*ny),ncts,replace=TRUE,prob=matrix(olay==i,nx,ny)*popvec)
        ncts <- length(idx)
        newx <- xls[idx] + runif(ncts,-0.5*xdiv,0.5*xdiv)
        newy <- yls[idx] + runif(ncts,-0.5*ydiv,0.5*ydiv)
        iow <- inside.owin(x=newx,y=newy,w=obj$window)
        test <- !all(iow)
        pidx <- which(!iow) # pidx==problem index: which of the new points are not in the observation window
        while(test){
            newx[pidx] <- xls[idx][pidx] + runif(length(pidx),-0.5*xdiv,0.5*xdiv)
            newy[pidx] <- yls[idx][pidx] + runif(length(pidx),-0.5*ydiv,0.5*ydiv)
            iow[pidx] <- inside.owin(x=newx[pidx],y=newy[pidx],w=obj$window)
            test <- !all(iow)
            pidx <- which(!iow)
        } 
        x <- c(x,newx)
        y <- c(y,newy)
        rctssub <- rcts
        rctssub[is.na(rcts)] <- 0
        t <- c(t,rep(obj$t,rctssub))
    }
    ord <- order(t) # put them in time order (for no particular reason)
    t <- t[ord]
    x <- x[ord]
    y <- y[ord]
    pointpat <- ppp(x=x,y=y,window=obj$window)   
    xyt <- stppp(pointpat,t=t,tlim=obj$tlim)    
    #attr(xyt,"cellInsideList") <- lapply(cellInside,function(M){matrix(M,nx,ny)}) no longer needed
    attr(xyt,"overlay") <- olay
    if (any(is.na(rcounts))){
        cat("... Missing data in regional counts, computing polygon holes...\n")
        ss <- apply(rcounts,1,function(x){any(is.na(x))})
        owlss <- as.owinlist(obj$spdf,dmin=dmin,check=check,subset=ss)
        owl <- rep(list(obj$window),dim(rcounts)[2])
        if(requireNamespace("gpclib",quietly=TRUE)){
            spatstat.options(gpclib=TRUE)
        }
        else{
            stop("This function requires the gpclib library to be installed. Please contact the lgcp package maintaner for further information.")
        }
        ind <- apply(rcounts,2,function(x){any(is.na(x))})
        for(i in 1:dim(rcounts)[2]){
            if(!ind[i]){
                next
            }
            regidx <- which(is.na(rcounts[,i]))
            for(j in 1:length(regidx)){
                owl[[i]] <- setminus.owin(owl[[i]],owlss[[regidx[j]]])
            }
        }
        attr(xyt,"owinlist") <- owl
    }
    return(xyt)
}

##' as.owinlist function
##'
##' Generic function for creating lists of owin objects
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method as.owinlist
##' @export

as.owinlist <- function(obj,...){
    UseMethod("as.owinlist")
}



##' as.owinlist.SpatialPolygonsDataFrame function
##'
##' A function to create a list of owin objects from a SpatialPolygonsDataFrame
##'
##' @method as.owinlist SpatialPolygonsDataFrame
##' @param obj a SpatialPolygonsDataFrame object
##' @param dmin purpose is to simplify the SpatialPolygons. A numeric value giving the smallest permissible length of an edge. See ? simplify.owin
##' @param check whether or not to use spatstat functions to check the validity of SpatialPolygons objects
##' @param subset logical vector. Subset of regions to extract and conver to owin objects. By default, all regions are extracted.
##' @param ... additional arguments
##' @return a list of owin objects corresponding to the constituent Polygons objects
##' @export

as.owinlist.SpatialPolygonsDataFrame <- function(obj,dmin=0,check=TRUE,subset=rep(TRUE,length(obj)),...){
    proj <- obj@proj4string
    nreg <- length(obj) # number of regions
    
    spatstat.options(checkpolygons = check)
    
    owinlist <- list()
    for (i in 1:nreg){
        if(!subset[i]){
            next
        }
        pol <- list()
        pol[[1]] <- obj@polygons[[i]]
        owinlist[[i]] <- as(SpatialPolygons(pol,proj4string=proj),"owin")
        if (dmin>0){
            owinlist[[i]] <- simplify.owin(owinlist[[i]],dmin=dmin)
        }
    }
    
    if(!check){
        spatstat.options(checkpolygons = !check)
    }
        
    class(owinlist) <- c("owinlist")
    return(owinlist)
}



##' as.owinlist.stapp function
##'
##' A function to create a list of owin objects from a stapp
##'
##' @method as.owinlist stapp
##' @param obj an stapp object
##' @param dmin purpose is to simplify the SpatialPolygons. A numeric value giving the smallest permissible length of an edge. See ? simplify.owin
##' @param check whether or not to use spatstat functions to check the validity of SpatialPolygons objects
##' @param ... additional arguments
##' @return a list of owin objects corresponding to the constituent Polygons objects
##' @export

as.owinlist.stapp <- function(obj,dmin=0,check=TRUE,...){
    return(as.owinlist.SpatialPolygonsDataFrame(obj=obj$spdf,dmin=dmin,check=check,...))
}


##' expectation.aggregatedPredict function
##'
##' \bold{This function requires data to have been dumped to disk}: see \code{?dump2dir} and \code{?setoutput}. This function computes the 
##' Monte Carlo Average of a function where data from a run of \code{lgcpPredict} has been dumped to disk. It returns the results averaged over the
##' regions defined in the original stapp object.
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
##' @method expectation aggregatedPredict
##' @param obj an object of class lgcpPredict
##' @param fun a function accepting a single argument that returns a numeric vector, matrix or array object
##' @param maxit Not used in ordinary circumstances. Defines subset of samples over which to compute expectation. Expectation is computed using information from iterations 1:maxit, where 1 is the first non-burn in iteration dumped to disk. 
##' @param ... additional arguments  
##' @return the expectated value of that quantity
##' @seealso \link{lgcpPredict}, \link{dump2dir}, \link{setoutput}
##' @export
#
#expectation.aggregatedPredict <- function(obj,fun,maxit=NULL,...){
#    if (is.null(obj$gridfunction$dirname)){
#        stop("dump2dir not specified, MCMC output must have be dumped to disk to use this function.  See ?dump2dir.")
#    }
#    fn <- paste(obj$gridfunction$dirname,"simout.nc",sep="")
#    ncdata <- open.ncdf(fn)
#    if (!is.null(maxit)){
#        if (maxit<2 | maxit>ncdata$dim$iter$len){
#            stop(paste("maxit must be between 2 and",ncdata$dim$iter$len))
#        }
#        endloop <- maxit
#    }
#    else{
#        endloop <- ncdata$dim$iter$len
#    }
#
#    gridinwindowlist <- lapply(lg$owinlist,function(W){gridInWindow(lg$mcens,lg$ncens,W)})
#    gridinwindowlist[[22]] <- matrix(FALSE,64,64)    
#    
#    Y <- lgcpgrid(get.var.ncdf(nc=ncdata, varid=ncdata$var[[1]], start=c(1,1,1,1), count=c(-1,-1,-1,1))) # first "simulation"
#    result <- lapply(gridinwindowlist,function(inW){lapply(Y$grid,function(y){fun(y[inW])})})
#
#    browser()    
#    
#    pb <- txtProgressBar(min=1,max=endloop,style=3)
#    setTxtProgressBar(pb,1)
#    for (i in 2:endloop){
#        setTxtProgressBar(pb,i)
#        Y <- lgcpgrid(get.var.ncdf(nc=ncdata, varid=ncdata$var[[1]], start=c(1,1,1,i), count=c(-1,-1,-1,1)))
#        tmp <- lapply(gridinwindowlist,function(inW){lapply(Y$grid,function(y){fun(y[inW])})})
#        for (j in 1:length(result)){
#            result[[j]] <- add.list(result[[j]],tmp[[j]])
#        }
#    }
#    close(pb)
#    for (i in 1:length(result)){
#        for (j in 1:length(result[[i]])){
#            result[[i]][[j]] <- result[[i]][[j]] / endloop
#        }
#    }
#    return(result)
#}

