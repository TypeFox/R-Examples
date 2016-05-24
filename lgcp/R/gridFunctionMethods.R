###
# Generic functions for the class gridFunction    
###

##' GFinitialise function
##'
##' Generic function defining the the initialisation step for the \code{gridFunction} class of objects. 
##' The function is called invisibly within \code{MALAlgcp} and facilitates the dumping of data to disk
##'
##' @param F an object    
##' @param ... additional arguments  
##' @return method GFinitialise
##' @seealso \link{setoutput}, \link{GFupdate}, \link{GFfinalise}, \link{GFreturnvalue} 
##' @export

GFinitialise <- function(F,...){
    UseMethod("GFinitialise")
}



##' GFupdate function
##'
##' Generic function defining the the update step for the \code{gridFunction} class of objects. 
##' The function is called invisibly within \code{MALAlgcp} and facilitates the dumping of data to disk
##'
##' @param F an object    
##' @param ... additional arguments  
##' @return method GFupdate
##' @seealso \link{setoutput}, \link{GFinitialise}, \link{GFfinalise}, \link{GFreturnvalue}
##' @export

GFupdate <- function(F,...){
    UseMethod("GFupdate")
}



##' GFfinalise function
##'
##' Generic function defining the the finalisation step for the \code{gridFunction} class of objects. 
##' The function is called invisibly within \code{MALAlgcp} and facilitates the dumping of data to disk
##'
##' @param F an object    
##' @param ... additional arguments  
##' @return method GFfinalise
##' @seealso \link{setoutput}, \link{GFinitialise}, \link{GFupdate}, \link{GFreturnvalue}
##' @export

GFfinalise <- function(F,...){
    UseMethod("GFfinalise")
}



##' GFreturnvalue function
##'
##' Generic function defining the the returned value for the \code{gridFunction} class of objects. 
##' The function is called invisibly within \code{MALAlgcp} and facilitates the dumping of data to disk
##'
##' @param F an object    
##' @param ... additional arguments  
##' @return method GFreturnvalue
##' @seealso \link{setoutput}, \link{GFinitialise}, \link{GFupdate}, \link{GFfinalise}
##' @export

GFreturnvalue <- function(F,...){
    UseMethod("GFreturnvalue")
}



###
# Functions to facilitate dumping of mcmc output to a directory
###


##' nullFunction function
##'
##' This is a null function and performs no action.
##'
##' @return object of class nullFunction
##' @seealso \link{setoutput}, \link{GFinitialise}, \link{GFupdate}, \link{GFfinalise}, \link{GFreturnvalue}
##' @export

nullFunction <- function(){
    obj <- "NULL"
    class(obj) <- c("nullFunction","gridFunction")
    return(obj)
}

##' dump2dir function
##'
##' This function, when set by the \code{gridfunction} argument of \link{setoutput}, in turn called by the argument 
##' \code{output.control} of \code{lgcpPredict} facilitates the dumping of data to disk. Data is dumped to a
##' netCDF file, \code{simout.nc}, stored in the directory specified by the user. If the directory does not exist, 
##' then it will be created. Since the requested data dumped to disk may be very large in a run of \code{lgcpPredict}, 
##' by default, the user is prompted as to whether to proceed with prediction, this can be turned off by setting 
##' the option \code{forceSave=TRUE} detailed here. To save space, or increase the number of simulations that can be
##' stored for a fixed disk space the option to only save the last time point is also available (\code{lastonly=TRUE}, 
##' which is the default setting).
##'
##' @param dirname character vector of length 1 containing the name of the directory to create
##' @param lastonly only save output from time T? (see ?lgcpPredict for definition of T)
##' @param forceSave option to override display of menu
##' @return object of class dump2dir
##' @seealso \link{setoutput}, \ \link{GFinitialise}, \link{GFupdate}, \link{GFfinalise}, \link{GFreturnvalue}
##' @export

dump2dir <- function(dirname,lastonly=TRUE,forceSave=FALSE){
    verifyclass(dirname,"character")
    if(length(dirname)>1){
        stop("dirname must have length 1")
    }
    obj <- list()
    if (substr(dirname,nchar(dirname),nchar(dirname))!="/"){
        dirname <- paste(dirname,"/",sep="")
    }
    obj$dirname <- path.expand(dirname)
    obj$lastonly <- lastonly
    obj$forceSave <- forceSave
    obj$runid <- as.numeric(Sys.time())
    dumpidx <- 0
    incridx <- function(){
        dumpidx <<- dumpidx + 1
    }
    retidx <- function(){
        return(dumpidx)
    }
    obj$incr <- incridx
    obj$ret <- retidx
    class(obj) <- c("dump2dir","gridFunction")
    return(obj)
}



##' print.dump2dir function
##'
##' Display function for dump2dir objects.
##'
##' @method print dump2dir
##' @param x an object of class dump2dir    
##' @param ... additional arguments 
##' @return nothing
##' @seealso \link{dump2dir}, 
##' @export

print.dump2dir <- function(x,...){
    cat("dump2dir object.\n")
    cat(paste("           Run ID: ",x$runid,"\n",sep=""))
    cat(paste("        Last only? ",x$lastonly,"\n",sep=""))
    cat(paste("        Directory: ",x$dirname,"\n",sep=""))
}



##' GFinitialise.nullFunction function
##'
##' This is a null function and performs no action.
##'
##' @method GFinitialise nullFunction
##' @param F an object of class dump2dir    
##' @param ... additional arguments 
##' @return nothing
##' @seealso \link{nullFunction}, \link{setoutput}, \link{GFinitialise}, \link{GFupdate}, \link{GFfinalise}, \link{GFreturnvalue}
##' @export

GFinitialise.nullFunction <- function(F,...){
    return(NULL)
}



##' GFinitialise.dump2dir function
##'
##' Creates a directory (if necessary) and allocates space for a netCDF dump.
##'
##' @method GFinitialise dump2dir
##' @param F an object of class dump2dir    
##' @param ... additional arguments 
##' @return creates initialisation file and folder
##' @seealso \link{dump2dir}, \link{setoutput}, \link{GFinitialise}, \link{GFupdate}, \link{GFfinalise}, \link{GFreturnvalue}
##' @export

GFinitialise.dump2dir <- function(F,...){
    # Create directory for dumping files
    dir.create(F$dirname)
    runid <- F$runid
    gridsize <- c(get("M",envir=parent.frame()),get("N",envir=parent.frame()),length(get("temporal.fitted",envir=parent.frame())))
    mLoop <- get("mcmcloop",envir=parent.frame())
    nsamp <- floor((mLoop$N-mLoop$burnin)/mLoop$thin)
    fn <- paste(F$dirname,"simulationinfo.RData",sep="")
    save(nsamp,gridsize,mLoop,runid,file=fn)
    # now prepare and open the netcdf file
    M <-get("M",envir=parent.frame())
    N <-get("N",envir=parent.frame())
    tlen <- length(get("temporal.fitted",envir=parent.frame()))
    nlevs <- get("nlevs",envir=parent.frame())
    MTmode <- get("MultiTypeMode",envir=parent.frame())
    fn <- paste(F$dirname,"simout.nc",sep="")  
    x <- ncdim_def( "X", "x coordinates", 1:M) 
    y <- ncdim_def( "Y", "y coordinates", 1:N)
    if (F$lastonly){
        t <- ncdim_def( "T", "time index", 1) 
    }
    else{
        t <- ncdim_def( "T", "time index", 1:tlen) 
    }   
    iter <- ncdim_def( "iter", "iteration", 1:nsamp)
    if(MTmode){
        proc <- ncdim_def( "P", "processes", 1:(nlevs+1)) # gives space for one common field and 'nlevs' other fields
        sout <- ncvar_def("simrun","none", list(x,y,t,proc,iter), missval=1.e30,prec="double")
    }
    else{ 
        sout <- ncvar_def("simrun","none", list(x,y,t,iter), missval=1.e30,prec="double")
    }
    ncdata <- nc_create(fn,sout) # allocates the disk space to be written
    nc_close(ncdata)
    cat(paste("Netcdf file: ",fn," created\n",sep=""))
}



##' GFupdate.nullFunction function
##'
##' This is a null function and performs no action.
##'
##' @method GFupdate nullFunction
##' @param F an object of class dump2dir    
##' @param ... additional arguments 
##' @return nothing
##' @seealso \link{nullFunction}, \link{setoutput}, \link{GFinitialise}, \link{GFupdate}, \link{GFfinalise}, \link{GFreturnvalue}
##' @export

GFupdate.nullFunction <- function(F,...){
    return(NULL)
}



##' GFupdate.dump2dir function
##'
##' This function gets the required information from \code{MALAlgcp} and writes the data to the netCDF file.
##'
##' @method GFupdate dump2dir
##' @param F an object    
##' @param ... additional arguments 
##' @return saves latent field
##' @seealso \link{dump2dir}, \link{setoutput}, \link{GFinitialise}, \link{GFupdate}, \link{GFfinalise}, \link{GFreturnvalue}
##' @export

GFupdate.dump2dir <- function(F,...){
    mLoop <- get("mcmcloop",envir=parent.frame())
    nsamp <- floor((mLoop$N-mLoop$burnin)/mLoop$thin)
    runid <- F$runid
    F$incr() 
    M <- get("M",envir=parent.frame())
    N <- get("N",envir=parent.frame())
    tfit <- get("temporal.fitted",envir=parent.frame())
    tlen <- length(tfit)
    ncdata <- nc_open(paste(F$dirname,"simout.nc",sep=""),write=TRUE)
    
    if(get("SpatialOnlyMode",envir=parent.frame())){
        if(!get("SpatialPlusParameters",envir=parent.frame())){
            Y <- list(get("oldtags",envir=parent.frame())$Y[1:M,1:N])
        }
        else{
            if(get("SpatialPlusParameters",envir=parent.frame())&!get("MultiTypeMode",envir=parent.frame())){
                Y <- list(get("GP",envir=parent.frame())$Y[1:M,1:N])
            }
            else if(get("SpatialPlusParameters",envir=parent.frame())&get("MultiTypeMode",envir=parent.frame())){
                Y <- list(get("Y",envir=parent.frame())[1:M,1:N,])
                nfields <- dim(Y[[length(Y)]])[3]
            }
            else{
                stop("unidentified MCMC method in GFupdate.dump2dir")            
            }
        }
    }
    else{ # in spatiotemporal mode
        if(get("SpatioTemporalPlusParameters",envir=parent.frame())){
            Y <- get("GP",envir=parent.frame())$Y
            Y <- lapply(Y,function(x){x[1:M,1:N]})
        }
        else{
            Y <- get("oldtags",envir=parent.frame())$Y
            Y <- lapply(Y,function(x){x[1:M,1:N]})
        }
    }

    if (F$lastonly){
        if(get("SpatialPlusParameters",envir=parent.frame())&get("MultiTypeMode",envir=parent.frame())){
            ncvar_put(ncdata,ncdata$var[[1]],Y[[length(Y)]],start=c(1,1,1,1,F$ret()),count=c(M,N,1,nfields,1))
        }
        else{       
            ncvar_put(ncdata,ncdata$var[[1]],Y[[length(Y)]],start=c(1,1,1,F$ret()),count=c(M,N,1,1))
        }
    }
    else{
        for (i in 1:tlen){ 
            ncvar_put(ncdata,ncdata$var[[1]],Y[[i]],start=c(1,1,i,F$ret()),count=c(M,N,1,1))
        }
    }
    nc_sync(ncdata)
    nc_close(ncdata) 
}



##' GFfinalise.nullFunction function
##'
##' This is a null function and performs no action. 
##'
##' @method GFfinalise nullFunction
##' @param F an object of class dump2dir    
##' @param ... additional arguments 
##' @return nothing
##' @seealso \link{nullFunction}, \link{setoutput}, \link{GFinitialise}, \link{GFupdate}, \link{GFfinalise}, \link{GFreturnvalue}
##' @export

GFfinalise.nullFunction <- function(F,...){
    return(NULL)
}



##' GFfinalise.dump2dir function
##'
##' This function finalises the dumping of data to a netCDF file.
##'
##' @method GFfinalise dump2dir
##' @param F an object    
##' @param ... additional arguments 
##' @return nothing
##' @seealso \link{dump2dir}, \link{setoutput}, \link{GFinitialise}, \link{GFupdate}, \link{GFfinalise}, \link{GFreturnvalue}
##' @export

GFfinalise.dump2dir <- function(F,...){
    return(NULL)
}



##' GFreturnvalue.nullFunction function
##'
##' This is a null function and performs no action. 
##'
##' @method GFreturnvalue nullFunction
##' @param F an object of class dump2dir    
##' @param ... additional arguments 
##' @return nothing
##' @seealso \link{nullFunction}, \link{setoutput}, \link{GFinitialise}, \link{GFupdate}, \link{GFfinalise}, \link{GFreturnvalue}
##' @export

GFreturnvalue.nullFunction <- function(F,...){
    return(NULL)
}



##' GFreturnvalue.dump2dir function
##'
##' This function returns the name of the directory the netCDF file was written to.
##'
##' @method GFreturnvalue dump2dir
##' @param F an object    
##' @param ... additional arguments 
##' @return display where files have been written to
##' @seealso \link{dump2dir}, \link{setoutput}, \link{GFinitialise}, \link{GFupdate}, \link{GFfinalise}, \link{GFreturnvalue}
##' @export

GFreturnvalue.dump2dir <- function(F,...){
    cat(paste("Files written to ",F$dirname,"\n",sep=""))
    return(F)
}
