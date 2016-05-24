##' cruts2raster function
##'
##' A function to convert Climatic Research Unit Time-Series in NetCDF format to raster format. 
##'
##' Data can be obtained from \url{http://catalogue.ceda.ac.uk/uuid/ac4ecbd554d0dd52a9b575d9666dc42d}
##'
##' @param ncfile character string giving name and location of the CRUTS time series NetCDF file (if the file you downloaded is zipped, then you will need to extract it) 
##' @param timeRange vector of length 2 giving the start and end dates in the first and second place. Dates are converted using the function ymd, please refer to the help for this funciton for details on appropriate formats.
##' @param poly an optional SpatialPolygonsDataFrame on which to crop the raster to
##' @param offset time offset for CRU TS data
##' @param type can be either 'brick' or 'stack' (thde default), this argument specifies what sort of raster object to return.
##' @return a raster stack or brick containing the raw data
##' @examples
##' \dontrun{crutsimport(ncfile="my_cruts_file.nc",timeRange=c("2000-01-01","2001-01-01"))}
##' @export

cruts2raster <- function(ncfile,timeRange=NULL,poly=NULL,offset="1900-01-01",type="stack"){

    nc <- nc_open(ncfile)
    lon <- nc$dim$lon$vals
    lat <- nc$dim$lat$vals
    time <- nc$dim$time$vals
    cstime <- time[1] + cumsum(diff(time))
    d <- sapply(nc$dim,function(x){x$len}) # dimension

    starttime <- as.Date(offset)
    times <- ymd(starttime+cstime)

    if(!is.null(poly)){
        poly <- spTransform(poly,CRS("+init=epsg:4326"))    
    }

    M <- length(lon)
    N <- length(lat)
    dx <- diff(lon[1:2])
    dy <- diff(lat[1:2])
    firstslice <- extractNetCDF(nc)
    firstslice <- raster(t(firstslice[,N:1]), xmn=lon[1]-dx/2, xmx=lon[M]+dx/2, ymn=lat[1]-dy/2, ymx=lat[N]+dy/2, crs=CRS("+init=epsg:4326"))


    if(is.null(timeRange)){
        tmax <- max(time,na.rm=TRUE)
        tmin <- tmax - min(10,length(time))
        warning("Returning last",tmax-tmin,"time points.\n",immediate.=TRUE)
    }
    else{
        tmin <- min(which(times>=ymd(timeRange[1])))    
        tmax <- max(which(times<=ymd(timeRange[2])))
    }

    rlist <- NULL
    pb <- txtProgressBar(min=tmin,max=tmax,style=3)
    for(i in tmin:tmax){

        lay <- extractNetCDF(nc,start=c(1,1,i))
        lay <- raster(t(lay[,N:1]), xmn=lon[1]-dx/2, xmx=lon[M]+dx/2, ymn=lat[1]-dy/2, ymx=lat[N]+dy/2, crs=CRS("+init=epsg:4326"))
        if(is.null(poly)){         
            rlist <- c(rlist,lay)
        }
        else{
            rlist <- c(rlist,crop(lay,poly))
        }
        setTxtProgressBar(pb,i) 
    }
    close(pb)
    nc_close(nc)
    names(rlist) <- as.character(times[tmin:tmax])

    if(type=="stack"){
        br <- stack(rlist)
        attr(br,"times") <- times[tmin:tmax]
    }
    else if(type=="brick"){
        br <- brick(rlist)
        attr(br,"times") <- times[tmin:tmax]
    }
    else{
        stop("type must be one of 'stack' or 'brick'")
    }
    

    return(br)
}

##' cruts2poly function
##'
##' A function to convert Climatic Research Unit Time-Series in NetCDF format to polygonal format, averaging over each region in question.
##'
##' Data can be obtained from \url{http://catalogue.ceda.ac.uk/uuid/ac4ecbd554d0dd52a9b575d9666dc42d}
##'
##' @param ncfile character string giving name and location of the CRUTS time series NetCDF file (if the file you downloaded is zipped, then you will need to extract it)  
##' @param poly a SpatialPolygonsDataFrame on which to average the variable in question 
##' @param timeRange vector of length 2 giving the start and end dates in the first and second place. Dates are converted using the function ymd, please refer to the help for this funciton for details on appropriate formats.
##' @param offset time offset for CRU TS data
##' @param na.rm logical, whether to ignore NA's in averaging, default is FALSE (to be consistent with other R functions in other packages), but option TRUE should probably be used on most occasions
##' @return a polygon with the averaged climate variable
##' @examples
##' \dontrun{crutsimport(ncfile="my_cruts_file.nc",timeRange=c("2000-01-01","2001-01-01"))}
##' @export

cruts2poly <- function(ncfile,poly,timeRange=NULL,offset="1900-01-01",na.rm=FALSE){

    polyp4 <- proj4string(poly)

    cat("Extracting data ...\n")
    br <- cruts2raster(ncfile=ncfile,timeRange=timeRange,poly=poly,offset=offset,type="brick")

    if(!is.null(poly)){
        poly <- spTransform(poly,CRS("+init=epsg:4326"))    
    }
    cat("Performing overlay operation ...\n")    
    ext <- raster::extract(br,poly,fun=mean,na.rm=na.rm)

    times <- attr(br,"times")

    colnames(ext) <- paste("Y",year(times),"M",month(times),sep="")

    poly@data <- as.data.frame(ext)
    poly <- spTransform(poly,CRS(polyp4)) #transform back
    attr(poly,"times") <- times
    return(poly)
}


##' extractNetCDF function
##'
##' A function to extract data from CRU TS NetCDF files. A wrapper function for get.var.ncdf.
##'
##' @param nc an object inheriting class ncdf
##' @param start the start index. A vector of indices indicating where to start reading the passed values (beginning at 1).  The length of this vector must equal the number of dimensions the variable has.  If not specified, reading starts at the beginning of the file (1,1,1,...).
##' @param count A vector of integers indicating the count of values to read along each dimension.  The length of this vector must equal the number of dimensions the variable has. If not specified and the variable does NOT have an unlimited dimension, the entire variable is read.  As a special case, the value '-1' indicates that all entries along that dimension should be read. By default this extracts data for the  first time point.
##' @return an array or matrix with the requested data
##' @examples
##' \dontrun{extractNetCDF(dat=dat)}
##' \dontrun{extractNetCDF(dat=dat,start=c(1,1,1,1),count=c(1,2,3,1))}
##' \dontrun{extractNetCDF(dat=dat,start=c(1,1,1,1),count=c(-1,1,1,1))}
##' @export


extractNetCDF <- function(nc,start=NULL,count=NULL){

    datadim <- nc$var[[1]]$varsize
    
    if(is.null(start)){
        start <- c(1,1,1)
    }

    if(is.null(count)){
        count <- c(-1,-1,1)           
    }
    
    x <- ncvar_get(nc=nc, varid=nc$var[[1]], start=start, count=count)

    return(x)
}




##' getAnomaly function
##'
##' A function to extract anomalies from the Climatic Research Unit Time-Series dataset
##'
##' Data can be obtained from \url{http://catalogue.ceda.ac.uk/uuid/ac4ecbd554d0dd52a9b575d9666dc42d}
##'
##' @param ncfile character string giving name and location of the CRUTS time series NetCDF file (if the file you downloaded is zipped, then you will need to extract it) 
##' @param poly an optional SpatialPolygonsDataFrame on which to compute the average anomalies if NULL (the default) a raster brick will be returned
##' @param timeRange vector of length 2 giving the start and end dates in the first and second place. Dates are converted using the function ymd, please refer to the help for this funciton for details on appropriate formats. 
##' @param offset time offset for CRU TS data
##' @param na.rm logical, whether to ignore NA's in averaging, default is FALSE (to be consistent with other R functions in other packages), but option TRUE should probably be used on most occasions 
##' @return a raster or polygon with the raw or spatially averaged anomalies
##' @export

getAnomaly <- function(ncfile,poly=NULL,timeRange=NULL,offset="1900-01-01",na.rm=FALSE){

    br <- cruts2raster(ncfile=ncfile,timeRange=timeRange,poly=poly,offset=offset,type="brick")

    arr <- raster::as.array(br)

    times <- attr(br,"times")
    ym <- cbind(y=year(times),m=month(times))

    month_mean <- list()
    month_sd <- list()
    for(i in 1:12){
        if(!any(ym[,2]==i)){
            next
        }
        else{
            a <- arr[,,ym[,2]==i]
            month_mean[[i]] <- apply(a,c(1,2),mean,na.rm=TRUE)
            month_sd[[i]] <- apply(a,c(1,2),sd,na.rm=TRUE)
        }
    }

    ans_arr <- array(NA,dim=dim(arr))
    for(i in 1:dim(arr)[3]){
        mnth <- ym[i,2]
        ans_arr[,,i] <- (arr[,,i] - month_mean[[mnth]]) / month_sd[[mnth]]
    }

    bra <- br

    raster::values(bra) <- raster::values(brick(ans_arr))

    if(is.null(poly)){
        return(bra)
    }
    else{  
        polyp4 <- proj4string(poly)    
        poly <- spTransform(poly,CRS("+init=epsg:4326"))    
        
        cat("Performing overlay operation ...\n")    
        ext <- raster::extract(bra,poly,fun=mean,na.rm=na.rm)

        colnames(ext) <- paste("Y",year(times),"M",month(times),sep="")

        poly@data <- as.data.frame(ext)
        poly <- spTransform(poly,CRS(polyp4))
        attr(poly,"times") <- times

        return(poly)
    }
}