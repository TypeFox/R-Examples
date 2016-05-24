##' aggCovInfo function
##'
##' Generic function for aggregation of covariate information.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method aggCovInfo
##' @export

aggCovInfo <- function(obj,...){
    UseMethod("aggCovInfo")
}


##' aggCovInfo.Majority function
##'
##' Aggregation via majority.
##'
##' @method aggCovInfo Majority
##' @param obj an Majority object
##' @param regwts regional (areal) weighting vector
##' @param ... additional arguments
##' @return The most popular cell type.
##' @export

aggCovInfo.Majority <- function(obj,regwts,...){
    if(all(is.na(obj))){
        return(NA)    
    }
    else{
        idx <- which(regwts==max(regwts,na.rm=TRUE))
        if(length(idx)>1){
            idx <- sample(idx,1)
        }
        return(obj[idx])
    }
}



##' aggCovInfo.ArealWeightedMean function
##'
##' Aggregation via weighted mean.
##'
##' @method aggCovInfo ArealWeightedMean
##' @param obj an ArealWeightedMean object
##' @param regwts regional (areal) weighting vector
##' @param ... additional arguments
##' @return Areal weighted mean.
##' @export

aggCovInfo.ArealWeightedMean <- function(obj,regwts,...){
    x <- regwts*obj/sum(regwts,na.rm=TRUE) 
    if(all(is.na(x))){
        return(NA)
    }
    else{
        return(sum(x,na.rm=TRUE))
    }    
    return()
}


##' aggCovInfo.ArealWeightedSum function
##'
##' Aggregation via weighted sum. Use to sum up population counts in regions. 
##'
##' @method aggCovInfo ArealWeightedSum
##' @param obj an ArealWeightedSum object
##' @param regwts regional (areal) weighting vector
##' @param ... additional arguments
##' @return Areal weighted Sum.
##' @export

aggCovInfo.ArealWeightedSum <- function(obj,regwts,...){
    x <- (regwts/attr(regwts,"polyareas"))*obj 
    if(all(is.na(x))){
        return(NA)
    }
    else{
        return(sum(x,na.rm=TRUE))
    }    
    return()
}



##' aggregateCovariateInfo function
##'
##' A function called by cov.interp.fft to allocate and perform interpolation of covariate infomation onto the FFT grid
##'
##' @param cellidx the index of the cell
##' @param cidx index of covariate, no longer used 
##' @param gidx grid index
##' @param df the data frame containing the covariate information
##' @param fftovl an overlay of the fft grid onto the SpatialPolygonsDataFrame or SpatialPixelsDataFrame objects
##' @param classes vector of class attributes of the dataframe
##' @param polyareas polygon areas of the SpatialPolygonsDataFrame or SpatialPixelsDataFrame objects
##' @return the interpolated covariate information onto the FFT grid
##' @export

aggregateCovariateInfo <- function(cellidx,cidx,gidx,df,fftovl,classes,polyareas){
    olinfo <- fftovl$info[fftovl$info$grididx==gidx[cellidx],]
    dat <- df[olinfo$polyidx,,drop=FALSE]
    sapply(1:ncol(dat),function(i){class(dat[[i]]) <<- classes[,i]}) 
    regwts <- olinfo$area
    attr(regwts,"polyareas") <- polyareas[olinfo$polyidx]
    return(sapply(dat,aggCovInfo,regwts=regwts))    
}

##' cov.interp.fft function
##'
##' A function to interpolate covariate values onto the fft grid, ready for analysis
##'
##' @param formula an object of class formula (or one that can be coerced to that class) starting with X ~ (eg X~var1+var2 *NOT for example* Y~var1+var2): a symbolic description of the model to be fitted.
##' @param W an owin observation window
##' @param regionalcovariates an optional SpatialPolygonsDataFrame
##' @param pixelcovariates an optional SpatialPixelsDataFrame
##' @param mcens x-coordinates of output grid centroids (not fft grid centroids ie *not* the extended grid)
##' @param ncens y-coordinates of output grid centroids (not fft grid centroids ie *not* the extended grid)
##' @param cellInside a 0-1 matrix indicating which computational cells are inside the observation window
##' @param overl an overlay of the computational grid onto the SpatialPolygonsDataFrame or SpatialPixelsDataFrame.
##' @return The interpolated design matrix, ready for analysis
##' @export

cov.interp.fft <- function(formula,W,regionalcovariates=NULL,pixelcovariates=NULL,mcens,ncens,cellInside,overl=NULL){

    if(!is.null(overl)){
        if(!(isTRUE(all.equal(mcens,overl$mcens)) & isTRUE(all.equal(ncens,overl$ncens)))){
            stop("Differing FFT grids ... check this is the correct overlay")
        }    
    }

    varn <- variablesinformula(formula)[-1] #independent variable names    
    
    if(!is.null(regionalcovariates)){
        if(any(is.na(getinterp(regionalcovariates@data)))){
            stop("No interpolation method specified for one or more regional covariates, see ?assigninterp and ?guessinterp. Note that this should also be supplied for pixelcovariates, if applicable.")
        }
    }
    if(!is.null(pixelcovariates)){
        if(any(is.na(getinterp(pixelcovariates@data)))){
            stop("No interpolation method specified for one or more pixel covariates, see ?assigninterp and ?guessinterp")
        }
    }

    M <- length(mcens)
    N <- length(ncens)

    if(is.null(regionalcovariates) & is.null(pixelcovariates)){

        testf <- as.formula(X~1)
        attr(testf, ".Environment") <- .GlobalEnv        

        if(identical(formula,testf)){
            Zmat <- matrix(cellInside,M*N,1)
            colnames(Zmat) <- "(Intercept)"
            attr(Zmat,"data.frame") <- data.frame(X=rep(1,sum(cellInside)))
            attr(Zmat,"cellInside") <- cellInside
            attr(Zmat,"anymiss") <- NULL
            attr(Zmat,"M") <- M
            attr(Zmat,"N") <- N
            attr(Zmat,"mcens") <- mcens
            attr(Zmat,"ncens") <- ncens
            attr(Zmat,"polygonOverlay") <- NULL
            attr(Zmat,"pixelOverlay") <- NULL
            attr(Zmat,"FORM") <- formula
            attr(Zmat,"fftpoly") <- NULL
            return(Zmat)
        }        
        
        stop("Must have either regional or pixel covariates.")
    }

    spoly <- as(W,"SpatialPolygons")    
    
    if(!is.null(overl)){
        fftpoly <- overl$fftpoly
    }
    else{
        fftpoly <- grid2spoly(mcens,ncens) # fft grid cells
    }   
    
   
    Zmat <- 0
    dmat <- as.data.frame(matrix(NA,length(fftpoly),length(varn)))
    clsvec <- c()
    polygonoverlay <- NULL
    pixeloverlay <- NULL
    if(!is.null(regionalcovariates)){
        cat("aggregating regional covariate information ...\n")
        s <- Sys.time()        
        polyareas <- sapply(1:length(regionalcovariates),function(ii){sum(sapply(slot(regionalcovariates[ii,],"polygons"),slot,"area"))})
        if(!is.null(overl)){
            cat("loading polygon overlay ...\n")
            polyol <- overl$polyol
        }
        else{
            cat("performing polygon overlay operations ...\n")
            polyol <- gOverlay(fftpoly,regionalcovariates)
        }         
        cidx <- match(varn,names(regionalcovariates))
        cidx <- cidx[!is.na(cidx)]       
        gidx <- sort(unique(polyol$info$grididx))        
        dtemp <- as.data.frame(matrix(NA,length(gidx),length(cidx)))        
        regionalcovariates <- regionalcovariates@data[,cidx,drop=FALSE] # drop=FALSE here prevents loss of interpolation information
        classes <- sapply(regionalcovariates,class) 
        cat("interpolating ...\n")
        sapply(1:length(gidx),function(i){dtemp[i,1:length(cidx)] <<- aggregateCovariateInfo(cellidx=i,cidx=cidx,gidx=gidx,df=regionalcovariates,fftovl=polyol,classes=classes,polyareas=polyareas)})
        dmat[gidx,1:length(cidx)] <- dtemp
        e <- Sys.time()
        cat("Time Taken: ",difftime(e, s, units = "secs"),"\n")
        clsvec <- c(clsvec,sapply(clearinterp(regionalcovariates),class))
        polygonoverlay <- polyol
    }
    if(!is.null(pixelcovariates)){
        cat("aggregating pixel covariate information ...\n")
        s <- Sys.time()
        cat("converting SpatialPixelsDataFrame to SpatialPolygonsDataFrame ...\n")
        oldclasses <- sapply(pixelcovariates@data,class)
        polyareas <- rep(prod(pixelcovariates@grid@cellsize),length(pixelcovariates))
        pixelcovariates <- as(pixelcovariates,"SpatialPolygonsDataFrame")
        sapply(1:ncol(pixelcovariates),function(i){class(pixelcovariates[[i]]) <<- oldclasses[,i]}) 
        if(!is.null(overl)){
            cat("loading polygon overlay ...\n")
            polyol <- overl$pixol
        }
        else{
            cat("performing polygon overlay operations ...\n")
            polyol <- gOverlay(fftpoly,pixelcovariates)   
        }
        cidx <- match(varn,names(pixelcovariates))
        cidx <- cidx[!is.na(cidx)]       
        gidx <- sort(unique(polyol$info$grididx))        
        dtemp <- as.data.frame(matrix(NA,length(gidx),length(cidx)))        
        pixelcovariates <- pixelcovariates@data[,cidx,drop=FALSE] # drop=FALSE here prevents loss of interpolation information
        classes <- sapply(pixelcovariates,class) 
        cat("interpolating ...\n")
        sapply(1:length(gidx),function(i){dtemp[i,1:length(cidx)] <<- aggregateCovariateInfo(cellidx=i,cidx=cidx,gidx=gidx,df=pixelcovariates,fftovl=polyol,classes=classes,polyareas=polyareas)})
        dmat[gidx,(length(varn)-length(cidx)+1):length(varn)] <- dtemp
        e <- Sys.time()
        cat("Time Taken: ",difftime(e, s, units = "secs"),"\n")
        clsvec <- c(clsvec,sapply(clearinterp(pixelcovariates),class))
        pixeloverlay <- polyol
    }
    cellin <- which(as.logical(cellInside))
    dmat <- dmat[cellin,,drop=FALSE]
    anymiss <- apply(dmat,1,function(x){any(is.na(x))})
    anymisscopy <- anymiss
    if(any(anymiss)){
        warning(paste("There is(are)",sum(anymiss),"missing value(s) in the covariate data. Imputing with median of non-missing values."),immediate.=TRUE)
        for(iii in 1:ncol(dmat)){
            if(any(is.na(dmat[,iii]))){
                dmat[,iii][is.na(dmat[,iii])] <- median(dmat[,iii],na.rm=TRUE)
            }
        }
        anymiss <- apply(dmat,1,function(x){any(is.na(x))})        
    }
    
    sapply(1:length(clsvec),function(i){dmat[[i]] <<- eval(call(paste("as.",clsvec[i],sep=""),dmat[[i]]))})
    names(dmat) <- varn
    dmat$X <- rep(1,sum(cellInside))            
    dmat <- as.data.frame(dmat)
    DM <- dmat
    dmat <- model.matrix(formula,data=dmat)
    
    missingind <- rep(0,dim(dmat)[2])
    
    Zmat <- matrix(0,M*N,dim(dmat)[2])
    colnames(Zmat) <- colnames(dmat)
    rownames(dmat) <- NULL
    idx <- as.logical(cellInside)
    idx[as.logical(cellInside)] <- !anymiss
    Zmat[idx,] <- dmat # [!anymiss] because NAs are removed by model.matrix
    if(any(anymisscopy)){
        missingind[idx] <- as.numeric(anymisscopy)
    }
    attr(Zmat,"data.frame") <- DM
    attr(Zmat,"cellInside") <- cellInside
    attr(Zmat,"anymiss") <- anymiss
    attr(Zmat,"mcens") <- mcens
    attr(Zmat,"ncens") <- ncens
    attr(Zmat,"M") <- M
    attr(Zmat,"N") <- N
    attr(Zmat,"polygonOverlay") <- polygonoverlay
    attr(Zmat,"pixelOverlay") <- pixeloverlay
    attr(Zmat,"FORM") <- formula
    attr(Zmat,"fftpoly") <- fftpoly
    attr(Zmat,"missingind") <- missingind
    
    return(Zmat)  
}

##' getpolyol function
##'
##' A function to perform polygon/polygon overlay operations and form the computational grid, on which inference will eventually take place.
##' For details and examples of using this fucntion, please see the package vignette "Bayesian_lgcp"
##'
##' @param data an object of class ppp or SpatialPolygonsDataFrame, containing the event counts, i.e. the dataset that will eventually be analysed 
##' @param regionalcovariates an object of class SpatialPolygonsDataFrame containng regionally measured covariate information 
##' @param pixelcovariates X an object of class SpatialPixelsDataFrame containng regionally measured covariate information
##' @param cellwidth the chosen cell width 
##' @param ext the amount by which to extend the observation window in forming the FFT grid, default is 2. In the case that the point pattern has long range spatial correlation, this may need to be increased.
##' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former, the default, includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window. 
##' @return an object of class lgcppolyol, which can then be fed into the function getZmat.
##' @seealso \link{minimum.contrast}, \link{minimum.contrast.spatiotemporal}, \link{chooseCellwidth}, \link{guessinterp}, \link{getZmat},
##' \link{addTemporalCovariates}, \link{lgcpPrior}, \link{lgcpInits}, \link{CovFunction}
##' \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}, \link{lgcpPredictSpatioTemporalPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars}
##' @export

getpolyol <- function(data,regionalcovariates=NULL,pixelcovariates=NULL,cellwidth,ext=2,inclusion="touching"){
    if(inherits(data,"SpatialPolygonsDataFrame")){
        spatstat.options(checkpolygons=FALSE)
        W <- as(gUnaryUnion(data),"owin")
        spatstat.options(checkpolygons=TRUE)
        sd <- ppp(window=W)         
    }
    else{    
        sd <- data
    }
    ow <- selectObsWindow(sd,cellwidth) 
	sd <- ow$xyt
	M <- ow$M # note for this function, M and N are powers of 2 
	N <- ow$N
	study.region <- sd$window
    gridobj <- genFFTgrid(study.region=study.region,M=M,N=N,ext=ext,inclusion=inclusion)
	del1 <- gridobj$del1
	del2 <- gridobj$del2
	Mext <- gridobj$Mext
	Next <- gridobj$Next
	mcens <- gridobj$mcens
	ncens <- gridobj$ncens
	cellarea <- gridobj$cellarea
	cellInside <- gridobj$cellInside
	
	W <-study.region
    mcens <- mcens[1:M]
    ncens <- ncens[1:N]
    cellInside <- cellInside[1:M,1:N]
    M <- length(mcens)
    N <- length(ncens)
    
    cat("Computing FFT grid ...\n") 
    fftpoly <- grid2spoly(mcens,ncens) # fft grid cells
    
    polyol <- NULL
    pixol <- NULL
    
    if(!is.null(regionalcovariates)){     
        cat("performing polygon overlay operations ...\n")
        polyol <- gOverlay(fftpoly,regionalcovariates) 
    }
    if(!is.null(pixelcovariates)){
        cat("converting SpatialPixelsDataFrame to SpatialPolygonsDataFrame ...\n")
        #oldclasses <- sapply(pixelcovariates@data,class)
        pixelcovariates <- as(pixelcovariates,"SpatialPolygonsDataFrame")
        #sapply(1:ncol(pixelcovariates),function(i){class(pixelcovariates[[i]]) <<- oldclasses[,i]}) 
        cat("performing polygon overlay operations ...\n")
        pixol <- gOverlay(fftpoly,pixelcovariates)
    }
    
    ans <- list()
    ans$gridobj <- gridobj
    ans$fftpoly <- fftpoly
    ans$polyol <- polyol
    ans$pixol <- pixol
    ans$mcens <- mcens
    ans$ncens <- ncens
    ans$cellwidth <- cellwidth
    ans$ext <- ext
    ans$inclusion <- inclusion
    
    class(ans) <- c("lgcppolyol","list")
    
    return(ans) 
}	

##' clearinterp function
##'
##' A function to remove the interpolation methods from a data frame.
##'
##' @param df a data frame
##' @return removes the interpolation methods
##' @export

clearinterp <- function(df){
    for(i in 1:length(df)){
        while(any(class(df[[i]])[1]==interptypes())){
            class(df[[i]]) <- class(df[[i]])[-1]
        }
    }
    return(df)
}

##' getinterp function
##'
##' A function to get the interpolation methods from a data frame\cr
##'
##' The three types of interpolation method employed in the package lgcp are:\cr
##' \enumerate{
##'    \item 'Majority' The interpolated value corresponds to the value of the covariate occupying 
##'        the largest area of the computational cell.
##'    \item 'ArealWeightedMean' The interpolated value corresponds to the mean of all covariate 
##'        values contributing to the computational cell weighted by their respective areas.
##'    \item 'ArealWeightedSum' The interpolated value is the sum of all contributing covariates 
##'        weighed by the proportion of area with respect to the covariate polygons. For example, 
##'        suppose region A has the same area as a computational grid cell and has 500 inhabitants. 
##'        If that region occupies half of a computational grid cell, then this interpolation type assigns 
##'        250 inhabitants from A to the computational grid cell.
##' }
##'
##' @param df a data frame
##' @return the interpolation methods
##' @export

getinterp <- function(df){
    int <- unlist(sapply(df,function(x){class(x)[1]}))
    int <- sapply(int,function(x){ifelse(!(any(x=="Majority")|any(x=="ArealWeightedMean")|any(x=="ArealWeightedSum")),NA,x)})
    return(int)
}

##' interptypes function
##'
##' A function to return the types of covariate interpolation available\cr
##'
##' The three types of interpolation method employed in the package lgcp are:\cr
##' \enumerate{
##'    \item 'Majority' The interpolated value corresponds to the value of the covariate occupying 
##'        the largest area of the computational cell.
##'    \item 'ArealWeightedMean' The interpolated value corresponds to the mean of all covariate 
##'        values contributing to the computational cell weighted by their respective areas.
##'    \item 'ArealWeightedSum' The interpolated value is the sum of all contributing covariates 
##'        weighed by the proportion of area with respect to the covariate polygons. For example, 
##'        suppose region A has the same area as a computational grid cell and has 500 inhabitants. 
##'        If that region occupies half of a computational grid cell, then this interpolation type assigns 
##'        250 inhabitants from A to the computational grid cell.
##' }
##'
##' @return character string of available interpolation types
##' @export

interptypes <- function(){
    return(c("Majority", "ArealWeightedMean","ArealWeightedSum"))
}

##' assigninterp function
##'
##' A function to assign an interpolation type to a variable in a data frame.\cr
##'
##' The three types of interpolation method employed in the package lgcp are:\cr
##' \enumerate{
##'    \item 'Majority' The interpolated value corresponds to the value of the covariate occupying 
##'        the largest area of the computational cell.
##'    \item 'ArealWeightedMean' The interpolated value corresponds to the mean of all covariate 
##'        values contributing to the computational cell weighted by their respective areas.
##'    \item 'ArealWeightedSum' The interpolated value is the sum of all contributing covariates 
##'        weighed by the proportion of area with respect to the covariate polygons. For example, 
##'        suppose region A has the same area as a computational grid cell and has 500 inhabitants. 
##'        If that region occupies half of a computational grid cell, then this interpolation type assigns 
##'        250 inhabitants from A to the computational grid cell.
##' }
##'
##' @param df a data frame
##' @param vars character vector giving name of variables
##' @param value an interpolation type, posssible options are given by typing interptypes(), see ?interptypes
##' @return assigns an interpolation type to a variable
##' @seealso \link{minimum.contrast}, \link{minimum.contrast.spatiotemporal}, \link{chooseCellwidth}, \link{getpolyol}, \link{guessinterp}, \link{getZmat},
##' \link{addTemporalCovariates}, \link{lgcpPrior}, \link{lgcpInits}, \link{CovFunction}
##' \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}, \link{lgcpPredictSpatioTemporalPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars}
##' @examples
##' \dontrun{spdf a SpatialPolygonsDataFrame}
##' \dontrun{spdf@@data <- assigninterp(df=spdf@@data,vars="pop",value="ArealWeightedSum")}
##' @export

assigninterp <- function(df,vars,value){
    df <- df # create a local copy
    if(!any(value==interptypes())){
        stop("invalid interptype, type 'interptypes()' to see possible list.")
    }
    vidx <- match(vars,names(df))
    if(any(is.na(vidx))){
        stop(paste("Variable(s)",paste(vars[which(is.na(vidx))],collapse=", ")," are not in the specified data frame"))
    }
    for (i in 1:length(vidx)){
        while(any(class(df[[vidx[i]]])[1]==interptypes())){ # strip off any existing interptypes
            class(df[[vidx[i]]]) <- class(df[[vidx[i]]])[-1]
        }
        class(df[[vidx[i]]]) <- c(value,class(df[[vidx[i]]]))
    }
    return(df)
}

##' guessinterp function
##'
##' A function to guess provisional interpolational methods to variables in a data frame. Numeric variables are assigned
##' interpolation by areal weighted mean (see below); factor, character and other types of variable are assigned
##' interpolation by majority vote (see below). Not that the interpolation type ArealWeightedSum is not assigned automatically.

##'
##' The three types of interpolation method employed in the package lgcp are:\cr
##' \enumerate{
##'    \item 'Majority' The interpolated value corresponds to the value of the covariate occupying 
##'        the largest area of the computational cell.
##'    \item 'ArealWeightedMean' The interpolated value corresponds to the mean of all covariate 
##'        values contributing to the computational cell weighted by their respective areas.
##'    \item 'ArealWeightedSum' The interpolated value is the sum of all contributing covariates 
##'        weighed by the proportion of area with respect to the covariate polygons. For example, 
##'        suppose region A has the same area as a computational grid cell and has 500 inhabitants. 
##'        If that region occupies half of a computational grid cell, then this interpolation type assigns 
##'        250 inhabitants from A to the computational grid cell.
##' }
##'
##'
##' @param df a data frame
##' @return the data frame, but with attributes describing the interpolation method for each variable
##' @seealso \link{minimum.contrast}, \link{minimum.contrast.spatiotemporal}, \link{chooseCellwidth}, \link{getpolyol}, \link{getZmat},
##' \link{addTemporalCovariates}, \link{lgcpPrior}, \link{lgcpInits}, \link{CovFunction}
##' \link{lgcpPredictSpatialPlusPars}, \link{lgcpPredictAggregateSpatialPlusPars}, \link{lgcpPredictSpatioTemporalPlusPars}, 
##' \link{lgcpPredictMultitypeSpatialPlusPars}
##' @examples 
##' \dontrun{spdf a SpatialPolygonsDataFrame}
##' \dontrun{spdf@@data <- guessinterp(spdf@@data)} 
##' @export

guessinterp <- function(df){
    df <- clearinterp(df)
    cl <- lapply(df,class)
    for(i in 1:length(cl)){
        if(any(cl[[i]]=="factor"|cl[[i]]=="character")){
            class(df[[i]]) <- c("Majority",class(df[[i]])) # interpolation via majority vote
        }
        else if(any(cl[[i]]=="numeric")){
            class(df[[i]]) <- c("ArealWeightedMean",class(df[[i]])) # example of this would be rainfall, so take area weighted mean ... weights are area(cell intersect region)/area(region)
        }
        else{
            class(df[[i]]) <- c("Majority",class(df[[i]])) # interpolation via majority vote e.g integer-valued variables
        }
        cat(names(df)[i],"interpolation via",class(df[[i]])[1],"\n")
    }
    return(df)
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
    m <- length(xgrid)
    n <- length(ygrid)
    return(SpatialPixels(SpatialPoints(cbind(rep(xgrid,n),rep(ygrid,each=m))),proj4string=proj4string))
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
    m <- length(xgrid)
    n <- length(ygrid)
    return(SpatialPoints(cbind(rep(xgrid,n),rep(ygrid,each=m)),proj4string=proj4string))
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
  spts <- SpatialPixels(SpatialPoints(cbind(rep(xgrid,n),rep(ygrid,each=m))),proj4string=proj4string)
  sps <- as(spts,"SpatialPolygons")
  spdf <- SpatialPolygonsDataFrame(sps,data=data.frame(grid=1:(m*n)),match.ID=FALSE)
  return(spdf)
}

