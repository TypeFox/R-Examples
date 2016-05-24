.TOCsp <- function(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, P=NA, Q=NA, progress=FALSE) {

if(!is.null(nthres) & !is.null(thres)) stop("Enter nthres OR thres as input to define thresholds, not both at the same time")

# extract cell values from the boolean and index maps 
boolval <- getValues(boolean)
indval <- getValues(index)

# extract cell values from the mask map if given
if(!is.null(mask)) mask <- getValues(mask)

# calculate population 
# mask out nodata cells in the index and boolean vectors if a mask vector is given
if(!is.null(mask)){
  mask[mask == NAval] <- NA
  boolval <- boolval*mask
}
# extract total number of cells with ones and zeros in the boolean vector
boolvals <- boolval[!is.na(boolval)]
ones.bool <- sum(as.bit(boolvals))
zeros.bool <- length(boolvals) - ones.bool

validPixels <- ones.bool + zeros.bool

population <- validPixels * res(index)[1] * res(index)[2]
if(!is.na(P) & !is.na(Q)){
  population <- P + Q
}

# extract map units for plotting purposes with plot.TOC
units <- paste("square", strsplit(strsplit(CRSargs(crs(index)), "+units=")[[1]][2], " ")[[1]][1])

tocd <- .TOCnosp(indval, boolval, mask=mask, nthres=nthres, thres=thres, NAval=NAval, progress=progress, 
                 ones.bool=ones.bool, zeros.bool=zeros.bool, population=population, units=units)

return(tocd) 
}
