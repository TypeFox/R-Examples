
setMethod("show", "pixmap",
function(object){
    cat("Pixmap image\n")
    cat("  Type          :", class(object), "\n")
    cat("  Size          :", paste(object@size, collapse="x"), "\n")
    cat("  Resolution    :", paste(object@cellres, collapse="x"), "\n")
    cat("  Bounding box  :", object@bbox, "\n")
    if(is(object, "pixmapIndexed"))
        cat("  Nr. of colors :",
            length(unique(as(object@index, "vector"))), "of", 
            length(object@col), "\n")
    cat("\n")
})


setMethod("plot", "pixmap",
function(x, y, xlab="", ylab="", axes=FALSE, asp=1, ...){
    x = as(x, "pixmapIndexed")
    X <- seq(x@bbox[1], x@bbox[3], by=x@cellres[1])
    Y <- seq(x@bbox[2], x@bbox[4], by=x@cellres[2])
    
    image(x=X, y=Y, z=t(x@index[nrow(x@index):1,,drop=FALSE]), col=x@col,
          xlab=xlab, ylab=ylab, axes=axes, asp=asp, ...)
})
          


###**********************************************************


pixmap <- function(data=NULL, nrow=dim(data)[1],
                   ncol=dim(data)[2],
                   bbox=NULL, bbcent=FALSE, cellres=NULL)
{
    cellres <- rep(cellres, length=2)
    if(is.null(bbox)){
        if(is.null(cellres))
            cellres <- c(1,1)
        
        if(is.null(nrow)){
            if(!is.null(ncol))
                nrow <- ceiling(length(data)/ncol)
            else
                stop("Too few dimension attributes (nrow, ncol, bbox)\n")
        }
        else if(is.null(ncol))
            ncol <- ceiling(length(data)/nrow)
        
        if(bbcent)
            bbox <- c(1,1,cellres[1]*ncol, cellres[2]*nrow)
        else
            bbox <- c(0,0,cellres[1]*ncol, cellres[2]*nrow)
    }
    else{
        if(is.null(cellres)){
            if(is.null(nrow)){
                if(!is.null(ncol))
                    nrow <- ceiling(length(data)/ncol)
                else
                    stop("Too few dimension attributes (nrow, ncol, bbox)\n")
            }
            else if(is.null(ncol))
                ncol <- ceiling(length(data)/nrow)

            cellres = .getCellres(bbox, bbcent, c(nrow, ncol))
        }
        else{
            if(bbcent){
                ncol <- (bbox[3]-bbox[1])/cellres[1]+1
                nrow <- (bbox[4]-bbox[2])/cellres[2]+1
            }
            else{
                ncol <- (bbox[3]-bbox[1])/cellres[1]
                nrow <- (bbox[4]-bbox[2])/cellres[2]
            }
        }
    }
    
    new("pixmap", size=as(c(nrow, ncol),"integer"),
        cellres=cellres, bbox=bbox, bbcent=bbcent)
}


pixmapGrey = function(data, ...)
{
    z = new("pixmapGrey", pixmap(data, ...))
    
    datamax <- max(data)
    datamin <- min(data)
    data <- as.numeric(data)
    if(datamax>1 || datamin<0)
        data <- (data - datamin)/(datamax-datamin)

    z@grey = matrix(data, nrow=z@size[1], ncol=z@size[2])
    z
}

pixmapRGB = function(data, ...)
{
    z = new("pixmapRGB", pixmap(data, ...))
    
    datamax <- max(data)
    datamin <- min(data)
    data <- as.numeric(data)
    if(datamax>1 || datamin<0)
        data <- (data - datamin)/(datamax-datamin)

    data = array(data, dim=c(z@size[1], z@size[2], 3))

    z@red = matrix(data[,,1], nrow=z@size[1], ncol=z@size[2])
    z@green = matrix(data[,,2], nrow=z@size[1], ncol=z@size[2])
    z@blue = matrix(data[,,3], nrow=z@size[1], ncol=z@size[2])
    z
}

pixmapIndexed = function(data, col=NULL, ...)
{
    z = new("pixmapIndexed", pixmap(data, ...))
    data <- as(data, "integer")
    datamin <- min(data)
    if(datamin<=0)
        data <- data - datamin + 1
    datamax <- max(data)
    
    z@index =  matrix(data, nrow=z@size[1], ncol=z@size[2])
    if(is.null(col))
        col <- heat.colors(datamax)
    else{
        if(is(col,"function"))
            col <- col(datamax)
        else {
            if(length(col) < datamax){
                warning("number of of colors smaller than number of data values, recycling\n")
                col <- rep(col, length=datamax)
            }
        }
    }
    z@col = col
    z
}

###**********************************************************

setAs("pixmapGrey", "pixmapRGB",
function(from, to){
    z = new(to, as(from, "pixmap"))
    z@red = from@grey
    z@green = from@grey
    z@blue = from@grey
    z@channels = c("red", "green", "blue")
    z
})

setAs("pixmapRGB", "pixmapGrey",
function(from, to){
    addChannels(from)
})

setAs("pixmapRGB", "pixmapIndexed",
function(from, to){
    z = new(to, as(from, "pixmap"))
    x = rgb(from@red,from@green,from@blue)
    col <- unique(x)
    x <- match(x, col)
    z@index <- matrix(x, nrow=z@size[1], ncol=z@size[2])
    z@col = col
    z
})

setAs("pixmapGrey", "pixmapIndexed",
function(from, to){
    z = new(to, as(from, "pixmap"))
    x = grey(from@grey)
    col <- unique(x)
    x <- match(x, col)
    z@index <- matrix(x, nrow=z@size[1], ncol=z@size[2])
    z@col = col
    z
})

setAs("pixmapIndexed", "pixmapRGB",
function(from, to){
    z = new(to, as(from, "pixmap"))
    x <- col2rgb(from@col[from@index])/255
    z@red <- matrix(x["red",], nrow=z@size[1], ncol=z@size[2])
    z@green <- matrix(x["green",], nrow=z@size[1], ncol=z@size[2])
    z@blue <- matrix(x["blue",], nrow=z@size[1], ncol=z@size[2])
    z@channels = c("red", "green", "blue")
    z
})


## the fallbacks: convert to RGB and then to target

setAs("ANY", "pixmapGrey",
function(from, to){
    as(as(from, "pixmapRGB"), to)
})

setAs("ANY", "pixmapIndexed",
function(from, to){
    as(as(from, "pixmapRGB"), to)
})

###**********************************************************

setGeneric("addChannels",
function(object, coef=NULL) standardGeneric("addChannels"))


## coercion from RGB to Grey calls addChannels, hence be careful when
## using as() methods (danger of infinite loops).
setMethod("addChannels", "pixmapRGB",
function(object, coef=NULL){
    if(is.null(coef)) coef = c(0.30, 0.59, 0.11)
    z = new("pixmapGrey", object)
    z@grey = coef[1] * object@red + coef[2] * object@green +
        coef[3] * object@blue
    z@channels = "grey"
    z
})
         
setGeneric("getChannels",
function(object, colors="all") standardGeneric("getChannels"))

setMethod("getChannels", "pixmapChannels",
function(object, colors="all"){

    for(k in 1:length(colors))
        colors[k] = match.arg(colors[k], c("all", object@channels))
    if(any(colors=="all")) colors = object@channels
    colors = unique(colors)
    if(length(colors)>1){
        z = array(0, dim=c(object@size, length(colors)))
        dimnames(z) = list(NULL, NULL, colors)
        for(k in colors){
            z[,,k] = slot(object, k)
        }
    }
    else{
        z = slot(object, colors)
    }
    z
})
          
                          
###**********************************************************

setMethod("[", "pixmap",
function(x, i, j, ..., drop=FALSE){
    if(missing(j))
        j = TRUE
    if(missing(i)) 
        i = TRUE
    osize = x@size
    if(is(x, "pixmapIndexed")){
        x@index = x@index[i,j,drop=FALSE]
        x@size = dim(x@index)
    }
    else if(is(x, "pixmapChannels")){
        for(k in x@channels)
            slot(x, k) = slot(x, k)[i,j,drop=FALSE]
        x@size = dim(slot(x, k))
    }
    else
        stop(paste("Cannot subset objects of class", class(x)))
    
    
    ## now we re-calculate bounding box and cellres
    bbox = numeric(4)
    if(x@bbcent){
        b = seq(x@bbox[1], x@bbox[3], length=osize[2])
        bbox[c(1,3)] = range(b[j])
        b = seq(x@bbox[2], x@bbox[4], length=osize[1])
        bbox[c(2,4)] = range(b[i])
    }
    else{
        b = seq(x@bbox[1], x@bbox[3]-x@cellres[1], length=osize[2])
        bbox[1] = min(b[j])
        bbox[3] = max(b[j]) + x@cellres[1]
        b = seq(x@bbox[2], x@bbox[4]-x@cellres[2], length=osize[1])
        bbox[2] = min(b[i])
        bbox[4] = max(b[i]) + x@cellres[2]
    }
    
    x@bbox = bbox
    x@cellres <- .getCellres(bbox, x@bbcent, x@size)
    x
})

.getCellres = function(bbox, bbcent, size)
{
    if(bbcent)
        cellres = c((bbox[3]-bbox[1])/(size[2]-1),
                    (bbox[4]-bbox[2])/(size[1]-1))
    else
        cellres = c((bbox[3]-bbox[1])/size[2],
                    (bbox[4]-bbox[2])/size[1])

    cellres
}
          


      
      
      
