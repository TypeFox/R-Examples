# Purpose        : Count the number of vector objects for some GridTopology 
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl); 
# Contributions  : Edzer Pebesma (edzer.pebesma@uni-muenster.de);
# Dev Status     : Alpha
# Note           : Compare with aggregate in utils, sp, and spacetime;


## get a summary of an object for a list of lines:
count.GridTopology <- function(x, vectL, ...){
    if(any(class(x)=="GridTopology")){    
    # rasterize each line separately:
    sg <- SpatialGridDataFrame(x, proj4string = vectL[[1]]@proj4string, data=data.frame(observed=rep(NA, x@cells.dim[1]*x@cells.dim[2]), observed.sd=rep(NA, x@cells.dim[1]*x@cells.dim[2])))
    xv <- NULL
    for(i in 1:length(vectL)){
     vv <- vectL[[i]]
     vv$x <- rep(1, nrow(vv))
     xv[[i]]  <- vect2rast(vv, fname="x", cell.size=sg@grid@cellsize[1], bbox=sg@bbox, ...)@data
    }
    # bind all rasters:
    xv <- do.call(cbind, xv)
    
    sg$observed <- rowSums(xv, na.rm=T, dims=1)/length(xv)
    sg$observed.sd <- ifelse(sg$observed==0|sg$observed==1, 0, -sg$observed*log2(sg$observed)-(1-sg$observed)*log2(1-sg$observed))  ## information entropy (H) of a Bernoulli trial 
    out = new("SpatialVectorsSimulations", realizations = vectL, summaries = sg)
    return(out)
    } else {
      stop("Object of class 'GridTopology' required")
    }
}

## end of script;