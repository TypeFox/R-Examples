`ensemble.analogue.object` <- function(
    ref.location, future.stack, current.stack, name="reference1", 
    method="mahal", an=10000, probs=c(0.025, 0.975), weights=NULL, z=2)
{
    target.values <- as.numeric(raster::extract(future.stack, ref.location))
    names(target.values) <- names(future.stack)
    if ((method %in% c("mahal", "quantile", "sd")) == F) {method <- "none"} 
    if (method == "mahal") {
        a <- dismo::randomPoints(current.stack, n=an, p=NULL, excludep=F)
        a <- data.frame(a)
        background.data <- raster::extract(current.stack, a)
        background.data <- data.frame(background.data)
        TrainValid <- complete.cases(background.data)
        background.data <- background.data[TrainValid,]
        cov.mahal <- cov(background.data)
        out <- list(name=name, ref.location=ref.location, stack.name=future.stack@title, 
            method=method, target.values=target.values, cov.mahal=cov.mahal)
        return(out)
    }
    if (method == "quantile") {
        lower.interval <- data.frame(t(quantile(current.stack, probs[1])))
        upper.interval <- data.frame(t(quantile(current.stack, probs[2])))
        norm.values <- as.numeric(upper.interval - lower.interval)
    }
    if (method == "sd") {
        norm.values <- raster::cellStats(current.stack, stat="sd")
    }
    if (method == "none") {
        norm.values <- rep(1, length=length(names(current.stack)))
    }
    names(norm.values) <- names(current.stack)
# problem if some of the norm values are zero
    zero.norm.values <- which(norm.values == 0)
    if(length(zero.norm.values) > 0) {
        cat(paste("WARNING: some of the normalizing values were zero", "\n\n", sep=""))
        print(names(zero.norm.values))
        cat(paste("\n", "respective values were now set to one", "\n", sep=""))
        norm.values[names(norm.values) %in% names(zero.norm.values)] <- 1
    }
#    
    if (is.null(weights)==T || all.equal(names(weights), names(current.stack))==F) {
        paste("WARNING: length of weights is different from number of variables in current RasterStack", "\n", sep="") 
        weight.values <- rep(1, raster::nlayers(current.stack))
        names(weight.values) <- names(current.stack)
    }else{
        weight.values <- weights
    }
    weight.values <- weight.values / sum(weight.values)
    out <- list(name=name, ref.location=ref.location, stack.name=future.stack@title, 
        method=method, target.values=target.values, norm.values=norm.values, weights=weight.values, z=z)
    return(out)
}


`ensemble.analogue` <- function(
    x=NULL, analogue.object=NULL, analogues=1,
    RASTER.object.name=analogue.object$name, RASTER.stack.name=x@title,
    RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    KML.out=T, KML.blur=10, KML.maxpixels = 100000, 
    limits=c(1, 5, 20, 50), limit.colours=c('red', 'orange', 'blue', 'grey')
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
    if(inherits(x, "RasterStack") == F) {stop("x is not a RasterStack object")}
    if (is.null(analogue.object) == T) {stop("value for parameter analogue.object is missing (hint: use the ensemble.analogue.object function)")}
# 
    predict.analogue <- function(object=analogue.object, newdata=newdata) {
        method <- object$method
        if (method == "mahal") {
            centroid <- object$target.values
            cov.mahal <- object$cov.mahal
            p <- mahalanobis(newdata, center=centroid, cov=cov.mahal)
            p <- as.numeric(p)
            return(p)
        }else{
            targetdata <- object$target.values
	    normdata <- object$norm.values
	    weightdata <- object$weights
            z <- object$z
            out <- newdata
            for (i in 1:ncol(out)) {
	        out[,i] <- as.numeric(out[,i]) - as.numeric(targetdata[as.numeric(na.omit(match(names(out)[i], names(targetdata))))])
	        out[,i] <- abs(out[,i])
	        out[,i] <- as.numeric(out[,i]) / as.numeric(normdata[as.numeric(na.omit(match(names(out)[i], names(normdata))))])
	        out[,i] <- (out[,i]) ^ z
	        out[,i] <- as.numeric(out[,i]) * as.numeric(weightdata[as.numeric(na.omit(match(names(out)[i], names(weightdata))))])
	    }
	    p <- rowSums(out)
	    z2 <- 1/z
	    p <- p ^ z2
	    return(p)
        }
    }
  
# avoid problems with non-existing directories and prepare for output
    dir.create("ensembles", showWarnings = F)
    dir.create("ensembles/analogue", showWarnings = F)
    if (KML.out == T) {
      dir.create("kml", showWarnings = F)
      dir.create("kml/analogue", showWarnings = F)
    }
    if(length(x@title) == 0) {x@title <- "stack1"}
    stack.title <- RASTER.stack.name
    if (gsub(".", "_", stack.title, fixed=T) != stack.title) {cat(paste("\n", "WARNING: title of stack (", stack.title, ") contains '.'", "\n\n", sep = ""))}
    rasterfull <- paste("ensembles/analogue/", RASTER.object.name, "_", stack.title , "_analogue", sep="")
    kmlfull <- paste("kml/analogue/", RASTER.object.name, "_", stack.title , "_analogue", sep="")
  
#
# predict
    tryCatch(analogue.raster <- raster::predict(object=x, model=analogue.object, fun=predict.analogue, na.rm=TRUE, 
                                           filename=rasterfull, progress='text', overwrite=T, format=RASTER.format),
           error= function(err) {print(paste("prediction of analogue raster failed"))},
           silent=F)
    analogue.raster <- round(analogue.raster, digits=8)
    raster::setMinMax(analogue.raster)
    print(analogue.raster)

#
# avoid possible problems with saving of names of the raster layers
    raster::writeRaster(analogue.raster, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(RASTER.object.name, "_", stack.title , "_", analogue.object$method, "_z", analogue.object$z, "_analogue", sep="")
    raster::writeRaster(working.raster, filename=rasterfull, progress='text', overwrite=T, format=RASTER.format)
#  
    limits.data <- data.frame(limits)
    limits.data <- data.frame(cbind(limits, limits))
    names(limits.data) <- c("threshold", "count")
    freqs <- raster::freq(analogue.raster, digits=8)
    for (i in 1:length(limits)) {
        j <- 1
        while (sum(freqs[1:j, 2]) < limits[i]) {j <- j+1}
        limits.data[i, 1] <- freqs[j, 1]
        limits.data[i, 2] <- sum(freqs[1:j, 2])
    }
    cat(paste("\n", "suggested breaks in colour scheme", "\n", sep=""))	
    print(limits.data)
#
    if (KML.out == T) {
        breaks1 <- c(raster::minValue(analogue.raster), limits.data[,1])
        raster::KML(working.raster, filename=kmlfull, overwrite=T, blur=KML.blur, col=limit.colours, breaks=breaks1)
    }
#
# get locations
    j <- 1
    while (sum(freqs[1:j, 2]) < analogues) {j <- j+1}
    threshold <- freqs[j, 1]
    cat(paste("\n", "Threshold (n =", j, "): ", threshold, "\n", sep=""))	
    index1 <- which(analogue.raster[,] <= threshold)
    pres1 <- raster::xyFromCell(analogue.raster, index1)
    vars <- length(names(x))
    output1 <- data.frame(array(dim=c(length(index1), 5+vars)))
    output2 <- data.frame(array(dim=c(1, 5+vars)))
    names(output1) <- names(output2) <- c("model", "method", "lon", "lat", "distance", names(x))
    output1[, 1] <- rep(analogue.object$stack.name, nrow(output1))
    output2[1, 1] <- analogue.object$stack.name
    if (analogue.object$method == "mahal") {
        method1 <- "mahal"
    }else{
        method1 <- paste(analogue.object$method, "_z_", analogue.object$z, sep="") 
    }
    output1[, 2] <- rep(method1, nrow(output1))
    output1[, c(3:4)] <- pres1
    point.data1 <- raster::extract(x, pres1)
    output1[, c(6:(5+vars))] <- point.data1
    point.data2 <- raster::extract(analogue.raster, pres1)
    output1[, 5] <- point.data2
    output1 <- output1[order(output1[,"distance"], decreasing=F),]
    output2[1, 1] <- analogue.object$stack.name
    output2[1, 2] <- "target"
    output2[1 ,3] <- as.numeric(analogue.object$ref.location[,1])
    output2[1 ,4] <- as.numeric(analogue.object$ref.location[,2])
    output2[, c(6:(5+vars))] <- analogue.object$target.values      
    output3 <- rbind(output2, output1)
    cat(paste("\n", "analogue raster provided in folder: ", getwd(), "//ensembles//analogue", "\n\n", sep=""))
    return(output3)
}

