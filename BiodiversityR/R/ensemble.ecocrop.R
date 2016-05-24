`ensemble.ecocrop.object` <- function(
    temp.thresholds, rain.thresholds, name="crop01", 
    temp.multiply=10, annual.temps=TRUE, transform=1)
{
    temps <- as.numeric(temp.thresholds[order(temp.thresholds)])
    temps <- temps*temp.multiply
    names(temps) <- c("tminabs", "tminopt", "tmaxopt", "tmaxabs")
    rains <- rain.thresholds[order(rain.thresholds)]
    names(rains) <- c("pminabs", "pminopt", "pmaxopt", "pmaxabs")
    ecocrop.object <- list(name=name, temp.thresholds=temps, rain.thresholds=rains, annual.temps=annual.temps, transform=transform)
    return(ecocrop.object)
}


`ensemble.ecocrop` <- function(
    x=NULL, ecocrop.object=NULL, 
    RASTER.object.name=ecocrop.object$name, RASTER.stack.name = x@title,
    RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    KML.out=TRUE, KML.blur=10, KML.maxpixels=100000 
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
    if(inherits(x, "RasterStack") == F) {stop("x is not a RasterStack object")}
    names(x)[which(names(x)=="bio01")] <- "bio1"
    names(x)[which(names(x)=="bio05")] <- "bio5"
    names(x)[which(names(x)=="bio06")] <- "bio6"
    if (is.null(ecocrop.object) == T) {stop("value for parameter ecocrop.object is missing (hint: use the ensemble.ecocrop.object function)")}
    vars <- names(x)
    if (any(vars == "bio12") == FALSE) {stop("Bioclimatic variable 'bio12' not provided with data")}
    if (ecocrop.object$annual.temps == F) {
        if (any(vars == "bio5") == FALSE) {stop("Bioclimatic variable 'bio5' not provided with data")}
        if (any(vars == "bio6") == FALSE) {stop("Bioclimatic variable 'bio6' not provided with data")}
    }else{
        if (any(vars == "bio1") == FALSE) {stop("Bioclimatic variable 'bio1' not provided with data")}
    }
# 
    predict.ecocrop <- function(object=ecocrop.object, newdata=newdata) {
        tminopt <- object$temp.thresholds["tminopt"]
        tminabs <- object$temp.thresholds["tminabs"]
        tmaxopt <- object$temp.thresholds["tmaxopt"]
        tmaxabs <- object$temp.thresholds["tmaxabs"]
        pminopt <- object$rain.thresholds["pminopt"]
        pminabs <- object$rain.thresholds["pminabs"]
        pmaxopt <- object$rain.thresholds["pmaxopt"]
        pmaxabs <- object$rain.thresholds["pmaxabs"]
        annual.temps <- object$annual.temps
        z <- object$transform

        result <- array(nrow(newdata))
        for (i in 1:nrow(newdata)) {
            datai <- newdata[i,,drop=F]
            P <- datai[, "bio12"]
            if (is.na(P) == T) {
                PS <- NA
            }else{
                PS1 <- PS2 <- PS3 <- 0 
                if (P > pminabs && P < pminopt) {PS1 <- (P-pminabs)/(pminopt-pminabs)}
                if (P >= pminopt && P <= pmaxopt) {PS2 <- 1.0}
                if (P > pmaxopt && P < pmaxabs) {PS3 <- (pmaxabs-P)/(pmaxabs-pmaxopt)}
                PS <- max(c(PS1, PS2, PS3))
            }
            if (annual.temps == F) {
                TMI <- datai[, "bio6"]
            }else{
                TMI <- datai[, "bio1"]
            }
            if (is.na(TMI) == T) {
                TMIS <- NA
            }else{
                TMI1 <- TMI2 <- 0
                if (TMI >= tminopt) {TMI1 <- 1.0}
                if (TMI > tminabs && TMI < tminopt) {TMI2 <- (TMI-tminabs)/(tminopt-tminabs)}
                TMIS <- max(c(TMI1, TMI2))
            }
            if (annual.temps == F) {
                TMX <- datai[, "bio5"]
            }else{
                TMX <- datai[, "bio1"]
            }
            if (is.na(TMX) == T) {
               TMXS <- NA
            }else{
                TMX1 <- TMX2 <- 0
                if (TMX <= tmaxopt) {TMX1 <- 1.0}
                if (TMX > tmaxopt && TMX < tmaxabs) {TMX2 <- (tmaxabs-TMX)/(tmaxabs-tmaxopt)}
                TMXS <- max(c(TMX1, TMX2))
            }
            SFINAL <- min(PS, TMIS, TMXS)
            result[i] <- SFINAL
        }
        p <- as.numeric(result)
        p <- p^z
        return(p)
    }
  
# avoid problems with non-existing directories and prepare for output
    dir.create("ensembles", showWarnings = F)
    dir.create("ensembles/ecocrop", showWarnings = F)
    if (KML.out == T) {
      dir.create("kml", showWarnings = F)
      dir.create("kml/ecocrop", showWarnings = F)
    }
    if(length(x@title) == 0) {x@title <- "stack1"}
    stack.title <- RASTER.stack.name
    rasterfull <- paste("ensembles/ecocrop/", RASTER.object.name, "_", stack.title , "_ecocrop", sep="")
    kmlfull <- paste("kml/ecocrop/", RASTER.object.name, "_", stack.title , "_ecocrop", sep="")
  
#
# predict
    tryCatch(ecocrop.raster <- raster::predict(object=x, model=ecocrop.object, fun=predict.ecocrop, na.rm=TRUE, 
                                           filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format),
           error= function(err) {print(paste("prediction of novel zones failed"))},
           silent=F)
    ecocrop.raster <- trunc(1000*ecocrop.raster)
    cat(paste("\n", "raster layer created (probabilities multiplied by 1000)", "\n", sep = ""))
    raster::setMinMax(ecocrop.raster)
    print(ecocrop.raster)

#
# avoid possible problems with saving of names of the raster layers
    raster::writeRaster(ecocrop.raster, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(RASTER.object.name, "_", stack.title , "_ecocrop", sep="")
    raster::writeRaster(working.raster, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#  
    if (KML.out == T) {
        seq1 <- seq(from = 0, to = 499, length.out = 10)
        seq2 <- seq(from = 500, to = 999, length.out = 10)
        raster::KML(working.raster, filename=kmlfull, col = c("grey", grDevices::rainbow(n = 9, start = 0, end = 1/6), grDevices::rainbow(n = 9, start = 3/6, end = 4/6), "green"), colNA = 0, 
            blur=KML.blur, maxpixels=KML.maxpixels, overwrite=T, breaks = c(seq1, seq2, 1001))
    }
  
    cat(paste("\n", "ecocrop raster provided in folder: ", getwd(), "//ensembles//ecocrop", "\n", sep=""))
    return(ecocrop.raster)
}
