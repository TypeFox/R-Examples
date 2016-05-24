`ensemble.mean` <- function(
    RASTER.species.name="Species001", RASTER.stack.name="base",
    positive.filters=c("grd","_ENSEMBLE_"), negative.filters=c("xml"), 
    RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    KML.out=FALSE, KML.maxpixels=100000, KML.blur=10,
    p=NULL, a=NULL,
    pt=NULL, at=NULL,
    threshold=-1,
    threshold.method="spec_sens", threshold.sensitivity=0.9, threshold.PresenceAbsence=FALSE
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (threshold < 0) {
        if (is.null(p)==T || is.null(a)==T) {stop(paste("Please provide locations p and a to calculate thresholds", "\n", sep = ""))}
    }
    retest <- F
    if (is.null(pt)==F && is.null(at)==F) {
        if(identical(pt, p) == F || identical(at, a) == F)  {retest <- T}
    }

# new methods for calculating thresholds
    threshold2 <- function(eval, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, 
            threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=Pres, Abs=Abs) {
        if (threshold.PresenceAbsence == T){        
            Pres2 <- cbind(rep(1, length(Pres)), Pres)
            Abs2 <- cbind(rep(0, length(Abs)), Abs)
            data1 <- rbind(Pres2, Abs2)
            data2 <- cbind(seq(1:nrow(data1)), data1)
            auc.value <- PresenceAbsence::auc(data2, st.dev=F)
            cat(paste("\n", "AUC from PresenceAbsence package (also used to calculate threshold): ", auc.value, "\n", sep = ""))
            if (threshold.method=="kappa") {threshold.method <- "MaxKappa"}
            if (threshold.method=="spec_sens") {threshold.method <- "MaxSens+Spec"}
            if (threshold.method=="prevalence") {threshold.method <- "ObsPrev"}
            if (threshold.method=="equal_sens_spec") {threshold.method <- "Sens=Spec"}
            if (threshold.method=="sensitivity") {threshold.method <- "ReqSens"}
            req.sens <- threshold.sensitivity
            if (threshold.method=="no_omission") {
                threshold.method <- "ReqSens"
                req.sens <- 1.0
            }
            result <- PresenceAbsence::optimal.thresholds(data2, threshold=seq(from=0, to=1, by=0.005), req.sens=req.sens)
            result2 <- as.numeric(result[, 2])
            names(result2) <- result[, 1]
            if (threshold.method == "threshold.min") {
                t1 <- result2[["MaxSens+Spec"]]
                t2 <- result2[["Sens=Spec"]]            
                t3 <- result2[["ObsPrev"]]
                thresholds <- as.numeric(c(t1, t2, t3))
                thresholds <- thresholds[thresholds > 0]
                return(min(thresholds))
            }
            if (threshold.method == "threshold.mean") {
                t1 <- result2[["MaxSens+Spec"]]
                t2 <- result2[["Sens=Spec"]]            
                t3 <- result2[["ObsPrev"]]
                thresholds <- as.numeric(c(t1, t2, t3))
                thresholds <- thresholds[thresholds > 0]
                return(mean(thresholds))
            }
            return(as.numeric(result2[[threshold.method]]))
        }else{
            result <- dismo::threshold(eval, sensitivity=threshold.sensitivity)        
            if (threshold.method == "threshold.min") {
                t1 <- result[["spec_sens"]]
                t2 <- result[["equal_sens_spec"]]            
                t3 <- result[["prevalence"]]
                thresholds <- as.numeric(c(t1, t2, t3))
                thresholds <- thresholds[thresholds > 0]
                return(min(thresholds))
            }
            if (threshold.method == "threshold.mean") {
                t1 <- result[["spec_sens"]]
                t2 <- result[["equal_sens_spec"]]            
                t3 <- result[["prevalence"]]
                thresholds <- as.numeric(c(t1, t2, t3))
                thresholds <- thresholds[thresholds > 0]
                return(mean(thresholds))
            }
            return(result[[threshold.method]])
        }
    }
#

# avoid problems with non-existing directories
    dir.create("ensembles", showWarnings = F)
    dir.create("ensembles/count", showWarnings = F)
    dir.create("ensembles/presence", showWarnings = F)
    if(KML.out == T) {
        dir.create("kml", showWarnings = F)
        dir.create("kml/count", showWarnings = F)
        dir.create("kml/presence", showWarnings = F)
    }

# get ensemble files
# basic assumption is that different ensemble files are named as species_ENSEMBLE_1, species_ENSEMBLE_2, ... i.e. output from ensemble.batch

    species_focus <- RASTER.species.name
    if (gsub(".", "_", RASTER.species.name, fixed=T) != RASTER.species.name) {cat(paste("\n", "WARNING: species name (", RASTER.species.name, ") contains '.'", "\n\n", sep = ""))}
    ensemble.files <- list.files(path=paste(getwd(), "//ensembles", sep=""), pattern=species_focus, full.names=TRUE)
    if (length(ensemble.files) < 1) {
        cat(paste("\n", "NOTE: not meaningful to provide means as there are no raster files for this species:", RASTER.species.name, "\n", sep = ""))
        return(NULL)
    }
    RASTER.stack.name2 <- RASTER.stack.name
    if (gsub(".", "_", RASTER.stack.name, fixed=T) != RASTER.stack.name) {cat(paste("\n", "WARNING: title of stack (", RASTER.stack.name, ") contains '.'", "\n\n", sep = ""))}
    if (RASTER.stack.name != "") {
        ensemble.files <- ensemble.files[grepl(pattern=RASTER.stack.name, x=ensemble.files)]
        RASTER.stack.name2 <- paste("_", RASTER.stack.name, sep="")
        if (length(ensemble.files) < 1) {
            cat(paste("\n", "NOTE: not meaningful to provide means as there are no raster files for this stack:", RASTER.stack.name, "\n", sep = ""))
            return(NULL)
        }
    }
    for (i in 1:length(positive.filters)) {
        ensemble.files <- ensemble.files[grepl(pattern=positive.filters[i], x=ensemble.files)]
    }
    for (i in 1:length(negative.filters)) {
        ensemble.files <- ensemble.files[grepl(pattern=negative.filters[i], x=ensemble.files) == FALSE]
    }
    if (length(ensemble.files) < 2) {
        cat(paste("\n", "NOTE: not meaningful to provide means as there are fewer than 2 ensemble files", "\n", sep = ""))
        return(NULL)
    }
    cat(paste("\n", "Files used to create mean ensemble", "\n\n", sep = ""))
    print(ensemble.files)
    ensemble.stack <- raster::stack(ensemble.files)
    cat(paste("\n", "RasterStack used to create mean ensemble", "\n\n", sep = ""))
    print(ensemble.stack)
    ensemble.mean <- raster::mean(ensemble.stack)
#    if (is.na(crs(ensemble.mean)) == T) {
#        crs(ensemble.mean) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
#    }
    ensemble.mean <- trunc(1.0 * ensemble.mean)
#    raster::setMinMax(ensemble.mean)
    names(ensemble.mean) <- paste(species_focus, "_MEAN", RASTER.stack.name2, sep="")
    cat(paste("\n", "Mean ensemble (truncated)", "\n\n", sep = ""))
    print(ensemble.mean)
    filename1 <- paste(getwd(), "//ensembles//", species_focus, "_MEAN", RASTER.stack.name2, sep="")
    raster::writeRaster(x=ensemble.mean, filename=filename1, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
# avoid possible problems with saving of names of the raster layers
    raster::writeRaster(ensemble.mean, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(species_focus, "_MEAN", RASTER.stack.name2, sep="")
    raster::writeRaster(working.raster, filename=filename1, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
# standard deviation
    ensemble.sd <- raster::calc(ensemble.stack, fun=sd)
    ensemble.sd <- trunc(1.0 * ensemble.sd)
#    raster::setMinMax(ensemble.mean)
    names(ensemble.sd) <- paste(species_focus, "_SD", RASTER.stack.name2, sep="")
    cat(paste("\n", "Standard deviation for ensemble (truncated)", "\n\n", sep = ""))
    print(ensemble.sd)
    filename1 <- paste(getwd(), "//ensembles//", species_focus, "_SD", RASTER.stack.name2, sep="")
    raster::writeRaster(x=ensemble.sd, filename=filename1, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(ensemble.sd, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(species_focus, "_SD", RASTER.stack.name2, sep="")
    raster::writeRaster(working.raster, filename=filename1, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
# thresholds apply to probabilities
    ensemble.mean <- ensemble.mean / 1000
#
    threshold.mean <- threshold
    if (threshold.mean < 0) {
        eval1 <- NULL
        cat(paste("\n", "Evaluation of created mean ensemble raster layer at locations p and a", "\n", sep = ""))
        if (ncol(p) == 3) {p <- p[p[,1]==species_focus, c(2:3)]}
        if (ncol(a) == 3) {a <- a[a[,1]==species_focus, c(2:3)]}
        pres_consensus <- raster::extract(ensemble.mean, p)
        abs_consensus <- raster::extract(ensemble.mean, a)
        eval1 <- dismo::evaluate(p=pres_consensus, a=abs_consensus)
        print(eval1)
#        threshold.mean <- threshold(eval1, sensitivity=threshold.sensitivity)[[threshold.method]]
#        threshold.mean <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity)
        threshold.mean <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=pres_consensus, Abs=abs_consensus)
        cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
        print(as.numeric(threshold.mean))
        if(retest == T) {
            eval1 <- NULL
            cat(paste("\n", "Evaluation of created mean ensemble raster layer at locations pt and at", "\n\n", sep = ""))
            if (ncol(pt) == 3) {pt <- pt[pt[,1]==species_focus, c(2:3)]}
            if (ncol(at) == 3) {at <- at[at[,1]==species_focus, c(2:3)]}
            pres_consensus <- raster::extract(ensemble.mean, pt)
            abs_consensus <- raster::extract(ensemble.mean, at)
            eval1 <- dismo::evaluate(p=pres_consensus, a=abs_consensus)
            print(eval1)
        }
    }
#
    if (KML.out == T) {
        seq1 <- seq(from = 0, to = threshold.mean, length.out = 10)
        seq2 <- seq(from = threshold.mean, to = 1, length.out = 11)
        filename2 <- paste(getwd(), "//kml//", species_focus, "_MEAN", RASTER.stack.name2, sep="")
        raster::KML(ensemble.mean, filename=filename2, col = c(grDevices::rainbow(n = 10, start = 0, end = 1/6), grDevices::rainbow(n = 10, start = 3/6, end = 4/6)), colNA = 0, 
            blur=KML.blur, maxpixels=KML.maxpixels, overwrite=TRUE, breaks = c(seq1, seq2))
        sd.max <- raster::cellStats(ensemble.sd, stat='max')
        seq1 <- seq(from = 0, to = sd.max, length.out = 20)
        filename2b <- paste(getwd(), "//kml//", species_focus, "_SD", RASTER.stack.name2, sep="")
        raster::KML(ensemble.sd, filename=filename2b, col = grDevices::rainbow(n = 19, start = 0, end = 4/6), colNA = 0, 
            blur=KML.blur, maxpixels=KML.maxpixels, overwrite=TRUE, breaks = seq1)
    }
#
# presence-absence maps based on the mean maps
    enspresence <- ensemble.mean >= threshold.mean
#    if (is.na(crs(enspresence)) == T) {crs(enspresence) <- crs(ensemble.mean)}
#    enspresence <- trunc(enspresence)
    raster::setMinMax(enspresence)
    names(enspresence) <- paste(species_focus, "_MEAN", RASTER.stack.name2, "_presence", sep="")
    filename3 <- paste(getwd(), "//ensembles//presence//", species_focus, "_MEAN", RASTER.stack.name2, sep="")
    raster::writeRaster(x=enspresence, filename=filename3, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(enspresence, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(species_focus, "_MEAN", RASTER.stack.name2, "_presence", sep="")
    raster::writeRaster(working.raster, filename=filename3, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
    if (KML.out == T) {
        filename4 <- paste(getwd(), "//kml//presence//", species_focus, "_MEAN", RASTER.stack.name2, sep="")
        raster::KML(enspresence, filename=filename4, col=c("grey", "green"),
            colNA=0, blur=KML.blur, maxpixels=KML.maxpixels, overwrite=TRUE)
    }
#
# count maps: counting the number of ensembles predicting presence

    presence.files <- list.files(path=paste(getwd(), "//ensembles//presence", sep=""), pattern=species_focus, full.names=TRUE)
    if (RASTER.stack.name != "") {
        presence.files <- presence.files[grepl(pattern=RASTER.stack.name, x=presence.files)]
    }
    for (i in 1:length(positive.filters)) {
        presence.files <- presence.files[grepl(pattern=positive.filters[i], x=presence.files)]
    }
    for (i in 1:length(negative.filters)) {
        presence.files <- presence.files[grepl(pattern=negative.filters[i], x=presence.files) == FALSE]
    }

    ensemble.stack <- raster::stack(presence.files)
    cat(paste("\n", "RasterStack (presence-absence) used to create mean ensemble (count)", "\n\n", sep = ""))
    print(ensemble.stack)
    ensemble.count <- sum(ensemble.stack)
#    if (is.na(crs(ensemble.count)) == T) {crs(ensemble.count) <- crs(ensemble.mean)}
#    ensemble.count <- trunc(ensemble.count)
    raster::setMinMax(ensemble.count)
    names(ensemble.count) <- paste(species_focus, "_MEAN", RASTER.stack.name2, "_count", sep="")
    filename5 <- paste(getwd(), "//ensembles//count//", species_focus, "_MEAN", RASTER.stack.name2, sep="")
    raster::writeRaster(x=ensemble.count, filename=filename5, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(ensemble.count, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(species_focus, "_MEAN", RASTER.stack.name2, "_count", sep="")
    raster::writeRaster(working.raster, filename=filename5, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
    if (KML.out == T) {
        filename6 <- paste(getwd(), "//kml//count//", species_focus, "_MEAN", RASTER.stack.name2, sep="")
        nmax <- length(presence.files)
        if (nmax > 3) {
            raster::KML(ensemble.count, filename=filename6, col=c("grey", grDevices::rainbow(n=(nmax-1), start=0, end=1/3), "blue"),
                colNA=0, blur=10, overwrite=TRUE, breaks=seq(from=-1, to=nmax, by=1))
        }else{
            raster::KML(ensemble.count, filename=filename6, col=c("grey", grDevices::rainbow(n=nmax, start=0, end=1/3)),
                colNA=0, blur=10, overwrite=TRUE, breaks=seq(from=-1, to=nmax, by=1))
        }
    }
    return(list(threshold=threshold.mean, call=match.call() ))
}

