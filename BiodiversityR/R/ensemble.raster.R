`ensemble.raster` <- function(
    xn=NULL, ext=NULL,
    models.list=NULL, 
    input.weights=models.list$output.weights,
    thresholds=models.list$thresholds,
    RASTER.species.name="Species001", RASTER.stack.name=xn@title, 
    RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    RASTER.models.overwrite=TRUE,
    KML.out=FALSE, KML.maxpixels=100000, KML.blur=10,
    evaluate=FALSE, SINK=FALSE,
    p=models.list$p, a=models.list$a,
    pt=models.list$pt, at=models.list$at
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (is.null(xn) == T) {stop("value for parameter xn is missing (RasterStack object)")}
    if (is.null(models.list) == T) {stop("provide 'models.list' as models will not be recalibrated and retested")}
    if (is.null(input.weights) == T) {input.weights <- models.list$output.weights}
    if (is.null(thresholds) == T) {stop("provide 'thresholds' as models will not be recalibrated and retested")}
    retest <- F
    if (evaluate == T) { 
        if (is.null(p)==T || is.null(a)==T) {
            cat(paste("\n", "NOTE: not possible to evaluate the models since locations p and a are not provided", "\n", sep = ""))
            evaluate <- F
        }
        if (is.null(pt)==F && is.null(at)==F) {
            if(identical(pt, p) == F || identical(at, a) == F)  {retest <- T}
        }
    }
    if (is.null(ext) == F) {
        if(length(xn@title) == 0) {xn@title <- "stack1"}
        title.old <- xn@title
        xn <- raster::crop(xn, y=ext, snap="in")
        xn@title <- title.old
    }

# create output file
    if (RASTER.species.name == "Species001") {
        RASTER.species.name <- models.list$species.name
    }
    dir.create("outputs", showWarnings = F)
    paste.file <- paste(getwd(), "/outputs/", RASTER.species.name, "_output.txt", sep="")
    OLD.SINK <- TRUE
    if (sink.number(type="output") == 0) {OLD.SINK <- F}
    if (SINK==T && OLD.SINK==F) {
        if (file.exists(paste.file) == F) {
            cat(paste("\n", "NOTE: results captured in file: ", paste.file, "\n", sep = ""))
        }else{
            cat(paste("\n", "NOTE: results appended in file: ", paste.file, "\n", sep = ""))
        }
        cat(paste("\n\n", "RESULTS (ensemble.raster function)", "\n", sep=""), file=paste.file, append=T)
        sink(file=paste.file, append=T)
        cat(paste(date(), "\n", sep=""))
        print(match.call())
    }

#
# check if all variables are present
    vars <- models.list$vars
    vars.xn <- names(xn)
    nv <- length(vars) 
    for (i in 1:nv) {
        if (any(vars.xn==vars[i]) == F) {stop("explanatory variable '", vars[i], "' not among grid layers of RasterStack xn \n", sep = "")}
    }
    nv <- length(vars.xn) 
    for (i in 1:nv) {
        if (any(vars==vars.xn[i]) == F) {
            cat(paste("\n", "NOTE: RasterStack layer '", vars.xn[i], "' was not calibrated as explanatory variable", "\n", sep = ""))
            xn <- raster::dropLayer(xn, which(names(xn) %in% c(vars.xn[i]) ))
        }
    }

#
# set minimum and maximum values for xn
    for (i in 1:raster::nlayers(xn)) {
        xn[[i]] <- raster::setMinMax(xn[[i]])
    }
    if(raster::projection(xn)=="NA") {
        raster::projection(xn) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    }
# declare categorical layers for xn
    factors <- models.list$factors
    categories <- NULL
    if(is.null(factors) == F) {
        for (i in 1:length(factors)) {
            j <- which(names(xn) == factors[i])
            xn[[j]] <- raster::as.factor(xn[[j]])
        }
        categories <- models.list$categories
    }
    dummy.vars <- models.list$dummy.vars 
#
    KML.blur <- trunc(KML.blur)
    if (KML.blur < 1) {KML.blur <- 1}
    if (is.null(input.weights)==F) {
        MAXENT <- max(c(input.weights["MAXENT"], -1), na.rm=T)
        GBM <- max(c(input.weights["GBM"], -1), na.rm=T)
        GBMSTEP <- max(c(input.weights["GBMSTEP"], -1), na.rm=T)
        RF <- max(c(input.weights["RF"], -1), na.rm=T)
        GLM <- max(c(input.weights["GLM"], -1), na.rm=T)
        GLMSTEP <- max(c(input.weights["GLMSTEP"], -1), na.rm=T)
        GAM <- max(c(input.weights["GAM"], -1), na.rm=T)
        GAMSTEP <- max(c(input.weights["GAMSTEP"], -1), na.rm=T)
        MGCV <- max(c(input.weights["MGCV"], -1), na.rm=T)
        MGCVFIX <- max(c(input.weights["MGCVFIX"], -1), na.rm=T)
        EARTH <- max(c(input.weights["EARTH"], -1), na.rm=T)
        RPART <- max(c(input.weights["RPART"], -1), na.rm=T)
        NNET <- max(c(input.weights["NNET"], -1), na.rm=T)
        FDA <- max(c(input.weights["FDA"], -1), na.rm=T)
        SVM <- max(c(input.weights["SVM"], -1), na.rm=T)
        SVME <- max(c(input.weights["SVME"], -1), na.rm=T)
        BIOCLIM <- max(c(input.weights["BIOCLIM"], -1), na.rm=T)
        DOMAIN <- max(c(input.weights["DOMAIN"], -1), na.rm=T)
        MAHAL<- max(c(input.weights["MAHAL"], -1), na.rm=T)
    }
    MAXENT.OLD <- GBM.OLD <- GBMSTEP.OLD <- RF.OLD <- GLM.OLD <- GLMSTEP.OLD <- GAM.OLD <- GAMSTEP.OLD <- MGCV.OLD <- NULL
    MGCVFIX.OLD <- EARTH.OLD <- RPART.OLD <- NNET.OLD <- FDA.OLD <- SVM.OLD <- SVME.OLD <- BIOCLIM.OLD <- DOMAIN.OLD <- MAHAL.OLD <- GEODIST.OLD <- NULL
# probit models, NULL if no probit model fitted
    MAXENT.PROBIT.OLD <- GBM.PROBIT.OLD <- GBMSTEP.PROBIT.OLD <- RF.PROBIT.OLD <- GLM.PROBIT.OLD <- GLMSTEP.PROBIT.OLD <- GAM.PROBIT.OLD <- GAMSTEP.PROBIT.OLD <- MGCV.PROBIT.OLD <- NULL
    MGCVFIX.PROBIT.OLD <- EARTH.PROBIT.OLD <- RPART.PROBIT.OLD <- NNET.PROBIT.OLD <- FDA.PROBIT.OLD <- SVM.PROBIT.OLD <- SVME.PROBIT.OLD <- BIOCLIM.PROBIT.OLD <- DOMAIN.PROBIT.OLD <- MAHAL.PROBIT.OLD <- NULL
    if (is.null(models.list) == F) {
        if (is.null(models.list$MAXENT) == F) {MAXENT.OLD <- models.list$MAXENT}
        if (is.null(models.list$GBM) == F) {GBM.OLD <- models.list$GBM}
        if (is.null(models.list$GBMSTEP) == F) {GBMSTEP.OLD <- models.list$GBMSTEP}
        if (is.null(models.list$RF) == F) {RF.OLD <- models.list$RF}
        if (is.null(models.list$GLM) == F) {GLM.OLD <- models.list$GLM}
        if (is.null(models.list$GLMSTEP) == F) {GLMSTEP.OLD <- models.list$GLMSTEP}
        if (is.null(models.list$GAM) == F) {GAM.OLD <- models.list$GAM}
        if (is.null(models.list$GAMSTEP) == F) {GAMSTEP.OLD <- models.list$GAMSTEP}
        if (is.null(models.list$MGCV) == F) {MGCV.OLD <- models.list$MGCV}
        if (is.null(models.list$MGCVFIX) == F) {MGCVFIX.OLD <- models.list$MGCVFIX}
        if (is.null(models.list$EARTH) == F) {EARTH.OLD <- models.list$EARTH}
        if (is.null(models.list$RPART) == F) {RPART.OLD <- models.list$RPART}
        if (is.null(models.list$NNET) == F) {NNET.OLD <- models.list$NNET}
        if (is.null(models.list$FDA) == F) {FDA.OLD <- models.list$FDA}
        if (is.null(models.list$SVM) == F) {SVM.OLD <- models.list$SVM}
        if (is.null(models.list$SVME) == F) {SVME.OLD <- models.list$SVME}
        if (is.null(models.list$BIOCLIM) == F) {BIOCLIM.OLD <- models.list$BIOCLIM}
        if (is.null(models.list$DOMAIN) == F) {DOMAIN.OLD <- models.list$DOMAIN}
        if (is.null(models.list$MAHAL) == F) {MAHAL.OLD <- models.list$MAHAL}
        if (is.null(models.list$GEODIST) == F) {GEODIST.OLD <- models.list$GEODIST}
# probit models
        if (is.null(models.list$MAXENT.PROBIT) == F) {MAXENT.PROBIT.OLD <- models.list$MAXENT.PROBIT}
        if (is.null(models.list$GBM.PROBIT) == F) {GBM.PROBIT.OLD <- models.list$GBM.PROBIT}
        if (is.null(models.list$GBMSTEP.PROBIT) == F) {GBMSTEP.PROBIT.OLD <- models.list$GBMSTEP.PROBIT}
        if (is.null(models.list$RF.PROBIT) == F) {RF.PROBIT.OLD <- models.list$RF.PROBIT}
        if (is.null(models.list$GLM.PROBIT) == F) {GLM.PROBIT.OLD <- models.list$GLM.PROBIT}
        if (is.null(models.list$GLMSTEP.PROBIT) == F) {GLMSTEP.PROBIT.OLD <- models.list$GLMSTEP.PROBIT}
        if (is.null(models.list$GAM.PROBIT) == F) {GAM.PROBIT.OLD <- models.list$GAM.PROBIT}
        if (is.null(models.list$GAMSTEP.PROBIT) == F) {GAMSTEP.PROBIT.OLD <- models.list$GAMSTEP.PROBIT}
        if (is.null(models.list$MGCV.PROBIT) == F) {MGCV.PROBIT.OLD <- models.list$MGCV.PROBIT}
        if (is.null(models.list$MGCVFIX.PROBIT) == F) {MGCVFIX.PROBIT.OLD <- models.list$MGCVFIX.PROBIT}
        if (is.null(models.list$EARTH.PROBIT) == F) {EARTH.PROBIT.OLD <- models.list$EARTH.PROBIT}
        if (is.null(models.list$RPART.PROBIT) == F) {RPART.PROBIT.OLD <- models.list$RPART.PROBIT}
        if (is.null(models.list$NNET.PROBIT) == F) {NNET.PROBIT.OLD <- models.list$NNET.PROBIT}
        if (is.null(models.list$FDA.PROBIT) == F) {FDA.PROBIT.OLD <- models.list$FDA.PROBIT}
        if (is.null(models.list$SVM.PROBIT) == F) {SVM.PROBIT.OLD <- models.list$SVM.PROBIT}
        if (is.null(models.list$SVME.PROBIT) == F) {SVME.PROBIT.OLD <- models.list$SVME.PROBIT}
        if (is.null(models.list$BIOCLIM.PROBIT) == F) {BIOCLIM.PROBIT.OLD <- models.list$BIOCLIM.PROBIT}
        if (is.null(models.list$DOMAIN.PROBIT) == F) {DOMAIN.PROBIT.OLD <- models.list$DOMAIN.PROBIT}
        if (is.null(models.list$MAHAL.PROBIT) == F) {MAHAL.PROBIT.OLD <- models.list$MAHAL.PROBIT}
    }
#
    if (MAXENT > 0) {
	jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
	if (!file.exists(jar)) {stop('maxent program is missing: ', jar, '\nPlease download it here: http://www.cs.princeton.edu/~schapire/maxent/')}
    }
    if (GBM > 0) {
        if (! requireNamespace("gbm")) {stop("Please install the gbm package")}
	requireNamespace("splines")
    }
    if (GBMSTEP > 0) {
#        if (! require(gbm)) {stop("Please install the gbm package")}
    }
    if (RF > 0) {
#        if (! require(randomForest)) {stop("Please install the randomForest package")}
    }
    if (GLMSTEP > 0) {
#        if (! require(MASS)) {stop("Please install the MASS package")}
    }
    if (GAM > 0  || GAMSTEP > 0) {
#        cat(paste("\n"))
#        try(detach(package:mgcv), silent=T)
#        suppressMessages(require(gam))
#        if (! require(gam)) {stop("Please install the gam package")}
    }
    if (MGCV > 0 || MGCVFIX > 0) {
#        cat(paste("\n"))
#        try(detach(package:gam), silent=T)
#        cat(paste("\n"))
#        options(warn=-1)
#        if (! require(mgcv)) {stop("Please install the mgcv package")}
#         get the probabilities from MGCV
            predict.mgcv <- function(object, newdata, type="response") {
                p <- mgcv::predict.gam(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
#        options(warn=0) 
    }
    if (EARTH > 0) {
#        if (! require(earth)) {stop("Please install the earth package")}
#         get the probabilities from earth
            predict.earth2 <- function(object, newdata, type="response") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
    }
    if (RPART > 0) {
#        if (! require(rpart)) {stop("Please install the rpart package")}
    }
    if (NNET > 0) {
#        if (! require(nnet)) {stop("Please install the nnet package")}
#         get the probabilities from nnet
            predict.nnet2 <- function(object, newdata, type="raw") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
    }
    if (FDA > 0) {
#        if (! require(mda)) {stop("Please install the mda package")}
    }
    if (SVM > 0) {
#        if (! require(kernlab)) {stop("Please install the kernlab package")}
    }
    if (SVME > 0) {
#        if (! require(e1071)) {stop("Please install the e1071 package")}
#         get the probabilities from svm
            predict.svme <- function(model, newdata, probability=T) {
                p <- predict(model, newdata, probability=T)
                return(attr(p, "probabilities")[,1])
             }
    }
    if (MAHAL > 0) {
        MAHAL.shape <- models.list$formulae$MAHAL.shape
#         get the probabilities from mahal
            predict.mahal <- function(model, newdata, MAHAL.shape) {
                p <- dismo::predict(object=model, x=newdata)
                p <- p - 1 - MAHAL.shape
                p <- abs(p)
                p <- MAHAL.shape / p
                return(p)
             }
    }
# 
    ws <- input.weights
    prediction.failures <- FALSE
#
# prepare for raster output
    dir.create("models", showWarnings = F)
    dir.create("ensembles", showWarnings = F)
    dir.create("ensembles/count", showWarnings = F)
    dir.create("ensembles/presence", showWarnings = F)
    stack.title <- RASTER.stack.name
    if (gsub(".", "_", stack.title, fixed=T) != stack.title) {cat(paste("\n", "WARNING: title of stack (", stack.title, ") contains '.'", "\n\n", sep = ""))}
#    stack.title <- xn@title
    if(KML.out == T) {
        dir.create("kml", showWarnings = F)
        dir.create("kml/count", showWarnings = F)
        dir.create("kml/presence", showWarnings = F)
    }
    rasterfull <- paste("ensembles/", RASTER.species.name, "_", stack.title , sep="")
    kmlfull <- paste("kml/", RASTER.species.name, "_", stack.title , sep="")
    raster.title <- paste(RASTER.species.name, "_", stack.title , sep="")
    rastercount <- paste("ensembles/count/", RASTER.species.name, "_", stack.title , sep="")
    kmlcount <- paste("kml/count/", RASTER.species.name, "_", stack.title , sep="")
    rasterpresence <- paste("ensembles/presence/", RASTER.species.name, "_", stack.title, sep="")
    kmlpresence <- paste("kml/presence/", RASTER.species.name, "_", stack.title, sep="")
    RASTER.species.orig <- RASTER.species.name
    if (RASTER.models.overwrite==T) {
        RASTER.species.name <- "working"
    }else{
        RASTER.species.name <- paste(RASTER.species.name, "_", stack.title, sep="")
    }
#
#
    cat(paste("\n", "Start of modelling for organism: ", RASTER.species.orig, "\n", sep = ""))
    cat(paste("Predictions for RasterStack: ", stack.title, "\n", sep = ""))
    ensemble.statistics <- NULL
    cat(paste("ensemble raster layers will be saved in folder ", getwd(), "/ensembles", "\n\n", sep = ""))
    statistics.names <- c("n.models", "ensemble.threshold", "ensemble.min", "ensemble.max", "count.min", "count.max") 
    ensemble.statistics <- numeric(6)
    names(ensemble.statistics) <- statistics.names
#
# sometimes still error warnings for minimum and maximum values of the layers
# set minimum and maximum values for xn
    for (i in 1:raster::nlayers(xn)) {
        xn[[i]] <- raster::setMinMax(xn[[i]])
    }
#
# since raster layers are scaled 0 - 1000, multiply the thresholds by 1000
    thresholds <- trunc(1000*thresholds)
    cat(paste("\n", "submodel thresholds used (scaled to raster output): ", "\n", sep = ""))
    print(thresholds)
#
# count models
    mc <- 0
#
# start raster layer creations
    if (ws["MAXENT"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Maximum entropy algorithm (package: dismo)\n", sep=""))
        # Put the file 'maxent.jar' in the 'java' folder of dismo
        # the file 'maxent.jar' can be obtained from from http://www.cs.princeton.edu/~schapire/maxent/.
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        results <- MAXENT.OLD
        pmaxent <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MAXENT", sep="")
        tryCatch(pmaxent <- raster::predict(object=results, x=xn, na.rm=TRUE, 
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("MAXENT prediction failed"))},
            silent=F)
        if (is.null(pmaxent) == F) {
            results2 <- MAXENT.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pmaxent, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "MAXENT"
                pmaxent <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pmaxent <- trunc(1000*pmaxent)
            raster::writeRaster(x=pmaxent, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pmaxent, p)/1000
                abs1 <- raster::extract(pmaxent, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pmaxent, pt)/1000
                abs1 <- raster::extract(pmaxent, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: MAXENT prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["MAXENT"] <- -1 
        }
    }
    if (ws["GBM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized boosted regression modeling (package: gbm) \n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- GBM.OLD
        pgbm <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GBM", sep="")
        tryCatch(pgbm <- raster::predict(object=xn, model=results, fun=gbm::predict.gbm, na.rm=TRUE, factors=categories,
                n.trees=results$n.trees, type="response", 
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("GBM prediction failed"))},
            silent=F)
        if (is.null(pgbm) == F) {
            results2 <- GBM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pgbm, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "GBM"
                pgbm <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pgbm <- trunc(1000*pgbm)
            raster::writeRaster(x=pgbm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pgbm, p)/1000
                abs1 <- raster::extract(pgbm, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pgbm, pt)/1000
                abs1 <- raster::extract(pgbm, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GBM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["GBM"] <- -1 
        }
    }
    if (ws["GBMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". gbm step algorithm (package: dismo)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- GBMSTEP.OLD
        pgbms <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GBMSTEP", sep="")
        tryCatch(pgbms <- raster::predict(object=xn, model=results, fun=gbm::predict.gbm, na.rm=TRUE, factors=categories,
                n.trees=results$n.trees, type="response", 
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("stepwise GBM prediction failed"))},
            silent=F)
        if (is.null(pgbms) == F) {
            results2 <- GBMSTEP.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pgbms, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "GBMSTEP"
                pgbms <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pgbms <- trunc(1000*pgbms)
            raster::writeRaster(x=pgbms, filename=fullname, progress='text', overwrite=TRUE)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pgbms, p)/1000
                abs1 <- raster::extract(pgbms, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pgbms, pt)/1000
                abs1 <- raster::extract(pgbms, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: stepwise GBM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["GBMSTEP"] <- -1 
        }
    }
    if (ws["RF"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Random forest algorithm (package: randomForest)\n", sep=""))
        results <- RF.OLD
        prf <- NULL
        fullname <- paste("models/", RASTER.species.name, "_RF", sep="")
        tryCatch(prf <- raster::predict(object=xn, model=results, na.rm=TRUE, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("random forest prediction failed"))},
            silent=F)
        if (is.null(prf) == F) {
            results2 <- RF.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=prf, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "RF"
                prf <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            prf <- trunc(1000*prf)
            raster::writeRaster(x=prf, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(prf, p)/1000
                abs1 <- raster::extract(prf, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(prf, pt)/1000
                abs1 <- raster::extract(prf, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: random forest prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["RF"] <- -1 
        }
    } 
    if (ws["GLM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Linear Model \n", sep=""))
        results <- GLM.OLD
        pglm <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GLM", sep="")
        tryCatch(pglm <- raster::predict(object=xn, model=results, na.rm=TRUE, type="response", factors=categories,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("GLM prediction failed"))},
            silent=F)
        if (is.null(pglm) == F) {
            results2 <- GLM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pglm, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "GLM"
                pglm <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pglm <- trunc(1000*pglm)
            raster::writeRaster(x=pglm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pglm, p)/1000
                abs1 <- raster::extract(pglm, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pglm, pt)/1000
                abs1 <- raster::extract(pglm, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GLM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["GLM"] <- -1 
        }
    }
    if (ws["GLMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Stepwise Generalized Linear Model \n", sep=""))
        results <- GLMSTEP.OLD
        pglms <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GLMSTEP", sep="")
        tryCatch(pglms <- raster::predict(object=xn, model=results, na.rm=TRUE, type="response", factors=categories,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("stepwise GLM prediction failed"))},
            silent=F)
        if (is.null(pglms) == F) {
            results2 <- GLMSTEP.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pglms, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "GLMSTEP"
                pglms <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pglms <- trunc(1000*pglms)
            raster::writeRaster(x=pglms, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pglms, p)/1000
                abs1 <- raster::extract(pglms, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pglms, pt)/1000
                abs1 <- raster::extract(pglms, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: stepwise GLM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["GLMSTEP"] <- -1 
        } 
    }
    if (ws["GAM"] > 0 || ws["GAMSTEP"] > 0) {
#        cat(paste("\n"))
#        try(detach(package:mgcv), silent=T)
#        suppressMessages(require(gam))
#        require(gam, quietly=T)
    }
    if (ws["GAM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Additive Model (package: gam)\n", sep=""))
        results <- GAM.OLD
        pgam <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GAM", sep="")
        tryCatch(pgam <- raster::predict(object=xn, model=results, fun=gam::predict.gam, na.rm=TRUE, type="response", factors=categories,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("GAM prediction (gam package) failed"))},
            silent=F)
        if (is.null(pgam) == F) {
            results2 <- GAM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pgam, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "GAM"
                pgam <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pgam <- trunc(1000*pgam)
            raster::writeRaster(x=pgam, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pgam, p)/1000
                abs1 <- raster::extract(pgam, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pgam, pt)/1000
                abs1 <- raster::extract(pgam, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GAM prediction (gam package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["GAM"] <- -1 
        } 
    }
    if (ws["GAMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Stepwise Generalized Additive Model (package: gam)\n", sep=""))
        results <- GAMSTEP.OLD
        pgams <- NULL
        fullname <- paste("models/", RASTER.species.name, "_GAMSTEP", sep="")
        tryCatch(pgams <- raster::predict(object=xn, model=results, fun=gam::predict.gam, type="response", na.rm=TRUE, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("stepwise GAM prediction (gam package) failed"))},
            silent=F)
        if (is.null(pgams) == F) {
            results2 <- GAMSTEP.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pgams, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "GAMSTEP"
                pgams <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pgams <- trunc(1000*pgams)
            raster::writeRaster(x=pgams, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pgams, p)/1000
                abs1 <- raster::extract(pgams, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pgams, pt)/1000
                abs1 <- raster::extract(pgams, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: stepwise GAM prediction (gam package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["GAMSTEP"] <- -1 
        } 
    }
    if (ws["MGCV"] > 0 || ws["MGCVFIX"] > 0) {
#        cat(paste("\n"))
#        try(detach(package:gam), silent=T)
#        options(warn=-1)
#        require(mgcv, quietly=T)
#        options(warn=0)
    }
    if (ws["MGCV"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Additive Model (package: mgcv)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- MGCV.OLD
        pmgcv <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MGCV", sep="")
        tryCatch(pmgcv <- raster::predict(object=xn, model=results, fun=predict.mgcv, na.rm=TRUE, type="response", factors=categories,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("GAM prediction (mgcv package) failed"))},
            silent=F)
        if (is.null(pmgcv) == F) {
            results2 <- MGCV.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pmgcv, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "MGCV"
                pmgcv <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pmgcv <- trunc(1000*pmgcv)
            raster::writeRaster(x=pmgcv, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pmgcv, p)/1000
                abs1 <- raster::extract(pmgcv, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pmgcv, pt)/1000
                abs1 <- raster::extract(pmgcv, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GAM prediction (mgcv package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["MGCV"] <- -1 
        } 
    }
    if (ws["MGCVFIX"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". GAM with fixed d.f. regression splines (package: mgcv)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- MGCVFIX.OLD
        pmgcvf <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MGCVFIX", sep="")
        tryCatch(pmgcvf <- raster::predict(object=xn, model=results, fun=predict.mgcv, na.rm=TRUE, type="response", factors=categories,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("MGCVFIX prediction (mgcv package) failed"))},
            silent=F)
        if (is.null(pmgcvf) == F) {
            results2 <- MGCVFIX.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pmgcvf, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "MGCVFIX"
                pmgcvf <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pmgcvf <- trunc(1000*pmgcvf)
            raster::writeRaster(x=pmgcvf, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pmgcvf, p)/1000
                abs1 <- raster::extract(pmgcvf, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pmgcvf, pt)/1000
                abs1 <- raster::extract(pmgcvf, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: GAM prediction (mgcv package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["MGCVFIX"] <- -1 
        } 
    }
    if (ws["EARTH"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Multivariate Adaptive Regression Splines (package: earth)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "NOTE: MARS (earth package) with factors may require explicit dummy variables", "\n", sep=""))
        } 
        results <- EARTH.OLD
        pearth <- NULL
        fullname <- paste("models/", RASTER.species.name, "_EARTH", sep="")
        tryCatch(pearth <- raster::predict(object=xn, model=results, fun=predict.earth2, na.rm=TRUE, type="response", factors=categories,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("MARS prediction (earth package) failed"))},
            silent=F)
        if (is.null(pearth) == F) {
            results2 <- EARTH.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pearth, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "EARTH"
                pearth <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pearth <- trunc(1000*pearth)
            raster::writeRaster(x=pearth, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pearth, p)/1000
                abs1 <- raster::extract(pearth, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pearth, pt)/1000
                abs1 <- raster::extract(pearth, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: MARS prediction (earth package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["EARTH"] <- -1 
        } 
    }
    if (ws["RPART"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Recursive Partitioning And Regression Trees (package: rpart)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- RPART.OLD
        prpart <- NULL
        fullname <- paste("models/", RASTER.species.name, "_RPART", sep="")
        tryCatch(prpart <- raster::predict(object=xn, model=results, na.rm=TRUE, type="prob", index=2, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("RPART prediction failed"))},
            silent=F)
        if (is.null(prpart) == F) {
            results2 <- RPART.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=prpart, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "RPART"
                prpart <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            prpart <- trunc(1000*prpart)
            raster::writeRaster(x=prpart, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(prpart, p)/1000
                abs1 <- raster::extract(prpart, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(prpart, pt)/1000
                abs1 <- raster::extract(prpart, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: RPART prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["RPART"] <- -1 
        } 
    }
    if (ws["NNET"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Artificial Neural Network (package: nnet)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- NNET.OLD
        pnnet <- NULL
        fullname <- paste("models/", RASTER.species.name, "_NNET", sep="")
        tryCatch(pnnet <- raster::predict(object=xn, model=results, fun=predict.nnet2, na.rm=TRUE, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("ANN prediction (nnet package) failed"))},
            silent=F)
        if (is.null(pnnet) == F) {
            results2 <- NNET.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pnnet, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "NNET"
                pnnet <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pnnet <- trunc(1000*pnnet)
            raster::writeRaster(x=pnnet, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pnnet, p)/1000
                abs1 <- raster::extract(pnnet, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pnnet, pt)/1000
                abs1 <- raster::extract(pnnet, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: ANN prediction (nnet package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["NNET"] <- -1 
        } 
    }
    if (ws["FDA"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Flexible Discriminant Analysis (package: mda)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "WARRNING: not certain whether the correct factor levels will be used", "\n", sep=""))
        } 
        results <- FDA.OLD
        pfda <- NULL
        fullname <- paste("models/", RASTER.species.name, "_FDA", sep="")
        tryCatch(pfda <- raster::predict(object=xn, model=results, na.rm=TRUE, type="posterior", index=2, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("FDA prediction failed"))},
            silent=F)
        if (is.null(pfda) == F) {
            results2 <- FDA.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pfda, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "FDA"
                pfda <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pfda <- trunc(1000*pfda)
            raster::writeRaster(x=pfda, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pfda, p)/1000
                abs1 <- raster::extract(pfda, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pfda, pt)/1000
                abs1 <- raster::extract(pfda, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: FDA prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["FDA"] <- -1 
        } 
    }

    if (ws["SVM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Support Vector Machines (package: kernlab)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "NOTE: SVM model with factors may require explicit dummy variables", "\n", sep=""))
        } 
        results <- SVM.OLD
        psvm <- NULL
        fullname <- paste("models/", RASTER.species.name, "_SVM", sep="")
        predict.svm2 <- as.function(kernlab::predict)
        tryCatch(psvm <- raster::predict(object=xn, model=results, fun=predict.svm2, na.rm=TRUE, type="probabilities", index=2, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("SVM prediction (kernlab package) failed"))},
            silent=F)
        if (is.null(psvm) == F) {
            results2 <- SVM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=psvm, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "SVM"
                psvm <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            psvm <- trunc(1000*psvm)
            raster::writeRaster(x=psvm, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(psvm, p)/1000
                abs1 <- raster::extract(psvm, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(psvm, pt)/1000
                abs1 <- raster::extract(psvm, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: SVM prediction (kernlab package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["SVM"] <- -1 
        } 
    }

    if (ws["SVME"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Support Vector Machines (package: e1071)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "NOTE: SVME model with factors may require explicit dummy variables", "\n", sep=""))
        }
        results <- SVME.OLD
        psvme <- NULL
        fullname <- paste("models/", RASTER.species.name, "_SVME", sep="")
        tryCatch(psvme <- raster::predict(object=xn, model=results, fun=predict.svme, na.rm=TRUE, factors=categories,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("SVM prediction (e1071 package) failed"))},
            warning= function(war) {print(paste("SVM prediction (e1071 package) failed"))},
            silent=F)
        if (is.null(psvme) == F) {
            results2 <- SVME.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=psvme, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "SVME"
                psvme <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            psvme <- trunc(1000*psvme)
            raster::writeRaster(x=psvme, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(psvme, p)/1000
                abs1 <- raster::extract(psvme, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(psvme, pt)/1000
                abs1 <- raster::extract(psvme, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: SVM prediction (e1071 package) failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["SVME"] <- -1 
        }
    }
    if (BIOCLIM > 0 || DOMAIN > 0 || MAHAL > 0) {
        if(is.null(factors) == F) {
            xn <- raster::dropLayer(xn, which(names(xn) %in% factors))               
        }
        if(is.null(dummy.vars) == F) {
            xn <- raster::dropLayer(xn, which(names(xn) %in% dummy.vars))               
        }
    }
    if (ws["BIOCLIM"] > 0) {  
        mc <- mc+1
        cat(paste("\n", mc, ". BIOCLIM algorithm (package: dismo)\n", sep=""))
        results <- BIOCLIM.OLD
        pbio <- NULL
        fullname <- paste("models/", RASTER.species.name, "_BIOCLIM", sep="")
        tryCatch(pbio <- dismo::predict(object=results, x=xn, na.rm=TRUE, 
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("BIOCLIM prediction failed"))},
            silent=F)
        if (is.null(pbio) == F) {
            results2 <- BIOCLIM.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pbio, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "BIOCLIM"
                pbio <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pbio <- trunc(1000*pbio)
            raster::writeRaster(x=pbio, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pbio, p)/1000
                abs1 <- raster::extract(pbio, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pbio, pt)/1000
                abs1 <- raster::extract(pbio, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: BIOCLIM prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["BIOCLIM"] <- -1 
        }
    }
    if (ws["DOMAIN"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". DOMAIN algorithm (package: dismo)\n", sep=""))
        results <- DOMAIN.OLD
        pdom <- NULL
        fullname <- paste("models/", RASTER.species.name, "_DOMAIN", sep="")
        tryCatch(pdom <- dismo::predict(object=results, x=xn, na.rm=TRUE, 
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("DOMAIN prediction failed"))},
            silent=F)
        if (is.null(pdom) == F) {
            results2 <- DOMAIN.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pdom, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "DOMAIN"
                pdom <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pdom <- trunc(1000*pdom)
            raster::writeRaster(x=pdom, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pdom, p)/1000
                abs1 <- raster::extract(pdom, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pdom, pt)/1000
                abs1 <- raster::extract(pdom, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: DOMAIN prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["DOMAIN"] <- -1 
        }
    }
    if (ws["MAHAL"] > 0) {  
        mc <- mc+1
        cat(paste("\n", mc, ". Mahalanobis algorithm (package: dismo)\n", sep=""))
        results <- MAHAL.OLD
        pmahal <- NULL
        fullname <- paste("models/", RASTER.species.name, "_MAHAL", sep="")
# not possible to use the predict.mahal function as raster::predict automatically reverts to dismo::predict for 'DistModel' objects
        tryCatch(pmahal <- dismo::predict(object=results, x=xn, na.rm=TRUE,
                filename=fullname, progress='text', overwrite=TRUE),
            error= function(err) {print(paste("Mahalanobis prediction failed"))},
            silent=F)
        if (is.null(pmahal) == F) {
            pmahal <- pmahal - 1 - MAHAL.shape
            pmahal <- abs(pmahal)
            pmahal <- MAHAL.shape / pmahal
            results2 <- MAHAL.PROBIT.OLD
            if (is.null(results2) == F) {
                cat(paste("Probit transformation", "\n", sep=""))
                fullname2 <- paste(fullname, "_step1", sep="")
                raster::writeRaster(x=pmahal, filename=fullname2, progress='text', overwrite=TRUE)
                explan.stack <- raster::stack(fullname2)
                names(explan.stack) <- "MAHAL"
                pmahal <- raster::predict(object=explan.stack, model=results2, na.rm=TRUE, type="response",
                    filename=fullname, progress='text', overwrite=TRUE)                
            }
            pmahal <- trunc(1000*pmahal)
            raster::writeRaster(x=pmahal, filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
            if(evaluate == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations p and a", "\n\n", sep = ""))
                pres1 <- raster::extract(pmahal, p)/1000
                abs1 <- raster::extract(pmahal, a)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
            if(retest == T) {
                eval1 <- pres1 <- abs1 <- NULL
                cat(paste("\n", "Evaluation at locations pt and at", "\n\n", sep = ""))
                pres1 <- raster::extract(pmahal, pt)/1000
                abs1 <- raster::extract(pmahal, at)/1000
                eval1 <- evaluate(p=pres1, a=abs1)
                print(eval1)
            }
        }else{
            cat(paste("\n", "WARNING: Mahalanobis prediction failed","\n\n", sep = ""))
            prediction.failures <- TRUE
            ws["MAHAL"] <- -1 
        }
    }
    if (prediction.failures == T) {
        cat(paste("\n", "WARNING: some predictions failed","\n", sep = ""))
        cat(paste("\n", "actual weights that were used were (-1 indicates failed predictions):","\n", sep = ""))
        print(ws)
        ws[which(ws==-1)] <- 0
    }
#
# create ensembles
    mc <- mc+1
    cat(paste("\n\n", mc, ". Ensemble algorithm\n", sep=""))
    ensemble.statistics["n.models"] <- sum(as.numeric(ws > 0))
    ensemble <- xn[[1]] == raster::NAvalue(xn[[1]])
    raster::setMinMax(ensemble)
    names(ensemble) <- raster.title
    raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
    enscount <- ensemble
    raster::setMinMax(enscount)
    names(enscount) <- paste(raster.title, "_count", sep="")
    raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    enspresence <- ensemble
    raster::setMinMax(enspresence)
    names(enspresence) <- paste(raster.title, "_presence", sep="")
    raster::writeRaster(x=enspresence, filename=rasterpresence, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    if (ws["MAXENT"] > 0) {
        ensemble <- ensemble + ws["MAXENT"] * pmaxent
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)       
        pmaxent <- pmaxent >= thresholds["MAXENT"]
        enscount <- enscount + pmaxent
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["GBM"] > 0) {
        ensemble <- ensemble + ws["GBM"] * pgbm
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        pgbm <- pgbm >= thresholds["GBM"]
        enscount <- enscount + pgbm
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["GBMSTEP"] > 0) {
        ensemble <- ensemble + ws["GBMSTEP"] * pgbms
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        pgbms <- pgbms >= thresholds["GBMSTEP"]
        enscount <- enscount + pgbms
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["RF"] > 0) {
        ensemble <- ensemble + ws["RF"] * prf
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        prf <- prf >= thresholds["RF"]
        enscount <- enscount + prf
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["GLM"] > 0) {
        ensemble <- ensemble + ws["GLM"] * pglm
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)    
        pglm <- pglm >= thresholds["GLM"]
        enscount <- enscount + pglm
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["GLMSTEP"] > 0) {
        ensemble <- ensemble + ws["GLMSTEP"] * pglms
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)      
        pglms <- pglms >= thresholds["GLMSTEP"]
        enscount <- enscount + pglms
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["GAM"] > 0) {
        ensemble <- ensemble + ws["GAM"] * pgam
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pgam <- pgam >= thresholds["GAM"]
        enscount <- enscount + pgam
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["GAMSTEP"] > 0) {
        ensemble <- ensemble + ws["GAMSTEP"] * pgams
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pgams <- pgams >= thresholds["GAMSTEP"]
        enscount <- enscount + pgams
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["MGCV"] > 0) {
        ensemble <- ensemble + ws["MGCV"] * pmgcv
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmgcv <- pmgcv >= thresholds["MGCV"]
        enscount <- enscount + pmgcv
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["MGCVFIX"] > 0) {
        ensemble <- ensemble + ws["MGCVFIX"] * pmgcvf
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmgcvf <- pmgcvf >= thresholds["MGCVFIX"]
        enscount <- enscount + pmgcvf
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["EARTH"] > 0) {
        ensemble <- ensemble + ws["EARTH"] * pearth
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pearth <- pearth >= thresholds["EARTH"]
        enscount <- enscount + pearth
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["RPART"] > 0) {
        ensemble <- ensemble + ws["RPART"] * prpart
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        prpart <- prpart >= thresholds["RPART"]
        enscount <- enscount + prpart
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["NNET"] > 0) {
        ensemble <- ensemble + ws["NNET"] * pnnet
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pnnet <- pnnet >= thresholds["NNET"]
        enscount <- enscount + pnnet
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["FDA"] > 0) {
        ensemble <- ensemble + ws["FDA"] * pfda
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pfda <- pfda >= thresholds["FDA"]
        enscount <- enscount + pfda
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["SVM"] > 0) {
        ensemble <- ensemble + ws["SVM"] * psvm
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        psvm <- psvm >= thresholds["SVM"]
        enscount <- enscount + psvm
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["SVME"] > 0) {
        ensemble <- ensemble + ws["SVME"] * psvme
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        psvme <- psvme >= thresholds["SVME"]
        enscount <- enscount + psvme
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["BIOCLIM"] > 0) {
        ensemble <- ensemble + ws["BIOCLIM"] * pbio
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pbio <- pbio >= thresholds["BIOCLIM"]
        enscount <- enscount + pbio
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["DOMAIN"] > 0) {
        ensemble <- ensemble + ws["DOMAIN"] * pdom
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pdom <- pdom >= thresholds["DOMAIN"]
        enscount <- enscount + pdom
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    if (ws["MAHAL"] > 0) {
        ensemble <- ensemble + ws["MAHAL"] * pmahal
        raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
        pmahal <- pmahal >= thresholds["MAHAL"]
        enscount <- enscount + pmahal
        raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
    }
    ensemble <- trunc(ensemble)
    raster::setMinMax(ensemble)
    ensemble.statistics["ensemble.min"] <- raster::minValue(ensemble)
    ensemble.statistics["ensemble.max"] <- raster::maxValue(ensemble)
#    names(ensemble) <- raster.title
#    raster::writeRaster(x=ensemble, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(ensemble, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- raster.title
    raster::writeRaster(working.raster, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
    if (KML.out == T) {
        thresholdx <- thresholds["ENSEMBLE"]
        seq1 <- seq(from = 0, to = thresholdx, length.out = 10)
        seq2 <- seq(from = thresholdx, to = 1000, length.out = 11)
        raster::KML(working.raster, filename=kmlfull, col = c(grDevices::rainbow(n = 10, start = 0, end = 1/6), grDevices::rainbow(n = 10, start = 3/6, end = 4/6)), colNA = 0, 
            blur=KML.blur, maxpixels=KML.maxpixels, overwrite=T, breaks = c(seq1, seq2))
    }
    raster::setMinMax(enscount)
    ensemble.statistics["count.min"] <- raster::minValue(enscount)
    ensemble.statistics["count.max"] <- raster::maxValue(enscount)
#    names(enscount) <- paste(raster.title, "_count", sep="")
#    raster::writeRaster(x=enscount, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(enscount, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(raster.title, "_count", sep="")
    raster::writeRaster(working.raster, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
    if (KML.out == T) {
        nmax <- sum(as.numeric(ws > 0))
        if (nmax > 3) {
            raster::KML(working.raster, filename=kmlcount, col=c("grey", grDevices::rainbow(n=(nmax-1), start=0, end=1/3), "blue"),
                colNA=0, blur=10, overwrite=T, breaks=seq(from=-1, to=nmax, by=1))
        }else{
            raster::KML(working.raster, filename=kmlcount, col=c("grey", grDevices::rainbow(n=nmax, start=0, end=1/3)),
                colNA=0, blur=10, overwrite=TRUE, breaks=seq(from=-1, to=nmax, by=1))
        }
    }
    ensemble.statistics["ensemble.threshold"] <- thresholds["ENSEMBLE"]
    enspresence <- ensemble >= thresholds["ENSEMBLE"]
    raster::setMinMax(enspresence)
#    names(enspresence) <- paste(raster.title, "_presence", sep="")
#    raster::writeRaster(x=enspresence, filename=rasterpresence, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(enspresence, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(raster.title, "_presence", sep="")
    raster::writeRaster(working.raster, filename=rasterpresence, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
    if (KML.out == T) {
        raster::KML(working.raster, filename=kmlpresence, col=c("grey", "green"),
            colNA=0, blur=KML.blur, maxpixels=KML.maxpixels, overwrite=T)
    }
    if(evaluate == T) {
        eval1 <- NULL
        cat(paste("\n", "Evaluation of created ensemble raster layer (", rasterfull, ") at locations p and a", "\n\n", sep = ""))
        pres_consensus <- raster::extract(ensemble, p)/1000
        abs_consensus <- raster::extract(ensemble, a)/1000
        eval1 <- evaluate(p=pres_consensus, a=abs_consensus)
        print(eval1)
    }
    if(retest == T) {
        eval1 <- NULL
        cat(paste("\n", "Evaluation of created ensemble raster layer (", rasterfull, ") at locations pt and at", "\n\n", sep = ""))
        pres_consensus <- raster::extract(ensemble, pt)/1000
        abs_consensus <- raster::extract(ensemble, at)/1000
        eval1 <- evaluate(p=pres_consensus, a=abs_consensus)
        print(eval1)
    }
    cat(paste("\n", "End of modelling for organism: ", RASTER.species.orig, "\n", sep = ""))
    cat(paste("Predictions were made for RasterStack: ", stack.title, "\n\n", sep = ""))
#
#  avoid possible problems with saving of names of the raster layers
    raster::writeRaster(ensemble, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- raster.title
    raster::writeRaster(working.raster, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
    raster::writeRaster(enscount, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(raster.title, "_count", sep="")
    raster::writeRaster(working.raster, filename=rastercount, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
    raster::writeRaster(enspresence, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(raster.title, "_presence", sep="")
    raster::writeRaster(working.raster, filename=rasterpresence, progress='text', overwrite=TRUE, format=RASTER.format, datatype="INT1U", NAflag=255)
#
    result <- list(ensemble.statistics=ensemble.statistics, call=match.call() )
    if (SINK==T && OLD.SINK==F) {sink(file=NULL, append=T)}
    return(result)
}

