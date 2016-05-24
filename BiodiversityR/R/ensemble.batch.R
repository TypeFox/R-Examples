`ensemble.batch` <- function(
    x=NULL, xn=c(x), ext=NULL, 
    species.presence=NULL, species.absence=NULL, 
    presence.min=20,
    an=1000, excludep=FALSE, CIRCLES.at=FALSE, CIRCLES.d=100000,
    k.splits=4, k.test=0, 
    n.ensembles=1, 
    SINK=FALSE,
    RASTER.format="raster", RASTER.datatype="INT2S", RASTER.NAflag=-32767,
    KML.out=FALSE, KML.maxpixels=100000, KML.blur=10,
    models.save=FALSE,
    threshold.method="spec_sens", threshold.sensitivity=0.9, threshold.PresenceAbsence=FALSE,
    ENSEMBLE.best=0, ENSEMBLE.min=0.7, ENSEMBLE.exponent=1,
    input.weights=NULL,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0,
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1, 
    PROBIT=FALSE, AUC.weights=TRUE,
    Yweights="BIOMOD", 
    layer.drops=NULL, factors=NULL, dummy.vars=NULL, 
    formulae.defaults=TRUE, maxit=100,
    MAXENT.a=NULL, MAXENT.an=10000, MAXENT.BackData=NULL, MAXENT.path=paste(getwd(), "/models/maxent", sep=""),
    GBM.formula=NULL, GBM.n.trees=2001,
    GBMSTEP.gbm.x=2:(1+raster::nlayers(x)), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100,
    RF.formula=NULL, RF.ntree=751, RF.mtry=floor(sqrt(raster::nlayers(x))), 
    GLM.formula=NULL, GLM.family=binomial(link="logit"), 
    GLMSTEP.steps=1000, STEP.formula=NULL, GLMSTEP.scope=NULL, GLMSTEP.k=2,
    GAM.formula=NULL, GAM.family=binomial(link="logit"), 
    GAMSTEP.steps=1000, GAMSTEP.scope=NULL, GAMSTEP.pos=1,
    MGCV.formula=NULL, MGCV.select=FALSE,
    MGCVFIX.formula=NULL, 
    EARTH.formula=NULL, EARTH.glm=list(family=binomial(link="logit"), maxit=maxit),
    RPART.formula=NULL, RPART.xval=50, 
    NNET.formula=NULL, NNET.size=8, NNET.decay=0.01,
    FDA.formula=NULL,
    SVM.formula=NULL, SVME.formula=NULL,
    MAHAL.shape=1  
)
{
    .BiodiversityR <- new.env()
    k.test <- as.integer(k.test)
    k.splits <- as.integer(k.splits)
    if (k.splits < 1) {
        cat(paste("\n", "NOTE: parameter k.splits was set to be smaller than 1", sep = ""))
        cat(paste("\n", "default value of 4 therefore set for parameter k.splits", sep = ""))
        k.splits <- 4
    }
    n.ensembles <- as.integer(n.ensembles)
    if (n.ensembles < 1) {n.ensembles <- 1}
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (is.null(xn) == T) {
        cat(paste("\n", "NOTE: new rasterStack assumed to be equal to the base rasterStack", sep = ""))
        xn <- x
    }
    xn <- c(xn)
# need to recalculate threshold for mean of ensembles
# therefore put x as first of new stacks
    if (n.ensembles > 1) {
        xn <- c(x, xn)
        i <- 1
        while (i < length(xn)) {
            i <- i+1
            if(identical(x, xn[[i]])) {xn[[i]] <- NULL}
        }
    }
    species.presence <- data.frame(species.presence)
    species.absence <- data.frame(species.absence)
    if (ncol(species.presence) < 2) {stop("species.presence expected to be 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns")}
    if (ncol(species.presence) == 2) {
        cat(paste("\n", "species.presence was expected to be 3-column data.frame with columns representing species, x (e.g., lon) and y (e.g., lat)", sep = ""))        
        cat(paste("\n", "only two columns were provided, it is therefore assumed that these reflect x and y coordinates for a single species", "\n\n", sep = ""))
        species.name <- rep("Species001", nrow(species.presence))
        species.presence <- cbind(species.name, species.presence)
        species.presence <- data.frame(species.presence)
        species.presence[,2] <- as.numeric(species.presence[,2])
        species.presence[,3] <- as.numeric(species.presence[,3])
    }
    if (ncol(species.presence) > 3) {
        cat(paste("\n", "species.presence was expected to be 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns", sep = ""))        
        cat(paste("\n", "only first three columns used", "\n\n", sep = ""))
        species.presence <- species.presence[,c(1:3)]
        species.presence[,2] <- as.numeric(species.presence[,2])
        species.presence[,3] <- as.numeric(species.presence[,3])
    }
    if (is.null(species.absence)==F && ncol(species.absence) < 2) {stop("species.absence expected to be a 2-column data.frame with x (e.g., lon) and y (e.g., lat),  or 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns")}
    if (is.null(species.absence)==F && ncol(species.absence)> 3) {
        cat(paste("\n", "species.absence was expected to be 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns", sep = ""))        
        cat(paste("\n", "only first three columns used", "\n\n", sep = ""))
        species.absence <- species.absence[,c(1:3)]
        species.absence[,2] <- as.numeric(species.absence[,2])
        species.absence[,3] <- as.numeric(species.absence[,3])
    }
    if (is.null(species.absence)==F && ncol(species.absence) == 2) {
        cat(paste("\n", "species.absence was expected to be 3-column data.frame with species, x (e.g., lon) and y (e.g., lat) columns", sep = ""))        
        cat(paste("\n", "only two columns were provided, it is therefore assumed that these reflect x and y coordinates for absence locations to be used for each species run", "\n\n", sep = ""))
        species.absence[,1] <- as.numeric(species.absence[,1])
        species.absence[,2] <- as.numeric(species.absence[,2])
        as <- species.absence
    }
# 
# process species by species
    species.names <- levels(droplevels(factor(species.presence[,1])))

    for (s in 1:length(species.names)) {
        focal.species <- species.names[s]

# check if species has required minimum number of presence points
        n.pres <- nrow(species.presence[species.presence[,1]==focal.species,])
        if (n.pres < presence.min) { 
            cat(paste("\n", "Species: ", focal.species, " only has ", n.pres, " presence locations", sep = ""))
            cat(paste("\n", "This species therefore not included in batch processing", "\n\n", sep = ""))

        }else{

# create output file
    if (s==1) {dir.create("outputs", showWarnings = F)}
    paste.file <- paste(getwd(), "/outputs/", focal.species, "_output.txt", sep="")
    OLD.SINK <- TRUE
    if (sink.number(type="output") == 0) {OLD.SINK <- F}
    if (SINK==T && OLD.SINK==F) {
        if (file.exists(paste.file) == F) {
            cat(paste("\n", "NOTE: results captured in file: ", paste.file, "\n", sep = ""))
        }else{
            cat(paste("\n", "NOTE: results appended in file: ", paste.file, "\n", sep = ""))
        }
        cat(paste("\n\n", "RESULTS (ensemble.batch function)", "\n", sep=""), file=paste.file, append=T)
        sink(file=paste.file, append=T)
        cat(paste(date(), "\n", sep=""))
        print(match.call())
    }

        cat(paste("\n", "Evaluations for species: ", focal.species, "\n", sep = ""))
        ps <- species.presence[species.presence[,1]==focal.species, c(2:3)]
        if (is.null(species.absence)==F && ncol(species.absence) == 3) {
            as <- species.absence[species.absence[,1]==focal.species, c(2:3)]
        }
        if (is.null(species.absence)==T) {
            if (excludep == T) {
                as <- dismo::randomPoints(x[[1]], n=an, p=ps, ext=ext, excludep=T)
            }else{
                as <- dismo::randomPoints(x[[1]], n=an, p=NULL, ext=ext, excludep=F)
            }
        }

    assign("ps", ps, envir=.BiodiversityR)
    assign("as", as, envir=.BiodiversityR)


# repeat the whole process for n.ensembles

    RASTER.species.name1 <- focal.species
    for (runs in 1:n.ensembles) {
        if (n.ensembles > 1) { 
            cat(paste("\n", focal.species, ": ENSEMBLE ", runs, "\n\n", sep = ""))
            RASTER.species.name1 <- paste(focal.species, "_ENSEMBLE_", runs, sep="")
        }

#1. first ensemble tests
    calibration.1 <- ensemble.test.splits(x=x, p=ps, a=as, ext=ext, k=k.splits, 
        CIRCLES.at=CIRCLES.at, CIRCLES.d=CIRCLES.d,
        ENSEMBLE.tune=T,
        ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min, 
        ENSEMBLE.exponent=ENSEMBLE.exponent,
        species.name = RASTER.species.name1,
        threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence,
        input.weights=input.weights,
        MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
        GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
        NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,
        PROBIT=PROBIT, VIF=T,
        Yweights=Yweights, 
        layer.drops=layer.drops, factors=factors, dummy.vars=dummy.vars,
        maxit=maxit,
        MAXENT.a=MAXENT.a, MAXENT.an=MAXENT.an, 
        MAXENT.BackData=MAXENT.BackData, MAXENT.path=MAXENT.path,
        GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees,
        GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
        GBMSTEP.step.size=GBMSTEP.step.size,
        RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
        GLM.formula=GLM.formula, GLM.family=GLM.family, 
        GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=STEP.formula, GLMSTEP.scope=GLMSTEP.scope, 
        GAM.formula=GAM.formula, GAM.family=GAM.family, 
        GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope, GAMSTEP.pos=GAMSTEP.pos,
        MGCV.formula=MGCV.formula, MGCV.select=MGCV.select,
        MGCVFIX.formula=MGCVFIX.formula, 
        EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm,
        RPART.formula=RPART.formula, RPART.xval=RPART.xval, 
        NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay,
        FDA.formula=FDA.formula, SVM.formula=SVM.formula, SVME.formula=SVME.formula,
        MAHAL.shape=MAHAL.shape)

#2. calibrate final model

#    xn.f <- eval(as.name(xn.focal))
    cat(paste("\n", "Final model calibrations for species: ", RASTER.species.name1,  "\n", sep = ""))
    cat(paste("\n", "Minimum input weight is 0.05", "\n", sep=""))

    if (AUC.weights == TRUE) {
        output.weights <- calibration.1$output.weights.AUC
    }else{
        output.weights <- calibration.1$output.weights
    }
    output.weights[output.weights < 0.05] <- 0
    print(output.weights)

    if (sum(output.weights) > 0) {

    calibration.2 <- ensemble.test(
        x=x, p=ps, a=as, ext=ext, k=k.test, pt=NULL, at=NULL,
        models.keep=TRUE, evaluations.keep=TRUE,
        PLOTS=F,
        models.save=models.save, species.name=RASTER.species.name1,
        AUC.weights=F, ENSEMBLE.tune=F,
        input.weights=output.weights,
        threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence,
        RASTER.format=RASTER.format,
        PROBIT=PROBIT, VIF=T,
        Yweights=Yweights, 
        layer.drops=layer.drops, factors=factors, dummy.vars=dummy.vars,
        maxit=maxit,
        MAXENT.BackData=MAXENT.BackData, MAXENT.path=MAXENT.path,
        GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees,
        GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
        GBMSTEP.step.size=GBMSTEP.step.size,
        RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
        GLM.formula=GLM.formula, GLM.family=GLM.family, 
        GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=STEP.formula, GLMSTEP.scope=GLMSTEP.scope, 
        GAM.formula=GAM.formula, GAM.family=GAM.family, 
        GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope, GAMSTEP.pos=GAMSTEP.pos,
        MGCV.formula=MGCV.formula, MGCV.select=MGCV.select,
        MGCVFIX.formula=MGCVFIX.formula, 
        EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm,
        RPART.formula=RPART.formula, RPART.xval=RPART.xval, 
        NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay,
        FDA.formula=FDA.formula, SVM.formula=SVM.formula, SVME.formula=SVME.formula,
        MAHAL.shape=MAHAL.shape) 


#3. predict for all the other rasters

    for (n in 1:length(xn)) {
        xn.f <- xn[[n]]
        if(length(xn.f@title) == 0) {xn.f@title <- paste("stack", n, sep="")}
        if (gsub(".", "_", xn.f@title, fixed=T) != xn.f@title) {cat(paste("\n", "WARNING: title of stack (", xn.f@title, ") contains '.'", "\n\n", sep = ""))}
        cat(paste("\n", "Predictions for species: ", RASTER.species.name1, " for rasterStack: ", xn.f@title, sep = ""))
        tryCatch(rasters2 <- ensemble.raster(xn=xn.f, ext=ext,
            models.list=calibration.2$models,            
            RASTER.species.name=RASTER.species.name1, 
            RASTER.format=RASTER.format, RASTER.datatype=RASTER.datatype, RASTER.NAflag=RASTER.NAflag,
            KML.out=KML.out, KML.maxpixels=KML.maxpixels, KML.blur=KML.blur),
                error= function(err) {print(paste("WARNING: prediction failed for stack: ", xn.f@title, sep=""))},
                silent=T)

        if(runs==n.ensembles && n.ensembles>1 && RASTER.format=="raster") {

# recalculate threshold for mean of predictions with calibration stack (xn[[1]])

            if (n == 1) {
                calibrate.mean <- NULL
                calibrate.mean <- ensemble.mean(RASTER.species.name=focal.species, RASTER.stack.name=xn.f@title,
                    positive.filters = c("grd", "_ENSEMBLE_"), negative.filters = c("xml"), 
                    RASTER.format=RASTER.format, RASTER.datatype=RASTER.datatype, RASTER.NAflag=RASTER.NAflag,
                    KML.out=KML.out, KML.maxpixels=KML.maxpixels, KML.blur=KML.blur,
                    p=ps, a=as,
                    pt = NULL, at = NULL,
                    threshold = -1,
                    threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence)
                cat(paste("\n", "threshold for mean suitability: ", calibrate.mean$threshold, "\n", sep = ""))
            }else{
                ensemble.mean(RASTER.species.name=focal.species, RASTER.stack.name=xn.f@title,
                    positive.filters = c("grd", "_ENSEMBLE_"), negative.filters = c("xml"), 
                    RASTER.format=RASTER.format, RASTER.datatype=RASTER.datatype, RASTER.NAflag=RASTER.NAflag,
                    KML.out=KML.out, KML.maxpixels=KML.maxpixels, KML.blur=KML.blur,
                    p=NULL, a=NULL,
                    pt = NULL, at = NULL,
                    threshold = calibrate.mean$threshold,
                    threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence)
            }
        }
    }


# sum output weights > 0 loop
    }

# n ensembles loop
    }

# if (sufficient presence locations) loop
    if (SINK==T && OLD.SINK==F) {sink(file=NULL, append=T)}
    }

# s (species) loop
    }

    result <- list(species=species.names, call=match.call())
    return(result)

}
