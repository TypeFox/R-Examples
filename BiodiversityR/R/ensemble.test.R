`ensemble.test` <- function(
    x=NULL, p=NULL, a=NULL, an=1000, excludep=FALSE, ext=NULL, 
    k=0, pt=NULL, at=NULL, CIRCLES.at=FALSE, CIRCLES.d=100000,
    TrainData=NULL, TestData=NULL,
    VIF=FALSE, COR=FALSE,
    SINK=FALSE, PLOTS=TRUE, 
    threshold.method="spec_sens", threshold.sensitivity=0.9, threshold.PresenceAbsence=FALSE,
    evaluations.keep=FALSE, 
    models.list=NULL, models.keep=FALSE, 
    models.save=FALSE, species.name="Species001",
    AUC.weights=TRUE, ENSEMBLE.tune=FALSE, 
    ENSEMBLE.best=0, ENSEMBLE.min=0.7, ENSEMBLE.exponent=1.0,
    input.weights=NULL, 
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0,
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1, 
    GEODIST=0, 
    PROBIT=FALSE,
    Yweights="BIOMOD", 
    layer.drops=NULL, factors=NULL, dummy.vars=NULL,
    formulae.defaults=TRUE, maxit=100,
    MAXENT.a=NULL, MAXENT.an=10000, MAXENT.BackData=NULL, 
    MAXENT.path=paste(getwd(), "/models/maxent_", species.name,  sep=""),
    GBM.formula=NULL, GBM.n.trees=2001, 
    GBMSTEP.gbm.x=2:(ncol(TrainData.vars)+1), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100, 
    RF.formula=NULL, RF.ntree=751, RF.mtry=floor(sqrt(ncol(TrainData.vars))),
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
    SVM.formula=NULL,
    SVME.formula=NULL, 
    MAHAL.shape=1,
    RASTER.format="raster"
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    k <- as.integer(k)
# check data
    if (is.null(TrainData) == T) {
        if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
        if(inherits(x,"RasterStack") == F) {stop("x is not a RasterStack object")}
        if(raster::projection(x)=="NA") {
            raster::projection(x) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        }
        if(is.null(p) == T) {stop("presence locations are missing (parameter p)")}
    }
# geoDist requires presence locations
    if ((GEODIST > 0) && (is.null(p) == T)) {stop("presence locations are missing for geoDist")}

#
    if(models.save==T) {
        models.keep <- TRUE
        dir.create("models", showWarnings = F)
    }

# create output file
    dir.create("outputs", showWarnings = F)
    paste.file <- paste(getwd(), "/outputs/", species.name, "_output.txt", sep="")
    OLD.SINK <- TRUE
    if (sink.number(type="output") == 0) {OLD.SINK <- F}
    if (SINK==T && OLD.SINK==F) {
        if (file.exists(paste.file) == F) {
            cat(paste("\n", "NOTE: results captured in file: ", paste.file, "\n", sep = ""))
        }else{
            cat(paste("\n", "NOTE: results appended in file: ", paste.file, "\n", sep = ""))
        }
        cat(paste(sep="\n\n", "RESULTS (ensemble.test function)"), file=paste.file, sep="\n", append=T)
        sink(file=paste.file, append=T)
        cat(paste(date()), sep="\n")
        print(match.call())
        cat(paste("   "), sep="\n")
    }

# check TrainData
    if (is.null(TrainData) == F) {
        TrainData <- data.frame(TrainData)
        if (colnames(TrainData)[1] !="pb") {stop("first column for TrainData should be 'pb' containing presence (1) and absence (0) data")}
        if ((is.null(x) == F) && (raster::nlayers(x) != (ncol(TrainData)-1))) {
            cat(paste("\n", "WARNING: different number of explanatory variables in rasterStack and TrainData", sep = ""))
        }
    }

# modify list of variables
# if TrainData is provided, then this data set takes precedence over raster x in the selection of variables

# 
# modify TrainData if layer.drops
    if (is.null(TrainData) == F) {
        if (is.null(layer.drops) == F) {
            vars <- names(TrainData)
            layer.drops <- as.character(layer.drops)
            factors <- as.character(factors)
            dummy.vars <- as.character(dummy.vars)
            nd <- length(layer.drops)
            for (i in 1:nd) {     
                if (any(vars==layer.drops[i]) == FALSE) {
                    cat(paste("\n", "WARNING: variable to exclude '", layer.drops[i], "' not among columns of TrainData", "\n", sep = ""))
                }else{
                    cat(paste("\n", "NOTE: variable '", layer.drops[i], "' will not be included as explanatory variable", "\n", sep = ""))
                    TrainData <- TrainData[, which(colnames(TrainData) != layer.drops[i])]
                    if (is.null(TestData) == F) {TestData <- TestData[, which(colnames(TestData) != layer.drops[i])]}
                    vars <- colnames(TrainData)
                    if (length(factors) > 0) {
                        factors <- factors[factors != layer.drops[i]]                        
                    }
                    if (length(dummy.vars) > 0) {
                        dummy.vars <- dummy.vars[dummy.vars != layer.drops[i]]
                    }
                }
            }
            if(length(layer.drops) == 0) {layer.drops <- NULL} 
            if(length(factors) == 0) {factors <- NULL}
            if(length(dummy.vars) == 0) {dummy.vars <- NULL}
        }
        if (is.null(factors) == F) {
            vars <- names(TrainData)
            factors <- as.character(factors)
            nf <- length(factors)
            old.factors <- factors
            for (i in 1:nf) {
                if (any(vars==old.factors[i]) == FALSE) {
                    cat(paste("\n", "WARNING: categorical variable '", old.factors[i], "' not among columns of TrainData", "\n", sep = ""))
                    factors <- factors[factors != old.factors[i]]
                }
            }
            if(length(factors) == 0) {factors <- NULL}
        }
        if (is.null(dummy.vars) == F) {
            vars <- names(TrainData)
            dummy.vars <- as.character(dummy.vars)
            nf <- length(dummy.vars)
            old.dummy.vars <- dummy.vars
            for (i in 1:nf) {
                if (any(vars==old.dummy.vars[i]) == FALSE) {
                    cat(paste("\n", "WARNING: dummy variable '", old.dummy.vars[i], "' not among columns of TrainData", "\n", sep = ""))
                    dummy.vars <- dummy.vars[dummy.vars != old.dummy.vars[i]]
                }
            }
            if(length(dummy.vars) == 0) {dummy.vars <- NULL}
        }
    }
# 

# modify RasterStack x only if this RasterStack was provided
    if (is.null(x) == F) {
        if (is.null(ext) == F) {
            if(length(x@title) == 0) {x@title <- "stack1"}
            title.old <- x@title
            x <- raster::crop(x, y=ext, snap="in")
            x@title <- title.old
        }
# same variables as TrainData in the rasterstack
        if (is.null(TrainData) == F) {
            vars <- names(TrainData)
            vars <- vars[which(vars!="pb")]
            x <- raster::subset(x, subset=vars, drop=FALSE)
        }
        if (is.null(TrainData) == T) {
            if (is.null(layer.drops) == F) {
                vars <- names(x)
                layer.drops <- as.character(layer.drops)
                factors <- as.character(factors)
                dummy.vars <- as.character(dummy.vars)
                nd <- length(layer.drops)
                for (i in 1:nd) {     
                    if (any(vars==layer.drops[i])==FALSE) {
                        cat(paste("\n", "WARNING: variable to exclude '", layer.drops[i], "' not among grid layers", "\n", sep = ""))
                    }else{
                        cat(paste("\n", "NOTE: variable '", layer.drops[i], "' will not be included as explanatory variable", "\n", sep = ""))
                        x <- raster::dropLayer(x, which(names(x) %in% c(layer.drops[i]) ))
                        vars <- names(x)
                        if (length(factors) > 0) {
                            factors <- factors[factors != layer.drops[i]]
                        }
                        if (length(dummy.vars) > 0) {
                            dummy.vars <- dummy.vars[dummy.vars != layer.drops[i]]
                        }
                    }
                }
                if(length(layer.drops) == 0) {layer.drops <- NULL}
                if(length(factors) == 0) {factors <- NULL}
                if(length(dummy.vars) == 0) {dummy.vars <- NULL}
            }
            if (is.null(factors) == F) {
                vars <- names(x)
                factors <- as.character(factors)
                nf <- length(factors)
                old.factors <- factors
                for (i in 1:nf) {
                    if (any(vars==old.factors[i])==FALSE) {
                        cat(paste("\n", "WARNING: categorical variable '", old.factors[i], "' not among grid layers", "\n", sep = ""))
                        factors <- factors[factors != old.factors[i]]
                    }
                }
                if(length(factors) == 0) {factors <- NULL}
            }
            if (is.null(dummy.vars) == F) {
                vars <- names(x)
                dummy.vars <- as.character(dummy.vars)
                nf <- length(dummy.vars)
                old.dummy.vars <- dummy.vars
                for (i in 1:nf) {
                    if (any(vars==old.dummy.vars[i]) == FALSE) {
                        cat(paste("\n", "WARNING: dummy variable '", old.dummy.vars[i], "' not among grid layers", "\n", sep = ""))
                        dummy.vars <- dummy.vars[dummy.vars != old.dummy.vars[i]]
                    }
                }
                if(length(dummy.vars) == 0) {dummy.vars <- NULL}
            }
        }
        # set minimum and maximum values
            for (i in 1:raster::nlayers(x)) {
                x[[i]] <- raster::setMinMax(x[[i]])
            }
        # declare factor layers
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                j <- which(names(x) == factors[i])
                x[[j]] <- raster::as.factor(x[[j]])
            }
        }
    }

# modify MAXENT.BackData if layer.drops
    if (is.null(MAXENT.BackData) == F) {
# same variables as TrainData
        if (is.null(TrainData) == F) {
            vars <- names(TrainData)
            vars <- vars[which(vars!="pb")]
            MAXENT.BackData <- MAXENT.BackData[, which(names(MAXENT.BackData) %in% vars)]
        }
        if (is.null(TrainData) == T) {
            if (is.null(layer.drops) == F) {
                vars <- names(MAXENT.BackData)
                layer.drops <- as.character(layer.drops)
                nd <- length(layer.drops)
                for (i in 1:nd) {     
                    if (any(vars==layer.drops[i])==FALSE) {
                        cat(paste("\n", "WARNING: variable to exclude '", layer.drops[i], "' not among columns of MAXENT.BackData", "\n", sep = ""))
                    }else{
                        MAXENT.BackData <- MAXENT.BackData[, which(names(MAXENT.BackData) != layer.drops[i])]
                    }
                }
                if(length(layer.drops) == 0) {layer.drops <- NULL}
            }
        }
    }
#
#
    if (is.null(input.weights) == F) {
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
        GEODIST <- max(c(input.weights["GEODIST"], -1), na.rm=T)
    }
    ws <- as.numeric(c(MAXENT, GBM, GBMSTEP, RF, GLM, GLMSTEP, GAM, GAMSTEP, MGCV, MGCVFIX, 
        EARTH, RPART, NNET, FDA, SVM, SVME, BIOCLIM, DOMAIN, MAHAL, GEODIST))
    names(ws) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX", 
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "GEODIST")
    ws <- ensemble.weights(weights=ws, exponent=1, best=0, min.weight=0)
#
    thresholds <- c(ws, NA)
    names(thresholds) <- c(names(ws), "ENSEMBLE")
#
#
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
# check formulae and packages
    if (ws["MAXENT"] > 0) {
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        if (!file.exists(jar)) {stop('maxent program is missing: ', jar, '\nPlease download it here: http://www.cs.princeton.edu/~schapire/maxent/')}
    }
    if (formulae.defaults == T) {
        if (is.null(TrainData) == T) {
            formulae <- ensemble.formulae(x, factors=factors, dummy.vars=dummy.vars)
        }else{
            formulae <- ensemble.formulae(TrainData, factors=factors, dummy.vars=dummy.vars)
        }
    }
    if (ws["GBM"] > 0) {
        if (! requireNamespace("gbm")) {stop("Please install the gbm package")}
	requireNamespace("splines")
        if (is.null(GBM.formula) == T && formulae.defaults == T) {GBM.formula <- formulae$GBM.formula}
        if (is.null(GBM.formula) == T) {stop("Please provide the GBM.formula (hint: use ensemble.formulae function)")}
        environment(GBM.formula) <- .BiodiversityR
    }
    if (ws["GBMSTEP"] > 0) {
#        if (! require(gbm)) {stop("Please install the gbm package")}
    }
    if (ws["RF"] > 0) {
#        if (! require(randomForest)) {stop("Please install the randomForest package")}
        if (is.null(RF.formula) == T && formulae.defaults == T) {RF.formula <- formulae$RF.formula}
        if (is.null(RF.formula) == T) {stop("Please provide the RF.formula (hint: use ensemble.formulae function)")}
        environment(RF.formula) <- .BiodiversityR
        if (identical(RF.ntree, trunc(RF.ntree/2)) == F) {RF.ntree <- RF.ntree + 1}
    }
    if (ws["GLM"] > 0) {
        if (is.null(GLM.formula) == T && formulae.defaults == T) {GLM.formula <- formulae$GLM.formula}
        if (is.null(GLM.formula) == T) {stop("Please provide the GLM.formula (hint: use ensemble.formulae function)")}
        environment(GLM.formula) <- .BiodiversityR
        assign("GLM.family", GLM.family, envir=.BiodiversityR)
    }
    if (ws["GLMSTEP"] > 0) {
#        if (! require(MASS)) {stop("Please install the MASS package")}
        if (is.null(STEP.formula) == T && formulae.defaults == T) {STEP.formula <- formulae$STEP.formula}
        if (is.null(GLMSTEP.scope) == T && formulae.defaults == T) {GLMSTEP.scope <- formulae$GLMSTEP.scope}
        if (is.null(STEP.formula) == T) {stop("Please provide the STEP.formula (hint: use ensemble.formulae function)")}
        if (is.null(GLMSTEP.scope) == T) {stop("Please provide the GLMSTEP.scope (hint: use ensemble.formulae function)")}
        environment(STEP.formula) <- .BiodiversityR
        assign("GLM.family", GLM.family, envir=.BiodiversityR)
    }
    if (ws["GAM"] > 0 || ws["GAMSTEP"] > 0) {
        cat(paste("\n\n"))
#        try(detach(package:mgcv), silent=T)
#        suppressMessages(require(gam))
#         if (! require("gam", quietly=T)) {stop("Please install the gam package")}
    }
    if (ws["GAM"] > 0) {
        if (is.null(GAM.formula) == T && formulae.defaults == T) {GAM.formula <- formulae$GAM.formula}
        if (is.null(GAM.formula) == T) {stop("Please provide the GAM.formula (hint: use ensemble.formulae function)")}
        environment(GAM.formula) <- .BiodiversityR
        assign("GAM.family", GAM.family, envir=.BiodiversityR)
    }
    if (ws["GAMSTEP"] > 0) {
        if (is.null(STEP.formula) == T && formulae.defaults == T) {STEP.formula <- formulae$STEP.formula}
        if (is.null(GAMSTEP.scope) == T && formulae.defaults == T) {GAMSTEP.scope <- formulae$GAMSTEP.scope}
        if (is.null(STEP.formula) == T) {stop("Please provide the STEP.formula (hint: use ensemble.formulae function)")}
        if (is.null(GAMSTEP.scope) == T) {stop("Please provide the GAMSTEP.scope (hint: use ensemble.formulae function)")}
        environment(STEP.formula) <- .BiodiversityR
        assign("GAM.family", GAM.family, envir=.BiodiversityR)
    }
    if (ws["MGCV"] > 0 || ws["MGCVFIX"] > 0) {
        cat(paste("\n\n"))
#        try(detach(package:gam), silent=T)
#        options(warn=-1)
#        if (! require(mgcv, quietly=T)) {stop("Please install the mgcv package")}
#         get the probabilities from MGCV
            predict.mgcv <- function(object, newdata, type="response") {
                p <- mgcv::predict.gam(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
#        options(warn=0)
    }
    if (ws["MGCV"] > 0) {
        if (is.null(MGCV.formula) == T && formulae.defaults == T) {MGCV.formula <- formulae$MGCV.formula}
        if (is.null(MGCV.formula) == T) {stop("Please provide the MGCV.formula (hint: use ensemble.formulae function)")}
        environment(MGCV.formula) <- .BiodiversityR
        assign("GAM.family", GAM.family, envir=.BiodiversityR)
    }
    if (ws["MGCVFIX"] > 0) {
        if (is.null(MGCVFIX.formula) == T && formulae.defaults == T) {MGCVFIX.formula <- formulae$MGCVFIX.formula}
        if (is.null(MGCVFIX.formula) == T) {stop("Please provide the MGCVFIX.formula (hint: use ensemble.formulae function)")}
        environment(MGCVFIX.formula) <- .BiodiversityR
        assign("GAM.family", GAM.family, envir=.BiodiversityR)
    }
    if (ws["EARTH"] > 0) {
#        if (! require(earth)) {stop("Please install the earth package")}
        if (is.null(EARTH.formula) == T && formulae.defaults == T) {EARTH.formula <- formulae$EARTH.formula}
        if (is.null(EARTH.formula) == T) {stop("Please provide the EARTH.formula (hint: use ensemble.formulae function)")}
        environment(EARTH.formula) <- .BiodiversityR
#         get the probabilities from earth
            predict.earth2 <- function(object, newdata, type="response") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
    }
    if (ws["RPART"] > 0) {
#        if (! require(rpart)) {stop("Please install the rpart package")}
        if (is.null(RPART.formula) == T && formulae.defaults == T) {RPART.formula <- formulae$RPART.formula}
        if (is.null(RPART.formula) == T) {stop("Please provide the RPART.formula (hint: use ensemble.formulae function)")}
        environment(RPART.formula) <- .BiodiversityR
    }
    if (ws["NNET"] > 0) {
#        if (! require(nnet)) {stop("Please install the nnet package")}
        if (is.null(NNET.formula) == T && formulae.defaults == T) {NNET.formula <- formulae$NNET.formula}
        if (is.null(NNET.formula) == T) {stop("Please provide the NNET.formula (hint: use ensemble.formulae function)")}
        environment(NNET.formula) <- .BiodiversityR
#         get the probabilities from nnet
            predict.nnet2 <- function(object, newdata, type="raw") {
                p <- predict(object=object, newdata=newdata, type=type)
                return(as.numeric(p))
            }
    }
    if (ws["FDA"] > 0) {
#        if (! require(mda)) {stop("Please install the mda package")}
        if (is.null(FDA.formula) == T && formulae.defaults == T) {FDA.formula <- formulae$FDA.formula}
        if (is.null(FDA.formula) == T) {stop("Please provide the FDA.formula (hint: use ensemble.formulae function)")}
        environment(FDA.formula) <- .BiodiversityR
    }
    if (ws["SVM"] > 0) {
#        if (! require(kernlab)) {stop("Please install the kernlab package")}
        if (is.null(SVM.formula) == T && formulae.defaults == T) {SVM.formula <- formulae$SVM.formula}
        if (is.null(SVM.formula) == T) {stop("Please provide the SVM.formula (hint: use ensemble.formulae function)")}
        environment(SVM.formula) <- .BiodiversityR
    }
    if (ws["SVME"] > 0) {
#        if (! require(e1071)) {stop("Please install the e1071 package")}
        if (is.null(SVME.formula) == T && formulae.defaults == T) {SVME.formula <- formulae$SVME.formula}
        if (is.null(SVME.formula) == T) {stop("Please provide the SVME.formula (hint: use ensemble.formulae function)")}
        environment(SVME.formula) <- .BiodiversityR
#         get the probabilities from svm
            predict.svme <- function(model, newdata) {
                p <- predict(model, newdata, probability=T)
                return(attr(p, "probabilities")[,1])
             }
    }
    if (ws["MAHAL"] > 0) {
#         get the probabilities from mahal
            predict.mahal <- function(model, newdata, MAHAL.shape) {
                p <- dismo::predict(object=model, x=newdata)
                p <- p - 1 - MAHAL.shape
                p <- abs(p)
                p <- MAHAL.shape / p
                return(p)
             }
    }

# create TrainData and TestData
    if (is.null(TrainData) == F) {
        if(any(is.na(TrainData))) {
            cat(paste("\n", "WARNING: sample units with missing data removed from calibration data","\n\n",sep = ""))
        }        
        TrainValid <- complete.cases(TrainData)
        TrainData <- TrainData[TrainValid,]
        if(is.null(TestData) == T) {
            TestData <- TrainData
            if (k > 1) {
                groupp <- dismo::kfold(TrainData, k=k, by=TrainData[,"pb"])
                TrainData.c <- TrainData[groupp != 1,]
                TestData <- TrainData[groupp == 1,]
                TrainData <- TrainData.c
            }
        } 
    }else{
        if (is.null(a)==T) {
            if (excludep == T) {
                a <- dismo::randomPoints(x[[1]], n=an, p=p, excludep=T)
            }else{
                a <- dismo::randomPoints(x[[1]], n=an, p=NULL, excludep=F)
            }        
        }
        if (is.null(pt)==T && is.null(TestData)) {pt <- p}
        if (k > 1 && identical(pt, p) == T) {
            groupp <- dismo::kfold(p, k=k)
            pc <- p[groupp != 1,]
            pt <- p[groupp == 1,]
            p <- pc
        }
        if (is.null(at)==T && is.null(TestData)) {at <- a}
        if (k > 1 && identical(at, a) == T) {
            groupa <- dismo::kfold(a, k=k)
            ac <- a[groupa != 1,]
            at <- a[groupa == 1,]
            a <- ac
        }
# check for spatial sorting bias (are the testing absences farther away than testing presences)
        if (is.null(p)==F && identical(pt, p)==F) {
            sb.bias <- dismo::ssb(p=pt, a=at, reference=p)
            sb.bias2 <- sb.bias[, 1]/sb.bias[, 2]
            cat(paste("\n", "Spatial sorting bias (dismo package, no bias=1, extreme bias=0): ",  sb.bias2, "\n", sep = ""))
        }
# attempt to reduce spatial bias by searching absence testing locations within circles around all known presences
        if (CIRCLES.at == T) {
            if (identical(pt, p) == T) {
                cat(paste("\n", "No search for testing absences in circular neighbourhoods since no separate testing presences", sep = ""))
                CIRCLES.at <- FALSE
            }else{
                cat(paste("\n", "Random selection of testing absences in circular neighbourhoods", sep = ""))
                pres_all <- rbind(pt, p)
                circles.calibrate <- dismo::circles(p=pres_all, lonlat=raster::isLonLat(x[[1]]), d=CIRCLES.d)
                circles.predicted <- predict(circles.calibrate, x[[1]])
                at <- dismo::randomPoints(circles.predicted, n=nrow(at), p=pres_all, excludep=T)
                sb.bias <- dismo::ssb(p=pt, a=at, reference=p)
                sb.bias2 <- sb.bias[, 1]/sb.bias[, 2]
                cat(paste("\n", "Spatial sorting bias with new testing absences (dismo package, no bias=1, extreme bias=0): ",  sb.bias2, "\n", sep = ""))
            }
        }
        TrainData <- dismo::prepareData(x, p, b=a, factors=factors, xy=FALSE)
        if(any(is.na(TrainData[TrainData[,"pb"]==1,]))) {
            cat(paste("\n", "WARNING: presence locations with missing data removed from calibration data","\n\n",sep = ""))
        }
        TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==1,])
        p <- p[TrainValid,]
        if(any(is.na(TrainData[TrainData[,"pb"]==0,]))) {
            cat(paste("\n", "WARNING: background locations with missing data removed from calibration data","\n\n",sep = ""))
        }
        TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==0,])
        a <- a[TrainValid,]
        TrainData <- dismo::prepareData(x, p, b=a, factors=factors, xy=FALSE)
    }
#
    if (is.null(TestData) == F) {
        TestData <- data.frame(TestData)
        if(any(is.na(TestData))) {
            cat(paste("\n", "WARNING: sample units with missing data removed from testing data","\n\n",sep = ""))
        }
        TestValid <- complete.cases(TestData)
        TestData <- TestData[TestValid,]
        if (all(colnames(TestData)!="pb") == T) {stop("one column needed of 'pb' with presence and absence for TestData")} 
    }else{
        TestData <- dismo::prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
        if(any(is.na(TestData[TestData[,"pb"]==1,]))) {
            cat(paste("\n", "WARNING: presence locations with missing data removed from evaluation data","\n\n",sep = ""))
        }
        TestValid <- complete.cases(TestData[TestData[,"pb"]==1,])
        pt <- pt[TestValid,]
        if(any(is.na(TestData[TestData[,"pb"]==0,]))) {
            cat(paste("\n", "WARNING: background locations with missing data removed from evaluation data","\n\n",sep = ""))
        }
        TestValid <- complete.cases(TestData[TestData[,"pb"]==0,])
        at <- at[TestValid,]
        TestData <- dismo::prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
    }
#
# check if TestData is different from TrainData
    no.tests <- FALSE
    if (identical(TrainData, TestData) == T) {no.tests <- TRUE}
#
# include all possible factor levels in TrainData (especially important if models are kept)
    if (is.null(factors)==F && is.null(x)==T) {
        for (i in 1:length(factors)) {
            if (identical(levels(droplevels(TrainData[,factors[i]])), levels(droplevels(TestData[,factors[i]])))==F) {
                cat(paste("\n", "WARNING: differences in factor levels between calibration and evaluation data (variable ", factors[i], ")", "\n", sep = ""))
                cat(paste("Same levels set for both data sets to avoid problems with some evaluations", "\n", sep = ""))
                cat(paste("However, some predictions may still fail", "\n", sep = ""))
                uniquelevels <- unique(c(levels(droplevels(TrainData[,factors[i]])), levels(droplevels(TestData[,factors[i]]))))
                levels(TrainData[,factors[i]]) <- uniquelevels
                levels(TestData[,factors[i]]) <- uniquelevels
            }
        }
    }
    if(is.null(factors)==F && is.null(x)==F) {
        if(models.keep==T) {
            categories <- as.list(factors)
            names(categories) <- factors
        }
        for (i in 1:length(factors)) {
            all.categories <- raster::freq(x[[which(names(x) == factors[i])]])[,1]
            all.categories <- all.categories[is.na(all.categories) == F]
            all.categories <- as.character(all.categories)
            if(models.keep==T) {
                categories[[as.name(factors[i])]] <- all.categories
            }
            train.categories <- levels(droplevels(TrainData[,factors[i]]))
            new.categories <- c(all.categories[is.na(match(all.categories, train.categories))])
            if (length(new.categories) > 0) {
                cat(paste("\n", "The following levels were initially not captured by TrainData for factor '", factors[i], "'\n", sep = ""))
                print(new.categories)
                if (is.null(x)==F && is.null(p)==F && is.null(a)==F) {
# step 1: search if suitable presence locations in TestData
                    for (j in 1:length(new.categories)) {
                        if (any(TestData[TestData[,"pb"]==1, factors[i]] == new.categories[j])) {    
                            cat(paste("Warning: level '", new.categories[j], "' available for presence location in Test Data", "\n", sep = ""))
                        }
                    }
# step 2: stratified background sample
                    strat1 <- raster::sampleStratified(x[[which(names(x) == factors[i])]], size=1, exp=1, na.rm=TRUE, xy=FALSE)
                    strat1 <- strat1[which(strat1[,2] %in% new.categories), 1]
                    xy1 <- raster::xyFromCell(x[[which(names(x) == factors[i])]], cell=strat1, spatial=FALSE)
                    a <- rbind(a, xy1)
                    TrainData <- dismo::prepareData(x, p, b=a, factors=factors, xy=FALSE)
                    TrainValid <- complete.cases(TrainData[TrainData[,"pb"]==0,])
                    a <- a[TrainValid,]
                    TrainData <- dismo::prepareData(x, p, b=a, factors=factors, xy=FALSE)
                    train.categories <- levels(droplevels(TrainData[,factors[i]]))
                    new.categories <- all.categories[is.na(match(all.categories, train.categories))]
                    if (length(new.categories) == 0) {
                        cat(paste("All levels have now been included as background data for TrainData for factor '", factors[i], "'\n", sep = ""))
                    }else{
                       cat(paste("The following levels were not captured by TrainData for factor '", factors[i], "'\n", sep = ""))
                       print(new.categories)
                       cat(paste("\n", "Attempt to include these levels was complicated by missing values in other layers", "\n", sep = ""))
                    }
                }
            }
# step 3: also modify test data, but only if no circular neighbourhood
            if (no.tests == F) {
                test.categories <- levels(droplevels(TestData[,factors[i]]))
                new.categories <- c(all.categories[is.na(match(all.categories, test.categories))])
                if (length(new.categories)>0  && CIRCLES.at==F) {
                    cat(paste("\n", "The following levels were initially not captured by TestData for factor '", factors[i], "'\n", sep = ""))
                    print(new.categories)
                    if (is.null(x)==F && is.null(pt)==F && is.null(at)==F) {
                        strat1 <- raster::sampleStratified(x[[which(names(x) == factors[i])]], size=1, exp=1, na.rm=TRUE, xy=FALSE)
                        strat1 <- strat1[which(strat1[,2] %in% new.categories), 1]
                        xy1 <- raster::xyFromCell(x[[which(names(x) == factors[i])]], cell=strat1, spatial=FALSE)
                        at <- rbind(at, xy1)
                        TestData <- dismo::prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
                        TestValid <- complete.cases(TestData[TestData[,"pb"]==0,])
                        at <- at[TestValid,]
                        TestData <- dismo::prepareData(x, p=pt, b=at, factors=factors, xy=FALSE)
                        test.categories <- levels(droplevels(TestData[,factors[i]]))
                        new.categories <- all.categories[is.na(match(all.categories, test.categories))]
                        if (length(new.categories) == 0) {
                            cat(paste("All levels have now been included as background data for TestData for factor '", factors[i], "'\n", sep = ""))
                        }else{
                            cat(paste("The following levels were not captured by TestData for factor '", factors[i], "'\n", sep = ""))
                            print(new.categories)
                            cat(paste("\n", "Attempt to include these levels was complicated by missing values in other layers", "\n", sep = ""))
                        }
                    }
                }
                if (length(new.categories)>0  && CIRCLES.at==T) {
                    cat(paste("\n", "Note that the following levels were not captured in the circular neighbourhood by TestData for factor '", factors[i], "'\n", sep = ""))
                    print(new.categories)
                }
            }
        }
    }
#
    if (sum(ws, na.rm=T) > 0) {
        cat(paste("\n", "Summary of Training data set used for calibrations (rows: ", nrow(TrainData),  ")\n", sep = ""))
        print(summary(TrainData))
        if (no.tests == F) {
            cat(paste("\n", "Summary of Testing data set used for evaluations (rows: ", nrow(TestData),  ")\n", sep = ""))
            print(summary(TestData))
        }else{
            cat(paste("\n", "(no tests with separate data set)", "\n", sep = ""))
        }
    }
#
    if(models.keep==T) {
        models <- list(MAXENT=NULL, GBM=NULL, GBMSTEP=NULL, RF=NULL, GLM=NULL, 
            GLMSTEP=NULL, GAM=NULL, GAMSTEP=NULL, MGCV=NULL, MGCVFIX=NULL, EARTH=NULL, RPART=NULL, 
            NNET=NULL, FDA=NULL, SVM=NULL, SVME=NULL, BIOCLIM=NULL, DOMAIN=NULL, MAHAL=NULL,
            formulae=NULL, TrainData=NULL, TestData=NULL, p=NULL, a=NULL, pt=NULL, at=NULL, MAXENT.BackData=NULL, 
            vars=NULL, factors=NULL, categories=NULL, dummy.vars=NULL, threshold.method=NULL, threshold.sensitivity=NULL, species.name=NULL)
        models$TrainData <- TrainData
        models$p <- p
        models$a <- a
        if (no.tests==F) {models$pt <- pt}
        if (no.tests==F) {models$at <- at}
        vars <- names(TrainData)
        vars <- vars[which(vars!="pb")]
        models$vars <- vars
        models$factors <- factors
        if(is.null(factors)==F) {models$categories <- categories}
        models$dummy.vars <- dummy.vars
        models$threshold.method <- threshold.method
        models$threshold.sensitivity <- threshold.sensitivity
        models$species.name <- species.name
        if (no.tests == F) {models$TestData <- TestData}
    }else{
        models <- NULL
    }
#
# Data frames for distance-based methods and SVME
    TrainData.vars <- TrainData[,colnames(TrainData) != "pb"]
    assign("TrainData.vars", TrainData.vars, envir=.BiodiversityR)
    TrainData.pres <- TrainData[TrainData[,"pb"]==1,]
    TrainData.pres <- TrainData.pres[,colnames(TrainData.pres) != "pb"]
    assign("TrainData.pres", TrainData.pres, envir=.BiodiversityR)
    var.names <- names(TrainData.pres)
    TestData.vars <- TestData[,colnames(TestData) != "pb"]
    assign("TestData.vars", TestData.vars, envir=.BiodiversityR)
#
# separate data set to calibrate MAXENT
# special case to use ensemble.test simply to provide data sets
    MAXENT2 <- 0
    if (sum(ws, na.rm=T) == 0) {MAXENT2 <- 1}
#
    if (ws["MAXENT"] > 0 || MAXENT2 > 0) {
        if (is.null(MAXENT.BackData) == T) {
            if (is.null(x) == T) {
                cat(paste("\n", "WARNING: not possible to create MAXENT.BackData as RasterStack x is missing", sep = "")) 
                cat(paste("\n", "MAXENT model will not be calibrated", "\n", sep = "")) 
                ws["MAXENT"] <- 0
            }else{
# default option of MAXENT is to exclude presence locations
                if (is.null(MAXENT.a) == T) {
                    MAXENT.a <- dismo::randomPoints(x[[1]], n=MAXENT.an, p=p, excludep=T)
                }
                MAXENT.BackData <- raster::extract(x, y=MAXENT.a)
            }
        }
        MAXENT.BackData <- data.frame(MAXENT.BackData)
        if (is.null(layer.drops) == F) {
            layer.drops <- as.character(layer.drops)
            nd <- length(layer.drops)
            for (i in 1:nd) {     
                MAXENT.BackData <- MAXENT.BackData[, which(colnames(MAXENT.BackData) != layer.drops[i])]
            }
            if(length(layer.drops) == 0) {layer.drops <- NULL}
        }
        TestValid <- complete.cases(MAXENT.BackData)
        MAXENT.BackData <- MAXENT.BackData[TestValid,]
        if(is.null(factors)==F) {
            for (i in 1:length(factors)) {
                j <- which(names(MAXENT.BackData) == factors[i])
                MAXENT.BackData[,j] <- as.factor(MAXENT.BackData[,j])
            }
        }
        if (sum(ws, na.rm=T) > 0) {
            cat(paste("\n", "Summary of Background data set used for calibration of MAXENT model (rows: ", nrow(MAXENT.BackData),  ")\n", sep = ""))
            print(summary(MAXENT.BackData))
        }
        if (all(names(MAXENT.BackData) %in% names(TrainData.vars)) == F) {
                cat(paste("\n", "WARNING: MAXENT.BackData has different (sequence of) colnames than TrainData", sep = ""))
                cat(paste("\n", "MAXENT model will not be calibrated", "\n", sep = "")) 
                ws["MAXENT"] <- 0
#                MAXENT.BackData <- NULL
        }else{
            MAXENT.TrainData <- rbind(TrainData.pres, MAXENT.BackData)
            cat(paste("\n", "Summary of Training data set used for calibration of MAXENT model (rows: ", nrow(MAXENT.TrainData),  ", presence locations: ", nrow(TrainData.pres), ")\n", sep = ""))
            print(summary(MAXENT.TrainData))
            MAXENT.pa <- c(rep(1, nrow(TrainData.pres)), rep(0, nrow(MAXENT.BackData)))
            assign("MAXENT.TrainData", MAXENT.TrainData, envir=.BiodiversityR)
            assign("MAXENT.pa", MAXENT.pa, envir=.BiodiversityR)            
        }
    }
#
#
    newVIF <- NULL
    if (VIF == T  && length(names(TrainData.vars)) > 1 ) {
#        if (! require(car)) {stop("Please install the car package")}
# only use background data
        TrainDataNum <- TrainData[TrainData[,"pb"]==0,]
        LM.formula <- ensemble.formulae(TrainData, factors=factors)$RF.formula
# create possible response
        TrainDataNum[,"pb"] <- mean(as.numeric(TrainDataNum[,2]))
        vifresult <- NULL
        tryCatch(vifresult <- car::vif(lm(formula=LM.formula, data=TrainDataNum)),
            error= function(err) {print(paste("VIF evaluation (package: car) failed"))},
                    silent=F)
        if (is.null(vifresult) == F) {
            cat(paste("\n", "Variance inflation (package: car)", "\n", sep = ""))        
            print(vifresult)
        }
        cat(paste("\n", "VIF directly calculated from linear model with focal numeric variable as response", "\n", sep = ""))
        TrainDataNum <- TrainDataNum[,names(TrainDataNum)!="pb"]
        varnames <- names(TrainDataNum)
        newVIF <- as.numeric(varnames)
        names(newVIF) <- varnames
        for (i in 1:length(varnames)) {
            response.name <- varnames[i]
            explan.names <- varnames[-i]
            if ((response.name  %in% factors) == F) {
                LM.formula <- as.formula(paste(response.name, "~", paste(explan.names, collapse="+"), sep=""))
                newVIF[i] <- summary(lm(formula=LM.formula, data=TrainDataNum))$r.squared
            }
        }
        newVIF <- 1/(1-newVIF)
        newVIF <- sort(newVIF, decreasing=T, na.last=T)
        print(newVIF)
    }
    if (COR == T) {
        TrainDataNum <- TrainData[, colnames(TrainData) != "pb"]
        if(is.null(factors)) {
            for (i in 1:length(factors)) {
                TrainDataNum <- TrainDataNum[, colnames(TrainDataNum) != factors[i]]            
            }
        }
        corresult <- cor(TrainDataNum)
        corresult <- round(100*corresult, digits=2)
        cat(paste("\n", "Correlation between numeric variables (as percentage)", "\n", sep = ""))        
        print(corresult)
    }
#
    modelresults <- data.frame(array(dim=c(nrow(TrainData), 20), 0))
    colnames(modelresults) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "ENSEMBLE")
    TrainData <- cbind(TrainData, modelresults)
    assign("TrainData", TrainData, envir=.BiodiversityR)
    modelresults <- data.frame(array(dim=c(nrow(TestData), 20), 0))
    colnames(modelresults) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "ENSEMBLE")
    TestData <- cbind(TestData, modelresults)
    assign("TestData", TestData, envir=.BiodiversityR)
    weights <- as.numeric(array(dim=19, 0))
    names(weights) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX", 
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL")
#
    if(evaluations.keep==T) {
        evaluations <- list(MAXENT.C=NULL, MAXENT.T=NULL, 
            GBM.trees=NULL, GBM.C=NULL, GBM.T=NULL, GBMSTEP.trees=NULL, GBMSTEP.C=NULL, GBMSTEP.T=NULL, 
            RF.C=NULL, RF.T=NULL, GLM.C=NULL, GLM.T=NULL, GLMS.C=NULL, GLMS.T=NULL, 
            GAM.C=NULL, GAM.T=NULL, GAMS.C=NULL, GAMS.T=NULL, MGCV.C=NULL, MGCV.T=NULL, MGCVF.C=NULL, MGCVF.T=NULL,
            EARTH.C=NULL, EARTH.T=NULL, RPART.C=NULL, RPART.T=NULL,
            NNET.C=NULL, NNET.T=NULL, FDA.C=NULL, FDA.T=NULL, SVM.C=NULL, SVM.T=NULL, SVME.C=NULL, SVME.T=NULL,
            BIOCLIM.C=NULL, BIOCLIM.T=NULL, DOMAIN.C=NULL, DOMAIN.T=NULL, MAHAL.C=NULL, MAHAL.T=NULL,
            GEODIST.C=NULL, GEODIST.T=NULL, ENSEMBLE.C=NULL, ENSEMBLE.T=NULL, STRATEGY.weights=NULL,
            TrainData=NULL, TestData=NULL, MAXENT.BackData=NULL, 
            factors=NULL, dummy.vars=NULL)
        evaluations$factors <- factors
        evaluations$dummy.vars <- dummy.vars
    }else{
        evaluations <- NULL
    }
#
    Yweights1 <- Yweights
    if (Yweights == "BIOMOD") {
        #have equal weight of presence vs. background
        Yweights1 <- numeric(length = nrow(TrainData))
        pres <- length(TrainData[, 1] == 1)
        abs <- length(TrainData[, 1] == 0)
        Yweights1[which(TrainData[, 1] == 1)] <- 1
        Yweights1[which(TrainData[, 1] == 0)] <- pres/abs
    }
    if (Yweights == "equal") {
        Yweights1 <- numeric(length = nrow(TrainData))
        Yweights1[] <- 1
    }
    Yweights <- Yweights1
    assign("Yweights", Yweights, envir=.BiodiversityR)
#
# prepare for calculation of deviance
    obs1 <- TrainData[, "pb"]
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
# count models
    mc <- 0
#
# Different modelling algorithms
#
    if (ws["MAXENT"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Maximum entropy algorithm (package: dismo)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        # Put the file 'maxent.jar' in the 'java' folder of dismo
        # the file 'maxent.jar' can be obtained from from http://www.cs.princeton.edu/~schapire/maxent/.
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        if(is.null(MAXENT.OLD) == T) {
            tryCatch(results <- dismo::maxent(x=MAXENT.TrainData, p=MAXENT.pa, factors=factors, path=MAXENT.path),
                error= function(err) {print(paste("MAXENT calibration failed"))},
                silent=T)
        }else{ 
            results <- MAXENT.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "MAXENT calibration tested with calibration data of other algorithms","\n\n", sep = ""))
            TrainData[,"MAXENT"] <- dismo::predict(object=results, x=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(MAXENT.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ MAXENT"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- MAXENT.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n", sep = ""))
                TrainData[,"MAXENT.step1"] <- TrainData[,"MAXENT"]
                TrainData[,"MAXENT"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "MAXENT"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MAXENT"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MAXENT"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
#            thresholds["MAXENT"] <- threshold(eval1, sensitivity=threshold.sensitivity)[[threshold.method]]
            thresholds["MAXENT"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["MAXENT"]))
            weights["MAXENT"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n", sep = ""))
                TestData[,"MAXENT"] <- dismo::predict(object=results, x=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"MAXENT.step1"] <- TestData[,"MAXENT"]
                    TestData[,"MAXENT"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"MAXENT"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MAXENT"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("MAXENT evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MAXENT"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="MAXENT", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: MAXENT evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MAXENT.C <- eval1
                evaluations$MAXENT.T <- eval2
            }
            if (models.keep==T) {
                models$MAXENT <- results
                models$MAXENT.PROBIT <- results2
            }
        }else{ cat(paste("\n", "WARNING: MAXENT calibration failed", "\n", "\n"))}
    }
    if (ws["GBM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized boosted regression modeling (package: gbm) \n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GBM.OLD) == T) {       
            tryCatch(results <- gbm::gbm(formula=GBM.formula, data=TrainData, weights=Yweights, distribution="bernoulli", 
                    interaction.depth=7, shrinkage=0.001, bag.fraction=0.5, train.fraction=1, 
                    n.trees=GBM.n.trees, verbose=F, cv.folds=5),
                error= function(err) {print(paste("GBM calibration failed"))},
                silent=T)
        }else{ 
            results <- GBM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"GBM"] <- gbm::predict.gbm(object=results, newdata=TrainData.vars, n.trees=results$n.trees, type="response")
            if (PROBIT == T) {
                if(is.null(GBM.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ GBM"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- GBM.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"GBM.step1"] <- TrainData[,"GBM"]
                TrainData[,"GBM"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "GBM"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GBM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GBM"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GBM"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["GBM"]))
            weights["GBM"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"GBM"] <- gbm::predict.gbm(object=results, newdata=TestData.vars, n.trees=results$n.trees, type="response")
                if (PROBIT == T) {
                    TestData[,"GBM.step1"] <- TestData[,"GBM"]
                    TestData[,"GBM"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"GBM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GBM"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GBM evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GBM"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="BRT (gbm)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: GBM evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GBM.trees <- results$n.trees 
                evaluations$GBM.C <- eval1
                evaluations$GBM.T <- eval2
            }
            if (models.keep==T) {
                models$GBM <- results
                models$GBM.PROBIT <- results2
                models$formulae$GBM.formula <- GBM.formula
            }
        }else{ cat(paste("\n", "WARNING: GBM calibration failed", "\n", "\n"))}
    }
    if (ws["GBMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". gbm step algorithm (package: dismo)\n", sep=""))
#        require(gbm, quietly=T)        
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GBMSTEP.OLD) == T) {       
            tryCatch(results <- dismo::gbm.step(data=TrainData, gbm.y=1, gbm.x=GBMSTEP.gbm.x, family="bernoulli",
                    site.weights=Yweights, tree.complexity=GBMSTEP.tree.complexity, learning.rate = GBMSTEP.learning.rate, 
                    bag.fraction=GBMSTEP.bag.fraction, step.size=GBMSTEP.step.size, verbose=F, silent=T, plot.main=F),
                error= function(err) {print(paste("stepwise GBM calibration failed"))},
                silent=T)
        }else{ 
            results <- GBMSTEP.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "stepwise GBM trees (target > 1000)", "\n", sep = ""))
            print(results$n.trees)
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"GBMSTEP"] <- gbm::predict.gbm(object=results, newdata=TrainData.vars, n.trees=results$n.trees, type="response")
            if (PROBIT == T) {
                if(is.null(GBMSTEP.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ GBMSTEP"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- GBMSTEP.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"GBMSTEP.step1"] <- TrainData[,"GBMSTEP"]
                TrainData[,"GBMSTEP"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "GBMSTEP"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GBMSTEP"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GBMSTEP"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GBMSTEP"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["GBMSTEP"]))
            weights["GBMSTEP"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"GBMSTEP"] <- gbm::predict.gbm(object=results, newdata=TestData.vars, n.trees=results$n.trees, type="response")
                if (PROBIT == T) {
                    TestData[,"GBMSTEP.step1"] <- TestData[,"GBMSTEP"]
                    TestData[,"GBMSTEP"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"GBMSTEP"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GBMSTEP"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("stepwise GBM evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GBMSTEP"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="BRT (gbm.step)", cex=2, adj=0, col.main="blue")
                   }
                }else{
                    cat(paste("\n", "WARNING: stepwise GBM evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GBMSTEP.trees <- results$n.trees 
                evaluations$GBMSTEP.C <- eval1
                evaluations$GBMSTEP.T <- eval2
            }
            if (models.keep==T) {
                models$GBMSTEP <- results
                models$GBMSTEP.PROBIT <- results2
            }
        }else{ cat(paste("\n", "WARNING: stepwise GBM calibration failed", "\n", "\n"))}
    }
    if (ws["RF"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Random forest algorithm (package: randomForest)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(RF.OLD) == T) {        
            tryCatch(results <- randomForest::randomForest(formula=RF.formula, ntree=RF.ntree, mtry=RF.mtry, data=TrainData, na.action=na.omit),
                error= function(err) {print(paste("random forest calibration failed"))},
                silent=T)
        }else{ 
            results <- RF.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"RF"] <- as.numeric(predict(object=results, newdata=TrainData.vars, type="response"))
            if (PROBIT == T) {
                if(is.null(RF.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ RF"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- RF.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"RF.step1"] <- TrainData[,"RF"]
                TrainData[,"RF"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "RF"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"RF"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"RF"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["RF"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["RF"]))
            weights["RF"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"RF"] <- as.numeric(predict(object=results, newdata=TestData.vars, type="response"))
                if (PROBIT == T) {
                    TestData[,"RF.step1"] <- TestData[,"RF"]
                    TestData[,"RF"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"RF"]
                TestAbs <- TestData[TestData[,"pb"]==0,"RF"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("RF evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["RF"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="RF", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: RF evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$RF.C <- eval1
                evaluations$RF.T <- eval2
            }
            if (models.keep==T) {
                models$RF <- results
                models$RF.PROBIT <- results2
                models$formulae$RF.formula <- RF.formula
            }
        }else{ cat(paste("\n", "WARNING: random forest calibration failed", "\n", "\n"))}
    }
    if (ws["GLM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Linear Model \n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GLM.OLD) == T) { 
            tryCatch(results <- glm(formula=GLM.formula, family=GLM.family, data=TrainData, weights=Yweights, control=glm.control(maxit=maxit)),
                error= function(err) {print(paste("GLM calibration failed"))},
                silent=T)
        }else{ 
            results <- GLM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"GLM"] <- predict.glm(object=results, newdata=TrainData.vars, type="response")
            if (PROBIT == T) {
                if(is.null(GLM.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ GLM"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- GLM.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"GLM.step1"] <- TrainData[,"GLM"]
                TrainData[,"GLM"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "GLM"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GLM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GLM"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GLM"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["GLM"]))
            weights["GLM"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"GLM"] <- predict.glm(object=results, newdata=TestData.vars, type="response")
                if (PROBIT == T) {
                    TestData[,"GLM.step1"] <- TestData[,"GLM"]
                    TestData[,"GLM"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"GLM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GLM"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GLM evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GLM"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="GLM", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: GLM evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GLM.C <- eval1
                evaluations$GLM.T <- eval2
            }
            if (models.keep==T) {
                models$GLM <- results
                models$GLM.PROBIT <- results2
                models$formulae$GLM.formula <- GLM.formula
            }
        }else{ cat(paste("\n", "WARNING: GLM calibration failed", "\n", "\n"))}
    }
    if (ws["GLMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Stepwise Generalized Linear Model \n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GLMSTEP.OLD) == T) {
            tryCatch(results <- glm(formula=STEP.formula, family=GLM.family, data=TrainData, weights=Yweights, control=glm.control(maxit=maxit)),
                error= function(err) {print(paste("first step of stepwise GLM calibration failed"))},
                silent=T)
            tryCatch(results2 <- MASS::stepAIC(results, scope=GLMSTEP.scope, direction="both", trace=F, steps=GLMSTEP.steps, k=GLMSTEP.k),
                error= function(err) {print(paste("stepwise GLM calibration failed"))},
                silent=T)
        }else{ 
            results2 <- GLMSTEP.OLD
        }
        if (is.null(results2) == F) {
            results <- results2
            results2 <- NULL
            cat(paste("\n", "stepwise GLM formula","\n\n",sep = ""))
            print(formula(results))
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"GLMSTEP"] <- predict.glm(object=results, newdata=TrainData.vars, type="response")
            if (PROBIT == T) {
                if(is.null(GLMSTEP.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ GLMSTEP"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- GLMSTEP.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"GLMSTEP.step1"] <- TrainData[,"GLMSTEP"]
                TrainData[,"GLMSTEP"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "GLMSTEP"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GLMSTEP"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GLMSTEP"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GLMSTEP"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["GLMSTEP"]))
            weights["GLMSTEP"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"GLMSTEP"] <- predict.glm(object=results, newdata=TestData.vars, type="response")
                if (PROBIT == T) {
                    TestData[,"GLMSTEP.step1"] <- TestData[,"GLMSTEP"]
                    TestData[,"GLMSTEP"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"GLMSTEP"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GLMSTEP"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("stepwise GLM evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GLMSTEP"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="STEPWISE GLM", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: stepwise GLM evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GLMS.C <- eval1
                evaluations$GLMS.T <- eval2
            }
            if (models.keep==T) {
                models$GLMSTEP <- results
                models$GLMSTEP.PROBIT <- results2
                models$formulae$STEP.formula <- STEP.formula
                models$formulae$GLMSTEP.scope <- GLMSTEP.scope
            }
        }else{ cat(paste("\n", "WARNING: stepwise GBM calibration failed", "\n", "\n"))}
    }
    if (ws["GAM"] > 0 || ws["GAMSTEP"] > 0) {
        cat(paste("\n\n"))
#        try(detach(package:mgcv), silent=T)
#        suppressMessages(require(gam))
#        require(gam, quietly=T)
    }
    if (ws["GAM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Additive Model (package: gam)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GAM.OLD) == T) {
            tryCatch(results <- gam::gam(formula=GAM.formula, family=GAM.family, data=TrainData, weights=Yweights, control=gam::gam.control(maxit=maxit, bf.maxit=50)),
                error= function(err) {print(paste("GAM calibration (gam package) failed"))},
                silent=F)
        }else{ 
            results <- GAM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"GAM"] <- gam::predict.gam(object=results, newdata=TrainData.vars, type="response")
            if (PROBIT == T) {
                if(is.null(GAM.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ GAM"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- GAM.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"GAM.step1"] <- TrainData[,"GAM"]
                TrainData[,"GAM"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "GAM"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GAM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GAM"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GAM"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["GAM"]))
            weights["GAM"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"GAM"] <- gam::predict.gam(object=results, newdata=TestData.vars, type="response")
                if (PROBIT == T) {
                    TestData[,"GAM.step1"] <- TestData[,"GAM"]
                    TestData[,"GAM"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"GAM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GAM"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GAM evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GAM"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="GAM (gam)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: GAM evaluation (gam package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GAM.C <- eval1
                evaluations$GAM.T <- eval2
            }
            if (models.keep==T) {
                models$GAM <- results
                models$GAM.PROBIT <- results2
                models$formulae$GAM.formula <- GAM.formula
            }
        }else{ cat(paste("\n", "WARNING: GAM calibration (gam package) failed", "\n", "\n"))}
    }
    if (ws["GAMSTEP"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Stepwise Generalized Additive Model (package: gam)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GAMSTEP.OLD) == T) {
            tryCatch(results <- gam::gam(formula=STEP.formula, family=GAM.family, data=TrainData, weights=Yweights, control=gam::gam.control(maxit=maxit, bf.maxit=50)), 
                error= function(err) {print(paste("first step of stepwise GAM calibration (gam package) failed"))},
                silent=T)
            assign("TrainData", TrainData, pos=GAMSTEP.pos)
            assign("GAM.family", GAM.family, pos=GAMSTEP.pos)
            assign("maxit", maxit, pos=GAMSTEP.pos)   
            tryCatch(results2 <- gam::step.gam(results, scope=GAMSTEP.scope, direction="both", trace=F, steps=GAMSTEP.steps), 
                error= function(err) {print(paste("stepwise GAM calibration (gam package) failed"))},
                silent=T)
            remove(TrainData, pos=GAMSTEP.pos)
            remove(GAM.family, pos=GAMSTEP.pos)
            remove(maxit, pos=GAMSTEP.pos)
        }else{ 
            results2 <- GAMSTEP.OLD
        }
        if (is.null(results2) == F) {
            results <- results2
            results2 <- NULL
            cat(paste("\n", "stepwise GAM formula (gam package)","\n\n",sep = ""))
            print(formula(results))
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"GAMSTEP"] <- gam::predict.gam(object=results, newdata=TrainData.vars, type="response")
            if (PROBIT == T) {
                if(is.null(GAMSTEP.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ GAMSTEP"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- GAMSTEP.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"GAMSTEP.step1"] <- TrainData[,"GAMSTEP"]
                TrainData[,"GAMSTEP"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "GAMSTEP"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"GAMSTEP"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"GAMSTEP"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["GAMSTEP"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["GAMSTEP"]))
            weights["GAMSTEP"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n","\n", sep = ""))
                TestData[,"GAMSTEP"] <- gam::predict.gam(object=results, newdata=TestData.vars, type="response")
                if (PROBIT == T) {
                    TestData[,"GAMSTEP.step1"] <- TestData[,"GAMSTEP"]
                    TestData[,"GAMSTEP"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"GAMSTEP"]
                TestAbs <- TestData[TestData[,"pb"]==0,"GAMSTEP"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("stepwise GAM evaluation (gam package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["GAMSTEP"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="STEPWISE GAM (gam)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: stepwise GAM evaluation (gam package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GAMS.C <- eval1
                evaluations$GAMS.T <- eval2
            }
            if (models.keep==T) {
                models$GAMSTEP <- results
                models$GAMSTEP.PROBIT <- results2
                models$formulae$STEP.formula <- STEP.formula
                models$formulae$GAMSTEP.scope <- GAMSTEP.scope
            }
        }else{ cat(paste("\n", "WARNING: stepwise GAM calibration (gam package) failed", "\n", "\n"))}
    }
    if (ws["MGCV"] > 0 || ws["MGCVFIX"] > 0) {
        cat(paste("\n\n"))
#        try(detach(package:gam), silent=T)
#        options(warn=-1)
#        require(mgcv, quietly=T)
#        options(warn=0)
    }
    if (ws["MGCV"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Generalized Additive Model (package: mgcv)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(MGCV.OLD) == T) {
            tryCatch(results <- mgcv::gam(formula=MGCV.formula, family=GAM.family, data=TrainData, weights=Yweights, 
                        select=MGCV.select, control=mgcv::gam.control(maxit=maxit)),
                error= function(err) {print(paste("GAM calibration (mgcv package) failed"))},
                silent=T)
        }else{ 
            results <- MGCV.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"MGCV"] <- predict.mgcv(object=results, newdata=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(MGCV.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ MGCV"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- MGCV.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"MGCV.step1"] <- TrainData[,"MGCV"]
                TrainData[,"MGCV"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "MGCV"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MGCV"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MGCV"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["MGCV"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["MGCV"]))
            weights["MGCV"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"MGCV"] <- predict.mgcv(object=results, newdata=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"MGCV.step1"] <- TestData[,"MGCV"]
                    TestData[,"MGCV"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"MGCV"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MGCV"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GAM evaluation (mgcv package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MGCV"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="GAM (mgcv)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: GAM evaluation (mgcv package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MGCV.C <- eval1
                evaluations$MGCV.T <- eval2
            }
            if (models.keep==T) {
                models$MGCV <- results
                models$MGCV.PROBIT <- results2
                models$formulae$MGCV.formula <- MGCV.formula
            }
        }else{ cat(paste("\n", "WARNING: GAM calibration (mgcv package) failed", "\n", "\n"))}
    }
    if (ws["MGCVFIX"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". GAM with fixed d.f. regression splines (package: mgcv)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(MGCVFIX.OLD) == T) {
            tryCatch(results <- mgcv::gam(formula=MGCVFIX.formula, family=GAM.family, data=TrainData, weights=Yweights, select=FALSE, control=mgcv::gam.control(maxit=maxit)),
                error= function(err) {print(paste("GAM calibration with fixed d.f. regression splines (mgcv package) failed"))},
                silent=T)
        }else{ 
            results <- MGCVFIX.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"MGCVFIX"] <- predict.mgcv(object=results, newdata=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(MGCVFIX.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ MGCVFIX"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- MGCVFIX.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"MGCVFIX.step1"] <- TrainData[,"MGCVFIX"]
                TrainData[,"MGCVFIX"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "MGCVFIX"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MGCVFIX"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MGCVFIX"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["MGCVFIX"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["MGCVFIX"]))
            weights["MGCVFIX"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"MGCVFIX"] <- predict.mgcv(object=results, newdata=TestData)
                if (PROBIT == T) {
                    TestData[,"MGCVFIX.step1"] <- TestData[,"MGCVFIX"]
                    TestData[,"MGCVFIX"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"MGCVFIX"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MGCVFIX"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("GAMFIX evaluation (mgcv package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MGCVFIX"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="GAM (mgcv)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: GAM with fixed d.f. regression splines evaluation (mgcv package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MGCVF.C <- eval1
                evaluations$MGCVF.T <- eval2
            }
            if (models.keep==T) {
                models$MGCVFIX <- results
                models$MGCVFIX.PROBIT <- results2
                models$formulae$MGCVFIX.formula <- MGCVFIX.formula
            }
        }else{ cat(paste("\n", "WARNING: MGCVFIX calibration (mgcv package) failed", "\n", "\n"))}
    }
    if (ws["EARTH"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Multivariate Adaptive Regression Splines (package: earth)\n", sep=""))
        if (!is.null(factors)) {
            cat(paste("\n", "NOTE: MARS evaluation (earth package) with factors probably requires dummy variables", "\n", sep=""))
        } 
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(EARTH.OLD) == T) {
            tryCatch(results <- earth::earth(formula=EARTH.formula, glm=EARTH.glm, data=TrainData, degree=2),
                error= function(err) {print(paste("MARS calibration (earth package) failed"))},
                silent=T)
        }else{ 
            results <- EARTH.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"EARTH"] <- predict.earth2(object=results, newdata=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(EARTH.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ EARTH"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- EARTH.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"EARTH.step1"] <- TrainData[,"EARTH"]
                TrainData[,"EARTH"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "EARTH"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"EARTH"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"EARTH"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["EARTH"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["EARTH"]))
            weights["EARTH"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"EARTH"] <- predict.earth2(object=results, newdata=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"EARTH.step1"] <- TestData[,"EARTH"]
                    TestData[,"EARTH"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"EARTH"]
                TestAbs <- TestData[TestData[,"pb"]==0,"EARTH"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("MARS evaluation (earth package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["EARTH"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="MARS (earth)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: MARS evaluation (earth package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$EARTH.C <- eval1
                evaluations$EARTH.T <- eval2
            }
            if (models.keep==T) {
                models$EARTH <- results
                models$EARTH.PROBIT <- results2
                models$formulae$EARTH.formula <- EARTH.formula
            }
        }else{ cat(paste("\n", "WARNING: MARS calibration (earth package) failed", "\n", "\n"))}
    }
    if (ws["RPART"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Recursive Partitioning And Regression Trees (package: rpart)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(RPART.OLD) == T) {
            tryCatch(results <- rpart::rpart(formula=RPART.formula, data=TrainData, weights=Yweights,
                    control=rpart::rpart.control(xval=RPART.xval, minbucket=5, minsplit=5, cp=0.001, maxdepth=25)),
                error= function(err) {print(paste("RPART calibration failed"))},
                silent=T)
        }else{ 
            results <- RPART.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"RPART"] <- predict(object=results, newdata=TrainData.vars, type="prob")[,2]
            if (PROBIT == T) {
                if(is.null(RPART.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ RPART"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- RPART.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"RPART.step1"] <- TrainData[,"RPART"]
                TrainData[,"RPART"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "RPART"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"RPART"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"RPART"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["RPART"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["RPART"]))
            weights["RPART"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"RPART"] <- predict(object=results, newdata=TestData.vars, type="prob")[,2]
                if (PROBIT == T) {
                    TestData[,"RPART.step1"] <- TestData[,"RPART"]
                    TestData[,"RPART"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"RPART"]
                TestAbs <- TestData[TestData[,"pb"]==0,"RPART"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("RPART evaluation failed"))},
                    silent=F)
                if (is.null(TestPres) == F && is.null(TestAbs) == F) {
                    eval2 <-  dismo::evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                    weights["RPART"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="RPART", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: RPART evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$RPART.C <- eval1
                evaluations$RPART.T <- eval2
            }
            if (models.keep==T) {
                models$RPART <- results
                models$RPART.PROBIT <- results2
                models$formulae$RPART.formula <- RPART.formula
            }
        }else{ cat(paste("\n", "WARNING: RPART calibration failed", "\n", "\n"))}
    }
    if (ws["NNET"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Artificial Neural Network (package: nnet)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(NNET.OLD) == T) {
            tryCatch(results <- nnet::nnet(formula=NNET.formula, size=NNET.size, decay=NNET.decay, data=TrainData, weights=Yweights, 
                    rang=0.1, maxit=maxit, trace=F),
                error= function(err) {print(paste("ANN calibration (nnet package) failed"))},
                silent=T)
        }else{ 
            results <- NNET.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"NNET"] <- predict.nnet2(object=results, newdata=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(NNET.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ NNET"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- NNET.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"NNET.step1"] <- TrainData[,"NNET"]
                TrainData[,"NNET"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "NNET"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"NNET"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"NNET"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["NNET"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["NNET"]))
            weights["NNET"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"NNET"] <- predict.nnet2(object=results, newdata=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"NNET.step1"] <- TestData[,"NNET"]
                    TestData[,"NNET"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"NNET"]
                TestAbs <- TestData[TestData[,"pb"]==0,"NNET"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("ANN evaluation (nnet package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["NNET"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="NNET", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: ANN evaluation (nnet package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$NNET.C <- eval1
                evaluations$NNET.T <- eval2
            }
            if (models.keep==T) {
                models$NNET <- results
                models$NNET.PROBIT <- results2
                models$formulae$NNET.formula <- NNET.formula
            }
        }else{ cat(paste("\n", "WARNING: ANN calibration (nnet package) failed", "\n", "\n"))}
    }
    if (ws["FDA"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Flexible Discriminant Analysis (package: mda)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(FDA.OLD) == T) {
            tryCatch(results <- mda::fda(formula=FDA.formula, method=mda::mars, data=TrainData, weights=Yweights),
                error= function(err) {print(paste("FDA calibration failed"))},
                silent=T)
        }else{ 
            results <- FDA.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"FDA"] <- predict(object=results, newdata=TrainData.vars, type="posterior")[,2]
            if (PROBIT == T) {
                if(is.null(FDA.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ FDA"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- FDA.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"FDA.step1"] <- TrainData[,"FDA"]
                TrainData[,"FDA"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "FDA"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"FDA"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"FDA"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["FDA"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["FDA"]))
            weights["FDA"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"FDA"] <- predict(object=results, newdata=TestData.vars, type="posterior")[,2]
                if (PROBIT == T) {
                    TestData[,"FDA.step1"] <- TestData[,"FDA"]
                    TestData[,"FDA"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"FDA"]
                TestAbs <- TestData[TestData[,"pb"]==0,"FDA"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("FDA evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["FDA"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="FDA", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: FDA evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$FDA.C <- eval1
                evaluations$FDA.T <- eval2
            }
            if (models.keep==T) {
                models$FDA <- results
                models$FDA.PROBIT <- results2
                models$formulae$FDA.formula <- FDA.formula
            }
        }else{ cat(paste("\n", "WARNING: FDA calibration failed", "\n", "\n"))}
    }
    if (ws["SVM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Support Vector Machines (package: kernlab)\n", sep=""))
        cat(paste("\n\n"))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(SVM.OLD) == T) {
            tryCatch(results <- kernlab::ksvm(SVM.formula, data=TrainData, type="C-svc", prob.model=T),
                error= function(err) {print(paste("SVM calibration failed"))},
                silent=T)
        }else{ 
            results <- SVM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"SVM"] <- kernlab::predict(object=results, newdata=TrainData.vars, type="probabilities")[,2]
            if (PROBIT == T) {
                if(is.null(SVM.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ SVM"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- SVM.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"SVM.step1"] <- TrainData[,"SVM"]
                TrainData[,"SVM"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "SVM"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"SVM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"SVM"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["SVM"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["SVM"]))
            weights["SVM"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"SVM"] <- kernlab::predict(object=results, newdata=TestData.vars, type="probabilities")[,2]
                if (PROBIT == T) {
                    TestData[,"SVM.step1"] <- TestData[,"SVM"]
                    TestData[,"SVM"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"SVM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"SVM"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("SVM evaluation (kernlab package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["SVM"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="SVM (kernlab package)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: SVM evaluation (kernlab package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$SVM.C <- eval1
                evaluations$SVM.T <- eval2
            }
            if (models.keep==T) {
                models$SVM <- results
                models$SVM.PROBIT <- results2
                models$formulae$SVM.formula <- SVM.formula
            }
        }else{ cat(paste("\n", "WARNING: SVM calibration failed", "\n", "\n"))}
    }
    if (ws["SVME"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Support Vector Machines (package: e1071)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(SVME.OLD) == T) {
            tryCatch(results <- e1071::svm(SVME.formula, data=TrainData, type="C-classification", kernel="polynomial", degree=3, probability=TRUE),
                error= function(err) {print(paste("SVM calibration (e1071 package) failed"))},
                silent=T)
        }else{ 
            results <- SVME.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"SVME"] <- predict.svme(model=results, newdata=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(SVME.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ SVME"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- SVME.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"SVME.step1"] <- TrainData[,"SVME"]
                TrainData[,"SVME"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "SVME"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"SVME"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"SVME"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["SVME"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["SVME"]))
            weights["SVME"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"SVME"] <- predict.svme(model=results, newdata=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"SVME.step1"] <- TestData[,"SVME"]
                    TestData[,"SVME"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"SVME"]
                TestAbs <- TestData[TestData[,"pb"]==0,"SVME"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("SVM evaluation (e1071 package) failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["SVME"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="SVM (e1017 package)", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: SVM evaluation (e1071 package) failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$SVME.C <- eval1
                evaluations$SVME.T <- eval2
            }
            if (models.keep==T) {
                models$SVME <- results
                models$SVME.PROBIT <- results2
                models$formulae$SVME.formula <- SVME.formula
            }
        }else{ cat(paste("\n", "WARNING: SVM calibration (e1071 package) failed", "\n", "\n"))}
    }
    if (ws["BIOCLIM"] > 0 || ws["DOMAIN"] > 0 || ws["MAHAL"] > 0) {
        if(is.null(factors) == F) {
            for (i in 1:length(factors)) {
                TrainData.vars <- TrainData.vars[, which(colnames(TrainData.vars) != factors[i])]
                TrainData.pres <- TrainData.pres[, which(colnames(TrainData.pres) != factors[i])]
                TestData.vars <- TestData.vars[, which(colnames(TestData.vars) != factors[i])]
            }
        }
        if(is.null(dummy.vars) == F) {
            for (i in 1:length(dummy.vars)) {
                TrainData.vars <- TrainData.vars[, which(colnames(TrainData.vars) != dummy.vars[i])]
                TrainData.pres <- TrainData.pres[, which(colnames(TrainData.pres) != dummy.vars[i])]
                TestData.vars <- TestData.vars[, which(colnames(TestData.vars) != dummy.vars[i])]
            }
        }
    }
    if (ws["BIOCLIM"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". BIOCLIM algorithm (package: dismo)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(BIOCLIM.OLD) == T) {
            tryCatch(results <- dismo::bioclim(x=TrainData.pres),
                error= function(err) {print(paste("BIOCLIM calibration failed"))},
                silent=T)
        }else{ 
            results <- BIOCLIM.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"BIOCLIM"] <- dismo::predict(object=results, x=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(BIOCLIM.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ BIOCLIM"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- BIOCLIM.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"BIOCLIM.step1"] <- TrainData[,"BIOCLIM"]
                TrainData[,"BIOCLIM"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "BIOCLIM"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"BIOCLIM"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"BIOCLIM"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["BIOCLIM"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["BIOCLIM"]))
            weights["BIOCLIM"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"BIOCLIM"] <- dismo::predict(object=results, x=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"BIOCLIM.step1"] <- TestData[,"BIOCLIM"]
                    TestData[,"BIOCLIM"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"BIOCLIM"]
                TestAbs <- TestData[TestData[,"pb"]==0,"BIOCLIM"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("BIOCLIM evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["BIOCLIM"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="BIOCLIM", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: BIOCLIM evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$BIOCLIM.C <- eval1
                evaluations$BIOCLIM.T <- eval2
            }
            if (models.keep==T) {
                models$BIOCLIM <- results
                models$BIOCLIM.PROBIT <- results2
            }
        }else{ cat(paste("\n", "WARNING: BIOCLIM calibration failed", "\n", "\n"))}
    }
    if (ws["DOMAIN"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". DOMAIN algorithm (package: dismo)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(DOMAIN.OLD) == T) {
            tryCatch(results <- dismo::domain(x=TrainData.pres),
                error= function(err) {print(paste("DOMAIN calibration failed"))},
                silent=T)
        }else{ 
            results <- DOMAIN.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"DOMAIN"] <- dismo::predict(object=results, x=TrainData.vars)
            if (PROBIT == T) {
                if(is.null(DOMAIN.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ DOMAIN"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- DOMAIN.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"DOMAIN.step1"] <- TrainData[,"DOMAIN"]
                TrainData[,"DOMAIN"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "DOMAIN"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"DOMAIN"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"DOMAIN"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["DOMAIN"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["DOMAIN"]))
            weights["DOMAIN"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"DOMAIN"] <- dismo::predict(object=results, x=TestData.vars)
                if (PROBIT == T) {
                    TestData[,"DOMAIN.step1"] <- TestData[,"DOMAIN"]
                    TestData[,"DOMAIN"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"DOMAIN"]
                TestAbs <- TestData[TestData[,"pb"]==0,"DOMAIN"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("DOMAIN evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["DOMAIN"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="DOMAIN", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: DOMAIN evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$DOMAIN.C <- eval1
                evaluations$DOMAIN.T <- eval2
            }
            if (models.keep==T) {
                models$DOMAIN <- results
                models$DOMAIN.PROBIT <- results2
            }
        }else{ cat(paste("\n", "WARNING: DOMAIN calibration failed", "\n", "\n"))}
    }
    if (ws["MAHAL"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". Mahalanobis algorithm (package: dismo)\n", sep=""))
        eval1 <- eval2 <- results <- results2 <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(MAHAL.OLD) == T) {
            tryCatch(results <- dismo::mahal(x=TrainData.pres),
                error= function(err) {print(paste("Mahalanobis calibration failed"))},
                silent=T)
        }else{ 
            results <- MAHAL.OLD
        }
        if (is.null(results) == F) {
            cat(paste("\n", "Evaluation with calibration data","\n",sep = ""))
            TrainData[,"MAHAL"] <- predict.mahal(model=results, newdata=TrainData.vars, MAHAL.shape=MAHAL.shape)
            if (PROBIT == T) {
                if(is.null(MAHAL.PROBIT.OLD) == T) { 
                    probit.formula <- as.formula(paste("pb ~ MAHAL"))
                    results2 <- glm(probit.formula, family=binomial(link="probit"), data=TrainData, weights=Yweights, control=glm.control(maxit=maxit))
                }else{ 
                    results2 <- MAHAL.PROBIT.OLD
                }
                cat(paste("(Predictions transformed with probit link)","\n\n", sep = ""))
                TrainData[,"MAHAL.step1"] <- TrainData[,"MAHAL"]
                TrainData[,"MAHAL"] <- predict.glm(object=results2, newdata=TrainData, type="response")
            }
            pred1 <- TrainData[, "MAHAL"]
            pred1[pred1 == 0] <- 0.0000000001
            pred1[pred1 == 1] <- 0.9999999999
            cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
            TrainPres <- TrainData[TrainData[,"pb"]==1,"MAHAL"]
            TrainAbs <- TrainData[TrainData[,"pb"]==0,"MAHAL"]
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["MAHAL"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["MAHAL"]))
            weights["MAHAL"] <- max(c(eval1@auc, 0), na.rm=T)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n",sep = ""))
                TestData[,"MAHAL"] <- predict.mahal(model=results, newdata=TestData.vars, MAHAL.shape=MAHAL.shape)
                if (PROBIT == T) {
                    TestData[,"MAHAL.step1"] <- TestData[,"MAHAL"]
                    TestData[,"MAHAL"] <- predict.glm(object=results2, newdata=TestData, type="response")
                }
                TestPres <- TestData[TestData[,"pb"]==1,"MAHAL"]
                TestAbs <- TestData[TestData[,"pb"]==0,"MAHAL"]
                tryCatch(eval2 <- dismo::evaluate(p=TestPres, a=TestAbs),
                    error= function(err) {print(paste("MAHAL evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    weights["MAHAL"] <- max(c(eval2@auc, 0), na.rm=T)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="MAHAL", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: Mahalanobis evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$MAHAL.C <- eval1
                evaluations$MAHAL.T <- eval2
            }
            if (models.keep==T) {
                models$MAHAL <- results
                models$MAHAL.PROBIT <- results2
                models$formulae$MAHAL.shape <- MAHAL.shape
            }
        }else{ cat(paste("\n", "WARNING: Mahalanobis calibration failed", "\n", "\n"))}
    }
    if (ws["GEODIST"] > 0) {
        mc <- mc+1
        cat(paste("\n", mc, ". geoDist algorithm (package: dismo)\n", sep=""))
        eval1 <- eval2 <- results <- TrainPres <- TrainAbs <- TestPres <- TestAbs <- NULL
        if(is.null(GEODIST.OLD) == T) {
            tryCatch(results <- dismo::geoDist(p=p, lonlat=TRUE),
                error= function(err) {print(paste("GEODIST calibration failed"))},
                silent=F)
        }else{ 
            results <- GEODIST.OLD
        }
        if (is.null(results) == F) {
            NAmask <- raster::crop(x[[1]], x[[1]])
            fullname <- paste("models/", species.name, "_GEO", sep="")
            pgeo <- dismo::predict(object=results, x=NAmask, mask=TRUE, 
                 filename=fullname, progress='text', overwrite=TRUE, format=RASTER.format)
            cat(paste("\n", "Evaluation with calibration data","\n\n",sep = ""))
            pres_geo <- raster::extract(pgeo, p)
            abs_geo <- raster::extract(pgeo, a)
            eval1 <- dismo::evaluate(p=pres_geo, a=abs_geo, tr=quantile(pgeo, na.rm=T, probs=c(1:50/50), names=F))
            print(eval1)
            if (no.tests == F) {
                cat(paste("\n", "Evaluation with test data","\n\n","\n", sep = ""))
                prest_geo <- raster::extract(pgeo, pt)
                abst_geo <- raster::extract(pgeo, at)
                tryCatch(eval2 <- dismo::evaluate(p=prest_geo, a=abst_geo, tr=quantile(pgeo, na.rm=T, probs=c(1:50/50), names=F)),
                    error= function(err) {print(paste("GEODIST evaluation failed"))},
                    silent=F)
                if (is.null(eval2) == F) {
                    print(eval2)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="GEODIST", cex=2, adj=0, col.main="blue")
                    }
                }else{
                    cat(paste("\n", "WARNING: GEODIST evaluation failed","\n\n",sep = ""))
                }
            }
            if(evaluations.keep ==T) {
                evaluations$GEODIST.C <- eval1
                evaluations$GEODIST.T <- eval2
            }
# results of GEODIST calibration not kept as GEODIST is a null model
        }else{ cat(paste("\n", "WARNING: GEODIST calibration failed", "\n", "\n"))}
    }
#
    if (AUC.weights == T) {
        if(length(ENSEMBLE.exponent) > 1 || length(ENSEMBLE.best) > 1 || length(ENSEMBLE.min) > 1) {ENSEMBLE.tune <- TRUE}
    }
    if(ENSEMBLE.tune == F) {
        if (sum(ws, na.rm=T) > 0) {
            if (AUC.weights == T) {
                cat(paste("\n", "Weights based on AUC", "\n", sep = ""))
                print(weights)
                weights <- ensemble.weights(weights=weights, exponent=ENSEMBLE.exponent, best=ENSEMBLE.best, min.weight=ENSEMBLE.min)
                cat(paste("\n", "Weights based on parameters ENSEMBLE.min=", ENSEMBLE.min, ", ENSEMBLE.best=", ENSEMBLE.best, 
                    " and ENSEMBLE.exponent=", ENSEMBLE.exponent, "\n", sep = ""))
                print(weights)
                ws <- weights
                cat(paste("\n", "Minimum input weight is 0.05", "\n", sep=""))
                ws[ws < 0.05] <- 0
                ws <- ensemble.weights(weights=ws, exponent=1, best=0, min.weight=0)
                cat(paste("\n", "Weights for ensemble forecasting", "\n", sep = ""))
                print(ws)
            }else{
                cat(paste("\n", "Ensemble weights based directly on input weights scaled to sum up to 1", "\n", sep = ""))
                print(ws)
            }
        }
        if(evaluations.keep == T) {evaluations$ensemble.weights <- ws}
        if(models.keep==T) {models$output.weights <- ws}
    }else{
# use different strategies for calculating the ensemble model
# similar to using test data for calculating input AUC, use internal test data for calculating best ensemble
# recalculating AUC does not require much computing time - initial calculations kept to spot problems for specific algorithms
        strategy.results <- ensemble.strategy(TrainData=TrainData, TestData=TestData,
            ENSEMBLE.exponent=ENSEMBLE.exponent, ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min)
        ws <- strategy.results$weights
        if (sum(ws, na.rm=T) > 0) {
            cat(paste("\n", "Minimum input weight is 0.05", "\n", sep=""))
            ws[ws < 0.05] <- 0
            ws <- ensemble.weights(weights=ws, exponent=1, best=0, min.weight=0)
            cat(paste("\n", "Weights for ensemble forecasting", "\n", sep = ""))
            print(ws)
        }else{
            ENSEMBLE.tune <- FALSE
        }
        if(evaluations.keep == T) {evaluations$STRATEGY.weights <- ws}
        if(models.keep==T) {models$output.weights <- ws}
    }
# do not return ensemble tests for no models when tuning is not implemented 
    if((sum(ws > 0, na.rm=T) > 0) || (ENSEMBLE.tune == T)) {
        TrainData[,"ENSEMBLE"] <- ws["MAXENT"]*TrainData[,"MAXENT"] + ws["GBM"]*TrainData[,"GBM"] +
            ws["GBMSTEP"]*TrainData[,"GBMSTEP"] + ws["RF"]*TrainData[,"RF"] + ws["GLM"]*TrainData[,"GLM"] +
            ws["GLMSTEP"]*TrainData[,"GLMSTEP"] + ws["GAM"]*TrainData[,"GAM"] + ws["GAMSTEP"]*TrainData[,"GAMSTEP"] +
            ws["MGCV"]*TrainData[,"MGCV"] + ws["MGCVFIX"]*TrainData[,"MGCVFIX"] + ws["EARTH"]*TrainData[,"EARTH"] +
            ws["RPART"]*TrainData[,"RPART"] + ws["NNET"]*TrainData[,"NNET"] + ws["FDA"]*TrainData[,"FDA"] +
            ws["SVM"]*TrainData[,"SVM"] + ws["SVME"]*TrainData[,"SVME"] + ws["BIOCLIM"]*TrainData[,"BIOCLIM"] +
            ws["DOMAIN"]*TrainData[,"DOMAIN"] + ws["MAHAL"]*TrainData[,"MAHAL"]
        pred1 <- TrainData[, "ENSEMBLE"]
        pred1[pred1 == 0] <- 0.0000000001
        pred1[pred1 == 1] <- 0.9999999999
        if (no.tests == F) {
            TestData[,"ENSEMBLE"] <- ws["MAXENT"]*TestData[,"MAXENT"] + ws["GBM"]*TestData[,"GBM"] +
                ws["GBMSTEP"]*TestData[,"GBMSTEP"] + ws["RF"]*TestData[,"RF"] + ws["GLM"]*TestData[,"GLM"] +
                ws["GLMSTEP"]*TestData[,"GLMSTEP"] + ws["GAM"]*TestData[,"GAM"] + ws["GAMSTEP"]*TestData[,"GAMSTEP"] +
                ws["MGCV"]*TestData[,"MGCV"] + ws["MGCVFIX"]*TestData[,"MGCVFIX"] + ws["EARTH"]*TestData[,"EARTH"] +
                ws["RPART"]*TestData[,"RPART"] + ws["NNET"]*TestData[,"NNET"] + ws["FDA"]*TestData[,"FDA"] +
                ws["SVM"]*TestData[,"SVM"] + ws["SVME"]*TestData[,"SVME"] + ws["BIOCLIM"]*TestData[,"BIOCLIM"] +
                ws["DOMAIN"]*TestData[,"DOMAIN"] + ws["MAHAL"]*TestData[,"MAHAL"]
        }
        mc <- mc+1
        cat(paste("\n\n", mc, ". Ensemble algorithm\n", sep=""))
        eval1 <- eval2 <- NULL
        cat(paste("\n", "Ensemble evaluation with calibration data", "\n\n", sep = ""))
        cat(paste("Residual deviance (dismo package): ", dismo::calc.deviance(obs=obs1, pred=pred1, calc.mean=F), "\n\n", sep = ""))
        TrainPres <- as.numeric(TrainData[TrainData[,"pb"]==1,"ENSEMBLE"])
        TrainAbs <- as.numeric(TrainData[TrainData[,"pb"]==0,"ENSEMBLE"])
        if (sum(TrainPres, na.rm=T) <= 0 || sum(TrainAbs, na.rm=T) <= 0) {
            cat(paste("\n", "NOTE: not possible to evaluate the ensemble model since calibration probabilities not available", "\n", sep = ""))
        }else{
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            print(eval1)
            thresholds["ENSEMBLE"] <- threshold2(eval1, threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence, Pres=TrainPres, Abs=TrainAbs)
            cat(paste("\n", "Threshold (method: ", threshold.method, ") \n", sep = ""))
            print(as.numeric(thresholds["ENSEMBLE"]))
            if (models.keep == T) {models$thresholds <- thresholds}
            if (no.tests == F) {
                cat(paste("\n", "Ensemble evaluation with testing data", "\n\n", sep = ""))
                TestPres <- as.numeric(TestData[TestData[,"pb"]==1,"ENSEMBLE"])
                TestAbs <- as.numeric(TestData[TestData[,"pb"]==0,"ENSEMBLE"])
                if (sum(TestPres, na.rm=T) <= 0 || sum(TestAbs, na.rm=T) <= 0) {
                    cat(paste("\n", "NOTE: not possible to evaluate the ensemble model since evaluation probabilities not available", "\n", sep = ""))
                }else{
                    eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
                    print(eval2)
                    if(PLOTS==T) {
                        graphics::plot(eval2, "ROC") 
                        graphics::title(main="ENSEMBLE", cex=2, adj=0, col.main="blue")
                    }
                }
            }
            if(evaluations.keep==T) {
                evaluations$ENSEMBLE.C <- eval1
                evaluations$ENSEMBLE.T <- eval2
            }
        }
    }
    if(models.keep==T) {
        models$TrainData <- TrainData
        models$TestData <- TestData
        models$var.names <- var.names
        if (ws["MAXENT"] > 0 || MAXENT2 > 0) {models$MAXENT.BackData <- MAXENT.BackData}
    }
    if(models.keep==F && evaluations.keep==T) {
        evaluations$TrainData <- TrainData
        evaluations$TestData <- TestData
        evaluations$p <- p
        evaluations$pt <- pt
        evaluations$a <- a
        evaluations$at <- at
        evaluations$var.names <- var.names
        if (ws["MAXENT"] > 0 || MAXENT2 > 0) {evaluations$MAXENT.BackData <- MAXENT.BackData}
    }
    remove(Yweights, envir=.BiodiversityR)
    remove(TrainData, envir=.BiodiversityR)
    remove(TestData, envir=.BiodiversityR)
    remove(TrainData.vars, envir=.BiodiversityR)
    remove(TrainData.pres, envir=.BiodiversityR)
    remove(TestData.vars, envir=.BiodiversityR)
    if (models.save==T && models.keep==T && MAXENT2==0) {
        ensemble.models <- models
        save(ensemble.models, file=paste(getwd(), "/models/", models$species.name, "_models", sep=""), compress="xz")
    }
    if (models.keep == F) {models <- NULL}     
    result <- list(evaluations=evaluations, models=models, VIF=newVIF, call=match.call() )
    cat(paste("\n\n"))
    if (SINK==T && OLD.SINK==F) {sink(file=NULL, append=T)}
    return(result)
}

