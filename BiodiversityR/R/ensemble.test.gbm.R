`ensemble.test.gbm` <- function(
    x=NULL, p=NULL, a=NULL, an=1000, excludep=FALSE, ext=NULL, k=4, 
    TrainData=NULL,
    VIF=FALSE, COR=FALSE,
    SINK=FALSE, PLOTS=FALSE, 
    species.name="Species001",
    Yweights="BIOMOD", 
    layer.drops=NULL, factors=NULL,
    GBMSTEP.gbm.x=2:(ncol(TrainData.orig)), 
    complexity=c(3:6), learning=c(0.005, 0.002, 0.001), 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100
)
{
    .BiodiversityR <- new.env()
#   if (! require(dismo)) {stop("Please install the dismo package")}
#    if (! require(gbm)) {stop("Please install the gbm package")}
    k <- as.integer(k)
    if (k < 2) {
        cat(paste("\n", "NOTE: parameter k was set to be smaller than 2", sep = ""))
        cat(paste("\n", "default value of 5 therefore set for parameter k", "\n", sep = ""))
        k <- 5
    }
# check data
    if (is.null(TrainData) == T) {
        if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
        if(inherits(x,"RasterStack") == F) {stop("x is not a RasterStack object")}
        if(raster::projection(x)=="NA") {
            raster::projection(x) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        }
        if(is.null(p) == T) {stop("presence locations are missing (parameter p)")}
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
        cat(paste("\n\n", "RESULTS (ensemble.test.gbm function)", "\n", sep=""), file=paste.file, append=T)
        sink(file=paste.file, append=T)
        cat(paste(date(), "\n", sep=""))
        print(match.call())
    }

# check TrainData
    if (is.null(TrainData) == F) {
        TrainData <- data.frame(TrainData)
        if (colnames(TrainData)[1] !="pb") {stop("first column for TrainData should be 'pb' containing presence (1) and absence (0) data")}
        if ((is.null(x) == F) && (raster::nlayers(x) != (ncol(TrainData)-1))) {
            cat(paste("\n", "WARNING: different number of explanatory variables in rasterStack and TrainData", sep = ""))
        }
    }

# modify RasterStack x only if this RasterStack was provided
    if (is.null(x) == F) {
        if (is.null(ext) == F) {
            if(length(x@title) == 0) {x@title <- "stack1"}
            title.old <- x@title
            x <- raster::crop(x, y=ext, snap="in")
            x@title <- title.old
        }
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
                    if (is.null(factors) == F) {
                        factors <- factors[factors != layer.drops[i]]
                        if(length(factors) == 0) {factors <- NULL}
                    }
                    if (is.null(dummy.vars) == F) {
                        dummy.vars <- dummy.vars[dummy.vars != layer.drops[i]]
                        if(length(dummy.vars) == 0) {dummy.vars <- NULL}
                    }
                }
            }
        }
        if (is.null(factors) == F) {
            vars <- names(x)
            factors <- as.character(factors)
            nf <- length(factors)
            for (i in 1:nf) {
                if (any(vars==factors[i])==FALSE) {
                    cat(paste("\n", "WARNING: categorical variable '", factors[i], "' not among grid layers", "\n", sep = ""))
                    factors <- factors[factors != factors[i]]
                    if(length(factors) == 0) {factors <- NULL}
                }
            }
        }
        if (is.null(dummy.vars) == F) {
            vars <- names(x)
            dummy.vars <- as.character(dummy.vars)
            nf <- length(dummy.vars)
            for (i in 1:nf) {
                if (any(vars==dummy.vars[i])==FALSE) {
                    cat(paste("\n", "WARNING: dummy variable '", dummy.vars[i], "' not among grid layers", "\n", sep = ""))
                    dummy.vars <- dummy.vars[dummy.vars != dummy.vars[i]]
                    if(length(dummy.vars) == 0) {dummy.vars <- NULL}
                }
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
                if (any(vars==layer.drops[i])==FALSE) {
                    cat(paste("\n", "WARNING: variable to exclude '", layer.drops[i], "' not among columns of TrainData", "\n", sep = ""))
                }else{
                    cat(paste("\n", "NOTE: variable '", layer.drops[i], "' will not be included as explanatory variable", "\n", sep = ""))
                    TrainData <- TrainData[, which(colnames(TrainData) != layer.drops[i])]
                    if (is.null(TestData) == F) {TestData <- TestData[, which(colnames(TestData) != layer.drops[i])]}
                    vars <- colnames(TrainData)
                    if (is.null(factors) == F) {
                        factors <- factors[factors != layer.drops[i]]
                        if(length(factors) == 0) {factors <- NULL}
                    }
                    if (is.null(dummy.vars) == F) {
                        dummy.vars <- dummy.vars[dummy.vars != layer.drops[i]]
                        if(length(dummy.vars) == 0) {dummy.vars <- NULL}
                    }
                }
            }
        }
        if (is.null(factors) == F) {
            vars <- names(TrainData)
            factors <- as.character(factors)
            nf <- length(factors)
            for (i in 1:nf) {
                if (any(vars==factors[i])==FALSE) {
                     cat(paste("\n", "WARNING: categorical variable '", factors[i], "' not among columns of TrainData", "\n", sep = ""))
                    factors <- factors[factors != factors[i]]
                    if(length(factors) == 0) {factors <- NULL}
                }
            }
        }
        if (is.null(dummy.vars) == F) {
            vars <- names(TrainData)
            dummy.vars <- as.character(dummy.vars)
            nf <- length(dummy.vars)
            for (i in 1:nf) {
                if (any(vars==dummy.vars[i])==FALSE) {
                    cat(paste("\n", "WARNING: dummy variable '", dummy.vars[i], "' not among columns of TrainData", "\n", sep = ""))
                    dummy.vars <- dummy.vars[dummy.vars != dummy.vars[i]]
                    if(length(dummy.vars) == 0) {dummy.vars <- NULL}
                }
            }
        }
    }

# create TrainData and TestData
    if (is.null(TrainData) == F) {
        if(any(is.na(TrainData))) {
            cat(paste("\n", "WARNING: sample units with missing data removed from calibration data","\n\n",sep = ""))
        }
        TrainValid <- complete.cases(TrainData)
        TrainData <- TrainData[TrainValid,]
    }else{
        if (is.null(a)==T) {
            if (excludep == T) {
                a <- dismo::randomPoints(x[[1]], n=an, p=p, ext=ext, excludep=T)
            }else{
                a <- dismo::randomPoints(x[[1]], n=an, p=NULL, ext=ext, excludep=F)
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
    TrainData.orig <- TrainData
    assign("TrainData.orig", TrainData.orig, envir=.BiodiversityR)
#
    nc <- length(complexity)
    nl <- length(learning)
    nt <- nc*nl
    output <- array(NA, dim=c(nt, 2*k+3))
    colnames(output) <- c("tree.complexity","learning.rate", 1:k,"MEAN",1:k)
    output[,"tree.complexity"] <- rep(complexity, nl)
    output[,"learning.rate"] <- rep(learning, each=nc) 
#
    groupp <- dismo::kfold(TrainData, k=k, by=TrainData[,"pb"])
    for (i in 1:k){
        cat(paste("\n", "EVALUATION RUN: ", i, "\n\n", sep = ""))
        TrainData.c <- TrainData[groupp != i,]
        TestData.c <- TrainData[groupp == i,]
        assign("TrainData.c", TrainData.c, envir=.BiodiversityR)
        assign("TestData.c", TestData.c, envir=.BiodiversityR)
        for (j in 1:nt) {
            complex <- output[j,"tree.complexity"]
            lr <- output[j, "learning.rate"]
            cat(paste("\n", "complexity: ", complex, ", learning: ", lr, "\n", sep=""))
            tests <- ensemble.test(x=x,
                TrainData=TrainData.c, TestData=TestData.c,
                VIF=VIF, COR=COR,
                PLOTS=PLOTS, evaluations.keep=T,
                MAXENT=0, GBM=0, GBMSTEP=1, RF=0, GLM=0, GLMSTEP=0, 
                GAM=0, GAMSTEP=0, MGCV=0, MGCVFIX=0, EARTH=0, RPART=0, 
                NNET=0, FDA=0, SVM=0, SVME=0, BIOCLIM=0, DOMAIN=0, MAHAL=0, GEODIST=0,   
                Yweights=Yweights, factors=factors,   
                GBMSTEP.gbm.x=2:(1+raster::nlayers(x)), 
                GBMSTEP.bag.fraction=GBMSTEP.bag.fraction, 
                GBMSTEP.tree.complexity=complex, GBMSTEP.learning.rate=lr, 
                GBMSTEP.step.size=GBMSTEP.step.size)
            output[j,2+i] <- tests$GBMSTEP.T@auc
            output[j,k+3+i] <- tests$GBMSTEP.trees 
        }
    }
    output[,k+3] <- rowMeans(output[,3:(k+2)], na.rm=T)
    output <- output[order(output[,"MEAN"], decreasing=T),]
    output[,3:(k+3)] <- 100*output[,3:(k+3)]
    print(output)
    cat(paste("\n\n"))
    if (SINK==T  && OLD.SINK==F) {sink(file=NULL, append=T)}
    return(list(table=output, call=match.call() ))
}


