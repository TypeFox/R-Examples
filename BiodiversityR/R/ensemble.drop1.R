`ensemble.drop1` <- function(
    x=NULL, p=NULL, a=NULL, an=1000, excludep=FALSE, ext=NULL, 
    k=0, pt=NULL, at=NULL, CIRCLES.at=FALSE, CIRCLES.d=100000,
    TrainData=NULL, TestData=NULL,
    VIF=FALSE, COR=FALSE,
    SINK=FALSE, species.name="Species001",
    difference=FALSE,
    ENSEMBLE.best=0, ENSEMBLE.min=0.7, ENSEMBLE.exponent=1,
    input.weights=NULL,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0, 
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1, 
    PROBIT=FALSE,
    Yweights="BIOMOD", 
    layer.drops=NULL, factors=NULL, dummy.vars=NULL,
    maxit=100,
    MAXENT.a=NULL, MAXENT.an=10000, MAXENT.BackData=NULL, 
    MAXENT.path=paste(getwd(), "/models/maxent_", species.name,  sep=""), 
    GBM.n.trees=2001, 
    GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100, 
    RF.ntree=751,  
    GLM.family=binomial(link="logit"), 
    GLMSTEP.steps=1000, GLMSTEP.scope=NULL, GLMSTEP.k=2, 
    GAM.family=binomial(link="logit"), 
    GAMSTEP.steps=1000, GAMSTEP.scope=NULL, GAMSTEP.pos=1, 
    MGCV.select=FALSE,  
    EARTH.glm=list(family=binomial(link="logit"), maxit=maxit), 
    RPART.xval=50, 
    NNET.size=8, NNET.decay=0.01, 
    MAHAL.shape=1
)
{

    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
#
    if (is.null(layer.drops) == F) {
        layer.drops <- as.character(layer.drops)
        if (is.null(x)==F) {x <- raster::dropLayer(x, which(names(x) %in% layer.drops))}
        factors <- as.character(factors)
        dummy.vars <- as.character(dummy.vars)
        nd <- length(layer.drops)
        for (i in 1:nd) {
            if (is.null(factors) == F) {
                factors <- factors[factors != layer.drops[i]]
                if(length(factors) == 0) {factors <- NULL}
            }
            if (is.null(dummy.vars) == F) {
                dummy.vars <- dummy.vars[dummy.vars != layer.drops[i]]
                if(length(dummy.vars) == 0) {dummy.vars <- NULL}
            }
        }
        if(length(layer.drops) == 0) {layer.drops <- NULL}
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
        cat(paste("\n\n", "RESULTS (ensemble.drop1 function)", "\n", sep=""), file=paste.file, append=T)
        sink(file=paste.file, append=T)
        cat(paste(date(), "\n", sep=""))
        print(match.call())
    }

# estimate deviance
    loglik.calculation <- function(obs=NULL, preds=NULL) {
        preds[preds<=0] <- 0.0000000001
        preds[preds>=1] <- 0.9999999999
        out <- dismo::calc.deviance(obs=obs, pred=preds, calc.mean=F)
        return(out)
    }

#
# first fit with all variables

    cat(paste("\n\n", "RESULTS WITH ALL VARIABLES", "\n", sep=""))

    tests <- ensemble.test(x=x, ext=ext,
        p=p, a=a, an=an, excludep=excludep, 
        k=k, pt=pt, at=at, CIRCLES.at=CIRCLES.at, CIRCLES.d=CIRCLES.d,
        TrainData=TrainData, TestData=TestData,
        PLOTS=FALSE, evaluations.keep=T, models.keep=F,
        VIF=VIF, COR=COR,
        formulae.defaults=T, maxit=maxit,
        AUC.weights=T,
        ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min,
        ENSEMBLE.exponent=ENSEMBLE.exponent,
        input.weights=input.weights,
        MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
        GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
        NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,
        GEODIST=0, 
        PROBIT=PROBIT,
        Yweights=Yweights, 
        layer.drops=NULL, factors=factors, dummy.vars=dummy.vars,
        MAXENT.BackData=MAXENT.BackData, MAXENT.path=MAXENT.path,
        GBM.formula=NULL, GBM.n.trees=GBM.n.trees,
        GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
        GBMSTEP.step.size=GBMSTEP.step.size,
        RF.formula=NULL, RF.ntree=RF.ntree, 
        GLM.formula=NULL, GLM.family=GLM.family, 
        GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=NULL, GLMSTEP.scope=NULL, 
        GAM.formula=NULL, GAM.family=GAM.family, 
        GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=NULL, GAMSTEP.pos=GAMSTEP.pos,
        MGCV.formula=NULL, MGCV.select=MGCV.select,
        MGCVFIX.formula=NULL, 
        EARTH.formula=NULL, EARTH.glm=EARTH.glm,
        RPART.formula=NULL, RPART.xval=RPART.xval, 
        NNET.formula=NULL, NNET.size=NNET.size, NNET.decay=NNET.decay,
        FDA.formula=NULL, SVM.formula=NULL, SVME.formula=NULL,
        MAHAL.shape=MAHAL.shape)

# use output to get names of the variables 
    var.names <- tests$evaluations$var.names
    nv <- length(var.names)

    output.C <- array(NA, dim=c(20, nv+1))
    rownames(output.C) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "ENSEMBLE")
    colnames(output.C) <- c("all_vars", paste("without_", var.names, sep=""))

    output.T <- array(NA, dim=c(20, nv+1))
    rownames(output.T) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "ENSEMBLE")
    colnames(output.T) <- c("all_vars", paste("without_", var.names, sep=""))

    output.LLC <- array(NA, dim=c(20, nv+1))
    rownames(output.LLC) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "ENSEMBLE")
    colnames(output.LLC) <- c("all_vars", paste("without_", var.names, sep=""))

    output.LLT <- array(NA, dim=c(20, nv+1))
    rownames(output.LLT) <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "ENSEMBLE")
    colnames(output.LLT) <- c("all_vars", paste("without_", var.names, sep=""))

    if(is.null(tests$evaluations$MAXENT.C)==F) {output.C["MAXENT",1] <- tests$evaluations$MAXENT.C@auc}
    if(is.null(tests$evaluations$MAXENT.T)==F) {output.T["MAXENT",1] <- tests$evaluations$MAXENT.T@auc}
    if(sum(tests$evaluations$TrainData$MAXENT) > 0) {output.LLC["MAXENT",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAXENT)}
    if(sum(tests$evaluations$TestData$MAXENT) > 0) {output.LLT["MAXENT",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAXENT)}

    if(is.null(tests$evaluations$GBM.C)==F) {output.C["GBM",1] <- tests$evaluations$GBM.C@auc} 
    if(is.null(tests$evaluations$GBM.T)==F) {output.T["GBM",1] <- tests$evaluations$GBM.T@auc} 
    if(sum(tests$evaluations$TrainData$GBM) > 0) {output.LLC["GBM",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GBM)}
    if(sum(tests$evaluations$TestData$GBM) > 0) {output.LLT["GBM",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GBM)}

    if(is.null(tests$evaluations$GBMSTEP.C)==F) {output.C["GBMSTEP",1] <- tests$evaluations$GBMSTEP.C@auc} 
    if(is.null(tests$evaluations$GBMSTEP.T)==F) {output.T["GBMSTEP",1] <- tests$evaluations$GBMSTEP.T@auc} 
    if(sum(tests$evaluations$TrainData$GBMSTEP) > 0) {output.LLC["MAXENT",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAXENT)}
    if(sum(tests$evaluations$TestData$GBMSTEP) > 0) {output.LLT["MAXENT",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAXENT)}

    if(is.null(tests$evaluations$RF.C)==F) {output.C["RF",1] <- tests$evaluations$RF.C@auc}
    if(is.null(tests$evaluations$RF.T)==F) {output.T["RF",1] <- tests$evaluations$RF.T@auc}
    if(sum(tests$evaluations$TrainData$RF) > 0) {output.LLC["RF",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$RF)}
    if(sum(tests$evaluations$TestData$RF) > 0) {output.LLT["RF",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$RF)}

    if(is.null(tests$evaluations$GLM.C)==F) {output.C["GLM",1] <- tests$evaluations$GLM.C@auc} 
    if(is.null(tests$evaluations$GLM.T)==F) {output.T["GLM",1] <- tests$evaluations$GLM.T@auc} 
    if(sum(tests$evaluations$TrainData$GLM) > 0) {output.LLC["GLM",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GLM)}
    if(sum(tests$evaluations$TestData$GLM) > 0) {output.LLT["GLM",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GLM)}

    if(is.null(tests$evaluations$GLMS.C)==F) {output.C["GLMSTEP",1] <- tests$evaluations$GLMS.C@auc}
    if(is.null(tests$evaluations$GLMS.T)==F) {output.T["GLMSTEP",1] <- tests$evaluations$GLMS.T@auc}
    if(sum(tests$evaluations$TrainData$GLMSTEP) > 0) {output.LLC["GLMSTEP",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GLMSTEP)}
    if(sum(tests$evaluations$TestData$GLMSTEP) > 0) {output.LLT["GLMSTEP",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GLMSTEP)}

    if(is.null(tests$evaluations$GAM.C)==F) {output.C["GAM",1] <- tests$evaluations$GAM.C@auc} 
    if(is.null(tests$evaluations$GAM.T)==F) {output.T["GAM",1] <- tests$evaluations$GAM.T@auc} 
    if(sum(tests$evaluations$TrainData$GAM) > 0) {output.LLC["GAM",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GAM)}
    if(sum(tests$evaluations$TestData$GAM) > 0) {output.LLT["GAM",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GAM)}

    if(is.null(tests$evaluations$GAMS.C)==F) {output.C["GAMSTEP",1] <- tests$evaluations$GAMS.C@auc}
    if(is.null(tests$evaluations$GAMS.T)==F) {output.T["GAMSTEP",1] <- tests$evaluations$GAMS.T@auc}
    if(sum(tests$evaluations$TrainData$GAMSTEP) > 0) {output.LLC["GAMSTEP",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GAMSTEP)}
    if(sum(tests$evaluations$TestData$GAMSTEP) > 0) {output.LLT["GAMSTEP",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GAMSTEP)}

    if(is.null(tests$evaluations$MGCV.C)==F) {output.C["MGCV",1] <- tests$evaluations$MGCV.C@auc} 
    if(is.null(tests$evaluations$MGCV.T)==F) {output.T["MGCV",1] <- tests$evaluations$MGCV.T@auc} 
    if(sum(tests$evaluations$TrainData$MGCV) > 0) {output.LLC["MGCV",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MGCV)}
    if(sum(tests$evaluations$TestData$MGCV) > 0) {output.LLT["MGCV",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MGCV)}

    if(is.null(tests$evaluations$MGCVF.C)==F) {output.C["MGCVFIX",1] <- tests$evaluations$MGCVF.C@auc} 
    if(is.null(tests$evaluations$MGCVF.T)==F) {output.T["MGCVFIX",1] <- tests$evaluations$MGCVF.T@auc} 
    if(sum(tests$evaluations$TrainData$MGCVFIX) > 0) {output.LLC["MGCVFIX",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MGCVFIX)}
    if(sum(tests$evaluations$TestData$MGCVFIX) > 0) {output.LLT["MGCVFIX",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MGCVFIX)}

    if(is.null(tests$evaluations$EARTH.C)==F) {output.C["EARTH",1] <- tests$evaluations$EARTH.C@auc} 
    if(is.null(tests$evaluations$EARTH.T)==F) {output.T["EARTH",1] <- tests$evaluations$EARTH.T@auc} 
    if(sum(tests$evaluations$TrainData$EARTH) > 0) {output.LLC["EARTH",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$EARTH)}
    if(sum(tests$evaluations$TestData$EARTH) > 0) {output.LLT["EARTH",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$EARTH)}

    if(is.null(tests$evaluations$RPART.C)==F) {output.C["RPART",1] <- tests$evaluations$RPART.C@auc}
    if(is.null(tests$evaluations$RPART.T)==F) {output.T["RPART",1] <- tests$evaluations$RPART.T@auc}
    if(sum(tests$evaluations$TrainData$RPART) > 0) {output.LLC["RPART",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$RPART)}
    if(sum(tests$evaluations$TestData$RPART) > 0) {output.LLT["RPART",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$RPART)}

    if(is.null(tests$evaluations$NNET.C)==F) {output.C["NNET",1] <- tests$evaluations$NNET.C@auc} 
    if(is.null(tests$evaluations$NNET.T)==F) {output.T["NNET",1] <- tests$evaluations$NNET.T@auc} 
    if(sum(tests$evaluations$TrainData$NNET) > 0) {output.LLC["NNET",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$NNET)}
    if(sum(tests$evaluations$TestData$NNET) > 0) {output.LLT["NNET",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$NNET)}

    if(is.null(tests$evaluations$FDA.C)==F) {output.C["FDA",1] <- tests$evaluations$FDA.C@auc}
    if(is.null(tests$evaluations$FDA.T)==F) {output.T["FDA",1] <- tests$evaluations$FDA.T@auc}
    if(sum(tests$evaluations$TrainData$FDA) > 0) {output.LLC["FDA",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$FDA)}
    if(sum(tests$evaluations$TestData$FDA) > 0) {output.LLT["FDA",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$FDA)}

    if(is.null(tests$evaluations$SVM.C)==F) {output.C["SVM",1] <- tests$evaluations$SVM.C@auc}
    if(is.null(tests$evaluations$SVM.T)==F) {output.T["SVM",1] <- tests$evaluations$SVM.T@auc}
    if(sum(tests$evaluations$TrainData$SVM) > 0) {output.LLC["SVM",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$SVM)}
    if(sum(tests$evaluations$TestData$SVM) > 0) {output.LLT["SVM",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$SVM)}

    if(is.null(tests$evaluations$SVME.C)==F) {output.C["SVME",1] <- tests$evaluations$SVME.C@auc}
    if(is.null(tests$evaluations$SVME.T)==F) {output.T["SVME",1] <- tests$evaluations$SVME.T@auc}
    if(sum(tests$evaluations$TrainData$SVME) > 0) {output.LLC["SVME",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$SVME)}
    if(sum(tests$evaluations$TestData$SVME) > 0) {output.LLT["SVME",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$SVME)}

    if(is.null(tests$evaluations$BIOCLIM.C)==F) {output.C["BIOCLIM",1] <- tests$evaluations$BIOCLIM.C@auc}
    if(is.null(tests$evaluations$BIOCLIM.T)==F) {output.T["BIOCLIM",1] <- tests$evaluations$BIOCLIM.T@auc}
    if(sum(tests$evaluations$TrainData$BIOCLIM) > 0) {output.LLC["BIOCLIM",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$BIOCLIM)}
    if(sum(tests$evaluations$TestData$BIOCLIM) > 0) {output.LLT["BIOCLIM",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$BIOCLIM)}

    if(is.null(tests$evaluations$DOMAIN.C)==F) {output.C["DOMAIN",1] <- tests$evaluations$DOMAIN.C@auc}
    if(is.null(tests$evaluations$DOMAIN.T)==F) {output.T["DOMAIN",1] <- tests$evaluations$DOMAIN.T@auc}
    if(sum(tests$evaluations$TrainData$DOMAIN) > 0) {output.LLC["DOMAIN",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$DOMAIN)}
    if(sum(tests$evaluations$TestData$DOMAIN) > 0) {output.LLT["DOMAIN",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$DOMAIN)}

    if(is.null(tests$evaluations$MAHAL.C)==F) {output.C["MAHAL",1] <- tests$evaluations$MAHAL.C@auc}
    if(is.null(tests$evaluations$MAHAL.T)==F) {output.T["MAHAL",1] <- tests$evaluations$MAHAL.T@auc}
    if(sum(tests$evaluations$TrainData$MAHAL) > 0) {output.LLC["MAHAL",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAHAL)}
    if(sum(tests$evaluations$TestData$MAHAL) > 0) {output.LLT["MAHAL",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAHAL)}

    if(is.null(tests$evaluations$ENSEMBLE.C)==F) {output.C["ENSEMBLE",1] <- tests$evaluations$ENSEMBLE.C@auc}
    if(is.null(tests$evaluations$ENSEMBLE.T)==F) {output.T["ENSEMBLE",1] <- tests$evaluations$ENSEMBLE.T@auc}
    if(sum(tests$evaluations$TrainData$ENSEMBLE) > 0) {output.LLC["ENSEMBLE",1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$ENSEMBLE)}
    if(sum(tests$evaluations$TestData$ENSEMBLE) > 0) {output.LLT["ENSEMBLE",1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$ENSEMBLE)}

# sequentially leave out the focal variable, then fit again
# the data sets are used from the "full" model

    var.names2 <- c("pb", var.names)
    TrainData1 <- tests$evaluations$TrainData
    TrainData1 <- TrainData1[, which(names(TrainData1) %in% var.names2)]
    TestData1 <- tests$evaluations$TestData
    TestData1 <- TestData1[, which(names(TestData1) %in% var.names2)]
    MAXENT.BackData1 <- tests$evaluations$MAXENT.BackData
    if (is.null(MAXENT.BackData1) == F) {MAXENT.BackData1 <- MAXENT.BackData1[, which(names(MAXENT.BackData1) %in% var.names)]}
    assign("MAXENT.BackData1", MAXENT.BackData1, envir=.BiodiversityR)


    for (i in 1:nv) {
        var.f <- var.names[i]
        cat(paste("\n", "Leaving out variable ", var.f, "\n\n", sep = ""))
        TrainData2 <- TrainData1[, which(names(TrainData1) != var.f)]
        TestData2 <- TestData1[, which(names(TestData1) != var.f)]
        MAXENT.BackData2 <- NULL
        if (is.null(MAXENT.BackData1) == F) {MAXENT.BackData2 <- MAXENT.BackData1[, which(names(MAXENT.BackData1) != var.f)]}
        assign("MAXENT.BackData2", MAXENT.BackData2, envir=.BiodiversityR)
        factors2 <- NULL
        if (is.null(factors) == F) {
            factors2 <- factors[which(factors != var.f)]
            if (length(factors2) == 0) {factors2 <- NULL}
        }
        dummy.vars2 <- NULL
        if (is.null(dummy.vars) == F) {
            dummy.vars2 <- dummy.vars[which(dummy.vars != var.f)]
            if (length(dummy.vars2) == 0) {dummy.vars2 <- NULL}
        }

    tests <- ensemble.test(x=NULL, ext=ext,
        TrainData=TrainData2, TestData=TestData2, 
        PLOTS=FALSE, evaluations.keep=T,  
        VIF=VIF, COR=COR,
        formulae.defaults=T, maxit=maxit,
        AUC.weights=T,
        ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min,
        ENSEMBLE.exponent=ENSEMBLE.exponent,
        input.weights=input.weights,
        MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
        GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
        NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,
        GEODIST=0, 
        PROBIT=PROBIT,
        Yweights=Yweights, 
        layer.drops=NULL, factors=factors2, dummy.vars=dummy.vars2,
        MAXENT.BackData=MAXENT.BackData2, MAXENT.path=MAXENT.path,
        GBM.formula=NULL, GBM.n.trees=GBM.n.trees,
        GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
        GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
        GBMSTEP.step.size=GBMSTEP.step.size,
        RF.formula=NULL, RF.ntree=RF.ntree, 
        GLM.formula=NULL, GLM.family=GLM.family, 
        GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=NULL, GLMSTEP.scope=NULL, 
        GAM.formula=NULL, GAM.family=GAM.family, 
        GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=NULL, GAMSTEP.pos=GAMSTEP.pos,
        MGCV.formula=NULL, MGCV.select=MGCV.select,
        MGCVFIX.formula=NULL, 
        EARTH.formula=NULL, EARTH.glm=EARTH.glm,
        RPART.formula=NULL, RPART.xval=RPART.xval, 
        NNET.formula=NULL, NNET.size=NNET.size, NNET.decay=NNET.decay,
        FDA.formula=NULL, SVM.formula=NULL, SVME.formula=NULL,
        MAHAL.shape=MAHAL.shape)

    if(is.null(tests$evaluations$MAXENT.C)==F) {output.C["MAXENT",i+1] <- tests$evaluations$MAXENT.C@auc}
    if(is.null(tests$evaluations$MAXENT.T)==F) {output.T["MAXENT",i+1] <- tests$evaluations$MAXENT.T@auc}
    if(sum(tests$evaluations$TrainData$MAXENT) > 0) {output.LLC["MAXENT",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAXENT)}
    if(sum(tests$evaluations$TestData$MAXENT) > 0) {output.LLT["MAXENT",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAXENT)}

    if(is.null(tests$evaluations$GBM.C)==F) {output.C["GBM",i+1] <- tests$evaluations$GBM.C@auc} 
    if(is.null(tests$evaluations$GBM.T)==F) {output.T["GBM",i+1] <- tests$evaluations$GBM.T@auc} 
    if(sum(tests$evaluations$TrainData$GBM) > 0) {output.LLC["GBM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GBM)}
    if(sum(tests$evaluations$TestData$GBM) > 0) {output.LLT["GBM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GBM)}

    if(is.null(tests$evaluations$GBMSTEP.C)==F) {output.C["GBMSTEP",i+1] <- tests$evaluations$GBMSTEP.C@auc} 
    if(is.null(tests$evaluations$GBMSTEP.T)==F) {output.T["GBMSTEP",i+1] <- tests$evaluations$GBMSTEP.T@auc} 
    if(sum(tests$evaluations$TrainData$GBMSTEP) > 0) {output.LLC["MAXENT",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAXENT)}
    if(sum(tests$evaluations$TestData$GBMSTEP) > 0) {output.LLT["MAXENT",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAXENT)}

    if(is.null(tests$evaluations$RF.C)==F) {output.C["RF",i+1] <- tests$evaluations$RF.C@auc}
    if(is.null(tests$evaluations$RF.T)==F) {output.T["RF",i+1] <- tests$evaluations$RF.T@auc}
    if(sum(tests$evaluations$TrainData$RF) > 0) {output.LLC["RF",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$RF)}
    if(sum(tests$evaluations$TestData$RF) > 0) {output.LLT["RF",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$RF)}

    if(is.null(tests$evaluations$GLM.C)==F) {output.C["GLM",i+1] <- tests$evaluations$GLM.C@auc} 
    if(is.null(tests$evaluations$GLM.T)==F) {output.T["GLM",i+1] <- tests$evaluations$GLM.T@auc} 
    if(sum(tests$evaluations$TrainData$GLM) > 0) {output.LLC["GLM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GLM)}
    if(sum(tests$evaluations$TestData$GLM) > 0) {output.LLT["GLM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GLM)}

    if(is.null(tests$evaluations$GLMS.C)==F) {output.C["GLMSTEP",i+1] <- tests$evaluations$GLMS.C@auc}
    if(is.null(tests$evaluations$GLMS.T)==F) {output.T["GLMSTEP",i+1] <- tests$evaluations$GLMS.T@auc}
    if(sum(tests$evaluations$TrainData$GLMSTEP) > 0) {output.LLC["GLMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GLMSTEP)}
    if(sum(tests$evaluations$TestData$GLMSTEP) > 0) {output.LLT["GLMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GLMSTEP)}

    if(is.null(tests$evaluations$GAM.C)==F) {output.C["GAM",i+1] <- tests$evaluations$GAM.C@auc} 
    if(is.null(tests$evaluations$GAM.T)==F) {output.T["GAM",i+1] <- tests$evaluations$GAM.T@auc} 
    if(sum(tests$evaluations$TrainData$GAM) > 0) {output.LLC["GAM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GAM)}
    if(sum(tests$evaluations$TestData$GAM) > 0) {output.LLT["GAM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GAM)}

    if(is.null(tests$evaluations$GAMS.C)==F) {output.C["GAMSTEP",i+1] <- tests$evaluations$GAMS.C@auc}
    if(is.null(tests$evaluations$GAMS.T)==F) {output.T["GAMSTEP",i+1] <- tests$evaluations$GAMS.T@auc}
    if(sum(tests$evaluations$TrainData$GAMSTEP) > 0) {output.LLC["GAMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$GAMSTEP)}
    if(sum(tests$evaluations$TestData$GAMSTEP) > 0) {output.LLT["GAMSTEP",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$GAMSTEP)}

    if(is.null(tests$evaluations$MGCV.C)==F) {output.C["MGCV",i+1] <- tests$evaluations$MGCV.C@auc} 
    if(is.null(tests$evaluations$MGCV.T)==F) {output.T["MGCV",i+1] <- tests$evaluations$MGCV.T@auc} 
    if(sum(tests$evaluations$TrainData$MGCV) > 0) {output.LLC["MGCV",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MGCV)}
    if(sum(tests$evaluations$TestData$MGCV) > 0) {output.LLT["MGCV",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MGCV)}

    if(is.null(tests$evaluations$MGCVF.C)==F) {output.C["MGCVFIX",i+1] <- tests$evaluations$MGCVF.C@auc} 
    if(is.null(tests$evaluations$MGCVF.T)==F) {output.T["MGCVFIX",i+1] <- tests$evaluations$MGCVF.T@auc} 
    if(sum(tests$evaluations$TrainData$MGCVFIX) > 0) {output.LLC["MGCVFIX",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MGCVFIX)}
    if(sum(tests$evaluations$TestData$MGCVFIX) > 0) {output.LLT["MGCVFIX",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MGCVFIX)}

    if(is.null(tests$evaluations$EARTH.C)==F) {output.C["EARTH",i+1] <- tests$evaluations$EARTH.C@auc} 
    if(is.null(tests$evaluations$EARTH.T)==F) {output.T["EARTH",i+1] <- tests$evaluations$EARTH.T@auc} 
    if(sum(tests$evaluations$TrainData$EARTH) > 0) {output.LLC["EARTH",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$EARTH)}
    if(sum(tests$evaluations$TestData$EARTH) > 0) {output.LLT["EARTH",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$EARTH)}

    if(is.null(tests$evaluations$RPART.C)==F) {output.C["RPART",i+1] <- tests$evaluations$RPART.C@auc}
    if(is.null(tests$evaluations$RPART.T)==F) {output.T["RPART",i+1] <- tests$evaluations$RPART.T@auc}
    if(sum(tests$evaluations$TrainData$RPART) > 0) {output.LLC["RPART",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$RPART)}
    if(sum(tests$evaluations$TestData$RPART) > 0) {output.LLT["RPART",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$RPART)}

    if(is.null(tests$evaluations$NNET.C)==F) {output.C["NNET",i+1] <- tests$evaluations$NNET.C@auc} 
    if(is.null(tests$evaluations$NNET.T)==F) {output.T["NNET",i+1] <- tests$evaluations$NNET.T@auc} 
    if(sum(tests$evaluations$TrainData$NNET) > 0) {output.LLC["NNET",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$NNET)}
    if(sum(tests$evaluations$TestData$NNET) > 0) {output.LLT["NNET",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$NNET)}

    if(is.null(tests$evaluations$FDA.C)==F) {output.C["FDA",i+1] <- tests$evaluations$FDA.C@auc}
    if(is.null(tests$evaluations$FDA.T)==F) {output.T["FDA",i+1] <- tests$evaluations$FDA.T@auc}
    if(sum(tests$evaluations$TrainData$FDA) > 0) {output.LLC["FDA",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$FDA)}
    if(sum(tests$evaluations$TestData$FDA) > 0) {output.LLT["FDA",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$FDA)}

    if(is.null(tests$evaluations$SVM.C)==F) {output.C["SVM",i+1] <- tests$evaluations$SVM.C@auc}
    if(is.null(tests$evaluations$SVM.T)==F) {output.T["SVM",i+1] <- tests$evaluations$SVM.T@auc}
    if(sum(tests$evaluations$TrainData$SVM) > 0) {output.LLC["SVM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$SVM)}
    if(sum(tests$evaluations$TestData$SVM) > 0) {output.LLT["SVM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$SVM)}

    if(is.null(tests$evaluations$SVME.C)==F) {output.C["SVME",i+1] <- tests$evaluations$SVME.C@auc}
    if(is.null(tests$evaluations$SVME.T)==F) {output.T["SVME",i+1] <- tests$evaluations$SVME.T@auc}
    if(sum(tests$evaluations$TrainData$SVME) > 0) {output.LLC["SVME",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$SVME)}
    if(sum(tests$evaluations$TestData$SVME) > 0) {output.LLT["SVME",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$SVME)}

    if(is.null(tests$evaluations$BIOCLIM.C)==F) {output.C["BIOCLIM",i+1] <- tests$evaluations$BIOCLIM.C@auc}
    if(is.null(tests$evaluations$BIOCLIM.T)==F) {output.T["BIOCLIM",i+1] <- tests$evaluations$BIOCLIM.T@auc}
    if(sum(tests$evaluations$TrainData$BIOCLIM) > 0) {output.LLC["BIOCLIM",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$BIOCLIM)}
    if(sum(tests$evaluations$TestData$BIOCLIM) > 0) {output.LLT["BIOCLIM",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$BIOCLIM)}

    if(is.null(tests$evaluations$DOMAIN.C)==F) {output.C["DOMAIN",i+1] <- tests$evaluations$DOMAIN.C@auc}
    if(is.null(tests$evaluations$DOMAIN.T)==F) {output.T["DOMAIN",i+1] <- tests$evaluations$DOMAIN.T@auc}
    if(sum(tests$evaluations$TrainData$DOMAIN) > 0) {output.LLC["DOMAIN",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$DOMAIN)}
    if(sum(tests$evaluations$TestData$DOMAIN) > 0) {output.LLT["DOMAIN",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$DOMAIN)}

    if(is.null(tests$evaluations$MAHAL.C)==F) {output.C["MAHAL",i+1] <- tests$evaluations$MAHAL.C@auc}
    if(is.null(tests$evaluations$MAHAL.T)==F) {output.T["MAHAL",i+1] <- tests$evaluations$MAHAL.T@auc}
    if(sum(tests$evaluations$TrainData$MAHAL) > 0) {output.LLC["MAHAL",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$MAHAL)}
    if(sum(tests$evaluations$TestData$MAHAL) > 0) {output.LLT["MAHAL",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$MAHAL)}

    if(is.null(tests$evaluations$ENSEMBLE.C)==F) {output.C["ENSEMBLE",i+1] <- tests$evaluations$ENSEMBLE.C@auc}
    if(is.null(tests$evaluations$ENSEMBLE.T)==F) {output.T["ENSEMBLE",i+1] <- tests$evaluations$ENSEMBLE.T@auc}
    if(sum(tests$evaluations$TrainData$ENSEMBLE) > 0) {output.LLC["ENSEMBLE",i+1] <- loglik.calculation(obs=tests$evaluations$TrainData$pb, preds=tests$evaluations$TrainData$ENSEMBLE)}
    if(sum(tests$evaluations$TestData$ENSEMBLE) > 0) {output.LLT["ENSEMBLE",i+1] <- loglik.calculation(obs=tests$evaluations$TestData$pb, preds=tests$evaluations$TestData$ENSEMBLE)}

    }

## possibly later: finally test with only the focal variable


    output.C <- 100*output.C
    if (difference == T) {
        for (i in 1:nv) {
            output.C[,i+1] <- output.C[,i+1] - output.C[,1]
        }
    }
    output.C <- output.C[order(output.C[,"all_vars"], decreasing=T),]
    cat(paste("\n", "AUC for testing data (as percentage)",  "\n\n", sep = ""))
    if (difference == T) {
        cat(paste("\n", "Note that positive differences indicate that the model without the variable", sep = ""))
        cat(paste("\n", "has higher AUC than the model with all the variables",  "\n\n", sep = ""))
    }
    print (output.C)
    cat(paste("\n\n"))

    output.T <- 100*output.T
    if (difference == T) {
        for (i in 1:nv) {
            output.T[,i+1] <- output.T[,i+1] - output.T[,1]
        }
    }
    output.T <- output.T[order(output.T[,"all_vars"], decreasing=T),]
    cat(paste("\n", "AUC for testing data (as percentage)",  "\n\n", sep = ""))
    if (difference == T) {
        cat(paste("\n", "Note that positive differences indicate that the model without the variable", sep = ""))
        cat(paste("\n", "has higher AUC than the model with all the variables",  "\n\n", sep = ""))
    }
    print (output.T)
    cat(paste("\n\n"))
#
# base null model on predictions of prevalence
    obs1 <- tests$evaluations$TrainData$pb
    prevalence <- sum(obs1) / length(obs1)
    preds1 <- obs1
    preds1[] <- prevalence
    null.dev <- loglik.calculation(obs=obs1, preds=preds1)
    cat(paste("\n", "Null deviance for calibration data (predicting prevalence): ", null.dev, "\n\n", sep=""))
#
    cat(paste("\n", "Residual deviance for calibration data",  "\n\n", sep = ""))
    output.LLC <- output.LLC[order(output.LLC[,"all_vars"], decreasing=F),]
    print (output.LLC)
    cat(paste("\n\n"))
#
    obs1 <- tests$evaluations$TestData$pb
    prevalence <- sum(obs1) / length(obs1)
    preds1 <- obs1
    preds1[] <- prevalence
    null.dev <- loglik.calculation(obs=obs1, preds=preds1)
    cat(paste("\n", "Null deviance for testing data (predicting prevalence): ", null.dev, "\n\n", sep=""))
#
    cat(paste("\n", "Residual deviance for testing data",  "\n\n", sep = ""))
    output.LLT <- output.LLT[order(output.LLT[,"all_vars"], decreasing=F),]
    print (output.LLT)
    cat(paste("\n\n"))

    if (SINK==T && OLD.SINK==F) {sink(file=NULL, append=T)}
    return(list(AUC.calibration=output.C, AUC.testing=output.T, 
        residual.deviance.calibration=output.LLC, residual.deviance.testing=output.LLT, call=match.call() ))
}

