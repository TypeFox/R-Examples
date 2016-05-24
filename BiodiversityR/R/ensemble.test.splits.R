`ensemble.test.splits` <- function(
    x=NULL, p=NULL, a=NULL, an=1000, CIRCLES.at=FALSE, CIRCLES.d=100000,
    excludep=FALSE, ext=NULL, k=4, 
    TrainData=NULL,
    VIF=FALSE, COR=FALSE,
    SINK=FALSE, PLOTS=FALSE, 
    data.keep=FALSE,
    species.name = "Species001",
    threshold.method="spec_sens", threshold.sensitivity=0.9, threshold.PresenceAbsence=FALSE,
    AUC.weights=TRUE, ENSEMBLE.tune=FALSE, 
    ENSEMBLE.best=0, ENSEMBLE.min=0.7, ENSEMBLE.exponent=1, 
    input.weights=NULL,
    MAXENT=1, GBM=1, GBMSTEP=1, RF=1, GLM=1, GLMSTEP=1, GAM=1, GAMSTEP=1, MGCV=1, MGCVFIX=0, 
    EARTH=1, RPART=1, NNET=1, FDA=1, SVM=1, SVME=1, BIOCLIM=1, DOMAIN=1, MAHAL=1, 
    PROBIT=FALSE,
    Yweights="BIOMOD", 
    layer.drops=NULL, factors=NULL, dummy.vars=NULL,
    formulae.defaults=TRUE, maxit=100,
    MAXENT.a=NULL, MAXENT.an=10000, MAXENT.BackData=NULL, 
    MAXENT.path=paste(getwd(), "/models/maxent_", species.name,  sep=""), 
    GBM.formula=NULL, GBM.n.trees=2001, 
    GBMSTEP.gbm.x=2:(ncol(TrainData1)), GBMSTEP.tree.complexity=5, GBMSTEP.learning.rate=0.005, 
    GBMSTEP.bag.fraction=0.5, GBMSTEP.step.size=100, 
    RF.formula=NULL, RF.ntree=751, RF.mtry=floor(sqrt(ncol(TrainData1)-1)), 
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
    MAHAL.shape=1 
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    k <- as.integer(k)
    if (k < 2) {
        cat(paste("\n", "NOTE: parameter k was set to be smaller than 2", sep = ""))
        cat(paste("\n", "default value of 4 therefore set for parameter k", "\n", sep = ""))
        k <- 4
    }
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
#
    output.rownames <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL", "ENSEMBLE")
    if(length(ENSEMBLE.exponent) > 1 || length(ENSEMBLE.best) > 1 || length(ENSEMBLE.min) > 1) {ENSEMBLE.tune <- TRUE}
    if(ENSEMBLE.tune == F) {
        output <- array(0, dim=c(length(output.rownames), k+1))
        rownames(output) <- output.rownames
        colnames(output) <- c(paste("T_", c(1:k), sep=""),"MEAN")
    }else{
        output <- array(0, dim=c(length(output.rownames), 2*k+2))
        rownames(output) <- output.rownames
        colnames(output) <- c(paste("T_", c(1:k), sep=""), "MEAN.T", paste("S_", c(1:k), sep=""), "MEAN")
    }

# keep data for final checks with suggested weights
    TestData.all <- vector("list", k)

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
        cat(paste("\n\n", "RESULTS (ensemble.test.splits function)", "\n", sep=""), file=paste.file, append=T)
        sink(file=paste.file, append=T)
        cat(paste(date(), "\n", sep=""))
        print(match.call())
    }

# 
# run ensemble.test first to obtain MAXENT.BackData and var.names
    tests <- ensemble.test(x=x, ext=ext,
        p=p, a=a, an=an, pt=NULL, at=NULL, excludep=excludep, k=0, 
        TrainData=TrainData, 
        VIF=F, COR=F,
        PLOTS=PLOTS, evaluations.keep=T, models.keep=F,
        AUC.weights=F, ENSEMBLE.tune=F,
        ENSEMBLE.exponent=1, ENSEMBLE.best=1, ENSEMBLE.min=0.7, 
        MAXENT=0, GBM=0, GBMSTEP=0, RF=0, GLM=0, GLMSTEP=0, 
        GAM=0, GAMSTEP=0, MGCV=0, MGCVFIX=0, EARTH=0, RPART=0, 
        NNET=0, FDA=0, SVM=0, SVME=0, BIOCLIM=0, DOMAIN=0, MAHAL=0,
        MAXENT.a=MAXENT.a, MAXENT.an=MAXENT.an, MAXENT.BackData=MAXENT.BackData,
        GEODIST=0,
        factors=factors)

    var.names <- tests$evaluations$var.names
    var.names2 <- c("pb", var.names)
    TrainData1 <- tests$evaluations$TrainData
    TrainData1 <- TrainData1[, which(names(TrainData1) %in% var.names2)]
    MAXENT.BackData1 <- tests$evaluations$MAXENT.BackData
    MAXENT.BackData2 <- NULL
    if (is.null(MAXENT.BackData1) == F) {MAXENT.BackData2 <- MAXENT.BackData1[, which(names(MAXENT.BackData1) %in% var.names)]}
    factors2 <- NULL
    if (is.null(factors) == F) {
        factors2 <- factors[which(factors %in% var.names)]
        if (length(factors2) == 0) {factors2 <- NULL}
    }
    dummy.vars2 <- NULL
    if (is.null(dummy.vars) == F) {
        dummy.vars2 <- dummy.vars[which(dummy.vars %in% var.names)]
        if (length(dummy.vars2) == 0) {dummy.vars2 <- NULL}
    }

# Different cross-validations if categorical variables or circular neighbourhoods (locations needed)
# If locations not needed, then process considerably faster

    if (length(factors2)>0  || CIRCLES.at==T) {
        p.all <- tests$evaluations$p
        a.all <- tests$evaluations$a
        groupp <- dismo::kfold(p.all, k=k)
        groupa <- dismo::kfold(a.all, k=k)
    }else{
        groupd <- dismo::kfold(TrainData1, k=k, by=TrainData1[,"pb"])
    }

# Start cross-validations
    
    for (i in 1:k){
        cat(paste(species.name, " K-FOLD CROSS-VALIDATION RUN: ", i, "\n", sep = ""))

        if (length(factors2)>0  || CIRCLES.at==T) {
            p1 <- p.all[groupp != i,]
            p2 <- p.all[groupp == i,]
            a1 <- a.all[groupa != i,]
            a2 <- a.all[groupa == i,]

            tests <- ensemble.test(x=x, ext=ext,
                TrainData=NULL, TestData=NULL,
                p=p1, a=a1, pt=p2, at=a2, CIRCLES.at=CIRCLES.at, CIRCLES.d=CIRCLES.d,
                VIF=VIF, COR=COR,
                threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence,
                PLOTS=PLOTS, evaluations.keep=T, models.keep=F,
                AUC.weights=AUC.weights, ENSEMBLE.tune=ENSEMBLE.tune,
                ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min, ENSEMBLE.exponent=ENSEMBLE.exponent, 
                MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
                GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
                NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,
                GEODIST=0, 
                PROBIT=PROBIT,  
                Yweights=Yweights, 
                factors=factors2, dummy.vars=dummy.vars2,
                maxit=maxit,
                MAXENT.BackData=MAXENT.BackData2, MAXENT.path=MAXENT.path,
                GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees, 
                GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
                GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
                GBMSTEP.step.size=GBMSTEP.step.size, 
                RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
                GLM.formula=GLM.formula, GLM.family=GLM.family, 
                GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=STEP.formula, 
                GLMSTEP.scope=GLMSTEP.scope, 
                GAM.formula=GAM.formula, GAM.family=GAM.family, 
                GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope, GAMSTEP.pos=GAMSTEP.pos, 
                MGCV.formula=MGCV.formula, MGCV.select=MGCV.select, 
                MGCVFIX.formula=MGCVFIX.formula, 
                EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm, 
                RPART.formula=RPART.formula, RPART.xval=RPART.xval, 
                NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay,
                FDA.formula=FDA.formula, 
                SVM.formula=SVM.formula, 
                SVME.formula=SVME.formula, 
                MAHAL.shape=MAHAL.shape)

# Different cross-validations if no categorical variables and no circular neighbourhoods

        }else{
            TrainData2 <- TrainData1[groupd != i,]
            TestData2 <- TrainData1[groupd == i,]

            tests <- ensemble.test(x=x, ext=ext, p=NULL, a=NULL, pt=NULL, at=NULL,
                TrainData=TrainData2, TestData=TestData2, CIRCLES.at=CIRCLES.at, CIRCLES.d=CIRCLES.d,
                VIF=VIF, COR=COR,
                threshold.method=threshold.method, threshold.sensitivity=threshold.sensitivity, threshold.PresenceAbsence=threshold.PresenceAbsence,
                PLOTS=PLOTS, evaluations.keep=T, models.keep=F,
                AUC.weights=AUC.weights, ENSEMBLE.tune=ENSEMBLE.tune,
                ENSEMBLE.best=ENSEMBLE.best, ENSEMBLE.min=ENSEMBLE.min, ENSEMBLE.exponent=ENSEMBLE.exponent,
                MAXENT=MAXENT, GBM=GBM, GBMSTEP=GBMSTEP, RF=RF, GLM=GLM, GLMSTEP=GLMSTEP, 
                GAM=GAM, GAMSTEP=GAMSTEP, MGCV=MGCV, MGCVFIX=MGCVFIX, EARTH=EARTH, RPART=RPART, 
                NNET=NNET, FDA=FDA, SVM=SVM, SVME=SVME, BIOCLIM=BIOCLIM, DOMAIN=DOMAIN, MAHAL=MAHAL,
                GEODIST=0, 
                PROBIT=PROBIT,  
                Yweights=Yweights, 
                factors=factors2, dummy.vars=dummy.vars2,
                maxit=maxit,
                MAXENT.BackData=MAXENT.BackData2, MAXENT.path=MAXENT.path,
                GBM.formula=GBM.formula, GBM.n.trees=GBM.n.trees, 
                GBMSTEP.gbm.x=GBMSTEP.gbm.x, GBMSTEP.tree.complexity=GBMSTEP.tree.complexity, 
                GBMSTEP.learning.rate=GBMSTEP.learning.rate, GBMSTEP.bag.fraction=GBMSTEP.bag.fraction,
                GBMSTEP.step.size=GBMSTEP.step.size, 
                RF.formula=RF.formula, RF.ntree=RF.ntree, RF.mtry=RF.mtry, 
                GLM.formula=GLM.formula, GLM.family=GLM.family, 
                GLMSTEP.k=GLMSTEP.k, GLMSTEP.steps=GLMSTEP.steps, STEP.formula=STEP.formula, 
                GLMSTEP.scope=GLMSTEP.scope, 
                GAM.formula=GAM.formula, GAM.family=GAM.family, 
                GAMSTEP.steps=GAMSTEP.steps, GAMSTEP.scope=GAMSTEP.scope, GAMSTEP.pos=GAMSTEP.pos, 
                MGCV.formula=MGCV.formula, MGCV.select=MGCV.select, 
                MGCVFIX.formula=MGCVFIX.formula, 
                EARTH.formula=EARTH.formula, EARTH.glm=EARTH.glm, 
                RPART.formula=RPART.formula, RPART.xval=RPART.xval, 
                NNET.formula=NNET.formula, NNET.size=NNET.size, NNET.decay=NNET.decay,
                FDA.formula=FDA.formula, 
                SVM.formula=SVM.formula, 
                SVME.formula=SVME.formula, 
                MAHAL.shape=MAHAL.shape)

        }

        if(is.null(tests$evaluations$MAXENT.T)==F) {output["MAXENT",i] <- tests$evaluations$MAXENT.T@auc}
        if(is.null(tests$evaluations$GBM.T)==F) {output["GBM",i] <- tests$evaluations$GBM.T@auc} 
        if(is.null(tests$evaluations$GBMSTEP.T)==F) {output["GBMSTEP",i] <- tests$evaluations$GBMSTEP.T@auc} 
        if(is.null(tests$evaluations$RF.T)==F) {output["RF",i] <- tests$evaluations$RF.T@auc}
        if(is.null(tests$evaluations$GLM.T)==F) {output["GLM",i] <- tests$evaluations$GLM.T@auc} 
        if(is.null(tests$evaluations$GLMS.T)==F) {output["GLMSTEP",i] <- tests$evaluations$GLMS.T@auc}
        if(is.null(tests$evaluations$GAM.T)==F) {output["GAM",i] <- tests$evaluations$GAM.T@auc} 
        if(is.null(tests$evaluations$GAMS.T)==F) {output["GAMSTEP",i] <- tests$evaluations$GAMS.T@auc}
        if(is.null(tests$evaluations$MGCV.T)==F) {output["MGCV",i] <- tests$evaluations$MGCV.T@auc} 
        if(is.null(tests$evaluations$MGCVF.T)==F) {output["MGCVFIX",i] <- tests$evaluations$MGCVF.T@auc} 
        if(is.null(tests$evaluations$EARTH.T)==F) {output["EARTH",i] <- tests$evaluations$EARTH.T@auc} 
        if(is.null(tests$evaluations$RPART.T)==F) {output["RPART",i] <- tests$evaluations$RPART.T@auc}
        if(is.null(tests$evaluations$NNET.T)==F) {output["NNET",i] <- tests$evaluations$NNET.T@auc} 
        if(is.null(tests$evaluations$FDA.T)==F) {output["FDA",i] <- tests$evaluations$FDA.T@auc}
        if(is.null(tests$evaluations$SVM.T)==F) {output["SVM",i] <- tests$evaluations$SVM.T@auc}
        if(is.null(tests$evaluations$SVME.T)==F) {output["SVME",i] <- tests$evaluations$SVME.T@auc}
        if(is.null(tests$evaluations$BIOCLIM.T)==F) {output["BIOCLIM",i] <- tests$evaluations$BIOCLIM.T@auc}
        if(is.null(tests$evaluations$DOMAIN.T)==F) {output["DOMAIN",i] <- tests$evaluations$DOMAIN.T@auc}
        if(is.null(tests$evaluations$MAHAL.T)==F) {output["MAHAL",i] <- tests$evaluations$MAHAL.T@auc}
        if(is.null(tests$evaluations$ENSEMBLE.T)==F) {output["ENSEMBLE",i] <- tests$evaluations$ENSEMBLE.T@auc}

        if(ENSEMBLE.tune == T) {
            output["MAXENT",k+1+i] <- tests$evaluations$STRATEGY.weights["MAXENT"]
            output["GBM",k+1+i] <- tests$evaluations$STRATEGY.weights["GBM"]
            output["GBMSTEP",k+1+i] <- tests$evaluations$STRATEGY.weights["GBMSTEP"]
            output["RF",k+1+i] <- tests$evaluations$STRATEGY.weights["RF"]
            output["GLM",k+1+i] <- tests$evaluations$STRATEGY.weights["GLM"]
            output["GLMSTEP",k+1+i] <- tests$evaluations$STRATEGY.weights["GLMSTEP"]
            output["GAM",k+1+i] <- tests$evaluations$STRATEGY.weights["GAM"]
            output["GAMSTEP",k+1+i] <- tests$evaluations$STRATEGY.weights["GAMSTEP"]
            output["MGCV",k+1+i] <- tests$evaluations$STRATEGY.weights["MGCV"]
            output["MGCVFIX",k+1+i] <- tests$evaluations$STRATEGY.weights["MGCVFIX"]
            output["EARTH",k+1+i] <- tests$evaluations$STRATEGY.weights["EARTH"]
            output["RPART",k+1+i] <- tests$evaluations$STRATEGY.weights["RPART"]
            output["NNET",k+1+i] <- tests$evaluations$STRATEGY.weights["NNET"]
            output["FDA",k+1+i] <- tests$evaluations$STRATEGY.weights["FDA"]
            output["SVM",k+1+i] <- tests$evaluations$STRATEGY.weights["SVM"]
            output["SVME",k+1+i] <- tests$evaluations$STRATEGY.weights["SVME"]
            output["BIOCLIM",k+1+i] <- tests$evaluations$STRATEGY.weights["BIOCLIM"]
            output["DOMAIN",k+1+i] <- tests$evaluations$STRATEGY.weights["DOMAIN"]
            output["MAHAL",k+1+i] <- tests$evaluations$STRATEGY.weights["MAHAL"]
        }

        TestData.all[[i]] <- tests$evaluations$TestData
    }

    output[,k+1] <- rowMeans(output[,c(1:k)], na.rm=T)
    output[is.na(output[,k+1]),(k+1)] <- 0
    if(ENSEMBLE.tune == T) {
        output[,2*k+2] <- rowMeans(output[,c((k+2):(2*k+1))], na.rm=T)
        output[is.na(output[,2*k+2]),(2*k+2)] <- 0
        output.weights.T <- output[,"MEAN.T"]
        output.weights.T <- output.weights.T[names(output.weights.T) != "ENSEMBLE"]
    }
#
    output.weights <- output[,"MEAN"]
    output.weights <- output.weights[names(output.weights) != "ENSEMBLE"]
    output.weights <- ensemble.weights(output.weights, exponent=1, best=0, min.weight=0)
    output <- output[order(output[,k+1], decreasing=T),]
    cat(paste("\n", "Results of ensemble.test.splits sorted by average AUC for tests T_1 to T_", k, "\n", sep = ""))
    if(ENSEMBLE.tune == T) {
        cat(paste("S_1 to S_", k, " show the weights for the ensemble model with best AUC", "\n", sep = ""))
        cat(paste("column MEAN shows the mean of these weights", "\n\n", sep = ""))
    }
    print(output)
#
    if(ENSEMBLE.tune == T) {
        cat(paste("\n", "suggested input weights for ensemble modelling (based on MEAN.T column)",  sep = ""))
        cat(paste("\n", "ENSEMBLE.exponent=2; ENSEMBLE.best=0; ENSEMBLE.min=", ENSEMBLE.min, "\n\n", sep = ""))
        output.weights.T <- ensemble.weights(weights=output.weights.T, exponent=2, best=0, min.weight=ENSEMBLE.min)
        print(output.weights.T)
    # test with suggested weights
        output3 <- numeric(length=k+1)
        names(output3)[1:k] <- paste("T_", c(1:k), sep="")
        names(output3)[k+1] <- c("MEAN.T")
        ws <- output.weights.T
        for (i in 1:k) {
            TestData <- TestData.all[[i]]
            TestData[,"ENSEMBLE"] <- ws["MAXENT"]*TestData[,"MAXENT"] + ws["GBM"]*TestData[,"GBM"] +
                ws["GBMSTEP"]*TestData[,"GBMSTEP"] + ws["RF"]*TestData[,"RF"] + ws["GLM"]*TestData[,"GLM"] +
                ws["GLMSTEP"]*TestData[,"GLMSTEP"] + ws["GAM"]*TestData[,"GAM"] + ws["GAMSTEP"]*TestData[,"GAMSTEP"] +
                ws["MGCV"]*TestData[,"MGCV"] + ws["MGCVFIX"]*TestData[,"MGCVFIX"] + ws["EARTH"]*TestData[,"EARTH"] +
                ws["RPART"]*TestData[,"RPART"] + ws["NNET"]*TestData[,"NNET"] + ws["FDA"]*TestData[,"FDA"] +
                ws["SVM"]*TestData[,"SVM"] + ws["SVME"]*TestData[,"SVME"] + ws["BIOCLIM"]*TestData[,"BIOCLIM"] +
                ws["DOMAIN"]*TestData[,"DOMAIN"] + ws["MAHAL"]*TestData[,"MAHAL"]
            eval1 <- eval2 <- NULL
            TestPres <- as.numeric(TestData[TestData[,"pb"]==1,"ENSEMBLE"])
            TestAbs <- as.numeric(TestData[TestData[,"pb"]==0,"ENSEMBLE"])
            eval1 <- dismo::evaluate(p=TestPres, a=TestAbs)
            output3[i] <- eval1@auc
        }
        output3[k+1] <- mean(output3[1:k])
        cat(paste("\n", "AUC for ensemble models based on suggested input weights",  "\n\n", sep = ""))
        print(output3)
    }
#
    cat(paste("\n", "suggested input weights for ensemble modelling (based on MEAN column)",  "\n\n", sep = ""))
    print(output.weights)
    cat(paste("\n", "Minimum input weight is 0.05", "\n", sep=""))
    output.weights[output.weights < 0.05] <- 0
    output.weights <- ensemble.weights(weights=output.weights, exponent=1, best=0, min.weight=0)
    cat(paste("\n", "Weights for ensemble forecasting", "\n", sep = ""))
    print(output.weights)

# test with suggested weights
    output2 <- numeric(length=k+1)
    names(output2)[1:k] <- paste("T_", c(1:k), sep="")
    names(output2)[k+1] <- c("MEAN.T")
    ws <- output.weights
    for (i in 1:k) {
        TestData <- TestData.all[[i]]
        TestData[,"ENSEMBLE"] <- ws["MAXENT"]*TestData[,"MAXENT"] + ws["GBM"]*TestData[,"GBM"] +
            ws["GBMSTEP"]*TestData[,"GBMSTEP"] + ws["RF"]*TestData[,"RF"] + ws["GLM"]*TestData[,"GLM"] +
            ws["GLMSTEP"]*TestData[,"GLMSTEP"] + ws["GAM"]*TestData[,"GAM"] + ws["GAMSTEP"]*TestData[,"GAMSTEP"] +
            ws["MGCV"]*TestData[,"MGCV"] + ws["MGCVFIX"]*TestData[,"MGCVFIX"] + ws["EARTH"]*TestData[,"EARTH"] +
            ws["RPART"]*TestData[,"RPART"] + ws["NNET"]*TestData[,"NNET"] + ws["FDA"]*TestData[,"FDA"] +
            ws["SVM"]*TestData[,"SVM"] + ws["SVME"]*TestData[,"SVME"] + ws["BIOCLIM"]*TestData[,"BIOCLIM"] +
            ws["DOMAIN"]*TestData[,"DOMAIN"] + ws["MAHAL"]*TestData[,"MAHAL"]
        eval1 <- eval2 <- NULL
        TestPres <- as.numeric(TestData[TestData[,"pb"]==1,"ENSEMBLE"])
        TestAbs <- as.numeric(TestData[TestData[,"pb"]==0,"ENSEMBLE"])
        eval1 <- dismo::evaluate(p=TestPres, a=TestAbs)
        output2[i] <- eval1@auc
    }
    output2[k+1] <- mean(output2[1:k])
    cat(paste("\n", "AUC for ensemble models based on suggested input weights",  "\n\n", sep = ""))
    print(output2)
    if(ENSEMBLE.tune == F) {
        output.weights.T <- output.weights
        output3 <- output2
    }
    if (SINK==T && OLD.SINK==F) {sink(file=NULL, append=T)}
    if (data.keep == F) {
        cat(paste("\n\n"))
        return(list(table=output, output.weights=output.weights, AUC.with.suggested.weights=output2, 
            output.weights.AUC=output.weights.T, AUC.with.suggested.weights2=output3, 
            call=match.call()))
    }else{
        cat(paste("\n", "(note that output data are integer values representing probabilities multiplied by 1000)",  "\n\n", sep = ""))
        return(list(table=output, output.weights=output.weights, AUC.with.suggested.weights=output2, 
            output.weights.AUC=output.weights.T, AUC.with.suggested.weights.T=output3, 
            data=TestData.all, call=match.call()))
    }
}

