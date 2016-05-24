`ensemble.strategy` <- function(
    TrainData=NULL, TestData=NULL,
    verbose=FALSE,
    ENSEMBLE.best=c(4:10), ENSEMBLE.min=c(0.7),
    ENSEMBLE.exponent=c(1, 2, 4, 6, 8) 
)
{
#    if (! require(dismo)) {stop("Please install the dismo package")}
#   input AUC
    modelnames <- c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL")
    weights <- numeric(length=length(modelnames))
    final.weights <- weights
    names(weights) <- modelnames
    bests <- length(ENSEMBLE.best)
    exponents <- length(ENSEMBLE.exponent)
    mins <- length(ENSEMBLE.min) 
#
# output for each cross-validation run
    output <- data.frame(array(dim=c(bests*exponents*mins, 7), NA))
    if (nrow(output) == 1) {cat(paste("\n", "NOTE: no alternatives available for choosing best strategy", "\n", sep=""))}
    colnames(output) <- c("ENSEMBLE.best", "ENSEMBLE.exponent", "ENSEMBLE.min", "model.C", "AUC.C", "model.T", "AUC.T")
    all.combinations <- expand.grid(ENSEMBLE.best, ENSEMBLE.exponent, ENSEMBLE.min)
    output[,c(1:3)] <- all.combinations
#
# recalculate AUC values 
    weights.cal <- c(0, weights)
    weights.eval <- c(0, weights)
    names(weights.cal) <- c("ENSEMBLE", modelnames)
    names(weights.eval) <- c("ENSEMBLE", modelnames)
    for (i in 1:length(weights)) {
        TrainPres <- TrainData[TrainData[,"pb"]==1, modelnames[i]]
        TrainAbs <- TrainData[TrainData[,"pb"]==0, modelnames[i]]
        if (sum(TrainPres, TrainAbs, na.rm=T) != 0) {
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            weights.cal[i+1] <- eval1@auc
        }
        TestPres <- TestData[TestData[,"pb"]==1, modelnames[i]]
        TestAbs <- TestData[TestData[,"pb"]==0, modelnames[i]]
        if (sum(TestPres, TestAbs, na.rm=T) != 0) {
            eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
            weights.eval[i+1] <- eval2@auc
        }
    }
    input.weights.c <- weights.cal[names(weights.cal) != "ENSEMBLE"]
    input.weights.e <- weights.eval[names(weights.eval) != "ENSEMBLE"]

    auc.target <- -1.0
    for (r in 1:nrow(output)) {
        if (verbose == T) {
            cat(paste("\n", "run ", r, ": best=", output[r, "ENSEMBLE.best"], ", exponent=", output[r, "ENSEMBLE.exponent"],  ", min=", output[r, "ENSEMBLE.min"] ,"\n", sep=""))
        }
#
# strategy based on evaluations
        ws <- ensemble.weights(input.weights.e, exponent=output[r, "ENSEMBLE.exponent"], best=output[r, "ENSEMBLE.best"], 
            min.weight=output[r, "ENSEMBLE.min"])
        if (verbose == T) {print(ws)}
        TrainData[,"ENSEMBLE"] <- ws["MAXENT"]*TrainData[,"MAXENT"] + ws["GBM"]*TrainData[,"GBM"] +
            ws["GBMSTEP"]*TrainData[,"GBMSTEP"] + ws["RF"]*TrainData[,"RF"] + ws["GLM"]*TrainData[,"GLM"] +
            ws["GLMSTEP"]*TrainData[,"GLMSTEP"] + ws["GAM"]*TrainData[,"GAM"] + ws["GAMSTEP"]*TrainData[,"GAMSTEP"] +
            ws["MGCV"]*TrainData[,"MGCV"] + ws["MGCVFIX"]*TrainData[,"MGCVFIX"] + ws["EARTH"]*TrainData[,"EARTH"] +
            ws["RPART"]*TrainData[,"RPART"] + ws["NNET"]*TrainData[,"NNET"] + ws["FDA"]*TrainData[,"FDA"] +
            ws["SVM"]*TrainData[,"SVM"] + ws["SVME"]*TrainData[,"SVME"] + ws["BIOCLIM"]*TrainData[,"BIOCLIM"] +
            ws["DOMAIN"]*TrainData[,"DOMAIN"] + ws["MAHAL"]*TrainData[,"MAHAL"]
#        TrainData[,"ENSEMBLE"] <- trunc(TrainData[,"ENSEMBLE"])
        TrainPres <- TrainData[TrainData[,"pb"]==1,"ENSEMBLE"]
        TrainAbs <- TrainData[TrainData[,"pb"]==0,"ENSEMBLE"]
        eval1 <- NULL
        if (sum(TrainPres, TrainAbs, na.rm=T) != 0) {
            eval1 <- dismo::evaluate(p=TrainPres, a=TrainAbs)
            if (verbose == T) {
                cat(paste("\n", "evaluation with train data", "\n", sep=""))
                print(eval1)
            }
            weights.cal["ENSEMBLE"] <- eval1@auc
            weights.cal <- weights.cal[order(weights.cal, decreasing=T)]
            output[r, "model.C"] <- names(weights.cal)[1]
            output[r, "AUC.C"] <- weights.cal[1]
        }
        TestData[,"ENSEMBLE"] <- ws["MAXENT"]*TestData[,"MAXENT"] + ws["GBM"]*TestData[,"GBM"] +
            ws["GBMSTEP"]*TestData[,"GBMSTEP"] + ws["RF"]*TestData[,"RF"] + ws["GLM"]*TestData[,"GLM"] +
            ws["GLMSTEP"]*TestData[,"GLMSTEP"] + ws["GAM"]*TestData[,"GAM"] + ws["GAMSTEP"]*TestData[,"GAMSTEP"] +
            ws["MGCV"]*TestData[,"MGCV"] + ws["MGCVFIX"]*TestData[,"MGCVFIX"] + ws["EARTH"]*TestData[,"EARTH"] +
            ws["RPART"]*TestData[,"RPART"] + ws["NNET"]*TestData[,"NNET"] + ws["FDA"]*TestData[,"FDA"] +
            ws["SVM"]*TestData[,"SVM"] + ws["SVME"]*TestData[,"SVME"] + ws["BIOCLIM"]*TestData[,"BIOCLIM"] +
            ws["DOMAIN"]*TestData[,"DOMAIN"] + ws["MAHAL"]*TestData[,"MAHAL"]
#        TestData[,"ENSEMBLE"] <- trunc(TestData[,"ENSEMBLE"])
        TestPres <- TestData[TestData[,"pb"]==1,"ENSEMBLE"]
        TestAbs <- TestData[TestData[,"pb"]==0,"ENSEMBLE"]
        eval2 <- NULL
        if (sum(TestPres, TestAbs, na.rm=T) != 0) {
            eval2 <- dismo::evaluate(p=TestPres, a=TestAbs)
            if (verbose == T) {
                cat(paste("\n", "evaluation with test data", "\n", sep=""))
                print(eval2)
            }
            weights.eval["ENSEMBLE"] <- eval2@auc
            weights.eval <- weights.eval[order(weights.eval, decreasing=T)]
            output[r, "model.T"] <- names(weights.eval)[1]
            output[r, "AUC.T"] <- weights.eval[1]
        }
        if (weights.eval[1] > auc.target) {
            auc.target <- weights.eval[1]
            weights.out <- ws
        }
    }
    output <- output[order(output[,"AUC.T"], decreasing=T), ]
    cat(paste("\n", "Ensemble tuning: best=", output[1, "ENSEMBLE.best"], ", exponent=", output[1, "ENSEMBLE.exponent"],  ", min=", output[1, "ENSEMBLE.min"] ,"\n", sep=""))
    cat(paste("\n", "Weights used for best strategy", "\n", sep = ""))
    print(weights.out)
    return(list(weights=weights.out, output=output))
}


