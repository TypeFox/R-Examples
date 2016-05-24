#############################################
### Standard Nonconformity Measure
#############################################
StandardMeasure <- function(obs,pred,error){
    if(!is.vector(obs) || !is.vector(pred) || !is.vector(error))
    stop("The function requires three numerical vectors as input")
    alpha <- abs(obs-pred) / error
    return(sort(alpha))
}

#############################################
### Extract Cross-validation Predictions from
### a caret model
#############################################
GetCVPreds <- function(model) {
    pred <- model$pred
    if (is.null(pred))
        stop("You must provide a model traind with cross-validatio. Cross-validation predictions are required to train the error model")
    bestTune <- model$bestTune
    for (name in names(bestTune)) {
        pred <- pred[pred[, name] == bestTune[, name], ]
    }
    pred <- pred[order(pred$rowIndex), ]
    return(pred)
}


#############################################
### Generate an Error Model
#############################################
ErrorModel <- function(PointPredictionModel,x.train,algorithm="svmRadial",...){
    predObsCV <- GetCVPreds(PointPredictionModel)
    error_model <- train(x.train, abs(predObsCV$pred - predObsCV$obs), algorithm, ...)
    return(error_model)
}


#############################################
### Create an Exponential Grid
#############################################
expGrid <- function (power.from, power.to, power.by, base)
{
    Grid <- c()
    for (i in seq(power.from, power.to, power.by)) {
        Grid <- append(Grid, base^i)
    }
    return(Grid)
}
