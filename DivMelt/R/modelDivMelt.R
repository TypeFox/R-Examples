modelDivMelt <- function( trngFile="", datFile="model.rda")
{
    if(trngFile == "") stop("Need training set argument to generate model")

    library ("glmnet")

    training <- read.csv (trngFile, header=T)

    if(ncol(training) != 13) stop("Incorrect number of columns in training set.  There should be 13 columns!")

    col.features <- 4:12
    col.classes <- 13

    names.features <- names (training)[col.features]

    training.features <- as.matrix (training[,col.features])
    training.classes <- training[,col.classes]

    model <- cv.glmnet (training.features, training.classes, alpha=1, family="binomial")
    pred.results <- cbind (training.classes, predict (model, training.features, s="lambda.min"))
    coefficients <- coef (model)

    tp <- sum (pred.results[pred.results[,1] == 1, 2] > 0)
    fp <- sum (pred.results[pred.results[,1] == 1, 2] <= 0)

    tn <- sum (pred.results[pred.results[,1] == 0, 2] <= 0)
    fn <- sum (pred.results[pred.results[,1] == 0, 2] > 0)

    save(model,file=datFile)
}
