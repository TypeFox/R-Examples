print.sgSEMp1 <- function(x){
    cat("--- An object of class sgSEMp1. A list of the following items ---\n")
    cat("--- Please see ?sgSEMp1 for details ---\n\n")
    cat("$data: The original data frame.\n")
    cat("$main predictor: The name of the main predictor used in the model.\n")
    cat("$response: The name of the response used in the model.\n")
    cat("$table: a table that gives the best models and the associated R-square, adjusted R-square and pvalues.\n")
    cat("$bestModels: A matrix of names of best models. Rows indicate predictors and columns indicate responses.\n")
    cat("$allModels: A three dimension array of all tried model objects. dim1 = predictors, dim2 = responses, dim3 = different models\n")
}
