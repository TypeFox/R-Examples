getvarnames <- 
# helper function for extracting variable names from model objects
function (model)
{
    if (is.null(model$terms))
        stop("model has no terms slot")
    string1 <- deparse(model$terms[[3L]])
    string2 <- unlist(strsplit(string1, split = NULL))
    string3 <- paste(string2[string2 != " "], collapse = "")
    predictors1 <- unlist(strsplit(string3, split = "+", fixed = TRUE))
    predictors2 <- unique(vapply(predictors1, cleanstring, character(1L)))
    response <- unlist(deparse(model$terms[[2L]]))
    list(response = response, predictors = predictors2)
}