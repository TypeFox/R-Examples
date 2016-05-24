sffs <-
function (data, method = c("lda", "knn", "rpart"), kvec = 5, 
    repet = 10) 
{
#    require("MASS")
#    require("class")
#    require("rpart")
    if (sum(is.na(data))> 0) 
        stop("This dataset has missing values, impute them before running this function.\n",call.=FALSE)
    if (!(method %in% c("lda", "knn", "rpart"))) {
        cat("The classifier entered is not supported by this function.\n")
        return(method)
    }
    n = dim(data)[1]
    p = dim(data)[2]
    grupos = data[, p]
    ngroups = dim(table(data[, p]))
    selected = rep(0, p)
    numselect = 0
    for (j in 1:repet) {
        indic <- rep(0, p - 1)
        correcto <- 0
        paso1 <- sfs1(data, indic, correcto, kvec, method)
        correcto <- paso1$accuracy
        indic <- paso1$indic
        i <- 2
        while (i <= (p - 1)) {
            paso2 <- sfs1(data, indic, correcto, kvec, method)
            if (paso2$accuracy > correcto) {
                correcto <- paso2$accuracy
                indic <- paso2$indic
                for (j in 1:(i - 1)) {
                  paso3 <- sbs1(data, indic, correcto, kvec, 
                    method)
                  correcto <- paso3$correcto
                  indic <- paso3$indic
                }
            }
            else {
                i <- p
            }
        }
        variables <- seq(1, (p - 1))
        variables <- variables[indic == 1]
        numselect = numselect + length(variables)
        selected[variables] = selected[variables] + 1
    }
    numselect = round(numselect/repet)
    fselect = order(selected, decreasing = TRUE)[1:numselect]
    cat("\nThe selected features are:\n")
    return(fselect)
}
