mod.cv.glmnet <- function(x, y, weights, offset = NULL, lambda = NULL, type.measure = c("mse", 
    "deviance", "class", "auc", "mae"), ..., nfolds = 10, foldid, 
    grouped = TRUE) {

    out <- glmnet::cv.glmnet(x=x,y=y, ..., nfolds = nfolds)
    val <- coef(out, s=out$lambda.min)[-1]
    return(val)
}
