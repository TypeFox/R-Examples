# rpart.R: plotmo methods for rpart objects

plotmo.type.rpart <- function(object, ...)
{
    # use same default as predict.rpart
    if(object$method == "class")
        "prob"
    else
        "vector"
}
plotmo.residtype.rpart <- function(object, ...)
{
    "usual"
}
plotmo.singles.rpart <- function(object, x, ...)
{
    # get all variables used in the tree
    varnames <- as.character(object$frame$var) # factor to character
    varnames <- unique(varnames[varnames != "<leaf>"])
    match(varnames, colnames(x))
}
plotmo.pairs.rpart <- function(object, x, ...)
{
    # we consider rpart variables paired if one is the direct
    # parent of the other in the tree.
    irow <- as.integer(row.names(object$frame))
    var.names <- character(length=max(irow))
    var.names[irow] <- as.character(object$frame$var) # factor to character
    ivar <- charmatch(var.names, colnames(x))
    # following is the same as var.names != "<leaf>" & var.names !=""
    is.split <- !is.na(ivar) & ivar > 0
    if(sum(is.split) == 0) # no splits? (intercept-only model)
        return(NULL)
    pairs <- NULL
    for(i in 1:length(ivar)) {
        if(is.split[i]) {
            left <- 2 * i
            if(left <= length(ivar) && is.split[left] && ivar[i] != ivar[left])
                pairs <- c(pairs, ivar[i], ivar[left])
            right <- left + 1
            if(right <= length(ivar) && is.split[right] && ivar[i] != ivar[right])
                pairs <- c(pairs, ivar[i], ivar[right])
        }
    }
    if(!is.null(pairs))
        pairs <- matrix(pairs, ncol=2, byrow=TRUE)
    pairs
}
plotmo.predict.rpart <- function(object, newdata, type, ..., TRACE)
{
    # do some hand holding to avoid obscure message from predict.rpart
    pmatch <- pmatch(object$method, c("anova", "class", "exp", "poisson"))
    if(pmatch == 2) { # class
        if(!pmatch(type, c("vector", "prob", "matrix", "class"), nomatch=0))
            stop0("predict.rpart does not support type=\"", type, "\"")
    } else if(!pmatch(type, c("vector", "matrix"), nomatch=0))
        stop0("predict.rpart does not support type=\"", type, "\" (for \"",
              object$method, "\" rpart objects)")
    plotmo.predict.default(object, newdata, type=type, ..., TRACE=TRACE)
}
