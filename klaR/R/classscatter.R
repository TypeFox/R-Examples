classscatter <- 
function(formula, data, method, col.correct = "black", 
    col.wrong = "red", gs = NULL, ...)
{
    variables <- dimnames(attributes(terms(formula, data=data))$factors)[[1]]
    response <- variables[1]
    discriminators <- variables[-1]
    if(any(discriminators == ".")) {
        exclude <- c(response, discriminators[discriminators != "."])
        discriminators <- colnames(data)[!is.element(colnames(data), exclude)]
    }
    x <- data[, discriminators] 
    grouping <- data[ , response]
    result <- do.call(method, list(x, grouping, ...))
    pr <- predict(result, x)
    if (is.list(pr)) pr <- pr$class
    nc <- pr != grouping
    color <- rep(col.correct, length(grouping))
    color[nc] <- col.wrong
    if (is.null(gs)) gs <- as.character(grouping)
    plot(x, col = color, pch = gs, ...)
    return(mean(nc))
}
