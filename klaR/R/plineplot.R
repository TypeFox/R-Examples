plineplot <- function(formula, data, method, x, col.wrong = "red", 
    ylim = c(0, 1), loo = FALSE, mfrow, ...)
{
    variables <- dimnames(attributes(terms(formula, data=data))$factors)[[1]]
    response <- variables[1]
    discriminators <- variables[-1]
    if(any(discriminators == ".")) {
        exclude <- c(response, discriminators[discriminators != "."])
        discriminators <- colnames(data)[!is.element(colnames(data), exclude)]
    }
    vars <- data[, discriminators] 
    grouping <- data[, response]
    if (loo){
        pr <- NULL
        n <- nrow(vars)
        for(i in 1:n){
            result <- do.call(method, list(vars[-i,], grouping[-i], ...))
            pr2 <- predict(result, vars[i,])
            if(is.list(pr2)) pr2 <- pr2$post  
            pr <- rbind(pr, pr2)
        }
    }
    else{
        result <- do.call(method, list(vars, grouping, ...))
        pr <- predict(result, vars)
        if(is.list(pr)) pr <- pr$post
    }
    wrong <- (factor(levels(grouping)[max.col(pr)])) != grouping
    if(length(x) == 1) x <- data[ , x]
    reihen <- order(x)
    pr <- pr[reihen,]
    x <- x[reihen]
    k <- ncol(pr)
    wrong <- wrong[reihen]
    if(missing(mfrow)) mfrow <- c(k, 1)
    opar <- par(mfrow = mfrow)
    on.exit(par(opar))
    for(i in 1:k){
        plot(x, pr[ , i], type = "b", ylab = levels(grouping)[i], ylim = ylim, ...)
        points(x[wrong], pr[wrong, i], col = col.wrong, ...)
    }
    return(mean(wrong))
}
