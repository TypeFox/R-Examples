makenewdata <- 
# create a dataframe to pass to predict in order to calculate fitted models
# on a section defined by xc.cond
function (xs, xc.cond)
{
    newdata <- cbind(xs, xc.cond[rep(1L, nrow(xs)), ])
    colnames(newdata) <- c(colnames(xs), colnames(xc.cond))
    rownames(newdata) <- NULL
    return(newdata)	
}
