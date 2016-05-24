relnorm <-
function(master, tofix, mask, method="MA", nperm=1000)
{
    # relative normalization of tofix to match master
    # based on values in mask and using regression
    # method OLS, MA, or SMA
    # mask should contain NA for values to include; all other values will be omitted

    results <- tofix

    master <- as.vector(as.matrix(master))
    tofix <- as.vector(as.matrix(tofix))

    if(missing(mask)) { 
        mask <- rep(NA, length(master))
    } else {
        mask <- as.vector(as.matrix(mask))
    }

    master.mask <- master[is.na(mask)]
    x.mask <- tofix[is.na(mask)]

    master.lm <- lmodel2(master.mask ~ x.mask, nperm=nperm)
    
    master.lm <- master.lm$regression.results[master.lm$regression.results[, "Method"] == method, ]
    names(master.lm) <- gsub("^ *", "", names(master.lm))
    
    x.transform <- master.lm$Slope * tofix + master.lm$Intercept
    
    x.transform[!is.na(mask)] <- NA


    # return the same structure as the input values
    if(class(results) == "SpatialGridDataFrame")
        results@data[,1] <- x.transform
    else if(is.data.frame(results))
        results <- data.frame(matrix(x.transform, nrow=nrow(results), ncol=ncol(results)))
    else if(is.matrix(results))
        results <- matrix(x.transform, nrow=nrow(results), ncol=ncol(results))
    else # return a vector 
        results <- x.transform

    list(regression.results = master.lm, newimage = results)
}

