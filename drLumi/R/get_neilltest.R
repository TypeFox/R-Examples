## Estimate of model fit based on neill.test
## The code is an adapatation of drc::neill.test but has been 
## changed in some sections because some errors occurs in the original 
## function with 'nls' class models. The drc::neill.test, version 2.3-96
get_neilltest <- function(object, x, method, grouping){
    noCluster <- floor(length(x)/2)
    if(is.null(grouping)){
        if (method == "finest") {
            lenx <- length(x)
            grouping <- floor((1 + 1:lenx)/2)
            grouping[lenx] <- grouping[lenx - 1]
        }
        if (method == "c-finest") {
            for (i in noCluster:(length(coef(object)) + 1)) {
                grouping <- cutree(hclust(dist(x)), k = i)
                if (all(tapply(x, grouping, length) > 1)) {
                    break
                }
            }
        }
        if (method == "percentiles") {
            cutVar <- c(-Inf, quantile(x, c(0.2, 0.4, 0.6, 0.8)), 
                    Inf)
            grouping <- cut(x, cutVar)
        }
    }

    


    M <- length(unique(grouping))
    lhs <- as.character(object$m$formula()[[2]])
    parameters <- names(object$m$getPars())
    allobj <- ls(object$m$getEnv())
    rhs <- allobj[-match(c(parameters,lhs),allobj)]
    ndf <- data.frame(get(lhs, object$m$getEnv()), get(rhs, object$m$getEnv()))
    names(ndf) <- c(lhs, rhs)
    norder <- order(ndf[,rhs], decreasing=TRUE)
    ndf <- ndf[norder,]

    N <- nrow(ndf)
    denDF <- N - M

    ## Checking the number of groups
    mes <- NULL
    if (denDF <= 0)  # (N <= M) 
    {
        # "error1.Too many groups in 'grouping'"
        mes <- "error1"
    }
    p <- N - df.residual(object)
    numDF <- M - p

    if (numDF <= 0)  # (M <= p)
    {
        # "error2. Too few groups in 'grouping'"
        mes <- "error2"
    }

    if(is.null(mes)){    
        ## Calculating the test statistic
        resVec <- residuals(object)
        resVec <- resVec[norder]
        resAver0 <- tapply(resVec, grouping, mean)
        resAver <- rep(resAver0, tapply(grouping, grouping, length))
    
        resDiff <- resVec - resAver
        FF <- (denDF/numDF)*(sum(resAver*resAver)/(sum(resDiff*resDiff)))
        p <- pf(FF ,numDF, denDF, lower.tail = FALSE)
        ans <- matrix(c(FF, p), 1, 2)
        colnames(ans) <- c("F","p")
    } else {
        ans <- mes
    }
    return(ans)
}
