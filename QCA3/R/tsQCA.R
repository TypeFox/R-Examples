## this is part of the QCA3 package
## by Ronggui HUANG (2010)

demean <- function(x, id, theta=1) {
    ## id is a variable in x,x is matrix or data.frame
    ## this  is for fixed effect QCA
    id <- as.factor(id)
    xm <- apply(x,2,function(y,z) tapply(y,z, mean), z=id)
    xdm<- x[] <- x-theta*xm[id,]
    res <- list(xm=xm, xdm=xdm)
    res
}


timeDiff <- function (x, id, time, na.action=na.omit) {
    k <- 1 ## The number of lags (in units of observations)
    id <- as.factor(id)
    time <- as.factor(time)

    singleLag <- function(x, ak) {
        ## x is a vector rather than matrix or data frame
        ## ak must > 0
        isNAtime <- c(rep(1, ak), diff(as.numeric(time),
                                       lag = ak)) != ak
        isNAid <- c(rep(1, ak), diff(as.numeric(id), lag = ak)) !=
            0
        isNA <- as.logical(isNAtime + isNAid)
        if (is.factor(x))
            levs <- levels(x)
        result <- c(rep(NA, ak), x[1:(length(x) - ak)])
        result[isNA] <- NA
        if (is.factor(x))
            result <- factor(result, labels = levs)
        structure(result, class = class(x), id = id, time=time)
    }

    col.ID <- sapply(x,is.numeric)
    xlag <- apply(x[,col.ID], 2, function(ii) singleLag(ii, k))
    ans <- x
    ans[,col.ID] <- x[,col.ID] - xlag
    ans <- match.fun(na.action)(ans)
    ans
}


## how to handle zero?
## how to handle don't care cases?
tsData_pooled <- function(mydata,FUN=thresholdssetter, ...){
    n_col <- NCOL(mydata)
    for (i in 1:n_col){
        if (is.numeric(mydata[,i])) {
            ## only handle the numeric variable
            mydata[,i] <- do.call(FUN,
                                  c(list(x=mydata[,i],print.table=FALSE),
                                    list(...))
                                  )
        }
    }
    mydata
}

tsData_tdiff <- function(mydata, id, time) {
    mydata <- timeDiff(mydata,id=id,time=time,na.action=na.omit)
    n_col <- NCOL(mydata)
    for (i in 1:n_col){
        if (is.numeric(mydata[,i])) {
            ## only handle the numeric variable
            mydata[,i] <- as.numeric(mydata[,i] > 0)
        }
    }
    mydata
}

