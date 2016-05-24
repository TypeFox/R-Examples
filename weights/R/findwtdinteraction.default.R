findwtdinteraction.default <- function(reg, across, by=NULL, at=NULL, acrosslevs=NULL, bylevs=NULL, atlevs=NULL, weight=NULL, dvname=NULL, acclevnames=NULL, bylevnames=NULL, atlevnames=NULL, stdzacross=FALSE, stdzby=FALSE, stdzat=FALSE, limitlevs=20, type="response"){
    df <- reg$model
    clsset <- sapply(lapply(df, class), function(x) x[1])
    acclass <- class(df[,across])[1]
    if(stdzacross==TRUE)
        df[,across] <- stdz(df[,across], weight)
    if(stdzby==TRUE & !is.null(by))
        df[,by] <- stdz(df[,by], weight)
    if(stdzat==TRUE & !is.null(at))
        df[,at] <- stdz(df[,at], weight)
    if(is.null(acrosslevs)){
        if(stdzacross==TRUE)
            acrosslevs <- c(-1,1)
        if(stdzacross==FALSE)
            acrosslevs <- sort(unique(df[,across]))
        if(length(acrosslevs)>limitlevs & is.numeric(acrosslevs))
            acrosslevs <- seq(min(df[,across], na.rm=TRUE), max(df[,across], na.rm=TRUE), (max(df[,across], na.rm=TRUE)-min(df[,across], na.rm=TRUE))/(limitlevs-1))
    }
    if(is.null(weight)){
        if(!is.null(reg$prior.weight))
            weight <- reg$prior.weight
        if(is.null(reg$prior.weight) & !is.null(reg$weights))
            weight <- reg$weights
        if(is.null(weight))
            weight <- rep(1, dim(df)[1])
    }
    if(length(weight)!=(dim(df)[1]))
        stop("Weight vector length must match the number of complete cases in the regression.")
    if(is.null(dvname))
        dvname <- names(df)[1]
    if(!is.null(by)){
        byclass <- class(df[,by])[1]
        if(is.null(bylevs)){
            if(stdzby==TRUE)
                bylevs <- c(-1,1)
            if(stdzby==FALSE)
                bylevs <- sort(unique(df[,by]))
            if(length(bylevs)>limitlevs & is.numeric(bylevs))
                bylevs <- seq(min(df[,by], na.rm=TRUE), max(df[,by], na.rm=TRUE), (max(df[,by], na.rm=TRUE)-min(df[,by], na.rm=TRUE))/(limitlevs-1))
        }
    }
    if(is.null(by)){
        by <- "All"
        df[,by] <- "All"
        bylevs <- "All"
        hasby <- FALSE
        stdzby <- FALSE
    }
    if(!is.null(at)){ 
        atclass <- class(df[,at])[1]
        if(is.null(atlevs)){
            if(stdzat==TRUE)
                atlevs <- c(-1,1)
            if(stdzat==FALSE)
                atlevs <- sort(unique(df[,at]))
            if(length(atlevs)>limitlevs & is.numeric(atlevs))
                atlevs <- seq(min(df[,at], na.rm=TRUE), max(df[,at], na.rm=TRUE), (max(df[,at], na.rm=TRUE)-min(df[,at], na.rm=TRUE))/(limitlevs-1))
        }
        hasat <- TRUE
        if(is.null(atlevnames)){
            if(stdzat==TRUE)
                atlevnames <- paste(atlevs, "SD", sep="")
            if(stdzat==FALSE)
                atlevnames <- paste(atlevs)
        }
    }
    if(is.null(at)){
        at <- "All"
        df[,at] <- "All"
        atlevs <- "All"
        atlevnames <- "All"
        hasat <- FALSE
        stdzat <- FALSE
    }
    if(is.null(bylevnames)){
        if(stdzby==TRUE)
            bylevnames <- paste(bylevs, "SD", sep="")
        if(stdzby==FALSE)
            bylevnames <- paste(bylevs)
    }
    if(is.null(acclevnames)){
        if(stdzacross==TRUE)
            acclevnames <- paste(acrosslevs, "SD", sep="")
        if(stdzacross==FALSE)
            acclevnames <- paste(acrosslevs)
    }
    levs <- acrosslevs
    ol <- acrosslevs
    lng <- length(acrosslevs)
    pd <- data.frame(df[1:lng,])
    for(i in 1:dim(pd)[2]){
        if(class(pd[,i])[1]=="numeric")
            pd[,i] <- rep(wtd.mean(df[,i], na.rm=TRUE), lng)
        if(class(pd[,i])[1]=="ordered")
            pd[,i] <- ordered(rep(names(cumsum(wtd.table(df[,i], weight)$sum.of.weights)/sum(wtd.table(df[,i], weight)$sum.of.weights))[cumsum(wtd.table(df[,i], weight)$sum.of.weights)/sum(wtd.table(df[,i], weight)$sum.of.weights)>=.5][1], lng), levels=levels(df[,i]))
        if(class(pd[,i])[1]=="factor")
            pd[,i] <- factor(rep(names(wtd.table(df[,i], weight)$sum.of.weights[wtd.table(df[,i], weight)$sum.of.weights==max(wtd.table(df[,i], weight)$sum.of.weights)]), lng), levels=levels(df[,i]))
        if(class(pd[,i])[1]=="logical")
            pd[,i] <- as.logical(rep(names(wtd.table(df[,i], weight)$sum.of.weights[wtd.table(df[,i], weight)$sum.of.weights==max(wtd.table(df[,i], weight)$sum.of.weights)]), lng))
    }    
    if(sum(clsset=="matrix")>0)
        stop(paste("Interactions Cannot Currently Be Resolved With Matrix Predictors, Please Insert Each Variable in", names(clsset)[clsset=="matrix"], "Separately in Regression Before Using This Tool")) # TRY TO MAKE THIS WORK EVENTUALLY
    names(pd) <- names(df)
    rownames(pd) <- 1:dim(pd)[1]
    pd[,across] <- acrosslevs
    out <- NULL
    out$Resp <- pd[1,!(colnames(pd) %in% c(dvname, at, across, by))]
    out$Meta <- list(dvname=dvname, across=across, by=by, at=at)
    out$Means <- as.list(1:length(atlevs))
    out$SEs <- as.list(1:length(atlevs))
    names(out$Means) <- names(out$SEs) <- atlevnames
    for(a in 1:length(atlevs)){
        pd[,at] <- atlevs[a]
        if(!is.null(df[,at]))
        class(pd[,at]) <- class(df[,at])
        bylist <- as.list(bylevs)
        out$Means[[a]] <- out$SEs[[a]] <- matrix(NA, length(bylevs), lng)
        for(i in 1:length(bylevs)){
            bylist[[i]] <- pd
            bylist[[i]][,by] <- factor(rep(bylevs[i], length(bylist[[i]][,by])), levels=bylevs)
            if(is.numeric(bylist[[i]][,by]))
                bylist[[i]] <- pd
            bylist[[i]][,by] <- rep(bylevs[i], length(bylist[[i]][,by]))
        }
        eachpred <- lapply(bylist, function(x) predict(reg, newdata=x, se.fit=TRUE, type=type))
        means <- t(sapply(eachpred, function(x) x$fit))
        out$Means[[a]] <- as.matrix(means)
        ses <- t(sapply(eachpred, function(x) x$se.fit))
        out$SEs[[a]] <- as.matrix(ses)
        try(rownames(out$Means[[a]]) <- rownames(out$SEs[[a]]) <- bylevnames)
        try(colnames(out$Means[[a]]) <- colnames(out$SEs[[a]]) <- acclevnames)
    }
    class(out) <- "interactpreds"
    out
}
