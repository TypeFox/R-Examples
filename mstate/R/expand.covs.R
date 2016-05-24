expand.covs <- function(data, ...) UseMethod("expand.covs")

expand.covs.default <-
    function (data, covs, append = TRUE, longnames = FALSE, event.types="failcode", ...)
{
    comp.risks <- unique(data[[event.types]])
    if(length(comp.risks)==1)
        stop("Function does not create type-specific covariates with only one event type")
    trans <- trans.comprisk(K=length(comp.risks), names=comp.risks)
    trans2 <- to.trans2(trans)
    K <- nrow(trans2)
    if (is.character(covs))
        form1 <- as.formula(paste("~ ", paste(covs, collapse = " + ")))
    else form1 <- as.formula(covs)
    mm4 <- NULL
    for (j in 1:length(covs)) {
        wh <- which(!is.na(data[[covs[j]]]))
        form1 <- as.formula(paste("~ ", covs[j]))
        mm <- model.matrix(form1, data = data)
        mm <- data.frame(mm)
        mm <- mm[, -1, drop = FALSE]
        if (!longnames) {
            nc <- ncol(mm)
            if (nc == 1)
                names(mm) <- covs[j]
            else names(mm) <- paste(covs[j], 1:nc, sep = "")
        }
        nms <- names(mm)
        ms <- data.frame(trans = data[[event.types]])
        ms$trans <- factor(ms$trans, levels=comp.risks)
        ms <- cbind(ms[wh, , drop = FALSE], mm)
        mm2 <- model.matrix(as.formula(paste("~ (", paste(nms,
                                                          collapse = " + "), "):trans")), data = ms)[, -1]
        mm3 <- matrix(NA, nrow(data), ncol(mm2))
        mm3[wh, ] <- mm2
        mm3 <- data.frame(mm3)
        nms <- as.vector(t(outer(nms, comp.risks, "paste", sep = ".")))
        names(mm3) <- nms
        if (j == 1)
            mm4 <- mm3
        else mm4 <- cbind(mm4, mm3)
    }
    if (!append)
        return(mm4)
    else {
        if (!all(is.na(match(names(data), nms))))
            warning("One or more names of appended data already in data!")
        mm4 <- cbind(data, mm4)
    }
    class(mm4) <- class(data)
    return(mm4)
}

expand.covs.msdata <- function(data, covs, append=TRUE, longnames=TRUE, ...)
{
    if (!inherits(data, "msdata"))
        stop("'data' must be an 'msdata' object")
    trans <- attr(data, "trans")
    data <- as.data.frame(data)
    trans2 <- to.trans2(trans)
    K <- nrow(trans2)
    if (is.character(covs)) form1 <- as.formula(paste("~ ",paste(covs,collapse=" + ")))
    else form1 <- as.formula(covs)
    # going to apply model.matrix, but NA's are not allowed, so have to deal with that
    mm4 <- NULL
    for (j in 1:length(covs)) {
        wh <- which(!is.na(data[[covs[j]]]))
        form1 <- as.formula(paste("~ ",covs[j]))
        mm <- model.matrix(form1,data=data)
        mm <- data.frame(mm)
        mm <- mm[,-1,drop=FALSE]
        if (!longnames) {
            nc <- ncol(mm)
            if (nc==1) names(mm) <- covs[j]
            else names(mm) <- paste(covs[j],1:nc,sep="")
        }
        nms <- names(mm)
        ms <- data.frame(trans=data[["trans"]])
        ms$trans <- factor(ms$trans)
        ms <- cbind(ms[wh,,drop=FALSE],mm)
        mm2 <- model.matrix(as.formula(paste("~ (",paste(nms,collapse=" + "),"):trans")),data=ms)[,-1]
        mm3 <- matrix(NA,nrow(data),ncol(mm2))
        mm3[wh,] <- mm2
        mm3 <- data.frame(mm3)
        nms <- as.vector(t(outer(nms,1:K,"paste",sep=".")))
        names(mm3) <- nms
        if (j==1) mm4 <- mm3 else mm4 <- cbind(mm4,mm3)
    }
    if (!append) return(mm4)
    else {
        if (!all(is.na(match(names(data),nms))))
            warning("One or more names of appended data already in data!")
        mm4 <- cbind(data,mm4)
    }
    attr(mm4, "trans") <- trans
    class(mm4) <- c("msdata", "data.frame")
    return(mm4)
}
