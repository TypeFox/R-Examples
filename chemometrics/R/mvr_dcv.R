mvr_dcv <-
function (formula, ncomp, data, subset, na.action, method = c("kernelpls", 
    "widekernelpls", "simpls", "oscorespls", "svdpc"), scale = FALSE, repl = 100, 
    sdfact = 2, segments0 = 4, segment0.type = c("random", "consecutive", 
        "interleaved"), length.seg0, segments = 10, segment.type = c("random", 
        "consecutive", "interleaved"), length.seg, trace = FALSE, 
    plot.opt = FALSE, selstrat = "hastie", ...) 
{
    error.bars <- function(x, upper, lower, width = 0.02, ...) {
        xlim <- range(x)
        barw <- diff(xlim) * width
        segments(x, upper, x, lower, ...)
        segments(x - barw, upper, x + barw, upper, ...)
        segments(x - barw, lower, x + barw, lower, ...)
        range(upper, lower)
    }
#    require(pls)
    mf <<- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 
        0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    method <- match.arg(method, c("kernelpls", "widekernelpls", "simpls", "oscorespls", 
        "svdpc", "model.frame"))
    if (method == "model.frame") 
        return(mf)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "numeric")
    if (is.matrix(Y)) {
        if (is.null(colnames(Y))) 
            colnames(Y) <- paste("Y", 1:dim(Y)[2], sep = "")
    }
    else {
        Y <- as.matrix(Y)
        colnames(Y) <- deparse(formula[[2]])
    }
    X <- delintercept(model.matrix(mt, mf))
    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    if (length(attr(mt, "term.labels")) == 1 && !is.null(colnames(mf[[attr(mt, 
        "term.labels")]]))) 
        colnames(X) <- sub(attr(mt, "term.labels"), "", colnames(X))
    if (missing(ncomp)) {
        ncomp <- min(nobj - 1, npred)
        ncompWarn <- FALSE
    }
    else {
        if (ncomp < 1 || ncomp > min(nobj - 1, npred)) 
            stop("Invalid number of components, ncomp")
        ncompWarn <- TRUE
    }
    Y <- as.matrix(Y)
    dy <- dim(Y)
    dx <- dim(X)
    dnX <- dimnames(X)
    dnY <- dimnames(Y)
    dimnames(X) <- dimnames(Y) <- NULL
    if (!is.logical(scale) || length(scale) != 1) 
        stop("'scale' must be 'TRUE' or 'FALSE'")
    ncomp <- min(ncomp, dx[1] - max(sapply(segments0, length)) - 
        1)
    pred <- array(dim = c(dx[1], dy[2], ncomp, repl))
    predopt <- array(dim = c(dx[1], dy[2], repl))
    optcomp <- matrix(NA, nrow = segments0, ncol = repl)
    for (i in 1:repl) {
        print(i)
        if (missing(length.seg0)) {
            segment0 <- cvsegments(dx[1], k = segments0, type = segment0.type)
        }
        else {
            segment0 <- cvsegments(dx[1], length.seg = length.seg0, 
                type = segment0.type)
        }
        if (trace) 
            cat(paste("Replication: ", i))
        for (n.seg0 in 1:length(segment0)) {
            if (trace) 
                cat(n.seg0, "")
            seg0 <- segment0[[n.seg0]]
            obsuse <- as.numeric(unlist(segment0[-n.seg0]))
            res <- mvr(Y ~ X, ncomp = ncomp, subset = obsuse, 
                na.action = na.action, method = method, scale = scale, 
                segments = segments, segment.type = segment.type, 
                trace = trace, validation = "CV")
            MSEPj <- matrix(NA, nrow = segments, ncol = ncomp)
            for (j in 1:segments) {
                predj <- res$vali$pred[res$vali$seg[[j]], , ]
                obsj <- obsuse[res$vali$seg[[j]]]
                resj <- predj - Y[obsj]
                MSEPj[j, ] <- apply(resj^2, 2, mean)
            }
            MSEPm <- apply(MSEPj, 2, mean)
            MSEPsd <- apply(MSEPj, 2, sd)/sqrt(segments)
            if (selstrat == "diffnext") {
                fvec <- (diff(MSEPm) + sdfact * MSEPsd[-1]) < 
                  0
                fvec <- c(TRUE, fvec)
                ind <- which.min(MSEPm)
                optcomp[n.seg0, i] <- max((1:ind)[fvec[1:ind]])
            }
            else if (selstrat == "hastie") {
                ind <- which.min(MSEPm)
                fvec <- (MSEPm < (MSEPm[ind] + sdfact * MSEPsd[ind]))
                optcomp[n.seg0, i] <- min((1:ind)[fvec[1:ind]])
            }
            else if (selstrat == "relchange") {
                ind <- which.min(MSEPm)
                MSEPsel <- MSEPm[1:ind]
                relchange <- (MSEPsel - MSEPm[ind])/max(MSEPsel) > 
                  0.001
                ind2 <- which.max((1:length(relchange))[relchange])
                MSEPm2 <- MSEPsel[1:ind2]
                MSEPsd2 <- MSEPsd[1:ind2]
                indm <- which.min(MSEPm2)
                fvec <- (MSEPm2 < (MSEPm2[indm] + sdfact * MSEPsd2[indm]))
                optcomp[n.seg0, i] <- min((1:indm)[fvec[1:indm]])
            }
            if (plot.opt) {
                plot(1:ncomp, MSEPm, xlab = "Component number", 
                  ylab = "Average MSEP", ylim = range(MSEPm, 
                    MSEPm + MSEPsd, MSEPm - MSEPsd), type = "b")
                error.bars(1:ncomp, MSEPm + MSEPsd, MSEPm - MSEPsd, 
                  width = 1/ncomp, col = 3)
                abline(v = optcomp[n.seg0, i], col = 2)
            }
            predopt[seg0, , i] <- predict(res, data, ncomp = optcomp[n.seg0, 
                i])[seg0, , ]
            pred[seg0, , , i] <- predict(res, data)[seg0, , ]
        }
    }
    resopt <- predopt - c(Y)
    MSEPopt <- mean(as.vector(resopt)^2)
    biasopt <- mean(resopt)
    SEPopt <- sqrt(sum((resopt - biasopt)^2)/(prod(dim(resopt)) - 
        1))
    sIQRopt <- IQR(resopt)/1.349
    sMADopt <- mad(resopt)
    objnames <- dnX[[1]]
    if (is.null(objnames)) 
        objnames <- dnY[[1]]
    yvarnames <- dnY[[2]]
    nCompnames <- paste(1:ncomp, "comps")
    nreplnames <- paste(1:repl, "repl")
    nsegnames <- paste(1:segments0, "segm")
    dimnames(pred) <- list(objnames, yvarnames, nCompnames, nreplnames)
    dimnames(predopt) <- list(objnames, yvarnames, nreplnames)
    dimnames(resopt) <- list(objnames, yvarnames, nreplnames)
    dimnames(optcomp) <- list(nsegnames, nreplnames)

    afinaldistr <- table(optcomp)/sum(table(optcomp))
    afinal <- as.numeric(names(which.max(afinaldistr)))
    residcomp <- pred-c(Y)
    biascomp <- apply(residcomp,3,mean)
    dimr <- dim(residcomp)
    biascomp1 <- array(biascomp,c(dimr[3],dimr[1],dimr[2],dimr[4]))
    biascomp1 <- aperm(biascomp1,c(2,3,1,4)) # permute array
    SEPfinal <- sqrt(apply((residcomp-biascomp1)^2,3,sum)/(prod(dim(resopt))-1))

    list(resopt = resopt, predopt = predopt, optcomp = optcomp, 
        pred = pred, SEPopt = SEPopt, sIQRopt = sIQRopt, sMADopt = sMADopt, 
        MSEPopt = MSEPopt,afinal=afinal,SEPfinal=SEPfinal)
}
