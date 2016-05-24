#################################################
## HELPER FUNCTIONS (ALL INTERNAL TO THE PACKAGE)
#################################################

##################
# LOGISTIC MODELS
##################

.nPL2 <- function(bottom, top, xmid, scal, s,  X){
    yfit <- 1/(1+10^((xmid-X)*scal))
    return(yfit)
}
.nPL3 <- function(bottom, top, xmid, scal, s,  X){
    yfit <- (top)/(1+10^((xmid-X)*scal))
    return(yfit)
}
.nPL4 <- function(bottom, top, xmid, scal, s,  X){
    yfit <- bottom+(top-bottom)/(1+10^((xmid-X)*scal))
    return(yfit)
}
.nPL5 <- function(bottom, top, xmid, scal, s,  X){
    yfit <- bottom+(top-bottom)/(1+10^((xmid-X)*scal))^s
    return(yfit)
}
.chooseModel <- function(npars){
  switch(as.character(npars),
        "2" = {nPL <- .nPL2},
        "3" = {nPL <- .nPL3},
        "4" = {nPL <- .nPL4},
        "5" = {nPL <- .nPL5}
  )
  return(nPL)
}

##################
# WEIGHT MODELS
##################

.wsqRes <- function(x, yobs, yfit, LPweight) {
    residuals <- (yobs - yfit)^2
    return((1/residuals)^(LPweight))
}
.sdWeight <- function(x, yobs, yfit, LPweight){
    v <- as.numeric(by(yobs, x, var, na.rm=TRUE))
    v <- ifelse(is.na(v), 1, v)
    return(1/rep(v, times=table(x)))
}
.generalWeight <- function(x, yobs, yfit, LPweight){
    return(1/yfit^LPweight)
}
.sce <- function(pars, x, yobs, .weight, LPweight, nPL){
    bottom <- pars[1]
    top <- pars[2]
    xmid <- pars[3]
    scal <- pars[4]
    s <- pars[5]
    yfit <- nPL(bottom, top, xmid, scal, s, x)
    residuals <- (yobs - yfit)^2
    w <- .weight(x, yobs, yfit, LPweight)
    return(sum(w*residuals))
}
.chooseWeight <- function(method){
  switch(method,
        res = {.weight <- .wsqRes},
        sdw = {.weight <- .sdWeight},
        gw = {.weight <- .generalWeight}
        )
    return(.weight)
}

##################
# PARS INIT
##################

.yTOz <- function(y){
    bottom <- as.numeric(quantile(y, .025, na.rm=TRUE))
    top <- as.numeric(quantile(y, .975, na.rm=TRUE))
    z <- (y - bottom)/(top - bottom)
    z[z<=0] <- 0.05; z[z>=1] <- 0.95
    z
}
.lmlz <- function(x, y){
    z <- .yTOz(y)
    z <- as.vector(by(z, x, mean, na.rm = TRUE))
    lz <- log(z/(1-z))
    lm(lz ~ unique(x))
}
.estimScal <- function(x, y){
    model <- .lmlz(x, y)
    return(coef(model)[2])
}
.estimMid <- function(x, y){
    model <- .lmlz(x, y)
    fit <- model$fitted.values
    predict(model, data.frame(x = median(fit, na.rm = TRUE)))  
}
.initPars <- function(x, y, npars){
    if(npars<5) s <- 1 else s <- 1
    if(npars<4) bottom <- 0 else bottom <- quantile(y, .05,na.rm=TRUE)
    if(npars<3) top <-1 else top <- quantile(y, .95,na.rm=TRUE)
    xmid <- (max(x)+min(x))/2 #.estimMid(x, y) #(max(x)+min(x))/2
    scal <- .estimScal(x, y)
    c(bottom, top, xmid, scal, s)
}

##################
# FITTING
##################

.getPars <- function(model){
    bottom <- model$estimate[1]
    top <- model$estimate[2]
    xmid<-model$estimate[3]
    scal <- model$estimate[4]
    s <- model$estimate[5]
    return(cbind.data.frame(bottom=bottom, top=top, xmid=xmid, scal=scal, s=s))
}
.fit <- function(x, y, npars, nPL, .sce, LPweight){
    best <- nlm(f=.sce, p=.initPars(x, y, npars), x=x, yobs=y,
        Weights=rep(1, length(x)), wcoef=LPweight, nPL)
    pars <- best$estimate
    return(nPL(pars[1], pars[2], pars[3], pars[4], pars[5], unique(x)))
}
.inflPoint <- function(pars){
    x = pars$xmid + (1/pars$scal)*log10(pars$s)
    y = pars$bottom + (pars$top - pars$bottom)*(pars$s/(pars$s+1))^pars$s
    return(cbind.data.frame(x=x, y=y))
}

##################
# PERFORMANCES
##################

.testAll <- function(.sce, x, y, .weight, LPweight, silent){
    if(!silent)
        message("Testing pars...")
    err <- sapply(1, function(p){
        test2 <- try(nlm(f=.sce, p=.initPars(x, y, 2), x=x, yobs=y, .weight, LPweight, .nPL2), silent=TRUE)
        test3 <- try(nlm(f=.sce, p=.initPars(x, y, 3), x=x, yobs=y, .weight, LPweight, .nPL3), silent=TRUE)
        test4 <- try(nlm(f=.sce, p=.initPars(x, y, 4), x=x, yobs=y, .weight, LPweight, .nPL4), silent=TRUE)
        test5 <- try(nlm(f=.sce, p=.initPars(x, y, 5), x=x, yobs=y, .weight, LPweight, .nPL5), silent=TRUE)
        scores <- sapply(list(test2, test3, test4, test5), function(t){
        if(class(t)!="try-error")
            return(t$minimum)
        else
            return(Inf) 
        })
        return(scores)
    })
  return( list(err=as.numeric(err), npars=which.min(err) + 1) )
}
.gof <- function(y, yfit, w){
    n <- length(y)
    S2y <- var(y)
    Rw <- 1 - sum(w^2)/((n-1)*S2y)
    return(Rw)
}
.getPerf <- function(y, yfit){
    w <- (y - yfit)^2
    lmtest <- summary(lm(y ~ yfit, weights=w))
    fstat <- lmtest$fstatistic
    p <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
    goodness <- .gof(y, yfit, w)
    n <- sum(w!=0)
    W <- n/((n-1)*sum(w))
    stdErr <- sqrt(W*sum(w*(yfit-y)^2))
    return(cbind.data.frame(goodness=goodness, stdErr=stdErr, p=p))
}

##################
# AREAS
##################

.AUC <- function(x, y){
    auc <- lapply(2:length(x), function(i){
        da <- x[i]-x[i-1]
        db <- y[i]-y[i-1]
        y[i]*da +1/2*db*da
      })
  return(do.call(sum, auc))
}
.Simpson <- function(x, y){
    dx <- mean(diff(x, lag = 1), na.rm = TRUE)
    n <- length(y)
    if(n%%2 != 0){
        x <- x[-n]
        y <- y[-n]
        n <- length(x)
        }
    f1 <- y[1]
    fn <- y[n]
    fy <- y[2:(n-1)]*rep(c(4, 2), (n-2)/2)
    return(dx/3*(f1 + sum(fy) + fn))
}


##################
# ESTIMATE RESP
##################

.invModel <- function(pars, target){
    return(pars$xmid - 1/pars$scal*
        log10(((pars$top - pars$bottom)/(target - pars$bottom))^(1/pars$s)-1))
}
.estimateRange <- function(target, stdErr, pars, B, useLog, conf.level){
    a1 <- (1-conf.level)/2
    a2 <- 1-(1-conf.level)/2
    #  xtarget = .invModel(pars, target)
    if(target<=pars$bottom || target>=pars$top){
        xlower <- xtarget <- xupper <- NA
    } else{
        Ytmp <- target + rnorm(B, 0, stdErr)
        if(any(Ytmp<=pars$bottom))
            Ytmp <- Ytmp[-which(Ytmp<=pars$bottom)]
        if(any(Ytmp>=pars$top))
            Ytmp <- Ytmp[-which(Ytmp>=pars$top)]
        Q <- quantile(Ytmp, probs=c(a1, a2), na.rm=T)
        estimates <- .invModel(pars, c(Q[1], target, Q[2]))
        if(useLog)
            estimates <- 10^estimates
        xlower <- min(estimates, na.rm=TRUE)
        xtarget <- estimates[2]
        xupper <- max(estimates, na.rm=TRUE)
        }
    return(as.numeric(c(xlower, xtarget, xupper)))
}
.confInt <- function(stdErr, yobs, yfit, newy){
    n <- length(yobs)
    ybar <- mean(yobs, na.rm = TRUE)
    t <- qt(.975, n-2)
    ci <- t*stdErr*sqrt((1/n+(newy - ybar)^2/sum((newy - ybar)^2)))
    lo <- newy - ci
    hi <- newy + ci
    return(list(lo = lo, hi = hi))
}

######################
# PLOT FUNCTIONS
######################
.plot <- function(object,...){
    x <- getX(object)
    y <- getY(object)
    gof <- format(getGoodness(object), digits=4, scientific = TRUE)
    plot(x, y, type = "n", bty = "n",...)
    #las = 1, cex.axis = 1.25, cex.lab = 1.5, 
}
.addPolygon <- function(object){
    newx <- getXcurve(object)
    newy <- getYcurve(object)
    bounds <- .confInt(getStdErr(object), getY(object), getFitValues(object), newy)
    xx <- c(newx, rev(newx))
    yy <- c(bounds$lo, rev(bounds$hi))
    polygon(xx, yy, border = NA, col = rgb(.8,.8,.8,.4))
}
.addEstim <- function(object, showEstim, unit, B, conf.level){
    stdErr <- getStdErr(object)
    estim <- .estimateRange(showEstim, stdErr, getPar(object)$params, B, object@useLog, conf.level)
    newx <- getXcurve(object)
    newy <- getYcurve(object)
    legend1 <- sprintf("IC%d : %s%s", showEstim*100, format(estim[2], scientific=TRUE, digits=2), unit)
    legend2 <- sprintf("[%s, %s]", format(estim[1], scientific=TRUE, digits=2), format(estim[3], scientific=TRUE, digits=2))
    legend(ifelse(newy[length(newy)]<newy[1], 'bottomleft', 'topleft'),
           legend = c(legend1, legend2), cex = 1.5, text.col = 'steelblue4', bty = 'n')
    
}
.addGOF <- function(object){
    gof <- format(getGoodness(object), digits=3, scientific = TRUE)
    newx <- getXcurve(object)
    newy <- getYcurve(object)
    legend(ifelse(newy[length(newy)]<newy[1], 'topright', 'bottomright'),
           legend = paste('Goodness of fit:', gof), bty = 'n', cex = 1.5)
}
.addPoints <- function(object, pcol, ...){
    x <- getX(object)
    y <- getY(object)
    points(x, y, col = pcol, pch = 19, ...)
    points(x, y, pch = 1)
}
.addCurve <- function(object, lcol, ...){
    x <- getXcurve(object)
    y <- getYcurve(object)
    lines(y ~ x, col=lcol, ...)
}
.SE <- function(x, y){
    .len <- function(x){ sum(!is.na(x)) }
    n <- by(y, x, .len)
    er <- by(y, x, sd, na.rm = TRUE)
    sEr <- as.vector(er/sqrt(n))
    sEr    
}
.addErr <- function(object, e, ...){
    x <- getX(object)
    y <- getY(object)
    my <- as.vector(by(y, x, mean, na.rm = TRUE))
    sEr <- .SE(x, y)
    points(unique(x), my, ...)
    e <- diff(range(x, na.rm = TRUE))/60
    segments(x0 = unique(x), y0 = my - sEr, y1 = my + sEr, ...)
    segments(x0 = unique(x) - e, x1 = unique(x) + e, y0 = my - sEr, ...)
    segments(x0 = unique(x) - e, x1 = unique(x) + e, y0 = my + sEr, ...)
}
.multiCurve <- function(modelList, showLegend, Cols,...){
    N <- length(modelList)
    Conc. <- do.call(c, lapply(modelList, function(tmp) getX(tmp) ))
    Resp <- do.call(c, lapply(modelList, function(tmp) getY(tmp) ))
#    args <- as.list( match.call() )

    if(is.null(Cols))
        Cols <- grey(seq_len(N)/(N+1))

    plot(range(Conc.),
        range(min(Resp, na.rm = TRUE), max(Resp, na.rm = TRUE) + .25),
        type = "n", bty = "n", ...)
    for(ii in seq_len(N)){
        tmp <- modelList[[ii]]
        .addErr(tmp, col = Cols[ii], ...)
        .addCurve(tmp, Cols[ii], ...)
    }

    if(showLegend)
        legend("top", legend = names(modelList), lwd = 2, pch = 19,
            col = Cols, ncol = length(modelList), bty = "n")

            #lwd = ifelse("lwd" %in% names(args), args$lwd, 2),
            #            pch = ifelse("pch" %in% names(args), args$pch, 19),        }
}

