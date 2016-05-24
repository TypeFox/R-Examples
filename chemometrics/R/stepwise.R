stepwise <-
function(formula,data, k, startM, maxTime=1800, direction="both", writeFile=FALSE, resname="stepres00", 
			maxsteps=500,...) {

# Stepwise regression, starting from the empty model, with scope to the full model
#
#require(pls)

    mf <<- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    y <- as.matrix(y)
    colnames(y) <- deparse(formula[[2]])
    x <- delintercept(model.matrix(mt, mf))
    x <- as.data.frame(x)

form <- function(a, namesX) {
    if(sum(a)==0) return(as.formula("y~1"))
    as.formula(paste("y~", paste(namesX[as.logical(a)], collapse="+"), sep=""))
}

    startTime <- as.integer(Sys.time())
    if(missing(k)) k <- log(nrow(x)) # BIC
    dat    <- data.frame(y, x)
    p      <- ncol(x)
    namesX <- names(x)

    bic      <- rep(Inf, maxsteps)
    usedTime <- rep(NA, maxsteps)
    models   <- matrix(NA, maxsteps, p)

    if(missing(startM)) startM <- rep(0, p)
    lmTmp  <- lm(y ~ 1, data=dat)
    lmTmp  <- update(lmTmp, form(startM, namesX))
    lmFull <- lm(form(rep(1,p), namesX), data=dat)

    lm0    <- lm(form(rep(0,p), namesX), data=dat)

    for(i in 1:maxsteps) {
        lmTmp <- step(object=lmTmp, scope=list(upper=lmFull, lower=lm0), direction=direction, trace=0, steps=1, k=k)
        
        usedTime[i] <- as.integer(Sys.time())-startTime
        bic[i]      <- AIC(lmTmp, k=k)
        models[i,]  <- as.numeric(is.element(el=namesX, set=names(lmTmp$coef)))
        if(i >= 2) {
            if(bic[i] >= bic[i-1]) break
        }
        if(as.integer(Sys.time())-startTime >= maxTime) break
        try(write.table(cbind(usedTime, bic, models), file=paste(resname, ".txt", sep=""), col.names=FALSE, row.names=FALSE))
    }
   if(any(is.na(models[,1]))) {
        tmp      <- which(is.na(models[,1]))[1]-1
        bic      <- bic[1:tmp]
        usedTime <- usedTime[1:tmp]
        models   <- models[1:tmp,]
    }
    if(length(bic) > 1) {
        if(bic[length(bic)] >= bic[length(bic)-1]) {
            tmp      <- length(bic)-1
            bic      <- bic[1:tmp]
            usedTime <- usedTime[1:tmp]
            models   <- models[1:tmp,]
        }
    }
    if (writeFile){
    write.table(cbind(usedTime, bic, models), file=paste(resname, ".txt", sep=""), col.names=FALSE, row.names=FALSE)
    }
    return(list(usedTime=usedTime, bic=bic, models=models))
}

