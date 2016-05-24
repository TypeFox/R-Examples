rbstpath <- function(x, y, rmstop=seq(40, 400, by=20), ctrl=bst_control(), del=1e-16, ...){
    fit <- vector("list", length(rmstop))
    for(jk in 1:length(rmstop)){
        if(jk==1) fk <- NULL
        else fk <- fit[[jk-1]]$yhat
        ctrl$mstop <- rmstop[jk]
        ctrl$fk <- fk
        fit[[jk]] <- rbst(x, y, ctrl = ctrl, del=del, ...)
    }
    fit
}

rbst <- function(x,y, cost=0.5, rfamily=c("tgaussian", "thuber", "thinge", "tbinom", "binomd", "texpo", "tpoisson"), ctrl = bst_control(), control.tree=list(maxdepth=1), learner=c("ls", "sm", "tree"), del=1e-10){
    call <- match.call()
    learner <- match.arg(learner)
    rfamily <- match.arg(rfamily)
    if(rfamily %in% c("thinge", "tbinom", "binomd", "texpo"))
        if(!all(names(table(y)) %in% c(1, -1)))
            stop("response variable must be 1/-1 for family ", rfamily, "\n")
    s <- ctrl$s
    sh <- ctrl$sh
    q <- ctrl$q
    qh <- ctrl$qh
    if(is.null(ctrl$s) && rfamily=="tgaussian")
        sa <- TRUE  ### adaptive s
    else sa <- FALSE
                                        #if(is.null(s) && is.null(q)) stop("s or q must be provided\n")
    if(!is.null(q)){
        if(q < 0 || q > 1)
            stop("proportion of outliers must between 0 and 1\n")
        if(rfamily=="thuber"){
            if(is.null(qh)) stop("qh must be provided\n")
            if(q <= qh) stop("q should be larger than qh\n")
        }
    }
    else if(!is.null(s)){
        if(rfamily=="tbinom")
            if(s > 0) stop("s must be non-negative for rf='tbinom'\n")
        if(rfamily=="binomd")
            if(s <= 0) stop("s must be positive for rf='binomd'\n")
        if(rfamily=="thuber"){
            if(is.null(sh)) stop("sh must be provided\n")
            if(s <= sh) stop("s should be larger than sh\n")
        }
    }
    fk <- ctrl$fk
    if(is.null(s)){
        if(rfamily=="tgaussian" || rfamily=="thuber"){
            stop("how to find s is not implemented\n")
        }
        else
            s <- switch(rfamily,       
                        "thinge"= -1,
                        "tbinom"= -log(3),
                        "binomd"= log(4),
                        "texpo"= log(0.5),
                        "tpoisson"= 5*mean(y))
    }
    famtype <- switch(rfamily,
                      "tgaussian"="tgaussianDC",
                      "thuber"="thuberDC",
                      "thinge"="thingeDC",
                      "tbinom"="tbinomDC",
                      "binomd"="binomdDC",
                      "texpo"="texpoDC",
                      "tpoisson"="tpoissonDC",
                      )
    ctrl$s <- s
    ctrl$sh <- sh
    iter <- ctrl$iter
    trace <- ctrl$trace
    if(trace) cat("\ngenerate initial values\n") 
    ly <- ifelse(y==1, 1-cost, cost)
### initiate values are important, best with nonrobust intercept models
### may need to upgrade for other nonrobust methods
    if(is.null(fk)){
        bsttype <- switch(rfamily,
                                        # "tgaussian"="gaussian",
                          "tgaussian"="huber",
                          "thuber"="huber",
                          "thinge"="hinge",
                          "tbinom"="binom",
                          "binomd"="binom",
                          "texpo"="expo",
                          "tpoisson"="poisson",
                          )
        RET <- bst(x, y, cost=cost, family=bsttype, ctrl = bst_control(mstop=1), control.tree=control.tree, learner=learner)
    }
    else {
        RET <- NULL
        RET$yhat <- fk
    }
    los <- loss(y, f=RET$yhat, cost, family = rfamily, s=ctrl$s, sh=ctrl$sh, fk=fk)
    d1 <- 10 
    k <- 1
    if(trace) {
        cat("\nrobust boosting ...\n")
        cat("\ninitial loss", mean(los), "\n")
    }
    los <- rep(NA, iter)
    while(d1 > del && k <= iter){
        ctrl$fk <- RET$yhat
        if(sa){
### fk is the previous fitted f
            if(is.null(ctrl$fk)) fk <- 0
            ctrl$s <- quantile(gaussloss(y, ctrl$fk), 0.5)   ### adaptive s,  test program
        }
        RET <- bst(x, y, cost=cost, family=famtype, ctrl = ctrl, control.tree=control.tree, learner=learner)
        los[k] <- mean(loss(y, f=RET$yhat, cost, family = rfamily, s=ctrl$s, sh=ctrl$sh, fk=NULL))
        d1 <- sum((RET$yhat - ctrl$fk)^2)/sum(ctrl$fk^2)
        if(trace) cat("\niteration", k, ": relative change of fk", d1, ", robust loss value", los[k], "\n") 
	if(k > 1){
            if(los[k] > los[k-1])
                k <- iter
        }
        k <- k + 1
    }
    RET$x <- x
    RET$y <- y
    RET$call <- call
    RET$cost <- cost
    RET$rfamily <- RET$family <- rfamily
    RET
}

"cv.rbst" <-
    function(x, y, K = 10, cost = 0.5, rfamily = c("tgaussian", "thuber", "thinge", "tbinom", "binomd", "texpo", "tpoisson"),
             learner = c("ls", "sm", "tree"), ctrl = bst_control(), type = c("loss", "error"), plot.it = TRUE, main = NULL, se = TRUE, n.cores=2, ...)
{
    call <- match.call()
    rfamily <- match.arg(rfamily)
    type <- match.arg(type)
    learner <- match.arg(learner)
    mstop <- ctrl$mstop
    nu <- ctrl$nu
    df <- ctrl$df
    twinboost <- ctrl$twinboost
    twintype <- ctrl$twintype
    trace <- ctrl$trace
    s <- ctrl$s
    sh <- ctrl$sh
    fk <- ctrl$fk
    ctrl.cv <- ctrl
    if(rfamily %in% c("tgaussian", "thuber", "tpoisson") && type =="error") stop("misclassification is Not applicable for rfamily ", rfamily, "\n")
    all.folds <- cv.folds(length(y), K)
    fraction <- 1:mstop
    registerDoParallel(cores=n.cores)
    i <- 1  ###needed to pass R CMD check with parallel code below
    residmat <- foreach(i=seq(K), .combine=cbind) %dopar% {
                                        #for(i in seq(K)) {
        omit <- all.folds[[i]]
        if(ctrl$twinboost){
            ctrl.cv$f.init <- ctrl$f.init[ - omit]
        }
	fit <- rbst(x[ - omit,,drop=FALSE  ], y[ - omit], cost = cost, rfamily = rfamily, learner = learner, ctrl = ctrl.cv, ...)
        predict(fit, newdata = x[omit,  ,drop=FALSE], newy=y[omit], mstop = mstop, type=type)

    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    object<-list(residmat = residmat, mstop = fraction, cv = cv, cv.error = cv.error, rfamily=rfamily)
    if(plot.it){
        if(type=="loss") ylab <- "Cross-validation loss values"
        else  if(type=="error") ylab <- "Cross-validation misclassification errors"
        plotCVbst(object,se=se, ylab=ylab, main=main)
    }
    invisible(object)
}

