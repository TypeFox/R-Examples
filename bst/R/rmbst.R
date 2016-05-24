rmbst <- function(x,y, cost=0.5, rfamily="thinge", ctrl = bst_control(), control.tree=list(maxdepth=1), learner=c("ls", "sm", "tree"), del=1e-10){
    call <- match.call()
    learner <- match.arg(learner)
    rfamily <- match.arg(rfamily)
    s <- ctrl$s
    if(!is.null(s)){
        if(s < 0) stop("s must be >= 0\n")
    }
    fk <- ctrl$fk
    if(is.null(s)){
        s <- switch(rfamily,       
                    "thinge"= 1)
    }
    famtype <- switch(rfamily,
                      "thinge"="thingeDC"
                      )
    ctrl$s <- s
    iter <- ctrl$iter
    trace <- ctrl$trace
    if(trace) cat("\ngenerate initial values\n") 
### initiate values are important, best with nonrobust intercept models
### may need to upgrade for other nonrobust methods
    if(is.null(fk)){
        bsttype <- switch(rfamily,
                          "thinge"="hinge2"
		          )
        RET <- mbst(x, y, cost=cost, family=bsttype, ctrl = bst_control(mstop=1), control.tree=control.tree, learner=learner)
    }
    else {
        RET <- NULL
        RET$yhat <- fk
    }
    los <- loss.mbst(y, f=RET$yhat, fk=fk, s=ctrl$s, k=RET$k, family = rfamily, cost=cost)
    d1 <- 10 
    k <- 1
    if(trace) {
        cat("\nrobust boosting ...\n")
        cat("\ninitial loss", mean(los), "\n")
    }
    los <- rep(NA, iter)
    while(d1 > del && k <= iter){
        ctrl$fk <- RET$yhat
        RET <- mbst(x, y, cost=cost, family=famtype, ctrl = ctrl, control.tree=control.tree, learner=learner)
	los[k] <- mean(loss.mbst(y, f=RET$yhat, fk=NULL, s=ctrl$s, k=RET$k, family = rfamily, cost=cost))
                                        #difference of convex function is linearly majorized, thus tmp1 - los[k] >= 0. cf Wang (2015)
	if(trace){
            tmp <- matrix(NA, nrow=length(y), ncol=RET$k)
            f <- RET$yhat; fk <- ctrl$fk
            for(j in 1:RET$k)
                tmp[,j] <- (y!=j)*(mapply(function(x) max(x, 0), f[,j]+1) - mapply(function(x) max(x, 0), fk[,j]-s)- (f[,j]-fk[,j])*(fk[,j] >= s))
            tmp1 <- sum(tmp)/length(y)
            cat("difference of convex function is linearly majorized, returning a non-negative number", tmp1-los[k], "\n")
        }
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

"cv.rmbst" <-
    function(x, y, balance=FALSE, K = 10, cost = NULL, rfamily = "thinge", learner = c("tree","ls", "sm"), ctrl = bst_control(), type = c("loss", "error"), plot.it = TRUE, main = NULL, se = TRUE, n.cores=2, ...)
{
    call <- match.call()
    rfamily <- match.arg(rfamily)
    learner <- match.arg(learner)
    type <- match.arg(type)
    mstop <- ctrl$mstop
    nu <- ctrl$nu
    df <- ctrl$df
    twinboost <- ctrl$twinboost
    trace <- ctrl$trace
    ctrl.cv <- ctrl
    if(balance)
        all.folds <- balanced.folds(y, K)
    else all.folds <- cv.folds(length(y), K)
    fraction <- seq(mstop)
    registerDoParallel(cores=n.cores)
    i <- 1  ###needed to pass R CMD check with parallel code below
    residmat <- foreach(i=seq(K), .combine=cbind) %dopar% {
        omit <- all.folds[[i]]
        if(ctrl$twinboost)
            ctrl.cv$f.init <- ctrl$f.init[ - omit, ]
        fit <- rmbst(x[ - omit,,drop=FALSE  ], y[ - omit], cost = cost, rfamily = rfamily, learner = learner, ctrl = ctrl.cv, ...)
	predict.mbst(fit, newdata = x[omit,  ,drop=FALSE], newy=y[ omit], mstop = mstop, type=type)
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    object<-list(residmat=residmat, mstop = fraction, cv = cv, cv.error = cv.error)
    if(plot.it){
        if(type=="loss") ylab <- "Cross-validation loss values"
        else  if(type=="error") ylab <- "Cross-validation misclassification errors"
        plotCVbst(object,se=se, ylab=ylab, main=main)
    }
    invisible(object)
}

