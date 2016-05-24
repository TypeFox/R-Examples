mhingebst_fit <- function(x, y, Fboost, yv, lq, b, learner, twinboost=FALSE, f.init=NULL, xselect.init=NULL, fixed.depth=TRUE, n.term.node=6, maxdepth=1, nu=0.1, df=4, inde=inde){
    p <- dim(x)[2]
    k <- length(table(y))
    ind <- coef <- rep(NA, k)
    if(learner=="tree")
        if(maxdepth == 1) xselect <- rep(NA, k)
        else xselect <- vector("list", k)
    mse.w <- coef.w <- matrix(NA,ncol=p, nrow=k)
    cor.w <- matrix(0, ncol=p, nrow=k)
    nob <- c(1:k)[-b]
    u <- -lq[,-b]*(sign(Fboost[,-b] - yv[,-b]) + 1)/2
    u1 <- lq[,b]*(sign(-apply(as.matrix(Fboost[,-b]), 1, sum) - yv[,b]) + 1)/2 ### sum-to-zero constraint
    
    u <- u + u1
    u <- as.matrix(u)
    tmp <- matrix(NA, nrow=dim(x)[1], ncol=k)
    tmp[,-b] <- u
    u <- tmp
    if(!twinboost) xselect.init <- 1:p
    ml.fit <- vector(mode = "list", length = k)
    pred.tr <- matrix(NA, nrow=dim(x)[1], ncol=k)
    coef0 <- matrix(NA, ncol(x), nrow=k)
    if(learner=="ls"){
        for (j in xselect.init){
            coef0[,j] <- 1/sum(x[,j]^2)*apply(as.matrix(x[,j] * u), 2, sum)
            for(i in nob)
                pred.tr[,i] <- x[,j] * coef0[i,j] 
            ss <- apply(pred.tr^2, 2, sum)
            if(twinboost){   
                a <- 2*colSums(u*pred.tr) - ss
                for(i in nob){
                    if(!all(pred.tr[,i] == 0))
                        cor.w[i,j] <- cov(f.init[,i],pred.tr[,i])/sqrt(sum(pred.tr[,i]^2))
                    mse.w[i,j] <- cor.w[i,j]^2 * a[i]
                }
            }
            else mse.w[,j] <- 2*colSums(u*pred.tr) - ss
        }
        for(i in nob)
            ind[i] <- which.max(mse.w[i,])
        for(i in nob){
            ml.fit[[i]] <- lm(u[,i] ~ x[, ind[i]] - 1)
            coef[i] <- coef(ml.fit[[i]])  ### this should be the same as above
        }
        xselect <- ind
    }
    else
        if(learner=="sm"){
            for(j in xselect.init){
###Twin L2 Boosting with genral weak learner, Buhlmann, page 8, step 4, in Twin boosting, improved feature selection andprediction
                for(i in nob){
                    if(length(unique(x[,j])) < 4)
                        res.fit <- lm(u[,i] ~ x[,j] - 1)
                    else res.fit <- smooth.spline(x=x[,j],y=u[,i],df=df)
                    pred.tr[,i] <- fitted(res.fit)
                }
                ss <- apply(pred.tr^2, 2, sum)
                if(twinboost){   
                    a <- 2*colSums(u*pred.tr) - ss
                    for(i in nob){
                        if(!all(pred.tr[,i] == 0))
                            cor.w[i,j] <- cov(f.init[,i],pred.tr[,i])/sqrt(sum(pred.tr[,i]^2))
                        mse.w[i,j] <- cor.w[i,j]^2 * a[i]
                    }
                }
        else mse.w[,j] <- 2*colSums(u*pred.tr) - ss
            }   
            for(i in nob)
                ind[i] <- which.max(mse.w[i,])
### this can be optimized with the previous results res.fit
            for(i in nob){
                if(length(unique(x[,ind[i]])) < 4)
                    ml.fit[[i]] <- lm(u[,i] ~ x[,ind[i]] - 1)
                else ml.fit[[i]] <- smooth.spline(x=x[,ind[i]],y=u[,i],df=df)
            }
            xselect <- ind
        }
    else
        if(learner=="tree"){
            cntrl <- rpart.control(maxdepth = maxdepth, #minsplit = nsample-1, #minbucket = 1,
                                   maxsurrogate = 0, maxcompete = 0, #usesurrogate=0
                                   cp = 0, xval = 0)
            cntrl0 <- rpart.control(maxdepth = 6, #minsplit = nsample-1, #minbucket = 1,
                                    maxsurrogate = 0, maxcompete = 0, #usesurrogate=0
                                    cp = 0, xval = 2)
            if(!twinboost){
                for(i in nob){
                    data.tr <- as.data.frame(cbind(u[,i],x)); 
                    colnames(data.tr) <- c("u",colnames(x))
                    if(fixed.depth)
                        ml.fit[[i]] <- rpart(u~.,data=data.tr,method="anova",control=cntrl)
                    else{
                        treefit <- rpart(u~.,data=data.tr,method="anova",control=cntrl0)
### find nsplit <= n.term.node
### (note: nsplit + 1 = number of terminal node = tree size) 
                        tmp <- which(treefit$cptable[,"nsplit"] <= n.term.node - 1)
### find the largest one among them
                        tmp1 <- tmp[length(tmp)]
### prune the tree to desired number of splits, which has the desirgd n.term.node 
                        ml.fit[[i]] <- prune(treefit, cp=treefit$cptable[,"CP"][which(treefit$cptable[,"nsplit"]==tmp1)])
                    }
                    labs <- rownames(ml.fit[[i]][["splits"]])
                    if(!is.null(labs))
                        xselect[[i]] <- which(colnames(x) %in% labs)
                }
            }
            else{
                for(j in 1:nrow(inde)){
                    for(i in nob){
                        if(maxdepth==1){
                            data.tr <- as.data.frame(cbind(u[,i],x[,j])); 
                            colnames(data.tr) <- c("u",colnames(x)[j])
                        }
              else{
                  warnings("Twin HingeBoost with base learner trees has not been fully tested\n")
                  data.tr <- as.data.frame(cbind(u[,i],x[,inde[j,]])); 
                  colnames(data.tr) <- c("u",colnames(x)[inde[j,]])
              }
###Twin L2 Boosting with genral weak learner, Buhlmann, page 127, step 4, in Twin boosting, improved feature selection and prediction, Statistics and Computing (2007) Volume: 20, Issue: 2, Pages: 119-138
                        if(fixed.depth)
                            res.fit <- rpart(u~.,data=data.tr,method="anova",control=cntrl)
                        else{
                            treefit <- rpart(u~.,data=data.tr,method="anova",control=cntrl0)
### find nsplit <= n.term.node
### (note: nsplit + 1 = number of terminal node = tree size) 
                            tmp <- which(treefit$cptable[,"nsplit"] <= n.term.node - 1) 
### find the largest one among them
                            tmp1 <- tmp[length(tmp)]
### prune the tree to desired number of splits, which has the desirgd n.term.node 
                            res.fit <- prune(treefit, cp=treefit$cptable[,"CP"][which(treefit$cptable[,"nsplit"]==tmp1)])
                        }
                        pred.tr[,i] <- predict(res.fit)
                    }  
                    ss <- apply(pred.tr^2, 2, sum)
                    a <- 2*colSums(u*pred.tr) - ss
                    for(i in nob){
                        if(!all(pred.tr[,i] == 0))
                            cor.w[i,j] <- cov(f.init[,i],pred.tr[,i])/sqrt(sum(pred.tr[,i]^2))
                        mse.w[i,j] <- cor.w[i,j]^2 * a[i]
                    }
                }
                for(i in nob)
                    ind[i] <- which.max(mse.w[i,])
### this can be optimized with the previous results res.fit
                if(maxdepth==1){
                    for(i in nob){
                        data.tr <- as.data.frame(cbind(u[,i],x[,ind[i]])); 
                        colnames(data.tr) <- c("u",colnames(x)[ind[i]])
                        treefit <- rpart(u~.,data=data.tr,method="anova",control=cntrl0)
                        if(fixed.depth)
                            ml.fit[[i]] <- rpart(u~.,data=data.tr,method="anova",control=cntrl)
                        else{
### find nsplit <= n.term.node
### (note: nsplit + 1 = number of terminal node = tree size) 
                            tmp <- which(treefit$cptable[,"nsplit"] <= n.term.node - 1) 
### find the largest one among them
                            tmp1 <- tmp[length(tmp)]
### prune the tree to desired number of splits, which has the desirgd n.term.node 
                            ml.fit[[i]] <- prune(treefit, cp=treefit$cptable[,"CP"][which(treefit$cptable[,"nsplit"]==tmp1)])
                        }
                    }
                    xselect[[i]] <- ind
                }
                else {
                    tmp1 <- NULL
                    for(i in nob){
                        data.tr <- as.data.frame(cbind(u[,i],x[,ind[i,]])); 
                        colnames(data.tr) <- c("u",colnames(x)[ind[i,]])
                        if(fixed.depth)
                            ml.fit[[i]] <- rpart(u~.,data=data.tr,method="anova",control=cntrl)
                        else{
                            treefit <- rpart(u~.,data=data.tr,method="anova",control=cntrl0)
### find nsplit <= n.term.node
### (note: nsplit + 1 = number of terminal node = tree size) 
                            tmp <- which(treefit$cptable[,"nsplit"] <= n.term.node - 1) 
### find the largest one among them
                            tmp1 <- tmp[length(tmp)]
### prune the tree to desired number of splits, which has the desirgd n.term.node 
                            ml.fit[[i]] <- prune(treefit, cp=treefit$cptable[,"CP"][which(treefit$cptable[,"nsplit"]==tmp1)])
                        } 
                        tmp <- ml.fit[[i]]$frame$var[ml.fit[[i]]$frame$var%in%colnames(x)]
                        tmp1 <- c(tmp1, tmp)
                    }
                    tmp1 <- unique(tmp1)
                    if(length(tmp1)!=0)
                        xselect[[i]] <- as.character(tmp1)
                }
            }          
        }
### update prediction
    for(i in nob){
        if(learner=="sm")
            Fboost[,i] <- Fboost[,i] + nu * fitted(ml.fit[[i]])
        else
            Fboost[,i] <- Fboost[,i] + nu * predict(ml.fit[[i]])
    }
    Fboost[,b] <- -apply(as.matrix(Fboost[,-b]), 1, sum) ### sum-to-zero
### empirical loss
    risk <- loss.mhingebst(y, Fboost, k)
    ensemble <- xselect
    return(list(b=b, Fboost=Fboost, ens=ml.fit, risk=risk, xselect=xselect, coef=coef))
} 

loss.mhingebst <- function(y, f, k, type=c("total","all"), cost=NULL){
    type <- match.arg(type)
    v <- matrix(rep(-1/(k-1), k*k), ncol=k)
    diag(v) <- 1
    yv <- v[y,]
    QQ <- matrix(rep(1, k*k), ncol=k)
    diag(QQ) <- 0
    lq <- QQ[y, ]
    los <- mapply(function(x) max(x, 0), f-yv)
    los <- matrix(los, byrow=FALSE, ncol=k)
    tmp <- lq * los 
    if(type=="total")
        return(sum(tmp)/length(y))
    else return(tmp/length(y))
}

#######################################################################################################################################################
mhingebst <- function(x,y, cost=NULL, family = c("hinge"), ctrl = bst_control(), control.tree=list(fixed.depth=TRUE, n.term.node=6, maxdepth=1), learner=c("ls", "sm", "tree")){
    call <- match.call()
    if(any(y < 1)) stop("y must > 0 \n")
    family <- match.arg(family)
    learner <- match.arg(learner)
    x <- as.matrix(x)
    if(learner == "tree" && is.null(colnames(x)))
        colnames(x) <- paste("x", 1:ncol(x), sep = "")
    mstop <- ctrl$mstop
    nu <- ctrl$nu
    twinboost <- ctrl$twinboost
    f.init <- ctrl$f.init
    xselect.init <- ctrl$xselect.init
    center <- ctrl$center
    trace <- ctrl$trace
    numsample <- ctrl$numsample
    df <- ctrl$df
    if(twinboost && (is.null(f.init) | is.null(xselect.init)))
        stop("Twin boosting requires initial function estimates and variable selected in the first round\n")
    nsample <- dim(x)[1]
    p <- dim(x)[2]
    if(learner == "tree" && p > 10 && twinboost && control.tree$maxdepth >= 2 && is.null(numsample))
        stop("for large p and degree >=2, random sample is suggested\n")   
    if(center){
        one <- rep(1,nrow(x))
        meanx <- drop(one %*% as.matrix(x))/length(y)
        x <- scale(x, meanx, FALSE) # centers x
    }
### multi-class coding, cf, Lee et al (2004), JASA, page 69. 
    k <- length(table(y))
    v <- matrix(rep(-1/(k-1), k*k), ncol=k)
    diag(v) <- 1
    yv <- v[y,]
    QQ <- matrix(rep(1, k*k), ncol=k)
    diag(QQ) <- 0
    lq <- QQ[y, ]
    oldx <- x; one <- rep(1,length(y))
    ens <- array(list(), c(mstop, k))
    Fboost <- offset <- pred.val <- 0
    Fboost <- matrix(Fboost, nrow=length(y), ncol=k)
    m <- 1
    baseclass <- risk <- rep(NA,mstop)
    sse <- minid <- rep(NA,ncol(x))
    coef <- matrix(NA, ncol=k, nrow=mstop)
    xselect <- vector("list", mstop)
    ind <- rep(NA, k)
    fixed.depth <- control.tree$fixed.depth
    if(is.null(fixed.depth)) fixed.depth <- TRUE
    maxdepth <- control.tree$maxdepth
    n.term.node <- control.tree$n.term.node

    if(learner=="tree" && twinboost){
        if(maxdepth==1){
            p1 <- ifelse(!twinboost, p,length(xselect.init))
            xselect.new <- xselect.init
            inde <- as.matrix(1:p1, ncol=1)
        }
    else if(maxdepth==2){
        if(missing(vimp.init)) vimp.init <- rep(1,length(xselect.init))
        if(p > 10){
            inde <- NULL
            for (i in 1:numsample)
                inde <- rbind(inde, sample(xselect.init,maxdepth,prob=vimp.init[vimp.init>0]))
        }  
        else
            inde <- t(combn(xselect.init,2)) #generate interactions
        xselect.new <- xselect.init[-(length(xselect.init))]
    }
    }
    while (m <= mstop){
        tmp <- rep(NA, k); res <- vector("list", k)
        for(i in 1:k){
	    res[[i]] <- mhingebst_fit(x=x, y=y, Fboost=Fboost, yv=yv, lq=lq, b=i, learner=learner, twinboost=twinboost, f.init=f.init, xselect.init=xselect.init, fixed.depth=fixed.depth, n.term.node=n.term.node, maxdepth=maxdepth, nu=nu, df=df, inde=inde)
            tmp[i] <- res[[i]]$risk
        }
        optb <- which.min(tmp)
        Fboost <- res[[optb]]$Fboost
        for(i in 1:k)
            ens[[m,i]] <- res[[optb]]$ens[i]
        risk[m] <- res[[optb]]$risk
        xselect[[m]] <- (res[[optb]]$xselect)
        coef[m,] <- res[[optb]]$coef
        baseclass[m] <- optb 
        if(trace){
            if(m %% 10==0) cat("m=", m, "  risk = ", risk[m], "\n")
        } 
        m <- m + 1
    }
    ensemble <- xselect
    xselect <- sort(unique(unlist(xselect)))
    xselect <- xselect[!is.na(xselect)]
    RET <- list(y=y,x=oldx, family = family, learner=learner, k=k, yhat=Fboost, offset=offset, ens=ens, control.tree=control.tree, risk=risk, ctrl = list(center=center, mstop=mstop,nu=nu, df=df), xselect=xselect, coef = coef, ensemble=ensemble, baseclass=baseclass)
    RET$call <- call
    class(RET) <- "mhingebst"
    return(RET)
}

predict.mhingebst <- function(object, newdata=NULL, newy=NULL, mstop=NULL, type=c("response", "class", "loss", "error"), ...){
    if(is.null(mstop))
        mstop <- object$ctrl$mstop
    else if(mstop > object$ctrl$mstop)
        stop("mstop must be equal or smaller than the one used for estimation ", object$ctrl$mstop)
                                        #  if((type=="loss" || type=="error") && (is.null(newdata) || is.null(newy)))
                                        #    stop("For estimation of loss or error, both newdata and newy are needed\n")
    type <- match.arg(type)
    one <- rep(1,nrow(object$x))
    x <- object$x
    y <- object$y
    if(is.null(newdata) && is.null(newy))
        ynow <- y
    else ynow <- newy
    if(!missing(newdata)){
        if(object$ctrl$center){
            meanx <- drop(one %*% as.matrix(x))/nrow(x)
            newdata <- scale(newdata, meanx, FALSE) # centers x
        }
    }
    learner <- object$learner
    ens <- object$ens
    k <- object$k
    nu <- object$ctrl$nu
    baseclass <- object$baseclass
    if(missing(newdata)) p <- dim(x)[1]
    else{ 
        if(!missing(newy))
            if(dim(newdata)[1] != length(newy))
                stop("Number of rows of newdata is different from length of newy\n")
        newdata <- as.matrix(newdata)
        p <- dim(newdata)[1]
    }
    risk <- rep(NA, mstop)
    lp <- matrix(object$offset, ncol=k, nrow=p)
    if (is.matrix(newdata)) newdata <- as.data.frame(newdata)
    for(m in 1:mstop){
        nob <- c(1:k)[-baseclass[m]]
        if(missing(newdata)){
            for(i in nob)
                if(learner=="tree")
                    lp[,i] <- lp[,i] + nu*predict(ens[[m,i]][[1]])
                else lp[,i] <- lp[,i] + nu*fitted(ens[[m,i]][[1]])
        }
        else{
            for(i in nob)
                if(learner=="tree")  
                    lp[,i] <- lp[,i] + nu*predict(ens[[m,i]][[1]], newdata = newdata)
                else if(learner=="sm"){
                    if(length(unique(x[,object$ensemble[[m]][i]])) < 4)
                    {      lp[,i] <- lp[,i] + nu * coef(ens[[m, i]][[1]])* newdata[, object$ensemble[[m]][i]]
                    }         else  
                        lp[,i] <- lp[,i] + nu * predict(ens[[m, i]][[1]], newdata[, object$ensemble[[m]][i]])$y
                }        else if(learner=="ls")
                    lp[,i] <- lp[,i] + nu * object$coef[m, i] * newdata[, object$ensemble[[m]][i]]
        }
        lp[,baseclass[m]] <- - rowSums(as.matrix(lp[,-baseclass[m]])) 
        if(type=="loss"){
            risk[m] <- loss.mhingebst(ynow, lp, k)
        }
        else if(type == "error"){
            tmp <- apply(lp, 1, which.max)
            risk[m] <- (mean(ynow != tmp))
        }
    }
    if(type=="loss" || type=="error")
        lp <- risk 
    else if(type == "class")
        lp <- apply(lp, 1, which.max)
    return(drop(lp))
}

"cv.folds" <-
    function(n, folds = 10)
    {
        split(sample(1:n), rep(1:folds, length = n))
    }
"cv.mhingebst" <-
    function(x, y, balance=FALSE, K = 10, cost = NULL, family = "hinge", learner = c("tree","ls", "sm"), ctrl = bst_control(), type = c("loss", "error"), plot.it = TRUE, main = NULL, se = TRUE, n.cores=2, ...)
    {
        call <- match.call()
        family <- match.arg(family)
        learner <- match.arg(learner)
        type <- match.arg(type)
        family <- match.arg(family)
        learner <- match.arg(learner)
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
            fit <- mhingebst(x[ - omit,,drop=FALSE  ], y[ - omit], cost = cost, family = family, learner = learner, ctrl = ctrl.cv, ...)
            predict.mhingebst(fit, newdata = x[omit,  ,drop=FALSE], newy=y[ omit], mstop = mstop, type=type)
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

"plotCVbst" <-
    function(cv.bst.object,se=TRUE,ylab=NULL, main=NULL, width=0.02, col="darkgrey", ...){
        mstop <- cv.bst.object$mstop
        cv <- cv.bst.object$cv
        cv.error <- cv.bst.object$cv.error
        family <- cv.bst.object$family
        if(is.null(ylab)){
            ylab <- "Cross-validation loss values"
        }
        plot(mstop, cv, type = "b", xlab = "Iteration", ylab= ylab, ylim = range(cv, cv + cv.error, cv - cv.error), main=main)
        if(se)
            error.bars(mstop, cv + cv.error, cv - cv.error,
                       width = width, col=col)
                                        #                 width = 1/length(fraction), col=col)
                                        #    detach(cv.bst.object)

        invisible()
    }

"error.bars" <-
    function(x, upper, lower, width=0.02, col="darkgrey", ...)
    {
        xlim <- range(x)
        barw <- diff(xlim) * width
        segments(x, upper, x, lower, col=col, ...)
        segments(x - barw, upper, x + barw, upper, col=col, ...)
        segments(x - barw, lower, x + barw, lower, col=col, ...)
        range(upper, lower)
    }

print.mhingebst <- function(x, ...) {

    cat("\n")
    cat("\t Models Fitted with Gradient Boosting\n")
    cat("\n")
    if (!is.null(x$call))
        cat("Call:\n", deparse(x$call), "\n\n", sep = "")
    show(x$family)
    cat("\n")
    if(!is.null(x$ctrl$twinboost))
        cat("Twin boosting", "\n")
    cat("Base learner: ", x$learner, "\n")
    cat("Number of boosting iterations: mstop =", x$ctrl$mstop, "\n")
    cat("Step size: ", x$ctrl$nu, "\n")
    cat("Offset: ", x$offset, "\n")
    cat("\n")
    if(x$learner=="ls"){
        cat("Coefficients: \n")
        cf <- coef(x)
        print(cf)
        cat("\n")
    }
    if(x$learner=="sm")
        cat("Degree of freedom used is: ", x$ctrl$df, "\n")
    invisible(x)
}

fpartial.mhingebst <- function (object, mstop=NULL, newdata=NULL)
{   
    if(is.null(mstop))
        mstop <- object$ctrl$mstop
    else if(mstop > object$ctrl$mstop)
        stop("mstop must be equal or smaller than the one used for estimation ", object$ctrl$mstop)
    one <- rep(1,nrow(object$x))
    x <- object$x
    if(is.null(newdata))
        newdata <- x
    if(!missing(newdata)){
        if(object$ctrl$center){
            meanx <- drop(one %*% as.matrix(x))/nrow(x)
            newdata <- scale(newdata, meanx, FALSE) # centers x
        }
    }
    ens <- object$ens
    k <- object$k
    nu <- object$ctrl$nu
    baseclass <- object$baseclass
    if(missing(newdata)) p <- dim(x)[1]
    else{
        newdata <- as.matrix(newdata)
        p <- dim(newdata)[1]
    }

    lp <- vector("list", k)
    for(i in 1:k)
        lp[[i]] <- matrix(0, ncol = NCOL(x), nrow = NROW(x))
    if (is.matrix(newdata)) newdata <- as.data.frame(newdata)
    for(m in 1:mstop){
        nob <- c(1:k)[-baseclass[m]]
        if(object$learner=="tree"){
            for(i in nob)
                xselect <- object$ensemble[[m]][i]
            lp[[i]][,xselect] <- lp[[i]][,xselect] + nu*predict(object$ens[[m,i]][[1]], newdata = newdata)
        }
        else if(object$learner=="sm"){
            for(i in nob)
                xselect <- object$ensemble[[m]][i]
            lp[[i]][,xselect] <- lp[[i]][,xselect] + nu * predict(object$ens[[m, i]][[1]], newdata[, object$ensemble[[m]][i]])$y
        }
        else if(object$learner=="ls"){
            for(i in nob)
                xselect <- object$ensemble[[m]][i]
            lp[[i]][,xselect] <- lp[[i]][,xselect] + nu * object$coef[m, i] * newdata[, object$ensemble[[m]][i]]
        }
        tmp <- 0
        for(i in nob)
            tmp <- tmp + lp[[i]]
        lp[[baseclass[m]]] <- -tmp
    }
    lp 
}

nsel <- function(object, mstop){
	#if(!class(object) %in% c("mhingebst", "mbst", "rmbst"))
	#stop("object is not a correct class\n")
    np <- rep(NA, mstop)
    tmp <- NULL
    for(i in 1:mstop){
        tmp <- c(tmp, unlist(object$ensemble[[i]]))
        np[i] <- length(unique(subset(tmp, !is.na(tmp))))
    }
    np
}

