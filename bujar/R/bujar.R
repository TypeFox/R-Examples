bujar <- function(y, cens, x, valdata = NULL, degree = 1, learner = "linear.regression", center=TRUE, mimpu = NULL, iter.bj = 20, max.cycle = 5, nu = 0.1, mstop = 50, twin = FALSE, mstop2 = 100, tuning = TRUE, cv = FALSE, nfold = 5, method = "corrected", vimpint = TRUE, gamma=3, lambda=NULL, whichlambda=NULL, lamb = 0, s = 0.5, nk = 4, wt.pow = 1, theta = NULL, rel.inf = FALSE, tol = .Machine$double.eps, n.cores=2, rng=123, trace = FALSE){
    call <- match.call()
    if(learner == "acosso")
        stop("learner = 'acosso' is no longer supported, see NEWS\n")
    if(!learner%in%c("linear.regression","mars","pspline","tree","acosso","enet", "enet2", "mnet","snet")) stop(sQuote("weak learner"), learner, " is not implemented")
    if(!is.null(valdata))
        if((dim(x)[2]+2) !=dim(valdata)[2])
            stop("check the dimension of x and valdata\n")
    if(learner=="acosso" && wt.pow < 0) stop(Quote("wt.pow should be > 0"))
    if(learner=="acosso" && cv) stop(Quote("if wt.pow is chosen by cross-validation, then BJ estimator is not stable, thus stop"))### check ?
    if(cv && nfold < 1)
        stop(sQuote("Number of CV folds"), " are less than 1")
    if(!all(unique(cens)%in%c(0,1)))
        stop(sQuote("censoring indicator"), " are not 0 or 1")
    if(!all(unique(valdata[,2]%in%c(0,1))))
        stop(sQuote("censoring indicator"), " are not 0 or 1")
    if(iter.bj < 2)
        stop(sQuote("iter.bj"), " should be greater than 1")
    if(max.cycle < 1)
        stop(sQuote("max.cycle"), " should be greater than 1")
    if(!learner %in% c("tree","mars","acosso") && degree !=1)
        stop("Not implemented for this degree with learner ",learner, "\n")
    if(learner=="tree" && degree > 1 && twin) stop(sQuote("learner"), learner, sQuote("degree"), degree, sQuote("twin"), twin, "Not implemented\n")
### match names to be used in other packages
    if(learner=="pspline") l2 <- "sm"
    else if(learner=="linear.regression") l2 <- "ls"
    else if(learner=="enet2") l2 <- "enet"
    else l2 <- learner
    xbar <- colMeans(x)                                    #componentwise linear least square
    f <- 0
    mse.bj.val <- nz.bj <- NA
    ynew <- y; ynew.m <- vim <- NULL; Fboost <- NA
    p <- ncol(x)
    sse <- rep(NA, p) #SSE for each univariate covariate
    res <- ystar <- matrix(NA, length(y), p) #residual matrix, one column for corresponding univariate covariate
    coef <- matrix(NA, 2, p)
    b <- matrix(NA,iter.bj+max.cycle,p)
    ybst <- matrix(NA, iter.bj+max.cycle, length(y)) #changed Apr 2, 2008
    ybstdiff <- rep(NA, iter.bj+max.cycle)       #changed Apr 2, 2008
    fnorm2 <- mseun <- ybstcon <- rep(NA, iter.bj+max.cycle)       #changed Apr 2, 2008
    mselect <- rep(NA,iter.bj+max.cycle)
    if(!tuning && learner%in%c("linear.regression","pspline","tree")){
        if(length(mstop) > 1){
            if(length(mstop) !=length(mselect))
                stop(sQuote("mstop must be one number or have the length iter.bj+max.cycle for a boosting learner"))
            else mstopRep <- mstop
        }
    else mstopRep <- rep(mstop, iter.bj+max.cycle)
    }
    else if(tuning && learner%in%c("linear.regression","pspline","tree")){
        if(length(mstop) > 1)
            stop(sQuote("mstop must be one number if tuning=TRUE for a boosting learner"))
    }    	
    tuningSwitch <- TRUE
    ydiff <- 100
    k <- 1; kt <- 1
    mse.bj <- pred.bj <- NA
    if(trace) cat("\nBJ with",learner,"\n")
    nm <- dim(x)
    n <- nm[1]
    if(trace){
        cat("\nNumber of observations:",n)
        cat("\nNumber of covariates:",nm[2],"\n")
    }
    one <- rep(1, n)
    normx <- rep(1,dim(x)[2])
    cycleperiod <- 0
    nonconv <- FALSE
    mse <- rep(NA,iter.bj+max.cycle)
    nz.bj.iter <- rep(NA,iter.bj+max.cycle)
                                        #mstop.mars <- mstop
    while (ydiff > tol && k <= iter.bj+max.cycle){
        oldydiff <- ydiff
                                        #when k=1, only for initializaiton purpose, find BJ estimator and predicted y-value                                      #determine ynew
### Initializtion of BJ estimator with different (3) methods
        if(is.null(mimpu) && k==1) ynew <- y #this is to say we begin with initial estimator of beta=0
        else if(mimpu==FALSE && k==1){
            for (i in 1:p){
                res.des <- try(bj(Surv(ynew,cens) ~ x[,i], link="identity",control=list(trace=FALSE)))
                ystar[,i] <- res.des$y.imputed #imputed y value
                res[,i] <- y - predict(res.des)
                sse[i] <- sum(res[,i][cens==1]^2) # SSE for uncensored
            }
            minid <- which.min(sse) #find the covariate with smallest loss
            ynew <- ystar[,minid] #construct the new 'response'
            cat("\nBJ step k=",k,"\n","\nInitial MSE for uncensored observations=", sse[minid]/sum(cens==1), "\n\n")
        }
        else{
            if(mimpu==TRUE && k==1){
                ynew <- bjboost.fit(cbind(y,cens),rep(0,length(y)))$y.imputed #imputing without covariates
                                        #        k <- iter.bj
            }
      else {
          ynew <- bjboost.fit(cbind(y,cens),Fboost)$y.imputed}
        }
### End of initialization with the imputed ynew replacing the original y values
        dat1 <- as.data.frame(cbind(ynew,x))
### Different methods for BJ estimation
        x <- as.matrix(x) 
        if(learner%in%c("linear.regression","pspline","tree")){
            if(!tuning) mselect.now <- mstopRep[k]  ### boosting iteration is fixed
            else if(k==1) mselect.now <- mstop ### initial value for boosting and will be updated if the tuning parameter is updated
        }
	else mselect.now <- NULL
	bstres <- bstfit(tuning, x, ynew, nu, mselect.now, mstop2, twin, center, interaction, degree, learner, l2, nfold, n.cores, cv, tuningSwitch, k, trace, gamma, lambda=lambda, lamb, whichlambda=whichlambda, method=method, rng)
        dat1.glm <- bstres$dat1.glm
	mselect.now <- mselect[k] <- bstres$mselect  ###update mselect.now for next BJ iteration, and store the value
### Compute predicted values and convergence criteria                                    
        predres <- predval(learner, twin, dat1.glm, b, k, x, s, mselect[k])
        Fboost <- predres$Fboost
	ybst[k,] <- Fboost
        beta0bj = predres$beta0bj
        betabj = predres$betabj
        b <- predres$b
        bdiff=predres$bdiff

	if(k>1){
            ydiff <- ybstdiff[k] <- max(abs(Fboost - ybst[k-1,]))
            ybstcon[k] <- sum((ybst[k,]-ybst[k-1,])^2)/sum(ybst[k-1,]^2)
        }
        mseun[k] <- mean((Fboost-ynew)^2) #mse
        if(k >1 && trace)
            cat("    k=",k,"   ybstdiff", ybstdiff[k],"  ybstcon", ybstcon[k],"\n")
                                        #cat("    k=",k,"   ybstdiff", ybstdiff[k],"  ybstcon", ybstcon[k]," mse", mse[k],"\n")
### Check convergence status
        if(!nonconv){
            if(k > 1)
                if((learner=="linear.regression" && bdiff <= tol)
		# || (learner %in% c("enet","enet2", "mnet","snet") && ybstcon[k] <= tol)
                   || (ybstcon[k] <= tol)){
                    contype <- 0
                    break
                }
            else if(k >= iter.bj) {
                cycleydiff <- NULL
                if(learner=="linear.regression") {
                    cycle.coef.bj <- NULL
                    firstb <- betabj
                }
                nonconv <- TRUE
                firstydiff <- ydiff
                first.ybstcon <- ybstcon[k]
                first.dat1.glm <- dat1.glm
                cycleb <- NULL
                tuningSwitch <- FALSE;
                                        #if(learner=="tree")  best.iter <- best.iter                           #use the best number of trees of the previous model
                if(learner=="mars") ynew <- dat1.glm$ynew
            }
        }
        else {
            if(learner=="linear.regression"){
                if(twin) coef.bj <- coef(dat1.glm)
                else coef.bj <- coef(dat1.glm, which = 1:length(variable.names(dat1.glm)), off2int=TRUE)
                cycle.coef.bj <- rbind(cycle.coef.bj,coef.bj)
            }
            cycleydiff <- c(cycleydiff,ydiff)
            cycleperiod <- cycleperiod + 1
            if(learner=="linear.regression"){
                if(twin) 
                    tmp <- (sum((firstb - coef.bj)^2) < tol || ydiff <= tol)
                else tmp <- (sum((firstb - coef.bj[-1])^2) < tol || ydiff <= tol)
                if(tmp){
                    contype <- 1
                    break
                }
                else if(cycleperiod >= max.cycle){
                    contype <- 2
                    break
                }
            }
            else {
                if(abs(ybstcon[k]-first.ybstcon) < tol){
                    contype <- 1
                    break
                }
        else if(cycleperiod >= max.cycle){
            contype <- 2
            break
        }
            }
        }
        k <- k + 1
    }
### END of BJ iteration loop
    if(trace)
        cat("\ncycle period is",cycleperiod,"\n")
    if(contype==2)
        dat1.glm <- first.dat1.glm                                    
    if(all(!is.na(Fboost))){
        tmpy <- y[cens==1]- Fboost[cens==1]
        tmpx <- (y[cens==1] + Fboost[cens==1])/2
        mse.bj <- mean(tmpy^2)
        if(learner=="linear.regression"){
      	    if(twin){
                beta0bj <- attr(coef(dat1.glm), "offset2int")
                names(beta0bj) <- "(Intercept)"
                betabj <- coef(dat1.glm)
            }
	    else{
                beta0bj <- coef(dat1.glm)[1] + dat1.glm$offset
                betabj <- coef(dat1.glm, which = 1:length(variable.names(dat1.glm)))[-1]
            }
        }
### prediction and MSE with censoring data in real applications
        if(!is.null(valdata)){
            if(learner=="linear.regression"){
                                        #pred.bj <- predict(dat1.glm, newdata=as.matrix(valdata)[,-(1:2)]) ### for twin, or bst function. needs modification if standardized variables
                pred.bj <- as.vector(beta0bj) + as.matrix(valdata)[,-(1:2)] %*% as.vector(betabj/normx) ### needs modification if standardized variables
                mse.bj.val <- mean((valdata[,1][valdata[,2]==1] - pred.bj[valdata[,2]==1])^2)
            } else if(learner %in% c("pspline", "mars"))
                pred.bj <- predict(dat1.glm, newdata=valdata[,-(1:2)])
            else if(learner=="tree"){
                if(!twin) pred.bj <- predict(dat1.glm, newdata=as.data.frame(valdata[,-(1:2)]),n.trees=dat1.glm$n.trees)
                else pred.bj <- predict(dat1.glm, newdata=valdata[,-(1:2)])
            }
            else if(learner=="enet")
                pred.bj <- predict(dat1.glm, newx=valdata[,-(1:2)], s=s, type="fit", mode="fraction")$fit
            else if(learner %in%c("enet2", "mnet", "snet"))
		pred.bj <- predict(dat1.glm, newx=as.matrix(valdata[,-(1:2)]), type="response", which=mselect[k])
            mse.bj.val <- mean((valdata[,1][valdata[,2]==1] - pred.bj[valdata[,2]==1])^2)
        }
### prediciton and MSE with simulations
        if(learner=="enet"){
            tmp <- predict(dat1.glm, type="coef", s=s, mode="fraction")$coef
            beta0.enet <- mean(ynew) - apply(x[,dat1.glm$allset], 2, mean) %*% tmp ### intercept
            beta.enet <- rep(0, p)
            beta.enet[dat1.glm$allset] <- tmp
        }
        else if(learner %in% c("enet2", "mnet", "snet"))
	    coef.ncv <- predict(dat1.glm, newx=x, type="coefficients", which=mselect[k])
        if(trace) {
            cat("mse.bj=",mse.bj,"\n","correlation of predicted and observed times in noncensoring training data is",cor(y[cens==1],Fboost[cens==1]),"\n\n")
            cat("mse of predicted times of validate data is\n")
            cat("mse.bj.val",mse.bj.val,"\n")
        }
        coef.bj <- NA
        if(learner=="linear.regression"){
            coef.bj <- c(beta0bj, betabj)
            coef.bj <- coef.bj/c(1,normx)
            if(twin) nz.bj <- sum(abs(coef(dat1.glm))>0) 
            else nz.bj <- sum(abs(coef(dat1.glm)[-1])>0) # -1 means intercept is not counted
            if(trace) {cat("Number of Non-zero coefficients with BJ boosting excluding but listing intercept is",nz.bj,"\n")
                print(coef.bj[abs(coef.bj)>0])
            }
        }
    }
                                        #ynew is the final imputed response used for boosting
    cycle.coef.diff <- NA
    if(exists("cycle.coef.bj") && !twin) #Todo with twin=TRUE
        cycle.coef.diff <- max(abs(scale(cycle.coef.bj, coef.bj, FALSE))) #compute max of absolute difference between coef in the cycle to the one claimed as the solution
    interactions=NULL
    d <- ncol(x)
                                        #construct a one-to-one corespondence table between input column and the ensemble index
    ind <- matrix(NA,ncol=2,nrow=d+d*(d-1)/2)
    kk <- 1
    for(i in 1:d)
        for(j in i:d){
            ind[kk,1] <- i; ind[kk,2] <- j
            kk <- kk + 1
        }
    if(rel.inf && learner=="tree" && degree > 1){
        vim <- summary(dat1.glm,plotit=FALSE,order=FALSE)[,2]
        interactions <- vim.interactions(dat1.glm,pred.data=x,x.pair=subset(ind,ind[,1]!=ind[,2]),learner="tree",verbose=FALSE)
    }
### Compute which variable is selected
    if(learner=="tree"){
        if(!twin){
            xselect <- summary(dat1.glm,order=FALSE,plotit=FALSE)[,2]
            xselect <- ifelse(xselect > 0, 1, 0)   #variable selected if xselect=1, o.w. 0.
        }
	    else{
                xselect <- rep(0,dim(x)[2])
                xselect[dat1.glm$xselect] <- 1 
            }
    }
    else if(learner=="mars"){
        vim <- evimp(update(dat1.glm),sqrt.=TRUE,trim=FALSE)[,c(1,4)]
                                        #  	  vim <- evimp(update.earth(dat1.glm),sqrt=TRUE,trim=FALSE)[,c(1,4)]
        vim <- vim[order(vim[,1]),]
        vim <- vim[,2]
        xselect <- ifelse(vim > 0, 1, 0)   #variable selected if xselect=1, o.w. 0.
    }
    else if(learner=="linear.regression")
        xselect <- ifelse(abs(coef.bj[-1]) > 0, 1, 0)   #without the intercept
    else if(learner=="enet"){
        tmp <- predict(dat1.glm, type="coef", s=s, mode="fraction")$coef
        xselect <- ifelse(abs(tmp) > 0, 1, 0) 
    }
    else if(learner %in% c("enet2", "mnet", "snet"))
        xselect <- ifelse(abs(coef.ncv[-1]) > 0, 1, 0)  ### the first element is intercept, thus removed
    else if(learner=="pspline"){
        if(!twin){
	    xselect <- rep(0,dim(x)[2])
            tmp <- unique(dat1.glm$xselect())-1
            xselect[tmp] <- 1
        }
	    else{
                xselect <- rep(0,dim(x)[2])
                xselect[dat1.glm$xselect] <- 1 
            }
    }
    else if(learner=="acosso"){
        if(dat1.glm$order==1)
            xselect <- ifelse(dat1.glm$theta > 0, 1, 0)
        else{
            ind <- gen.ind(p,learner="acosso") 
            xselect <- rep(0,dim(x)[2])
            tmp <- unique(as.vector(ind[dat1.glm$theta > 0,]))
            xselect[tmp] <- 1
        }}   
### Compute variable importance and interaction measures for MARS
    if(learner=="mars" && degree > 1 && vimpint){
        ind <- gen.ind(p) 
        interactions <- vim.interactions(dat1.glm,pred.data=x,x.pair=subset(ind,ind[,1]!=ind[,2]),learner="mars",verbose=FALSE)
    }
    if(learner=="tree" && !twin)
        vim <- summary(dat1.glm,plotit=FALSE,order=FALSE)[,2]
    mse.tr <- NULL
    if(learner=="enet" && tuning){
        mse.tr <- sum((ynew - predict(dat1.glm, x, s=s, mode="frac", type="fit")$fit)^2)
        b <- predict(dat1.glm, type="coef", s=s, mode="frac")$coef
        if(any(abs(b) > 0)){
            b <- which(abs(b) > 0)
            x0 <- as.matrix(x[,b])
            if(lamb==0) q <- dim(x0)[2]
            else {q <- sum(diag(x0 %*% solve(t(x0) %*% x0 + diag(lamb, nrow=dim(x0)[2])) %*% t(x0))) ### trace
            }
        }
        else q <- 0
        mse.tr <- mse.tr/(length(ynew) - q)^2 
    }
    else if(learner=="enet") 
        coef.bj <- c(beta0.enet, beta.enet)
    else if(learner %in%c("enet2", "mnet", "snet"))
        coef.bj <- coef.ncv 
    if(!is.null(vim)) vim <- 100*vim/sum(vim)
    RET <- list(x=x,y=y,cens=cens,ynew=ynew,res.fit=dat1.glm,learner=learner,degree=degree,mse=mse,nz.bj.iter=nz.bj.iter,mse.bj=mse.bj,mse.bj.val=mse.bj.val,nz.bj=nz.bj,mse.all=mseun[1:(k-1)],yhat=Fboost,ybstdiff=c(NA,ybstdiff[1:(k-1)]),ybstcon = ybstcon,coef.bj=coef.bj,pred.bj=pred.bj,cycleperiod=cycleperiod,cycle.coef.diff = cycle.coef.diff,nonconv=nonconv,fnorm2=fnorm2,vim=vim,interactions=interactions,mselect=mselect,contype=contype,xselect=xselect,lamb=lamb, s=s, mse.tr=mse.tr,valdata=valdata, twin=twin)
    RET$call <- match.call()
    class(RET) <- "bujar"
    return(RET)
}

gen.ind <- function(d,learner="tree"){
    ind <- matrix(NA,ncol=2,nrow=d+d*(d-1)/2)
    if(learner=="mars"){
                                        #construct a one-to-one corespondence table between input column and the ensemble index
        kk <- 1
        for(i in 1:d)
            for(j in i:d){
                ind[kk,1] <- i; ind[kk,2] <- j
                kk <- kk + 1
            }
    }
    else if(learner=="tree" || learner=="acosso"){       #cf: get.gram function in acosso.R 
        ind[1:d,] <- cbind(1:d,1:d)
        next.ind <- d+1
        for(i in 1:(d-1))
            for(j in ((i+1): d)){
                ind[next.ind,] <- cbind(i,j)
                next.ind <- next.ind + 1
            }
    }
    ind
}

convbujar <- function(x){
    ybstdiff <- x$ybstdiff
    ybstcon <- x$ybstcon
    mseun <- x$mse.all
    mse <- x$mse
    fnorm2 <- x$fnorm2
    plot(ybstcon, type="b",xlab="Buckley-James estimator iteration",ylab="Convergence criterion",ylim=c(0,0.01))
}

###compute the number of covariates selected based on the position of x
nxselect <- function(obj, varpos) sum(obj$xselect[varpos] == 1)

print.bujar <- function(x, ...) {

    cat("\n")
    cat("\t Models Fitted with Buckley-James Regression\n")
    cat("\n")
    if (!is.null(x$call))
        cat("Call:\n", deparse(x$call), "\n\n", sep = "")
    cat("\n")
    if(x$learner%in%c("linear.regression","mars","pspline","tree"))
        cat("Base learner: ", x$learner, "\n")
    else cat("Regression methods: ", x$learner, "\n") 
    cat("\n")
    if(x$learner=="linear.regression"){
        cat("Coefficients: \n")
        cf <- x$coef.bj
        print(cf)
        cat("\n")
    }
    invisible(x)
}

### methods: coefficients
coef.bujar <- function(object, ...) {
    if(!object$learner %in% c("linear.regression","pspline","enet", "enet2", "mnet", "snet"))
        stop("Coefficients Not implemented for learner ",object$learner,"\n")
    object$coef.bj
### check if coef.bj is computed for learner="pspline"
}

plot.bujar <- function(x, ...){
    if(!x$learner %in% c("mars", "pspline", "acosso"))
        plot(x$res.fit)
    else stop("Not implemented for learner ",x$learner,"\n")
}

predict.bujar <- function(object, newx=NULL, ...){
    if(is.null(newx)) return(object$yhat)
    if(dim(newx)[2]!=dim(object$x)[2]) stop("newx should have the same number of predictors as x\n")
    learner <- object$learner
    dat1.glm <- object$res.fit
    if(learner=="linear.regression")
        object$coef.bj[1] + as.matrix(newx) %*% as.vector(object$coef.bj[-1])
    else if(learner %in% c("pspline", "mars"))
        pred.bj <- predict(dat1.glm, newdata=newx)
    else if(learner=="tree"){
        twin <- object$twin
    	if(!twin) pred.bj <- predict(dat1.glm, newdata=as.data.frame(newx),n.trees=dat1.glm$n.trees)
        else pred.bj <- predict(dat1.glm, newdata=newx)
    }
    else if(learner=="enet")
        pred.bj <- predict(dat1.glm, newx=newx, s=object$s, type="fit", mode="fraction")$fit
    else if(learner %in%c("enet2", "mnet", "snet")){
	mselect <- object$mselect
        k <- length(mselect)
	pred.bj <- predict(dat1.glm, newx=as.matrix(newx), type="response", which=mselect[k])
    }
}

summary.bujar <- function(object, ...)
    summary(object$res.fit, ...)
