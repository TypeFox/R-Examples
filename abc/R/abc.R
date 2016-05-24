######################################################################
#
# abc.R
#
# copyright (c) 2011-05-30, Katalin Csillery, Olivier Francois and
# Michael GB Blum with some initial code from Mark Beaumont
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
# Part of the R/abc package
# Contains: abc, return.abc, print.abc, is.abc, summary.abc, hist.abc, plot.abc
#
######################################################################

abc <- function(target, param, sumstat, tol, method, hcorr = TRUE,
                transf = "none", logit.bounds = c(0,0), subset = NULL,
                kernel = "epanechnikov", numnet = 10, sizenet = 5,
                lambda = c(0.0001,0.001,0.01), trace = FALSE, maxit =
                500, ...){
    
    call <- match.call()
    
    ## general checks that the function is used correctly
    ## ###################################################
    
    if(missing(target)) stop("'target' is missing")
    if(missing(param)) stop("'param' is missing")
    if(missing(sumstat)) stop("'sumstat' is missing")
    if(!is.matrix(param) && !is.data.frame(param) && !is.vector(param)) stop("'param' has to be a matrix, data.frame or vector.")
    if(!is.matrix(sumstat) && !is.data.frame(sumstat) && !is.vector(sumstat)) stop("'sumstat' has to be a matrix, data.frame or vector.")
    if(missing(tol)) stop("'tol' is missing")
    if(missing(method)) stop("'method' is missing with no default")  
    if(!any(method == c("rejection", "loclinear", "neuralnet","ridge"))){
        stop("Method must be rejection, loclinear, or neuralnet or ridge")
    }
    if(method == "rejection") rejmethod <- TRUE
    else rejmethod <- FALSE
    
    if(!any(kernel == c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine"))){
        kernel <- "epanechnikov"
        warning("Kernel is incorrectly defined. Setting to default (Epanechnikov)")
    }

    if(is.data.frame(param)) param <- as.matrix(param)
    if(is.data.frame(sumstat)) sumstat <- as.matrix(sumstat)
    if(is.list(target)) target <- unlist(target)
    if(is.vector(sumstat)) sumstat <- matrix(sumstat, ncol=1)
    if(length(target)!=dim(sumstat)[2]) stop("Number of summary statistics in 'target' has to be the same as in 'sumstat'.")
    
    ## stop if zero var in sumstat
    ## #########################
    nss <- length(sumstat[1,])
    cond1 <- !any(as.logical(apply(sumstat, 2, function(x) length(unique(x))-1)))
    if(cond1) stop("Zero variance in the summary statistics.")
    
    
    ## transformations
    ## ################
    ltransf <- length(transf)
    if(is.vector(param)){
        numparam <- 1
        param <- matrix(param, ncol=1)
    }
    else numparam <- dim(param)[2]
    for (i in 1:ltransf){
        if(sum(transf[i] == c("none","log","logit")) == 0){
            stop("Transformations must be none, log, or logit.")
        }
        if(transf[i]=="logit"){
            if(logit.bounds[i,1] >= logit.bounds[i,2]){
                stop("Logit bounds are incorrect.")       
            }
        }
    }
    ## if no logit change logit.bounds to NULL
    ## if(any(transf) == "logit") logit.bounds <- NULL
    
    ## no transformation should be applied when rejmethod is true
    if(rejmethod){
        if(!all(transf == "none")){
            warning("No transformation is applied when the simple rejection is used.", call.=F)
        }
        transf[1:numparam] <- "none"
    }
    else{
        if(numparam != ltransf){
            if(length(transf) == 1){
                transf <- rep(transf[1], numparam)
                warning("All parameters are \"", transf[1], "\" transformed.", sep="", call.=F)
            }
            else stop("Number of parameters is not the same as number of transformations.", sep="", call.=F)
        }
    }
    
    ## parameter and/or sumstat values that are to be excluded
    ## #######################################################
    gwt <- rep(TRUE,length(sumstat[,1]))
    gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
    if(missing(subset)) subset <- rep(TRUE,length(sumstat[,1]))
    gwt <- as.logical(gwt*subset)
    
    ## extract names of parameters and statistics if given
    ## ###################################################
    if(!length(colnames(param))){
        warning("No parameter names are given, using P1, P2, ...")
        paramnames <- paste("P", 1:numparam, sep="")
    }
    else paramnames <- colnames(param)
    
    if(!length(colnames(sumstat))){
        warning("No summary statistics names are given, using S1, S2, ...")
        statnames <- paste("S", 1:nss, sep="")
    }
    else statnames <- colnames(sumstat)
    
    ## scale everything
    ## #################
    
    scaled.sumstat <- sumstat
    for(j in 1:nss){
        scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
    }
    
    for(j in 1:nss){
        target[j] <- normalise(target[j],sumstat[,j][gwt])
    }
    
    ## calculate euclidean distance
    ## ############################
    sum1 <- 0
    for(j in 1:nss){
        sum1 <- sum1 + (scaled.sumstat[,j]-target[j])^2
    }
    dist <- sqrt(sum1)
    
    ## includes the effect of gwt in the tolerance
    dist[!gwt] <- floor(max(dist[gwt])+10)
    
    ## wt1 defines the region we're interested in

    ceiling(length(dist)*tol)->nacc
    sort(dist)[nacc]->ds
    wt1 <- (dist <= ds)
    aux<-cumsum(wt1)
    wt1 <- wt1 & (aux<=nacc)
    if(kernel == "gaussian") 
    {
    	wt1 <- rep(TRUE, length(dist))
    }
    ## transform parameters
    ## ######################
    for (i in 1:numparam){
        if(transf[i] == "log"){
            if(min(param[,i]) <= 0){
                cat("log transform: values out of bounds - correcting...")
                x.tmp <- ifelse(param[,i] <= 0,max(param[,i]),param[,i])
                x.tmp.min <- min(x.tmp)
                param[,i] <- ifelse(param[,i] <= 0, x.tmp.min,param[,i])
            }
            param[,i] <- log(param[,i])
        }
        else if(transf[i] == "logit"){
            if(min(param[,i]) <= logit.bounds[i,1]){
                x.tmp <- ifelse(param[,i] <= logit.bounds[i,1],max(param[,i]),param[,i])
                x.tmp.min <- min(x.tmp)
                param[,i] <- ifelse(param[,i] <= logit.bounds[i,1], x.tmp.min,param[,i])
            }
            if(max(param[,i]) >= logit.bounds[i,2]){
                x.tmp <- ifelse(param[,i] >= logit.bounds[i,2],min(param[,i]),param[,i])
                x.tmp.max <- max(x.tmp)
                param[,i] <- ifelse(param[,i] >= logit.bounds[i,2], x.tmp.max,param[,i])
            }
            param[,i] <- (param[,i]-logit.bounds[i,1])/(logit.bounds[i,2]-logit.bounds[i,1])
            param[,i] <- log(param[,i]/(1-param[,i]))
        }
    } # end of parameter transformations
    
    ## select summary statistics in region
    ## ###################################
    ss <- sumstat[wt1,]
    unadj.values <- param[wt1,]
    
    ## if simple rejection or in the selected region there is no var in sumstat
    ## ########################################################################
    statvar <- as.logical(apply(cbind(sumstat[wt1,]), 2, function(x) length(unique(x))-1))
    cond2 <- !any(statvar)
    
    if(cond2 && !rejmethod)
        stop("Zero variance in the summary statistics in the selected region. Try: checking summary statistics, choosing larger tolerance, or rejection method.")
    
    if(rejmethod){
        if(cond2) warning("Zero variance in the summary statistics in the selected region. Check summary statistics, consider larger tolerance.")
        weights <- rep(1,length=sum(wt1))
        adj.values <- NULL
        residuals <- NULL
        lambda <- NULL
    }
    
    else{
        if(cond2) cat("Warning messages:\nStatistic(s)", statnames[!statvar], "has/have zero variance in the selected region.\nConsider using larger tolerance or the rejection method or discard this/these statistics, which might solve the collinearity problem in 'lsfit'.\n", sep=", ")
        
        ## weights
        if(kernel == "epanechnikov") weights <- 1 - (dist[wt1]/ds)^2
        if(kernel == "rectangular") weights <- dist[wt1]/ds
        if(kernel == "gaussian") weights <- 1/sqrt(2*pi)*exp(-0.5*(dist/(ds/2))^2)
        if(kernel == "triangular") weights <- 1 - abs(dist[wt1]/ds)
        if(kernel == "biweight") weights <- (1 - (dist[wt1]/ds)^2)^2
        if(kernel == "cosine") weights <- cos(pi/2*dist[wt1]/ds)
        
        ## regression correction
        ## ######################
        if(method == "loclinear"){
            
           fit1 <- lsfit(scaled.sumstat[wt1,],param[wt1,],wt=weights)
            #pred <- t(fit1$coeff) %*% c(1,target)
            pred <- t(structure(cbind(fit1$coefficients)[fit1$qr$pivot,], names=names(fit1$coefficients))) %*% c(1,target)
            
            pred <- matrix(pred, ncol=numparam, nrow=sum(wt1), byrow=TRUE)
            #residuals <- fit1$residuals
            residuals<-param[wt1,]-t(t(structure(cbind(fit1$coefficients)[fit1$qr$pivot,], names=names(fit1$coefficients))) %*%t(cbind(1,scaled.sumstat[wt1,])))
	    	residuals<-cbind(residuals)
	    	#rrr<<-residuals
	    	the_m<-apply(residuals,FUN=mean,2)
	    	residuals<-sapply(1:numparam,FUN=function(x){residuals[,x]-the_m[x]})
	    	pred <- sapply(1:numparam,FUN=function(x){pred[,x]+the_m[x]})

            
            sigma2<-apply(as.matrix(residuals),FUN=function(x){sum((x)^2 *weights)/sum(weights)},MARGIN=2)	
            aic<-sum(wt1)*sum(log(sigma2))+2*(nss+1)*numparam
            bic<-sum(wt1)*sum(log(sigma2))+log(sum(wt1))*(nss+1)*numparam
            
            if(hcorr == TRUE){
                 fit2 <- lsfit(scaled.sumstat[wt1,],log(residuals^2),wt=weights)
                auxaux<-t(structure(cbind(fit2$coefficients)[fit2$qr$pivot,], names=names(fit2$coefficients))) %*% c(1,target)
                pred.sd <- sqrt(exp(auxaux))
                pred.sd <- matrix(pred.sd, nrow=sum(wt1), ncol=numparam, byrow=T)
                pred.si<-t(t(structure(cbind(fit2$coefficients)[fit2$qr$pivot,], names=names(fit2$coefficients))) %*%t(cbind(1,scaled.sumstat[wt1,])))
                pred.si<-sqrt(exp(pred.si))
               	#residuals2<-log(residuals^2)-t(t(structure(cbind(fit2$coefficients)[fit2$qr$pivot,], names=names(fit2$coefficients))) %*%t(cbind(1,scaled.sumstat[wt1,])))
                #fv <- sqrt(exp(log(fit1$residuals^2) - fit2$residuals))
                #fv <- sqrt(exp(residuals2))
               # adj.values <- pred + (pred.sd*residuals) / fv # correction heteroscedasticity
               adj.values <- pred + (pred.sd*residuals) / pred.si
               residuals<-(pred.sd*residuals) / pred.si
            }
            else{
                adj.values <- pred + residuals
            }
            colnames(adj.values) <- colnames(unadj.values)
            lambda <- NULL
        }
        
        ## neural network regression
        ## ##########################
        if(method == "neuralnet"){
            linout <- TRUE
            
            ## normalise parameters
            param.mad <- c()
            for(i in 1:numparam){
                param.mad[i] <- mad(param[,i][gwt]) # save for renormalisation
                param[,i] <- normalise(param[,i],param[,i][gwt])
            }
            
            lambda <- sample(lambda, numnet, replace=T)
            fv <- array(dim=c(sum(wt1), numparam, numnet))
            pred <- matrix(nrow=numparam, ncol=numnet)
            for (i in 1:numnet){
                fit1 <- nnet(scaled.sumstat[wt1,], param[wt1,], weights = weights, decay = lambda[i],
                             size = sizenet, trace = trace, linout = linout, maxit = maxit, ...)
                cat(i)
                fv[,,i] <- fit1$fitted.values
                pred[,i] <- predict(fit1, data.frame(rbind(target)))
            }
            cat("\n")
            pred.med <- apply(pred, 1, median)
            pred.med <- matrix(pred.med, nrow=sum(wt1), ncol=numparam, byrow=T)
            fitted.values <- apply(fv, c(1,2), median)
            residuals <- param[wt1,] - fitted.values # median of fitted values par nnets for each accepted point and parameter
            
            if(hcorr == TRUE){
                pred2 <- matrix(nrow=numparam, ncol=numnet)
                fv2 <- array(dim=c(sum(wt1), numparam, numnet))
                for (i in 1:numnet){
                    fit2 <- nnet(scaled.sumstat[wt1,], log(residuals^2), weights = weights, decay = lambda[i], size = sizenet, trace = trace, linout = linout, ...)
                    cat(i)
                    fv2[,,i] <- fit2$fitted.values
                    pred2[,i] <- predict(fit2, data.frame(rbind(target)))
                }
                cat("\n")
                pred.sd <- sqrt(exp(apply(pred2, 1, median)))
                pred.sd <- matrix(pred.sd, nrow=sum(wt1), ncol=numparam, byrow=T)
                fv.sd <- sqrt(exp(apply(fv2, c(1,2), median)))
                adj.values <- pred.med + (pred.sd*residuals)/fv.sd # correction heteroscedasticity
           		residuals<-(pred.sd*residuals)/fv.sd
           }
            else{
                adj.values <- pred.med + residuals
            }
            colnames(adj.values) <- colnames(unadj.values)
            
            ## renormalise
            for(i in 1:numparam){
                ##       residuals[,i] <- residuals[,i]*param.mad[i] not much sense...
                adj.values[,i] <- adj.values[,i]*param.mad[i]
            }
            
        } # end of neuralnet
             if(method == "ridge"){
            #print("la")
            ## normalise parameters
            param.mad <- c()
            for(i in 1:numparam){
                param.mad[i] <- mad(param[,i][gwt]) # save for renormalisation
                param[,i] <- normalise(param[,i],param[,i][gwt])
            }
           # print("ici")
            numnet<-length(lambda)
           #lambda <- sample(lambda, numnet, replace=T)
            fv <- array(dim=c(sum(wt1), numparam, numnet))
            pred <- matrix(nrow=numparam, ncol=numnet)
            mataux<-sqrt(diag(weights))
            paramaux<-as.matrix(mataux%*%param[wt1,])
            
            #print(dim(paramaux))
            scaledaux<-mataux%*%scaled.sumstat[wt1,]
            #targetaux<-drop(mataux%*%target
            for(parcur in (1:numparam))
            {
            	#cat(parcur)
            	fit1<-lm.ridge(paramaux[,parcur]~scaledaux,lambda=lambda)
            	for (i in 1:numnet){
            #   fit1 <- nnet(scaled.sumstat[wt1,], param[wt1,], weights = weights, decay = lambda[i],
            #                 size = sizenet, trace = trace, linout = linout, maxit = maxit, ...)
            #    cat(i)
            	    fv[,parcur,i] <- drop(cbind(1,scaled.sumstat[wt1,])%*%(rbind(coef(fit1))[i,]))
                	pred[parcur,i] <- drop(c(1, target) %*%(rbind(coef(fit1))[i,]))
            					}
            }
           # cat("\n")
            pred.med <- apply(pred, 1, median)
            pred.med <- matrix(pred.med, nrow=sum(wt1), ncol=numparam, byrow=T)
            fitted.values <- apply(fv, c(1,2), median)
            residuals <- param[wt1,] - fitted.values # median of fitted values par nnets for each accepted point and parameter
            
            if(hcorr == TRUE){
                pred2 <- matrix(nrow=numparam, ncol=numnet)
                fv2 <- array(dim=c(sum(wt1), numparam, numnet))
                 for(parcur in (1:numparam))
            	{
            		lresidaux<-(mataux%*%(log(residuals[,parcur]^2)))
              		fit2<-lm.ridge(lresidaux~scaledaux,lambda=lambda)
                	for (i in 1:numnet){
                    	#cat(i)
                    	fv2[,parcur,i] <- drop(cbind(1,scaled.sumstat[wt1,])%*%(rbind(coef(fit2))[i,]))
						pred2[parcur,i] <- drop(c(1, target) %*%(rbind(coef(fit2))[i,]))         
									}
                }
                cat("\n")
                pred.sd <- sqrt(exp(apply(pred2, 1, median)))
                pred.sd <- matrix(pred.sd, nrow=sum(wt1), ncol=numparam, byrow=T)
                fv.sd <- sqrt(exp(apply(fv2, c(1,2), median)))
                adj.values <- pred.med + (pred.sd*residuals)/fv.sd # correction heteroscedasticity
            	residuals<- (pred.sd*residuals)/fv.sd
            }
            else{
                adj.values <- pred.med + residuals
            }
            colnames(adj.values) <- colnames(unadj.values)
            
            ## renormalise
            for(i in 1:numparam){
                ##       residuals[,i] <- residuals[,i]*param.mad[i] not much sense...
                adj.values[,i] <- adj.values[,i]*param.mad[i]
            }
            
        } # end of ridge
        
    } # end of else than rejmethod
    
    ## back transform parameter values
    ## ################################
    if(numparam == 1){
        unadj.values <- matrix(unadj.values, ncol=1)
        if(method!="rejection")
        {adj.values <- matrix(adj.values, ncol=1)
        residuals <- matrix(residuals, ncol=1)}
    }
    
    for (i in 1:numparam){
        if(transf[i] == "log"){
            unadj.values[,i] <- exp(unadj.values[,i])
            adj.values[,i] <- exp(adj.values[,i])
        }
        else if(transf[i] == "logit"){
            unadj.values[,i] <- exp(unadj.values[,i])/(1+exp(unadj.values[,i]))
            unadj.values[,i] <- unadj.values[,i]*(logit.bounds[i,2]-logit.bounds[i,1])+logit.bounds[i,1]
            adj.values[,i] <- exp(adj.values[,i])/(1+exp(adj.values[,i]))
            adj.values[,i] <- adj.values[,i]*(logit.bounds[i,2]-logit.bounds[i,1])+logit.bounds[i,1]
        }
    }
    
    abc.return(transf, logit.bounds, method, call, numparam, nss, paramnames, statnames,
               unadj.values, adj.values, ss, weights, residuals, dist, wt1, gwt, lambda, hcorr,aic,bic)
    
}


abc.return <- function(transf, logit.bounds, method, call, numparam, nss, paramnames, statnames,
                       unadj.values, adj.values, ss, weights, residuals, dist, wt1, gwt, lambda, hcorr,aic,bic){
    
    if(method == "rejection"){
        out <- list(unadj.values=unadj.values, ss=ss, dist=dist,
                    call=call, na.action=gwt, region=wt1, transf=transf, logit.bounds = logit.bounds,
                    method="rejection", numparam=numparam, numstat=nss, names=list(parameter.names=paramnames, statistics.names=statnames))
    }
    else if(method == "loclinear"){
        out <- list(adj.values=adj.values, unadj.values=unadj.values,
                    ss=ss, weights=weights, residuals=residuals, dist=dist,
                    call=call, na.action=gwt, region=wt1, transf=transf, logit.bounds = logit.bounds,
                    method="loclinear", hcorr = hcorr, numparam=numparam, numstat=nss,aic=aic,bic=bic,
                    names=list(parameter.names=paramnames, statistics.names=statnames))
    }
    else if(method == "ridge"){
        out <- list(adj.values=adj.values, unadj.values=unadj.values,
                    ss=ss, weights=weights, residuals=residuals, dist=dist,
                    call=call, na.action=gwt, region=wt1, transf=transf, logit.bounds = logit.bounds,
                    method="ridge", hcorr = hcorr, numparam=numparam, numstat=nss,
                    names=list(parameter.names=paramnames, statistics.names=statnames))
    }
    else if(method == "neuralnet"){
        out <- list(adj.values=adj.values, unadj.values=unadj.values,
                    ss=ss, weights=weights, residuals=residuals, dist=dist,
                    call=call, na.action=gwt, region=wt1, transf=transf, logit.bounds = logit.bounds,
                    method="neuralnet", hcorr = hcorr, lambda=lambda, numparam=numparam, numstat=nss,
                    names=list(parameter.names=paramnames, statistics.names=statnames))
    }
    
    class(out) <- "abc"
    out
    
}

is.abc <- function(x){
    if (inherits(x, "abc")) TRUE
    else FALSE    
}

print.abc <- function(x, ...){
    if (!inherits(x, "abc")) 
        stop("Use only with objects of class \"abc\".", call.=F)
    abc.out <- x
    cl <- abc.out$call
    cat("Call:\n")
    dput(cl, control=NULL)
    
    cat("Method:\n")
    ## rejection
    if (abc.out$method == "rejection")
        cat("Rejection\n\n")
    ## loclinear
    else if (abc.out$method == "loclinear"){
        if(abc.out$hcorr == T){
            cat("Local linear regression\n")
            cat("with correction for heteroscedasticity\n\n")
        }
        else cat("Local linear regression\n\n")
    }
    ## ridge
    else if (abc.out$method == "ridge"){
        if(abc.out$hcorr == T){
            cat("Ridge regression\n")
            cat("with correction for heteroscedasticity\n\n")
        }
        else cat("Ridge regression\n\n")
    }
    ## nnet
    else if (abc.out$method == "neuralnet"){
        if(abc.out$hcorr == T){
            cat("Non-linear regression via neural networks\n")
            cat("with correction for heteroscedasticity\n\n")
        }
        else cat("Non-linear regression via neural networks\n\n")
    }
    cat("Parameters:\n")
    cat(abc.out$names$parameter.names, sep=", ")
    cat("\n\n")
    cat("Statistics:\n")
    cat(abc.out$names$statistics.names, sep=", ")
    cat("\n\n")        
    cat("Total number of simulations", length(abc.out$na.action), "\n\n")
    cat("Number of accepted simulations: ", dim(abc.out$unadj.values)[1], "\n\n")
    invisible(abc.out)
}

summary.abc <- function(object, unadj = FALSE, intvl = .95, print = TRUE, digits = max(3, getOption("digits")-3), ...){
  
  if (!inherits(object, "abc")) 
    stop("Use only with objects of class \"abc\".", call.=F)
  abc.out <- object
  np <- abc.out$numparam
  cl <- abc.out$call
  parnames <- abc.out$names[[1]]
  npri <- dim(abc.out$unadj.values)[1]
  npost <- dim(abc.out$unadj.values)[1]
  myprobs <- c((1-intvl)/2, 0.5, 1-(1-intvl)/2)
  
  if(print){
    cat("Call: \n")
    dput(cl, control=NULL)
  }
  
  if(abc.out$method == "rejection" || unadj){
    if(print) cat(paste("Data:\n abc.out$unadj.values (",npost," posterior samples)\n\n", sep=""))
    res <- matrix(abc.out$unadj.values, ncol=np)
    mymin <- apply(res, 2, min)
    mymax <- apply(res, 2, max)
    mymean <- apply(res, 2, mean)
    mymode <- apply(res, 2, getmode, ...)
    quants <- apply(res, 2, quantile, probs = myprobs)
    sums <- rbind(mymin, quants[1,], quants[2,], mymean, mymode, quants[3,], mymax)
    dimnames(sums) <- list(c("Min.:", paste(c(myprobs[1])*100, "% Perc.:", sep=""),
                             "Median:", "Mean:", "Mode:",
                             paste(c(myprobs[3])*100, "% Perc.:", sep=""), "Max.:"),
                           parnames)
  }
  else{
    if(print){
      cat(paste("Data:\n abc.out$adj.values (",npost," posterior samples)\n", sep=""))
      cat(paste("Weights:\n abc.out$weights\n\n", sep=""))
    }
    res <- matrix(abc.out$adj.values, ncol=np)
    wt <- abc.out$weights
    mymin <- apply(res, 2, min)
    mymax <- apply(res, 2, max)
    mymean <- apply(res, 2, weighted.mean, w = wt)
    mymode <- apply(res, 2, getmode, wt/sum(wt), ...)
    quants <- apply(res, 2, function(x) rq(x~1, tau = myprobs, weights = wt)$coef)
    sums <- rbind(mymin, quants[1,], quants[2,], mymean, mymode, quants[3,], mymax)
    dimnames(sums) <- list(c("Min.:", paste("Weighted ", c(myprobs[1])*100, " % Perc.:", sep=""),
                             "Weighted Median:", "Weighted Mean:", "Weighted Mode:",
                             paste("Weighted ", c(myprobs[3])*100, " % Perc.:", sep=""), "Max.:"),
                           parnames)
  }
  class(sums) <- "table"
  if(print){
    if(length(digits) == 1) 
      print(round(sums, digits = digits), quote=FALSE)
    else
      print(apply(rbind(digits, sums), 2, function(a) round(a[-1], digits = a[1])), quote=FALSE)
    invisible(sums)
  }
  else sums
}

getmode <- function(x, weights = NULL, ...){

  ##  if(missing(bw) | missing(kernel) | missing(window) | missing(n)) warning("density.default() was used to calculate the posterior mode")
  d <- density(x, weights = weights)
  d$x[which(d$y == max(d$y))][1] 
}


hist.abc <- function(x, unadj = FALSE, true = NULL,
                     file = NULL, postscript = FALSE, onefile = TRUE, ask = !is.null(deviceIsInteractive()),
                     col.hist = "grey", col.true = "red", caption = NULL, ...){

  
  if (!inherits(x, "abc")) 
      stop("Use only with objects of class \"abc\".", call.=F)
  abc.out <- x
  np <- abc.out$numparam
  parnames <- abc.out$names[[1]]

  if(is.null(caption)) mycaption <- as.graphicsAnnot(parnames)
  else mycaption <- caption
  
  ## checks if true is given
  if(!is.null(true)){
      if(!is.vector(true)){
          stop("Supply true parameter value(s) as a vector.", call.=F)
      }
      if (length(true) != np){
          stop("Number of true values has to be the same as number of parameters in 'x'.", call.=F)
      }
      cond <- isTRUE(c(match(names(true), parnames), match(names(true), parnames)))
      if(cond) stop("Names do not match in 'x' and 'true'.", call.=F)
  }

  if(abc.out$method == "rejection") res <- abc.out$unadj.values
  else if (unadj){
    rej <- abc.out$unadj.values
    res <- abc.out$adj.values
  }
  else res <- abc.out$adj.values

  ## Devices
  ## ##########
  save.devAskNewPage <- devAskNewPage()
  if(!is.null(file)){
    file <- substitute(file)
    if(!postscript) pdf(file = paste(file, "pdf", sep="."), onefile=onefile)
    if(postscript) postscript(file = paste(file, "ps", sep="."), onefile=onefile)
  }
  else{
    if (ask && prod(par("mfcol")) < np) {
      devAskNewPage(TRUE)
    }
  }

  myhist <- list()
  for(i in 1:np){
    myxlim <- range(res[,i], true[i])
    myhist[[i]] <- hist(res[,i], xlab = parnames[i], main = paste("Posterior histogram of ", mycaption[i], sep=""),
                        col = col.hist, xlim = myxlim, ...)
    myhist[[i]]$xname <- parnames[i]
    if(!is.null(true)){
      abline(v = true[i], lwd = 2, col = col.true)
    }
  }

  if(!is.null(file)){
    dev.off()
  }
  else devAskNewPage(save.devAskNewPage)
  
  names(myhist) <- parnames
  invisible(myhist)

} # end of hist.abc


##########################################################################################

## plots: for each parameter
## 1. prior density
## 2. posterior + prior density
## 3. distances vs parameters
## 4. histogram of the residuals

##########################################################################################
plot.abc <- function(x, param, subsample = 1000, true = NULL,
                     file = NULL, postscript = FALSE, onefile = TRUE, ask = !is.null(deviceIsInteractive()), ...){

  if (!inherits(x, "abc")) 
    stop("Use only with objects of class \"abc\".", call.=F)
  
  abc.out <- x
  mymethod <- abc.out$method

  if(mymethod == "rejection")
    stop("Diagnostic plots can be displayed only when method is \"loclinear\", \"neuralnet\" or \"ridge\".", cal.=F)

  if(!is.matrix(param) && !is.data.frame(param) && !is.vector(param)) stop("'param' has to be a matrix, data.frame or vector.", call.=F)
  if(is.null(dim(param))) param <- matrix(param, ncol=1)
  if(is.data.frame(param)) param <- as.matrix(param)
  
  np <- abc.out$numparam
  numsim <- length(param)/np
  alldist <- log(abc.out$dist)
  myregion <- abc.out$region
  residuals <- abc.out$residuals
  parnames <- abc.out$names$parameter.names
  transf <- abc.out$transf
  logit.bounds <- abc.out$logit.bounds
  
  ##check if param is compatible with x
  cond <- isTRUE(c(match(colnames(param), parnames), match(parnames, colnames(param))))
  if(cond) stop("'abc.out' and 'param' are not compatible; paramater names are different.", call.=F)
  
  ## checks if true is given
  if(!is.null(true)){
      if(!is.vector(true)){
          stop("Supply true parameter value(s) as a vector.", call.=F)
      }
      if (length(true) != np){
          stop("Number of true values has to be the same as number of parameters in 'x'.", call.=F)
      }
      cond <- isTRUE(c(match(names(true), parnames), match(names(true), parnames)))
      if(cond) stop("Names do not match in 'x' and 'true'.", call.=F)
  }
  
  rej <- abc.out$unadj.values
  res <- abc.out$adj.values
  if(np == 1) rej <- matrix(rej, ncol=1)
  
  if(is.vector(param)){
    np.orig <- 1
    nsim <- length(param)
  }
  else if(is.matrix(param)){
    np.orig <- dim(param)[2]
    nsim <- dim(param)[1]

    myorder <- match(parnames, colnames(param))
    if(isTRUE(myorder-1:np)){
      param <- param[, myorder]
      warning("'param' is being re-ordered according to 'abc.out'...", call.=F, immediate.=T)
    }
  }
  
  ## check if param has the right dimensions
  if(np.orig != np){
    stop("The number parameters supplied in \"param\" is different from that in \"x\".", call.=F)
  }
  
  ## for (i in 1:np){
  ##   if(transf[i] == "log"){
  ##     if(min(param[,i]) <= 0){
  ##       warning("Correcting out of bounds values for \"log\" transformed parameters.", call.=F, immediate.=T)
  ##       x.tmp <- ifelse(param[,i] <= 0,max(param[,i]),param[,i])
  ##       x.tmp.min <- min(x.tmp)
  ##       param[,i] <- ifelse(param[,i] <= 0, x.tmp.min,param[,i])
  ##     }
  ##     param[,i] <- log(param[,i])
  ##     if(min(rej[,i]) <= 0){
  ##       warning("Correcting out of bounds values for \"log\" transformed parameters.", call.=F, immediate.=T)
  ##       x.tmp <- ifelse(rej[,i] <= 0,max(rej[,i]),rej[,i])
  ##       x.tmp.min <- min(x.tmp)
  ##       rej[,i] <- ifelse(rej[,i] <= 0, x.tmp.min,rej[,i])
  ##     }
  ##     rej[,i] <- log(rej[,i])
  ##     if(min(res[,i]) <= 0){
  ##       warning("Correcting out of bounds values for \"log\" transformed parameters.", call.=F, immediate.=T)
  ##       x.tmp <- ifelse(res[,i] <= 0,max(res[,i]),res[,i])
  ##       x.tmp.min <- min(x.tmp)
  ##       res[,i] <- ifelse(res[,i] <= 0, x.tmp.min,res[,i])
  ##     }
  ##     res[,i] <- log(res[,i])
  ##     if(!is.null(true)){
  ##       if(min(true[i]) <= 0){
  ##         x.tmp <- ifelse(true[i] <= 0,max(true[i]),true[i])
  ##         x.tmp.min <- min(x.tmp)
  ##         true[,i] <- ifelse(true[i] <= 0, x.tmp.min,true[i])
  ##       }
  ##       true[i] <- log(true[i])
  ##     }
  ##   }
  ##   else if(transf[i] == "logit"){
  ##     if(min(param[,i]) <= logit.bounds[i,1]){
  ##       x.tmp <- ifelse(param[,i] <= logit.bounds[i,1],max(param[,i]),param[,i])
  ##       x.tmp.min <- min(x.tmp)
  ##       param[,i] <- ifelse(param[,i] <= logit.bounds[i,1], x.tmp.min,param[,i])
  ##     }
  ##     if(max(param[,i]) >= logit.bounds[i,2]){
  ##       x.tmp <- ifelse(param[,i] >= logit.bounds[i,2],min(param[,i]),param[,i])
  ##       x.tmp.max <- max(x.tmp)
  ##       param[,i] <- ifelse(param[,i] >= logit.bounds[i,2], x.tmp.max,param[,i])
  ##     }
  ##     param[,i] <- (param[,i]-logit.bounds[i,1])/(logit.bounds[i,2]-logit.bounds[i,1])
  ##     param[,i] <- log(param[,i]/(1-param[,i]))

  ##     if(min(rej[,i]) <= logit.bounds[i,1]){
  ##       x.tmp <- ifelse(rej[,i] <= logit.bounds[i,1],max(rej[,i]),rej[,i])
  ##       x.tmp.min <- min(x.tmp)
  ##       rej[,i] <- ifelse(rej[,i] <= logit.bounds[i,1], x.tmp.min,rej[,i])
  ##     }
  ##     if(max(rej[,i]) >= logit.bounds[i,2]){
  ##       x.tmp <- ifelse(rej[,i] >= logit.bounds[i,2],min(rej[,i]),rej[,i])
  ##       x.tmp.max <- max(x.tmp)
  ##       rej[,i] <- ifelse(rej[,i] >= logit.bounds[i,2], x.tmp.max,rej[,i])
  ##     }
  ##     rej[,i] <- (rej[,i]-logit.bounds[i,1])/(logit.bounds[i,2]-logit.bounds[i,1])
  ##     rej[,i] <- log(rej[,i]/(1-rej[,i]))

  ##     if(min(res[,i]) <= logit.bounds[i,1]){
  ##       x.tmp <- ifelse(res[,i] <= logit.bounds[i,1],max(res[,i]),res[,i])
  ##       x.tmp.min <- min(x.tmp)
  ##       res[,i] <- ifelse(res[,i] <= logit.bounds[i,1], x.tmp.min,res[,i])
  ##     }
  ##     if(max(res[,i]) >= logit.bounds[i,2]){
  ##       x.tmp <- ifelse(res[,i] >= logit.bounds[i,2],min(res[,i]),res[,i])
  ##       x.tmp.max <- max(x.tmp)
  ##       res[,i] <- ifelse(res[,i] >= logit.bounds[i,2], x.tmp.max,res[,i])
  ##     }
  ##     res[,i] <- (res[,i]-logit.bounds[i,1])/(logit.bounds[i,2]-logit.bounds[i,1])
  ##     res[,i] <- log(res[,i]/(1-res[,i]))



  ##     if(!is.null(true)){
  ##       if(min(true[i]) <= logit.bounds[i,1]){
  ##         x.tmp <- ifelse(true[i] <= logit.bounds[i,1],max(true[i]),true[i])
  ##         x.tmp.min <- min(x.tmp)
  ##         true[i] <- ifelse(true[i] <= logit.bounds[i,1], x.tmp.min,true[i])
  ##       }
  ##       if(max(true[i]) >= logit.bounds[i,2]){
  ##         x.tmp <- ifelse(true[i] >= logit.bounds[i,2],min(true[i]),true[i])
  ##         x.tmp.max <- max(x.tmp)
  ##         true[i] <- ifelse(true[i] >= logit.bounds[i,2], x.tmp.max,true[i])
  ##       }
  ##       true[i] <- (true[i]-logit.bounds[i,1])/(logit.bounds[i,2]-logit.bounds[i,1])
  ##       true[i] <- log(true[i]/(1-true[i]))
  ##     }
  ##   }
  ## } # end of parameter transformations
  
  ## devices
  ## ########
  save.devAskNewPage <- devAskNewPage()
  if(!is.null(file)){
    file <- substitute(file)
    if(!postscript) pdf(file = paste(file, "pdf", sep="."), onefile=onefile)
    if(postscript) postscript(file = paste(file, "ps", sep="."), onefile=onefile)
  }
  else{
    if (ask && 1 < np) {
      devAskNewPage(TRUE)
    }
  }
  
  ## if param is too large, draw a random sample from it the size of
  ## which can be set by the user
  mysample <- sample(1:nsim, subsample)
  
  ## plots
  ## #####

  # prior, reg, rej, true
  ltys <- c(3, 1, 1, 3)
  lwds <- c(1, 2, 1, 2)
  cols <- c("black", "red", "black", "black")
  
  unadj <- TRUE
  for(i in 1:np){

    par(mfcol=c(2,2))    
    prior.d <- density(param[mysample,i])
    post.d <- density(res[,i])
    rej.d <- density(rej[,i])
    myxlim <- range(c(post.d$x, rej.d$x, true))
    myylim <- range(c(post.d$y, rej.d$y, prior.d$y))
    par(cex = 1, cex.main = 1.2, cex.lab = 1.1, lwd = 2)

    ##if(transf[i] == "none")
    myxlab <- parnames[i]
    ##else if(transf[i] == "log") myxlab <- paste("log(", parnames[i], ")", sep="")
    ##else if(transf[i] == "logit") myxlab <- paste("log(", parnames[i], "/(1-",parnames[i],"))", sep="")
    
    ## 1. prior
    ## ########
    plot(prior.d, main = "", col=cols[1], lty = ltys[1], lwd=lwds[1], xlab=myxlab, ...)
    title("Prior", sub = paste("N =", prior.d$n, "  Bandwidth =", formatC(prior.d$bw)))

    ## 2. posterior
    ## ############
    plot(post.d, col=cols[2], lty = ltys[2], lwd=lwds[2], xlim = myxlim, ylim = myylim, main = "", xlab=myxlab, ...)
    lines(prior.d, col=cols[1], lty = ltys[1], lwd=lwds[1], ...)
    lines(rej.d, col = cols[3], lty = ltys[3], lwd=lwds[3], ...)
    if(!is.null(true)){
      segments(x0 = true[i], x1 = true[i], y0 = myylim[1], y1 = myylim[2], col=cols[4], lty = ltys[4], lwd=lwds[4])
      mtext(text="True value", side = 1, line=2, at=true[i])
    }
      
    title(paste("Posterior with \"", mymethod, "\"", sep=""), sub=paste("N =", post.d$n, "  Bandwidth =", formatC(post.d$bw)))
    mtext("\"rejection\" and prior as reference", side=3, line = 0.5)
    
    ## 3. distances
    ## ############

    mypch <- 19
    myylim <- range(alldist[mysample])#*c(1,1.2)

    plot(param[mysample,i], alldist[mysample], col=cols[1],
         xlab=myxlab, ylab="log Euclidean distance", ylim=myylim, pch=mypch, ...)
    title("Euclidean distances")
    mtext(paste("N(All / plotted) = ", numsim, "/", subsample, sep=" "), side = 3, line = .5)
    points(param[myregion,i], alldist[myregion], col=cols[2], pch=mypch, ...)
    mypar <- par()
    if(!is.null(true)){
      segments(x0 = true[i], x1 = true[i], y0 = myylim[1], y1 = myylim[2], col=cols[4], lty = ltys[4], lwd=lwds[4])
      mtext(text="True value", side = 1, line=2, at=true[i])
    }
  
    ## 4. residuals
    ## ############

    if(mymethod=="loclinear") mymain <- "Residuals from lsfit()"
    if(mymethod=="neuralnet") mymain <- "Residuals from nnet()"
    if(mymethod=="ridge") mymain <- "Residuals from lm.ridge()"
    qqnorm(residuals, pch=mypch, main=mymain, sub="Normal Q-Q plot", xlab="Theoretical quantiles", ylab="Residuals",...)
    qqline(residuals)
    
  } # np

  par(mfcol=c(1,1))
  
  if(!is.null(file)){
    dev.off()
  }
  else devAskNewPage(save.devAskNewPage)

  invisible()
  
} # end of plot.abc
