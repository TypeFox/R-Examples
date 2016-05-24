nsccv <- function(x, y=NULL, proby=NULL, nfold = min(table(y)), folds = NULL, threshold =
        NULL, threshold.scale = NULL, survival.time=NULL, censoring.status=NULL, ngroup.survival=NULL,prior, object, ...)
{
        this.call <- match.call()

        argy <- y
        
#         if( !is.null(y) & !is.null(proby)){
#           stop("Must have at most one of y and  proby  present in the data object")
#         }

        if(is.null(y)){ y <- as.factor(apply(proby,1,which.is.max))}
        
        n <- length(y)

if(is.null(nfold) & is.null(survival.time)) {nfold <- min(table(y))}
if(is.null(nfold) & !is.null(survival.time)) {nfold <- 10}


 if(is.null(survival.time)){
        if(is.null(folds)) {
                folds <-balanced.folds(y)
        }
       }


        if(!is.null(survival.time)){
        if(is.null(folds)) {
                folds <- split(sample(1:n), rep(1:nfold, length = n))
        }
        }
         
nfold<- length(folds)

        if(missing(prior)) {
                if(missing(object))
                        prior <- table(y)/n
                else prior <- object$prior
        }
    
        if(missing(threshold)) {
                if(missing(object))
                        stop("Must either supply threshold argument, or an nsc object"
                                )
                else {
                        threshold <- object$threshold
                        threshold.scale <- object$threshold.scale
                        se.scale <- object$se.scale
                }
        }
       
        n.threshold <- length(threshold)        ### Set up the data structures
        yhat <- rep(list(y), n.threshold)
        names(yhat) <- paste(seq(n.threshold))
        yhat <- data.frame(yhat)
        n.class <- table(y)
        prob <- array(1, c(n, length(n.class), n.threshold))
        size <- double(n.threshold)
        hetero <-object$hetero
        cv.objects=vector("list",nfold)
        for(ii in 1:nfold) {
                cat("Fold", ii, ":")
                a <- nsc(x[,  - folds[[ii]]], y=argy[ - folds[[ii]]], x[, folds[[ii
                        ]], drop = FALSE], proby=proby[-folds[[ii]],],
                         threshold = threshold, threshold.scale
                         = threshold.scale, se.scale = se.scale, prior = prior,
                          hetero=hetero,
                        ..., remove.zeros = FALSE)
                size <- size + a$nonzero
                prob[folds[[ii]],  ,  ] <- a$prob
                yhat[folds[[ii]],  ] <- a$yhat
                cat("\n")
        cv.objects[[ii]]=a
        }
        if(missing(object))
                size <- round(size/nfold)
        else size <- object$nonzero
        error <- rep(NA, n.threshold)
        loglik <- error
        pvalue.survival <- error
        
        pvalue.survival.func <- function(group, survival.time, censoring.status,ngroup.survival){
            temp <- coxph(Surv(survival.time, censoring.status)~as.factor(group))
            loglik <- 2*(temp$loglik[2]-temp$loglik[1])
            return(1-pchisq(loglik, ngroup.survival-1))
          }
        
        if(!is.null(proby)){proby.temp <-proby}
        else if(!is.null(survival.time)){proby.temp <- pamr.surv.to.class2(survival.time,
                                       censoring.status, n.class=ngroup.survival)$prob
                                       }
        
        for(i in 1:n.threshold) {
      
                if(is.null(survival.time) & is.null(proby)){error[i] <- sum(yhat[, i] != y)/n}
                if(!is.null(survival.time)){
                    
                    temp <- c(yhat[,i],names(table(y)))
                    Yhat <- model.matrix( ~ factor(temp) - 1,
                                       data = list(y = temp))
                     Yhat <- Yhat[1:length(yhat[[ii]]),]
                     error[i] <- (length(yhat[,i])-sum(Yhat*proby.temp))/n
                  }
            
                
                if(is.null(survival.time)){
                  loglik[i] <- sum(log(prob[,  , i][cbind(seq(1, n), unclass(y))]))/                        n}
                
                if(!is.null(survival.time)){
                  pvalue.survival[i]<- pvalue.survival.func(yhat[,i], survival.time,censoring.status, ngroup.survival)
                }
        }

obj<- list(threshold=threshold, error=error, loglik=loglik,size=size, yhat=yhat,y=y,prob=prob,folds=folds, cv.objects=cv.objects, pvalue.survival=pvalue.survival,
                call = this.call)
        class(obj) <- "nsccv"
        obj
}

