estimPVal <- function(object,x,y,permute.n=10,per.covariate=FALSE,parallel=FALSE,multicore=FALSE,trace=FALSE,...) {
    if (is.null(object$scoremat) || length(object$predictors) > 1) {
        stop("Uncertainty can only be estimated for binary response models with score-based selection and linear components only")
    }
    
    if (object$standardize.linear) x <- scale(x,center=object$mean.linear,scale=object$sd.linear)

    permute.index <- matrix(NA,permute.n,nrow(x))    
    for (actual.permute in 1:permute.n) permute.index[actual.permute,] <- sample(nrow(x))

    unpen.index <- which(object$penalty.linear == 0)
    if (length(unpen.index) > 0) {
        pen.index <- (1:ncol(x))[-unpen.index]
    } else {
        pen.index <- 1:ncol(x)
    }
    stepno <- object$stepno
    penalty <- c(object$penalty.linear,object$penalty.linear[pen.index])
    family <- object$family

    if (stepno == 0) return(rep(NA,length(pen.index)))

    eval.permute <- function(actual.permute,...) {
        if (trace) cat("permutation",actual.permute,"\n")

        actual.x <- cbind(x,x[permute.index[actual.permute,],pen.index])
        
        permute.res <- GAMBoost(y=y,x.linear=actual.x,penalty.linear=penalty,family=family,
                                standardize.linear=FALSE,stepno=stepno,
                                calc.hat=FALSE,calc.se=FALSE,
                                criterion="score",return.score=TRUE,trace=FALSE,...)
    
        actual.score <- colMeans(permute.res$scoremat[,1:length(pen.index),drop=FALSE])
        null.score <- colMeans(permute.res$scoremat[,(length(pen.index)+1):ncol(permute.res$scoremat),drop=FALSE])

        if (per.covariate) {
            return(actual.score <= null.score)
        } else {
            return(unlist(lapply(actual.score,function(arg) mean(arg <= null.score))))
        }
    }

    done.parallel <- FALSE

    if (parallel) {
        if (!require(snowfall)) {
            warning("package 'snowfall' not found, i.e., parallelization cannot be performed")
        } else {
            snowfall::sfLibrary(GAMBoost)
            snowfall::sfExport("y","x","family","unpen.index","pen.index","stepno","penalty","trace","permute.index")
            permute.mat <- matrix(unlist(snowfall::sfClusterApplyLB(1:permute.n,eval.permute,...)),length(pen.index),permute.n)
            done.parallel <- TRUE            
        }
    } 

    if (!done.parallel & multicore) {
        if (!require(parallel)) {
            warning("package 'parallel' not found, i.e., parallelization cannot be performed using this package")
        } else {
            if (multicore > 1) {
                permute.mat <- matrix(unlist(mclapply(1:permute.n,eval.permute,mc.preschedule=FALSE,mc.cores=multicore,...)),length(pen.index),permute.n)
            } else {
                permute.mat <- matrix(unlist(mclapply(1:permute.n,eval.permute,mc.preschedule=FALSE,...)),length(pen.index),permute.n)
            }
            done.parallel <- TRUE
        }        
    }

    if (!done.parallel) {
        permute.mat <- matrix(unlist(lapply(1:permute.n,eval.permute)),length(pen.index),permute.n)
    }

    rowMeans(permute.mat)
}
