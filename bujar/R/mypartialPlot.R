#Partial dependence (plots) for boosting, derived from randomForest R package partialPlot.R. A similar implementation can be found in brt.functions.R for MART
#partialPlot <- function(x, ...) UseMethod("partialPlot")

#partialPlot.default <- function(x, ...)
#    stop("partial dependence plot not implemented for this class of objects.\n")
#note: this function is only intended to be useful when with interactions, i.e. degree=2 
partialPlot <-
    function (x, pred.data, x.var, which.class, w, learner="linear.regression",plot=TRUE, nseq=51, add=FALSE,
              n.pt = min(length(unique(pred.data[, xname])), nseq), rug = TRUE,
              xlab=deparse(substitute(x.var)), ylab="",
              main=paste("Partial Dependence on", deparse(substitute(x.var))),
              ...) 
{
    if(!learner%in%c("linear.regression","tensor.product"))
    stop("learner",learner,"is not implemented\n")
    classRF <- FALSE
    if(learner=="linear.regression"){
     d <- ncol(pred.data)
     p <- -0.5+0.5*sqrt(1+8*d)
    }
#    classRF <- x$type != "regression"
#    if (is.null(x$forest)) 
#        stop("The randomForest object must contain the forest.\n")
#    x.var <- substitute(x.var)
    xname <- if (is.character(x.var)) x.var else {
        if (is.name(x.var)) deparse(x.var) else {
            eval(x.var)
        }
    }
    xv <- pred.data[, xname]
    n <- nrow(pred.data)
    nreplace <- ifelse(learner=="tensor.product", floor(n*3/5), n)       # I had to do this because if replacing all, the predicted values are all NaN
    if (missing(w)) w <- rep(1, n)
    if (classRF) {
        if (missing(which.class)) {
            focus <- 1
        }
        else {
            focus <- charmatch(which.class, colnames(x$votes))
            if (is.na(focus)) 
                stop(which.class, "is not one of the class labels.")
        }
    }
    if (is.factor(xv) && !is.ordered(xv)) {
        x.pt <- levels(xv)
        y.pt <- numeric(length(x.pt))
        for (i in seq(along = x.pt)) {
            x.data <- pred.data
            x.data[, xname] <- factor(rep(x.pt[i], n), levels = x.pt)
            if(learner=="linear.regression"){
              x.data <- x.data[,1:p]   #main effects
              for(j in 1:(p-1))
               for (k in (j+1):p)
                x.data <- cbind(x.data,x.data[,j]*x.data[,k])
            } 
            if (classRF) {
                pr <- predict(x, x.data[1:nreplace,], type = "prob")
                y.pt[i] <- weighted.mean(log(ifelse(pr[, focus] > 0,
                                                    pr[, focus], 1)) -
                                         rowMeans(log(ifelse(pr > 0, pr, 1))),
                                         w[1:nreplace], na.rm=TRUE)
            } else 
               if(learner=="linear.regression") 
                 y.pt[i] <- weighted.mean(predict(x, cbind(1,x.data))[1:nreplace], w[1:nreplace], na.rm=TRUE)
               else y.pt[i] <- weighted.mean(predict(x, x.data)[1:nreplace], w[1:nreplace], na.rm=TRUE)
        }
        if (add) {
            points(1:length(x.pt), y.pt, type="h", lwd=2, ...)
        } else {
            if (plot) barplot(y.pt, width=rep(1, length(y.pt)), col="blue",
                              xlab = xlab, ylab = ylab, main=main,
                              names.arg=x.pt, ...)
        }
    } else {
        if (is.ordered(xv)) 
            xv <- as.numeric(xv)
        x.pt <- seq(min(xv), max(xv), length = n.pt)
        y.pt <- numeric(length(x.pt));
        for (i in seq(along = x.pt)) {
            x.data <- pred.data
            x.data[1:nreplace, xname] <- rep(x.pt[i], nreplace)
            if(learner=="linear.regression"){
              x.data <- x.data[,1:p]   #main effects
              for(j in 1:(p-1))
               for (k in (j+1):p)
                x.data <- cbind(x.data,x.data[,j]*x.data[,k])
            } 
            if (classRF) {
                pr <- predict(x, x.data, type = "prob")
                y.pt[i] <- weighted.mean(log(ifelse(pr[, focus] == 0, 1, pr[, focus]))
                                         - rowMeans(log(ifelse(pr == 0, 1, pr))),
                                         w, na.rm=TRUE)
            } else {
               if(learner=="linear.regression") 
                 y.pt[i] <- weighted.mean(predict(x, cbind(1,x.data))[1:nreplace], w[1:nreplace], na.rm=TRUE)
               else y.pt[i] <- weighted.mean(predict(x, x.data)[1:nreplace], w[1:nreplace], na.rm=TRUE)
   #cat("i=",i,"",y.pt[i],"\n")
            }
        }
        if (add) {
            lines(x.pt, y.pt, ...)
        } else {
            if (plot) plot(x.pt, y.pt, type = "l", xlab=xlab, ylab=ylab,
                           main = main, ...)
        }
        if (rug && plot) {
            if (n.pt > 10) {
                rug(quantile(xv, seq(0.1, 0.9, by = 0.1)), side = 1)
            } else {
                rug(unique(xv, side = 1))
            }
        }
    }
    invisible(list(x = x.pt, y = y.pt))
}
