nsc <-
  
  function(x, y = NULL, xtest = NULL, proby = NULL, ytest = NULL, prob.ytest = 
        NULL, threshold = NULL, n.threshold = 30, hetero = NULL, scale.sd = 
        TRUE, threshold.scale = NULL, se.scale = NULL, offset.percent = 50, 
        prior = table(y)/length(y), remove.zeros = TRUE, sign.contrast = "both",
           problem.type=c("class", "surv.km","surv.latent"))
{

# modified aug 2003 to add survival analysis facilities
#
#         problem.type can be "class", "surv.km", or "surv.latent"         
#
#         y= class variable => classification problem (problem.type="class")
#         proby= matrix of class probabilities => "soft classification"
#               from Kaplan-Meier estimate, for survival analysis
#         in this case, nsc computes probability-weighted centroids and
#             training error 

        this.call <- match.call()

        argy <- ytest
        if(is.null(ytest)) {
                argy <- y
        }
        if(!is.null(y) & !is.null(proby) & problem.type!="surv.latent") {
                stop("Can't specify both y and proby")
        }
        if(!is.null(ytest) & !is.null(prob.ytest)) {
                stop("Can't specify both ytest and prob.ytest")
        }
        if(is.null(y)) {
                y <- apply(proby, 1, which.is.max)
        }
        n.class <- table(y)
        if(min(n.class) == 1) {
                cat("Warning: a class contains only 1 sample")
        }
        if(is.null(xtest)) {
                xtest <- x
                ytest <- y
                prob.ytest <- proby
        }
        norm.cent <- NULL
        if(!is.null(hetero)) {
                norm.cent <- apply(x[, y == hetero], 1, mean)
                x <- abs(t(scale(t(x), center = norm.cent, scale = FALSE)))
                if(!missing(xtest)) {
                        xtest <- abs(t(scale(t(xtest), center = norm.cent, 
                                scale = FALSE)))
                }
        }
        n <- sum(n.class)
        ntest <- ncol(xtest)
        K <- length(prior)
        p <- nrow(x)
        if(is.null(proby)) {
                Y <- model.matrix( ~ factor(y) - 1, data = list(y = y))
        }
        if(!is.null(proby)) {
                Y <- proby
        }
        dimnames(Y) <- list(NULL, names(n.class))

        centroids <- scale(x %*% Y, FALSE, n.class)
        sd <- rep(1, p)
        if(scale.sd) {
                xdif <- x - centroids %*% t(Y)
                sd <- (xdif^2) %*% rep(1/(n - K), n)
                sd <- drop(sqrt(sd))
                offset <- quantile(sd, offset.percent/100)
                sd <- sd + offset
        }
        centroid.overall <- drop(x %*% rep(1/n, n))
        if(is.null(threshold.scale)) {
                threshold.scale <- rep(1, K)
                names(threshold.scale) <- names(n.class)
        }
### Now make an adjustment for the sample sizes in the "t" ratios

        if(is.null(se.scale))
                se.scale <- sqrt(1/n.class - 1/n)
        delta <- (centroids - centroid.overall)/sd
        delta <- scale(delta, FALSE, threshold.scale * se.scale)
        if(sign.contrast == "positive") {
                delta <- delta * (delta > 0)
        }
        if(sign.contrast == "negative") {
                delta <- delta * (delta < 0)
        }
#allows differential shrinkage
        if(!is.null(threshold)) {
                n.threshold <- length(threshold)
        }
        else {
                threshold <- seq(0, max(abs(delta)), length = n.threshold)
        }
        nonzero <- seq(n.threshold)
        errors <- threshold
        yhat <- as.list(seq(n.threshold))
        prob <- array(0, c(ntest, K, n.threshold))
        for(ii in 1:n.threshold) {
                cat(ii)
                delta.shrunk <- soft.shrink(delta, threshold[ii])
                delta.shrunk <- scale(delta.shrunk, FALSE, 1/(threshold.scale * 
                        se.scale))
                nonzero[ii] <- attr(delta.shrunk, "nonzero")
                posid <- drop(abs(delta.shrunk) %*% rep(1, K)) > 0
                dd <- diag.disc((xtest - centroid.overall)/sd, delta.shrunk, 
                        prior, weight = posid)
                yhat[[ii]] <- softmax(dd)
                dd <- safe.exp(dd)
                prob[,  , ii] <- dd/drop(dd %*% rep(1, K))
                if(!is.null(ytest)) {
                        errors[ii] <- sum(yhat[[ii]] != ytest)
                }
                if(!is.null(prob.ytest)) {

# use of temp below is to ensure that Yhat doesn;t drop a column
#  when no predictions are made to that class
                        temp <- c(yhat[[ii]], names(table(y)))
                        Yhat <- model.matrix( ~ factor(temp) - 1, data = list(y
                                 = temp))
                        Yhat <- Yhat[1:length(yhat[[ii]]),  ]
                     
                        errors[ii] <- length(yhat[[ii]]) - sum(Yhat * prob.ytest)
                }
               
        }
        thresh.names <- format(round(threshold, 3))
        names(yhat) <- thresh.names
        attr(yhat, "row.names") <- paste(seq(ntest))
        class(yhat) <- "data.frame"
        if(remove.zeros)
                n.threshold <- match(0, nonzero, n.threshold)
        dimnames(prob) <- list(paste(seq(ntest)), names(n.class), thresh.names)
        object <- list(y = argy, proby = prob.ytest, yhat = yhat[, seq(
                n.threshold)], prob = prob[,  , seq(n.threshold)], centroids = 
                centroids, centroid.overall = centroid.overall, sd = sd, 
                threshold = threshold[seq(n.threshold)], nonzero = nonzero[seq(
                n.threshold)], threshold.scale = threshold.scale, se.scale = 
                se.scale, scale.sd=scale.sd, call = this.call, hetero = hetero, norm.cent = 
                norm.cent, prior = prior, offset = offset, sign.contrast = 
                sign.contrast)
        if(!is.null(ytest) | !is.null(prob.ytest))
                object$errors <- errors[seq(n.threshold)]
        class(object) <- "nsc"
        object
}


safe.exp=function(x){
 xx=sign(x)*pmin(abs(x),500)
 return(exp(xx))
}
