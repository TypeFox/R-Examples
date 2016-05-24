###############################################################
# baseline survival/ hazard Breslow estimator
# function essentially based on gbm::basehaz.gbm
###############################################################
basesurv <- function (response, lp, times.eval = NULL, centered = FALSE)
{
    if (is.null(times.eval)) times.eval <- sort(unique(response[,1]))
    
    t.unique <- sort(unique(response[,1][response[,2] == 1]))
    alpha    <- length(t.unique)

    for (i in 1:length(t.unique)) {
        alpha[i] <- sum(response[,1][response[,2] == 1] == t.unique[i])/sum(exp(lp[response[,1] >=  t.unique[i]]))
    }

    obj   <- approx(t.unique, cumsum(alpha), yleft=0, xout = times.eval, rule=2)

    if (centered) obj$y <- obj$y * exp(mean(lp))
    obj$z <- exp(-obj$y)

    names(obj) <- c("times","cumBaseHaz","BaseSurv")
    return(obj)
}

###############################################
### wrapper for glmnet
###############################################

fit.glmnet <- function (response, x, cplx, ...) 
{
    #require(glmnet)
    res <- NULL
    tryerr <- try(res <- glmnet(y = response, x = data.matrix(x), lambda = cplx,  ...), silent=TRUE)

    if(class(tryerr) != 'try-error' && "coxnet" %in% class(res)) {
          res$linear.predictor  <- as.numeric(predict(res, newx=data.matrix(x), type="link"))
          res$response          <- response
    }
    res
}

complexity.glmnet <- function (response, x, full.data, ...) 
{
    #require(glmnet)
    lambda <- NULL
    tryerr <- try(cv <- cv.glmnet(y = response, x = data.matrix(x),  ...), silent=TRUE)
    
    if(class(tryerr) != 'try-error'){
      lambda <-cv$lambda.min
    }    
    lambda
}

predictProb.coxnet <- predictProb.glmnet <- function (object, response, x, times, complexity,  ...) 
{
    #require(glmnet)    
    lp       <- as.numeric(predict(object, newx=data.matrix(x),s=complexity, type="link"))
    basesurv <- basesurv(object$response,object$linear.predictor, sort(unique(times)))
    p        <- exp(exp(lp) %*% -t(basesurv$cumBaseHaz))
    
    if (NROW(p) != NROW(x) || NCOL(p) != length(times)) 
        stop("Prediction failed")
    p
}

PLL.coxnet <- function(object, newdata, newtime, newstatus, complexity, ...) 
{
   #require(glmnet)
   PLL <- glmnet::coxnet.deviance(pred = NULL, Surv(newtime,newstatus), x = data.matrix(newdata), offset = NULL, weights = NULL, beta = coef(object,s=complexity)) 
   PLL / -2
}


##################################################
###### classification aggregation functions   ####
##################################################

aggregation.misclass <- function (full.data = NULL, response, x, model, cplx = NULL, 
    type = c("apparent", "noinf"), fullsample.attr = NULL, ...) 
{
    data <- as.data.frame(x)
    data$response <- response

    if ("glmnet" %in% class(model)) {
        probs <- as.numeric(predict(model, newx = data.matrix(x), type="response", ...))
    }
    else if (class(model)[1] == "penfit") {
      probs <- predict(model, data = data, penalized = x, ...)
    }
    else if (class(model)[1]=="glm") {
        probs <- predict(model, newdata = data, type="response", ...)
    }
    else {
        probs <- predict(model, data = data, type = "response", 
            ...)
    }
    type <- match.arg(type)
    if (type == "apparent") {
        mr <- sum(abs(round(probs) - response))/length(response)
    }
    if (type == "noinf") {
        mr <- mean(abs((matrix(response, length(response), length(response), 
            byrow = TRUE) - round(probs))))
    }
    mr
}

aggregation.brier <- function (full.data = NULL, response, x, model, cplx = NULL, 
    type = c("apparent", "noinf"), fullsample.attr = NULL, ...) 
{
    data          <- as.data.frame(x)
    data$response <- response

    if ("glmnet" %in% class(model)) {
        probs <- as.numeric(predict(model, newx = data.matrix(x), type="response", ...))
    }
    else if (class(model)[1] == "penfit") {
      probs <- predict(model, data = data, penalized = x, ...)
    }
    else if (class(model)[1]=="glm") {
        probs <- predict(model, newdata = data, type="response", ...)
    }    
    else {
        probs <- predict(model, data = data, type = "response", 
            ...)
    }
    type <- match.arg(type)
    if (type == "apparent") {
        brier.score <- sum((probs - response)^2)/length(response)
    }
    if (type == "noinf") {
        brier.score <- mean((matrix(response, length(response), 
            length(response), byrow = TRUE) - probs)^2)
    }
    brier.score
}

aggregation.auc <- function (full.data = NULL, response, x, model, cplx = NULL, 
    type = c("apparent", "noinf"), fullsample.attr = NULL, ...) 
{
    data <- as.data.frame(x)
    data$response <- response
    if ("glmnet" %in% class(model)) {
        probs <- as.numeric(predict(model, newx = data.matrix(x), type="response", ...))
    }
    else if (class(model)[1] == "penfit") {
      probs <- predict(model, data = data, penalized = x, ...)
    }    
    else if (class(model)[1]=="glm") {
        probs <- predict(model, newdata = data, type="response", ...)
    }    
    else {
        probs <- predict(model, data = data, type = "response", 
            ...)
    }
    type <- match.arg(type)
    if (type == "apparent") {
        auc <- glmnet::auc(response,probs)
    }
    if (type == "noinf") {
        resp.mat <- matrix(response, length(response),  length(response), byrow = TRUE)
        auc      <- mean(apply(resp.mat, 1, function(d) glmnet::auc(d,probs)))
    }
    auc
}

########################
### plot pecs        ###
########################

Plot.peperr.curves <- function(x,at.risk=TRUE,allErrors=FALSE,  bootRuns=FALSE, bootQuants=TRUE, bootQuants.level=0.95, leg.cex=0.7, ...) {
  
  #require(peperr)

  if (bootRuns) bootQuants <- FALSE
  
  plot(x$attribute, x$null.model, type = "n", las=1,
       col = "blue", xlab = "Evaluation time points", 
       ylab = "Prediction error", main = "Prediction error curves", 
       ylim = c(0, max(perr(x), x$full.apparent, x$null.model) + 0.1))
  
  if (length(x$sample.error) > 1 & bootRuns==TRUE) {
    for (i in 1:(length(x$sample.error))) {
      lines(x$attribute, x$sample.error[[i]], type = "l", col = "light grey", lty = 1)
    }
  }

  if (length(x$sample.error) > 1 & bootQuants==TRUE) {
    boots  <- do.call("rbind",x$sample.error)
    quants <- apply(boots, 2, function(d) quantile(d, probs=c((1-bootQuants.level)/2,1 - (1-bootQuants.level)/2)))
    polygon(c(x$attribute,rev(x$attribute)),c(quants[1,],rev(quants[2,])), col="light grey", border="light grey")
  }
  
  if (allErrors==FALSE) {
     lines(x$attribute, x$null.model, type = "l", col = "blue", lwd = 2, lty = 1)
     lines(x$attribute, perr(x, "632p"), type = "l", col= "black", lty = 1, lwd = 2)
     lines(x$attribute, x$full.apparent, type = "l", col = "red", lty = 1, lwd = 2)
     if (bootRuns==TRUE) {
       legend(x = "topleft", col = c("blue", "black", "red", "light grey"), lwd=c(2,2,2,1), cex=leg.cex,
            lty = 1, legend = c("Null model", ".632+ estimate", "Full apparent", "Bootstrap samples"))
     } else {
       legend(x = "topleft", col = c("blue", "black", "red"), lwd=c(2,2,2), cex=leg.cex,
              lty = 1, legend = c("Null model", ".632+ estimate", "Full apparent"))
     }   
  }

  if (allErrors==TRUE) {
    lines(x$attribute, x$null.model, type = "l", col = "blue", lwd = 2, lty = 1)
    lines(x$attribute, perr(x, "632p"), type = "l", col= "black", lty = 1, lwd = 2)
    lines(x$attribute, perr(x, "632"), type = "l", col= "brown", lty = 1, lwd = 2)
    lines(x$attribute, perr(x, "NoInf"), type = "l", col= "green", lty = 1, lwd = 2)
    lines(x$attribute, perr(x, "resample"), type = "l", col= "dark grey", lty = 1, lwd = 2)
    lines(x$attribute, x$full.apparent, type = "l", col = "red", lty = 1, lwd = 2)
    if (bootRuns==TRUE) {
      legend(x = "topleft", ncol=2, col = c("blue", "black","brown","green","dark grey","red", "light grey"), lwd=c(2,2,2,2,2,2,1), cex=leg.cex,
           lty = 1, legend = c("Null model", ".632+ estimate",".632 estimate", "No Information","Out-of-bag average","Full apparent", "Bootstrap samples"))
    } else {
      legend(x = "topleft", ncol=2, col = c("blue", "black","brown","green","dark grey","red"), lwd=c(2,2,2,2,2,2), cex=leg.cex,
             lty = 1, legend = c("Null model", ".632+ estimate",".632 estimate", "No Information","Out-of-bag average","Full apparent"))
    }
  }
  
  if (at.risk) {
     tmpxaxp   <- par("xaxp")
     tmpusr    <- par("usr")
     at.loc    <- seq(tmpxaxp[1],tmpxaxp[2],length=tmpxaxp[3]+1)
     n.at.risk <- summary(survfit(x$response ~ 1),times=at.loc)$n.risk     

     text(x=at.loc, y=tmpusr[3], labels=n.at.risk, cex=0.8, pos=3)
     text(x=tmpxaxp[2]+(tmpusr[2]-tmpxaxp[2])/2, y=tmpusr[3], labels="at\nrisk", cex=0.8, pos=3)
  }  
}

