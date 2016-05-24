#' Model interpretation functions: Partial Dependence Plots
#'
#' \code{parDepPlot} creates partial dependence plots for binary (cross-validated) classification models and regression models. Currently only binary classification models estimated with the packages \code{randomForest} and \code{ada} are supported.  In addition \code{randomForest} regression models are supported.
#'
#' @param x.name the name of the predictor as a character string for which a partial dependence plot has to be created.
#' @param object can be a model or a list of cross- validated models. Currently only binary classification models built using the packages \code{randomForest} and \code{ada} are supported.
#' @param data a data frame containing the predictors for the model or a list of data frames for cross-validation with length equal to the number of models.
#' @param rm.outliers boolean, remove the outliers in x.name. Outliers are values that are smaller than max(Q1-fact*IQR,min) or greater than min(Q3+fact*IQR,max). Overridden if xlim is used.
#' @param fact factor to use in rm.outliers. The default is 1.5.
#' @param n.pt if x.name is a continuous predictor, the number of points that will be used to plot the curve.
#' @param robust if TRUE then the median is used to plot the central tendency (recommended when logit=FALSE). If FALSE the mean is used. 
#' @param ci boolean. Should a confidence interval based on quantiles be plotted? This only works if robust=TRUE.
#' @param u.quant Upper quantile for ci. This only works if ci=TRUE and robust=TRUE.
#' @param l.quant Lower quantile for ci. This only works if ci=TRUE and robust=TRUE.
#' @param xlab  label for the x-axis. Is determined automatically if NULL.
#' @param ylab	label for the y-axis.
#' @param main main title for the plot.
#' @param logit boolean. Should the y-axis be on a logit scale or not? If FALSE, it is recommended to set robust=TRUE. Only applicable for classifcation.
#' @param ylim The y limits of the plot
#' @param ...  other graphical parameters for \code{plot}.
#' @references The code in this function uses part of the code from the \code{partialPlot} function in \code{randomForest}. It is expanded and generalized to support cross-validation and other packages. 
#' @details For classification, the response variable in the model is always assumed to take on the values \{0,1\}. Resulting partial dependence plots always refer to class 1. Whenever strange results are obtained the user has three options. First set rm.outliers=TRUE. Second, if that doesn't help, set robust=TRUE. Finally, if that doesn't help, the user can also try setting ci=TRUE. Areas with larger confidence intervals typically indicate problem areas. These options help the user tease out the root of strange results and converge to better parameter values.
#' @examples
#' 
#' library(randomForest)
#' #Prepare data
#' data(iris)
#' iris <- iris[1:100,]
#' iris$Species <- as.factor(ifelse(factor(iris$Species)=="setosa",0,1))
#' 
#' #Cross-validated models
#' #Estimate 10 models and create 10 test sets
#' data <- list()
#' rf <- list()
#' for (i in 1:10) {
#'   ind <- sample(nrow(iris),50)
#'   rf[[i]] <- randomForest(Species~., iris[ind,])
#'   data[[i]] <- iris[-ind,]
#' }
#' 
#' 
#' parDepPlot(x.name="Petal.Width", object=rf, data=data)
#' 
#' #Single model
#' #Estimate a single model
#' ind <- sample(nrow(iris),50)
#' rf <- randomForest(Species~., iris[ind,])
#' parDepPlot(x.name="Petal.Width", object=rf, data=iris[-ind,])
#' 
#' @seealso \code{\link{variableImportance}}
#' @author Authors: Michel Ballings, and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
parDepPlot <-
    function (x.name,
              object, 
              data,
              rm.outliers=TRUE,
              fact=1.5,
              n.pt=50,
              robust=FALSE,
              ci=FALSE,
              u.quant=0.75,
              l.quant=0.25,
              xlab=substr(x.name,1,50), 
              ylab=NULL,
              main= if (any(class(object) %in% c("randomForest","ada"))) paste("Partial Dependence on",substr(x.name,1,20)) else paste("Cross-Validated Partial Dependence on",substr(x.name,1,10)),
              logit=TRUE,
              ylim=NULL,
              ...)
{


      
regression <- FALSE
if (is.null(ylab)==TRUE && robust==FALSE) {
        if (logit==TRUE)  ylab <- if (any(class(object) %in% c("randomForest","ada"))) bquote(paste('mean(0.5*logit(P'[1],'))')) else bquote(paste('CV mean(0.5*logit(P'[1],'))'))
        if (logit==FALSE) ylab <- if (any(class(object) %in% c("randomForest","ada"))) bquote(paste('mean(P'[1],')')) else bquote(paste('CV mean(P'[1],')'))
        #regression randomForest
        if((any(class(object[[1]])=="randomForest") || any(class(object)=="randomForest")) && c(object[[1]]$type,object$type)=="regression") {
          ylab <- if (any(class(object) %in% c("randomForest","ada"))) paste('mean(Prediction)') else paste('CV mean(Prediction)')
        }
} else if (is.null(ylab)==TRUE && robust==TRUE) {
        if (logit==TRUE)  ylab <- if (any(class(object) %in% c("randomForest","ada"))) bquote(paste('median(0.5*logit(P'[1],'))')) else bquote(paste('CV median(0.5*logit(P'[1],'))'))
        if (logit==FALSE) ylab <- if (any(class(object) %in% c("randomForest","ada"))) bquote(paste('median(P'[1],')')) else bquote(paste('CV median(P'[1],')'))
        #regression randomForest
        if((any(class(object[[1]])=="randomForest") || any(class(object)=="randomForest")) && c(object[[1]]$type,object$type)=="regression") {
          ylab <- if (any(class(object) %in% c("randomForest","ada"))) paste('median(Prediction)') else paste('CV median(Prediction)')
        }
}      

#only randomForest supports regression
if (any(class(object) == "randomForest")==TRUE) {
  if (c(object[[1]]$type,object$type)=="regression") regression <- TRUE
}
      
if (robust==FALSE) ci <- FALSE

#solve labeling issue by assigning labels before anyting else  
ylab <- ylab
main <- main

#if object has only one model (=object has a class) then assign object to a list. Same for data.
#If there are multiple models (=object has no class) then do nothing
if (any(class(object) %in% c("randomForest","ada"))) {
    object <- list(object)
    data <- list(data)
}

if (length(table(data[[1]][,x.name]))==1) stop("The variable is a constant and therefore a partial dependence plot cannot be created. Make sure to have at least two values.")

xy <- vector(mode="list", length=length(object))
xub <- vector(mode="list", length=length(object))
xlb <- vector(mode="list", length=length(object))

# x <- as.numeric(x)
# y <- as.numeric(y)

# if ('xlim' %in% names(args)) {
#   y <- y[ x >= min(args$xlim) & x <= max(args$xlim) ]  
#   x <- x[ x >= min(args$xlim) & x <= max(args$xlim) ]    
# }      

  

      
args <- list(...)
if ('xlim' %in% names(args)) {
  x <- seq(min(args$xlim), max(args$xlim), length = n.pt )
#   length = min(length(unique(data[[1]][,x.name])), n.pt)
} else {
#get min and max of predictor across all data sets for numeric predictors
if (!is.factor(data[[1]][,x.name]) && length(unique(data[[1]][,x.name]))>2 ) { 
    ranges <- matrix(NA,ncol=2,nrow=length(object))
    for (ii in 1:length(object)) {
      
      ranges[ii,] <- range(data[[ii]][,x.name])
      if (rm.outliers==TRUE) {
        IQRval <- fact*IQR(data[[ii]][,x.name])
        ranges[ii,] <- c(max(quantile(data[[ii]][,x.name],probs=0.25)-IQRval,min(data[[ii]][,x.name])),
                         min(quantile(data[[ii]][,x.name],probs=0.75)+IQRval,max(data[[ii]][,x.name])))
      }
    }
x <- seq(min(ranges[,1]), max(ranges[,2]), length = n.pt)
} else if (!is.factor(data[[1]][,x.name]) && length(table(data[[1]][,x.name]))==2 ) {
  x <- c(0,1)
  }
}


if (regression==TRUE) {

#start loop through models
for (ii in 1:length(object)) {
              
              predictor <- data[[ii]][,x.name]
              n <- nrow(data[[ii]])
            
              if (is.factor(predictor) ) { #if predictor is a factor
                  x <- levels(predictor)
                  y <- numeric(length(x))
                  ub <- numeric(length(x))
                  lb <- numeric(length(x))
                  for (i in seq(along = x)) {
                      data.after <- data[[ii]]
                      data.after[, x.name] <- factor(rep(x[i], n), levels = x)
          
                          pred <- predict(object[[ii]], data.after)
                          
                      
                      if (robust==TRUE) {
                          
                          
                          y[i] <- median(pred, na.rm=TRUE)
                          
                          if (ci==TRUE){

                          ub[i] <- quantile(pred, probs=  u.quant ,na.rm=TRUE)
                          lb[i] <- quantile(pred, probs=  l.quant , na.rm=TRUE)
 
                         
                          }
                      } else {
                        
                          y[i] <- mean(pred , na.rm=TRUE)
                        
                      }
                           
                  }
   
              } else { #if predictor is a numeric
                 

                  y <- numeric(length(x))
                  ub <- numeric(length(x))
                  lb <- numeric(length(x))
                  for (i in seq(along.with = x)) {
                      data.after <- data[[ii]]
                      data.after[, x.name] <- rep(x[i], n)
                      
                       pred <- predict(object[[ii]], data.after)
                      
                          
                      if (robust==TRUE) {
                     
                          y[i] <- median(pred, na.rm=TRUE)
                          
                          if (ci==TRUE){
                          
                          ub[i] <- quantile(pred, probs= u.quant , na.rm=TRUE)
                          lb[i] <- quantile(pred, probs= l.quant , na.rm=TRUE)
                          
                          }
                          
                      } else {
                        
                          y[i] <- mean(pred, na.rm=TRUE)
                        
                      }
          
                      
                    }
                  }
              
          #store x and y as a data.frame in one list
          xy[[ii]] <- data.frame(x,y)
          if (ci==TRUE){
          xub[[ii]] <- data.frame(x,ub)
          xlb[[ii]] <- data.frame(x,lb)
          }

#end loop through models
}  
  
} else if (regression==FALSE){

#start loop through models
for (ii in 1:length(object)) {
              
              predictor <- data[[ii]][,x.name]
              n <- nrow(data[[ii]])
            
              if (is.factor(predictor) ) { #if predictor is a factor
                  x <- levels(predictor)
                  y <- numeric(length(x))
                  ub <- numeric(length(x))
                  lb <- numeric(length(x))
                  for (i in seq(along = x)) {
                      data.after <- data[[ii]]
                      data.after[, x.name] <- factor(rep(x[i], n), levels = x)
          
                          if(any(class(object[[ii]]) %in% "randomForest")) pred <- predict(object[[ii]], data.after, type = "prob")
                          if(any(class(object[[ii]])=="ada")) pred <- predict(object[[ii]], data.after, type="probs")
                      
                      if (robust==TRUE) {
                          
                          if (logit==FALSE) {
                          y[i] <- median(pred[, 2], na.rm=TRUE)
                          } else {
                          y[i] <- median(log(ifelse(pred[, 2] > 0,
                                                              pred[, 2], .Machine$double.eps)) -
                                                   rowMeans(log(ifelse(pred > 0, pred, .Machine$double.eps))),
                                                   na.rm=TRUE)
                          }
                          if (ci==TRUE){
                          if (logit==FALSE) {
                          ub[i] <- quantile(pred[, 2], probs=  u.quant ,na.rm=TRUE)
                          lb[i] <- quantile(pred[, 2] , probs=  l.quant , na.rm=TRUE)
                          } else {
                          ub[i] <- quantile(log(ifelse(pred[, 2] > 0,
                                                    pred[, 2], .Machine$double.eps)) -
                                           rowMeans(log(ifelse(pred > 0, pred, .Machine$double.eps))), probs=  u.quant ,
                                         na.rm=TRUE)
                          
                          lb[i] <- quantile(log(ifelse(pred[, 2] > 0,
                                                       pred[, 2], .Machine$double.eps)) -
                                              rowMeans(log(ifelse(pred > 0, pred, .Machine$double.eps))), probs=  l.quant ,
                                            na.rm=TRUE)
                          }
                          }
                      } else {
                        if (logit==FALSE) {
                          y[i] <- mean(pred[, 2] , na.rm=TRUE)
                        } else {
                          y[i] <- mean(log(ifelse(pred[, 2] > 0,
                                                pred[, 2], .Machine$double.eps)) -
                                       rowMeans(log(ifelse(pred > 0, pred, .Machine$double.eps))),
                                       na.rm=TRUE)
                        }
                      }
                           
                  }
                
              } else { #if predictor is a numeric
                 

                  y <- numeric(length(x))
                  ub <- numeric(length(x))
                  lb <- numeric(length(x))
                  for (i in seq(along.with = x)) {
                      data.after <- data[[ii]]
                      data.after[, x.name] <- rep(x[i], n)
                      
                      if(any(class(object[[ii]]) %in% "randomForest"))  pred <- predict(object[[ii]], data.after, type = "prob")
                      if(any(class(object[[ii]])=="ada")) pred <- predict(object[[ii]], data.after, type="probs")
                          
                      if (robust==TRUE) {
                          if (logit==FALSE) {
                          y[i] <- median(pred[, 2], na.rm=TRUE)
                          } else {
                          y[i] <- median(log(ifelse(pred[, 2] == 0,
                                                .Machine$double.eps, pred[, 2]))
                                     - rowMeans(log(ifelse(pred == 0, .Machine$double.eps, pred))),
                                     na.rm=TRUE)
                          }
                          if (ci==TRUE){
                          if (logit==FALSE) {
                          ub[i] <- quantile(pred[, 2], probs= u.quant , na.rm=TRUE)
                          lb[i] <- quantile(pred[, 2], probs= l.quant , na.rm=TRUE)
                          } else {
                          ub[i] <- quantile(log(ifelse(pred[, 2] == 0,
                                                    .Machine$double.eps, pred[, 2]))
                                         - rowMeans(log(ifelse(pred == 0, .Machine$double.eps, pred))), probs= u.quant ,
                                         na.rm=TRUE)
                          
                          lb[i] <- quantile(log(ifelse(pred[, 2] == 0,
                                                    .Machine$double.eps, pred[, 2]))
                                         - rowMeans(log(ifelse(pred == 0, .Machine$double.eps, pred))), probs= l.quant ,
                                         na.rm=TRUE)
                          }
                          }
                          
                      } else {
                        if (logit==FALSE) {
                          y[i] <- mean(pred[, 2], na.rm=TRUE)
                        } else {
                          y[i] <- mean(log(ifelse(pred[, 2] == 0,
                                                              .Machine$double.eps, pred[, 2]))
                                                   - rowMeans(log(ifelse(pred == 0, .Machine$double.eps, pred))),
                                                    na.rm=TRUE)
                        }
                      }
          
                      
                    }
                  }
              
          #store x and y as a data.frame in one list
          xy[[ii]] <- data.frame(x,y)
          if (ci==TRUE){
          xub[[ii]] <- data.frame(x,ub)
          xlb[[ii]] <- data.frame(x,lb)
          }

#end loop through models
}
#if not regression
}
      
#merge all elements of the list by x 
#compute mean of y by x
options(warn=-1)
xy <- data.frame(Reduce(function(x, y) merge(x, y, by='x',all=TRUE),xy, accumulate=FALSE))
if (ci==TRUE){
  xub <- data.frame(Reduce(function(x, y) merge(x, y, by='x',all=TRUE),xub, accumulate=FALSE))
  xlb <- data.frame(Reduce(function(x, y) merge(x, y, by='x',all=TRUE),xlb, accumulate=FALSE))
}
x <- xy[,1]

y <- rowMeans(data.frame(xy[,2:ncol(xy)]),na.rm=TRUE)
if (ci==TRUE){
  ub <- rowMeans(data.frame(xub[,2:ncol(xub)]),na.rm=TRUE)
  lb <- rowMeans(data.frame(xlb[,2:ncol(xlb)]),na.rm=TRUE)
}

options(warn=0)



if (ci==FALSE){
  if (is.factor(predictor)==FALSE && length(unique(x)) > 2){
    plot(as.numeric(x), y, type = "l", xlab=xlab, ylab=ylab, ylim=ylim, main = main, ...)
  } else {
    plot(as.numeric(x), y, type = "l", xlab=xlab, ylab=ylab, ylim=ylim, main = main, xaxt="n",...)
    axis(1,at=as.numeric(x), labels=x)
  }
}
else if (ci==TRUE){
  if (is.factor(predictor)==FALSE && length(unique(x)) > 2){
    if (is.null(ylim)) ylim <- c(min(lb),max(ub))
    plot(as.numeric(x), y, type = "l", xlab=xlab, ylab=ylab, ylim=ylim, main = main, ...)
    polygon(c(as.numeric(x),rev(as.numeric(x))),c(lb,rev(ub)),col = "grey80", border = FALSE)
    lines(as.numeric(x), y, type = "l", xlab=xlab, ylab=ylab, main = main)
  } else {
    if (is.null(ylim)) ylim <- c(min(lb),max(ub))
    plot(as.numeric(x), y, type = "l", xlab=xlab, ylab=ylab, ylim=ylim, main = main, xaxt="n",...)
    polygon(c(as.numeric(x),rev(as.numeric(x))),c(lb,rev(ub)),col = "grey80", border = FALSE)
    lines(as.numeric(x), y, type = "l", xlab=xlab, ylab=ylab, main = main)
    axis(1,at=as.numeric(x), labels=x)
  }
}
      invisible(list(x = x, y = y))
}