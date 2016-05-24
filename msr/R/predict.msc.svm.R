predict.msc.svm <- function (object, newdata, ...) 
{ 
    x <- model.matrix(object, newdata)
    ms <- object
    if(length(ms$level[[ms$predictLevel]]$partitionSize) == 1 ){
      d <- matrix(1, nrow = nrow(x), ncol = 1)
    }
    else{
      if( is.null(ms$level[[ms$predictLevel]]$svm) ){
        df <- data.frame(pId = as.factor(ms$level[[ms$predictLevel]]$partition), ms$x)
        ms$level[[ms$predictLevel]]$svm <- svm(pId ~ ., df, probability=TRUE, cost=ms$cost, scale = FALSE)
        object <<- ms
      }
      d <- predict(ms$level[[ms$predictLevel]]$svm, x, probability=TRUE)
      d <- attr(d, "probabilities")
      d <- d[, order(as.numeric(colnames(d)))]
    }
    if(nrow(x) == 1){
      d <- matrix(d, nrow=1)
    }
    d 
}

