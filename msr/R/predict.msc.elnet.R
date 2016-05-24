predict.msc.elnet <- function (object, newdata, ...) 
{    
    ms = object$ms 
    msLevel = ms$level[[ms$predictLevel]]
    lms = object$elnet[[ms$predictLevel]]

    x <- model.matrix(ms, newdata)

    colnames(x) <- colnames(ms$x)
 
    #predict model for each crystals 
    f = matrix(ncol=length(msLevel$mins), nrow=nrow(x));
    for(i in 1:length(lms$lm)){
      f[,i] = predict(lms$lm[[i]], x)
    }
    
    #piecewise linear or blending of the models
    d <- predict(ms, x)
    if(object$blend){ 
      f <- f * d
      predicted <- rowSums(f)
    }
    else{
      predicted <- vector(mode = "double", length = nrow(f)) 
      for(i in 1:nrow(f)){
        ind = which.max(d[i,])
        predicted[i] = f[i, ind]
      }
    }
    predicted

}

