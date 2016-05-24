msc.nn.svm <- function(y, x, knn = 3*ncol(x), pLevelP = 0.2, pLevel, nLevels, 
                        cost = 1,  type = 1, smooth=FALSE, precompute = FALSE, 
                        eps=0.01 ) 
{
    obj <- msc.nn(y,x,knn, pLevelP, pLevel, nLevels, type = type, smooth =
        smooth, eps = eps)
    class(obj) <- paste(class(obj), ".svm", sep="")
    obj$cost <- cost
    if(precompute){
      nLevels <- obj$nLevels
      if(precompute){
        for( i in 2:nLevels ){
          df <- data.frame(cry = as.factor(obj$level[[i]]$partition), obj$x)
          obj$level[[i]]$svm <- svm(cry ~ ., df, probability=TRUE, cost = obj$cost, scale = FALSE)
        }
      }
    }
    obj
}

