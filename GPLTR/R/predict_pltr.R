predict_pltr <- function(xtree, xdata, Y.name, X.names, newdata, type = "response", family = 'binomial', thresshold = seq(0.1, 0.9, by = 0.1))
{
  if(!inherits(xtree, 'rpart')) stop('xtree have to be an rpart object!')
  
  pltr_fit <- tree2glm(xtree, xdata, Y.name, X.names, family)
  
  if(family == 'binomial'){
    pred <- predict.glm(pltr_fit, newdata = newdata, type = type)
    predict_glm <- sapply(thresshold, function(uw) return(as.numeric(pred > uw)) )
    ERR_PRED <- sapply(1:length(thresshold), function(u) mean(predict_glm[, u]!= xdata[, Y.name]))
  } else{
    predict_glm <- predict.glm(pltr_fit, newdata = newdata, type = type)
    ERR_PRED <- mean((newdata[, Y.name] - predict_glm)^2)
  }
  
  return(list(predict_glm = predict_glm, ERR_PRED = ERR_PRED))
}
    
   