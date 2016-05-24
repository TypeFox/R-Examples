R2CV <- function(model){
  if (!any(class(model)=="lm")) stop("model must be class ``lm'' or ``glm''")
  yhat <- fitted(model)
  ehat <- residuals(model)
  y <- yhat + ehat
  hi <- influence(model)$hat
  yhatcve <- (yhat-hi*y)/(1-hi)
  cor(y,yhatcve)^2
}
  
