diagnostics <-
function(model,residuals){
  
  plot(residuals$swr2, main="Standardized residuals 2 vs. Indices of Obs")  
  abline(a=2,b=0,col=4)
  abline(a=-2,b=0,col=4)
  
#  windows()
  qqnorm(residuals$swr2)
  qqline(residuals$swr2,col=4)
  
# windows()
  plot(log(fitted(model))/(1-log(fitted(model))),residuals$swr2, main="Standardized residuals 2  vs. Linear Predictor",xlab="Linear Predictor", ylab="Standardized Residuals")
  abline(a=0,b=0,col=4)
}
