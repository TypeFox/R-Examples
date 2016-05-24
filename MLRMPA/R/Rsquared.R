Rsquared <-
function(y.pre,y){
  SSR = sum((y.pre-mean(y))^2)
  SSE = sum((y-y.pre)^2)         
  R2 = SSR/(SSE+SSR)
  return
  R2
}
