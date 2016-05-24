Qsquared.ext <-
function(y.pre,y,ytr){
  numerator = sum((y-y.pre)^2)
  denominator= sum((y-mean(ytr))^2)
  Qsquared.ext=1-(numerator/denominator)
  return
  Qsquared.ext
}
