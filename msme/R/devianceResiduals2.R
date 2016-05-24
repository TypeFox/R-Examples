devianceResiduals2 <- function(y, b.hat, X, p, offset = 0) {
  y.hat <- predict(y, b.hat[1:p], X[,1:p], offset)
  scale <- predict_s(y, b.hat[-(1:p)], X[,-(1:p)])
  sign(y - y.hat) *  
    sqrt(2 * getDispersion(y, scale) * 
         (   jll2(y, y,     scale) -
             jll2(y, y.hat, scale)    )) 
}
