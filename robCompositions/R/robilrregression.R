robilrregression <- function(X,y){
# delivers appropriate inference for robust regression of y on a compositional matrix X
# PF, Aug 18, 2010
#
# classical regression
d <- data.frame(y=y,X=X)
lmcla <- robustbase::ltsReg(y~.,data=d)
lmcla.sum <- summary(lmcla)
# ilr regressions:
require(robCompositions)
ilr.sum <- lmcla.sum
for (j in 1:ncol(X)){
  Zj <- isomLR(cbind(X[,j],X[,-j]))
  dj <- data.frame(y=y,Z=Zj)
  res <- robustbase::ltsReg(y~.,data=dj)
  res.sum <- summary(res)
  if (j==1){
    ilr.sum$coefficients[1:2,] <- res.sum$coefficients[1:2,]
    ilr.sum$residuals <- res.sum$residuals
    ilr.sum$sigma <- res.sum$sigma
    ilr.sum$r.squared <- res.sum$r.squared
    ilr.sum$adj.r.squared <- res.sum$adj.r.squared
    ilr.sum$fstatistic <- res.sum$fstatistic
  }
  else{
    ilr.sum$coefficients[j+1,] <- res.sum$coefficients[2,]
  }
}
list(lm=lmcla.sum,ilr=ilr.sum)
}
