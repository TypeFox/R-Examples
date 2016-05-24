assign("rbf.cv",
  function(formula, data, eta, rho, n.neigh, func){
  z = extractFormula(formula, data, newdata=data)$z
  rbf.pred <- as.data.frame(matrix(NA,nrow= length(z), ncol=4))
  colnames(rbf.pred) <- c("x","y","var1.pred","var1.var")
  for(i in 1:length(z)){
  rbf.pred[i,3] <- rbf(formula, data[-i,], eta, rho, newdata=data[i,], n.neigh, func)[,3]
  }
RMSPE <-  sqrt(sum((rbf.pred$var1.pred-z)^2)/length(z))
RMSPE
}
)