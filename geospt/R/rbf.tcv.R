assign("rbf.tcv",
   function(formula, data, eta, rho, n.neigh, func){
   z = extractFormula(formula, data, newdata=data)$z
   s = coordinates(data)
   rbf.pred <- as.data.frame(matrix(NA,nrow= length(z), ncol=8))
   colnames(rbf.pred) <- c("var1.pred","var1.var","observed","residual","zscore","fold","x","y")
   pb <- txtProgressBar(min = 0, max = length(z), char = "=", style = 3)
   for(i in 1:(length(z))){
   rbf.pred[i,1] <- rbf(formula, data[-i,], eta, rho, newdata=data[i,], n.neigh, func)[,3]
   rbf.pred[i,6] <- i
   setTxtProgressBar(pb, i)
}             
close(pb)
rbf.pred[,3]<- z
rbf.pred[,7:8]<-s
rbf.pred[,4]<- rbf.pred[,3]-rbf.pred[,1]
rbf.pred
}
)