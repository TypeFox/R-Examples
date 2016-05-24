"calc.difference" <-
function(obs){
   n.obs <- length(obs)
   obs.mat1 <- 
matrix(rep(obs,each=n.obs),nrow=n.obs,ncol=n.obs,byrow=TRUE)
   obs.mat2 <- t(obs.mat1)
   diff.matrix <- (obs.mat1-obs.mat2)^2
return(diff.matrix)
}
