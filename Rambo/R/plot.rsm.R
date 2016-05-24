plot.rsm <-
function(x, frame = 0, ...){
  Y = x
  if(class(Y) != "rsm") stop("Entry must be from class rsm")
  
  vec <- rep(-Inf, max(Y$Klist))
  
  for(K in Y$Klist)
    vec[K] <- Y$output[[K]]$lower
  vec <- vec[-c(1:(Y$Klist[1]-1))] 
  
  if(frame == 1){
    par(mfrow=c(1,1))
    plot(Y$Klist, vec, xlab="Number of classes", ylab= "Lower bound criterion", main = "Lower bound", type="o", col="blue") 
    abline(v= Y$K_star, pch=22, col="red", lty=2)
  } 
  else if (frame == 2){
    par(mfrow=c(1,1))
    barplot(t(Y$output[[Y$K_star]]$alpha), beside=T, main="Repartition of clusters into subgraphs", ylab="Proportions", 
            legend.text=paste(rep("cluster", Y$K_star), 1:Y$K_star),
            names.arg = paste(rep("subgraph", Y$R), 1:Y$R),
            args.legend = list(x = "topright", bty="n"))
  } 
 # else if(frame == 3){
 #   par(mfrow=c(1,1))
 #   library(sna)
 #   gplot(Y$X , vertex.col = Y$output[[Y$K_star]]$Zcol)
 # }
  else{
    par(mfrow=c(2,1))
    
    plot(Y$Klist, vec, xlab="Number of classes", ylab= "Lower bound criterion", main = "Lower bound", type="o", col="blue") 
    abline(v= Y$K_star, pch=22, col="red", lty=2)
    
    barplot(t(Y$output[[Y$K_star]]$alpha), beside=T, main="Repartition of clusters into subgraphs", ylab="Proportions", xlab="",
            names.arg = paste(rep("subgraph", Y$R), 1:Y$R), legend.text=paste(rep("cluster", Y$K_star), 1:Y$K_star),
            args.legend = list(x = "topright", bty="n"))  }
}
