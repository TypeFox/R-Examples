# FUNCTION TO RESCALE LATENT COMPONENTS
rescaleB <- function(x,
                     Aoi,
                     Boi){
 pplusj <- colSums(x)/sum(x)

 piplus <- rowSums(x)/sum(x)
 # The margins of the unconditional proportions

 K <- ncol(Boi)
 J <- nrow(Boi)

 lbudget <- paste("LB", 1:K)
 # Labels of budgets

 pk <- piplus %*% Aoi 

 colnames(pk) <- lbudget
 # budget proportions

 aux <- matrix(NA,nrow=J,ncol=K)
 for(i in 1:K){
  aux[,i] <- Boi[,i] * pk[i]
 }

 Bres <- aux/pplusj

}
