compute_alphas <- function(Z, gs=gf, dgs=dgf, name=gnames)
{
   ng <- length(gs)
   p <- ncol(Z)
   alphas <- matrix(0,ng,p)
   
   if(length(name)!=length(gs)){
    name <- NULL
    for(i in 1:ng){ 
     name[i] <- paste("g",i)
    }
   }

   rownames(alphas) <- name
 
   if(is.null(colnames(Z))){
    cnam <- NULL
    for(j in 1:p){
     cnam[j] <- paste("IC",j)
    } 
    colnames(alphas) <- cnam  
   }else colnames(alphas) <- colnames(Z)
 
   for(j in 1:p){
    for(i in 1:ng){
     Eg <- mean(gs[[i]](Z[,j]))
     Eg2 <- mean(gs[[i]](Z[,j])^2)
     Egx <- mean(gs[[i]](Z[,j])*Z[,j])
     Edg <- mean(dgs[[i]](Z[,j]))
    if((Egx-Edg)==0){
     alphas[i] <- Inf
    }else alphas[i,j] <- (Eg2-Eg^2-Egx^2)/(Egx-Edg)^2
   } 
   }
   
   alphas
}


