trimmOutl <- function(resMethod, nsizes){
 UseMethod("trimmOutl")
}

trimmOutl.default <- function(resMethod, nsizes){
 discarded <- resMethod$discarded
 class(discarded) <- "trimmOutl"
 return(discarded)
}

trimmOutl.trimowa <- function(resMethod, nsizes){
 discarded <- list()
 for (i in 1 : nsizes){
  discarded[[i]] <- resMethod[[i]]$discarded
 } 
    
 class(discarded) <- "trimmOutl"
 return(discarded)
}

trimmOutl.hipamAnthropom <- function(resMethod, nsizes){
 discarded <- list()
 for (i in 1 : nsizes){ 
  aux <- table(resMethod[[i]]$clustering)
  aux <- as.numeric(aux)
  auxNoBig <- which(aux == 1 | aux == 2)
  discarded[[i]] <- rownames(unique(resMethod[[i]]$cases))[auxNoBig]
 }
  
 class(discarded) <- "trimmOutl"
 return(discarded)
}  

