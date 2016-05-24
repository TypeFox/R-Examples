get.neighbor <- function(map, regions)
{
   if(! inherits(map, "gra"))
      stop("Argument 'map' is not an object of class 'gra'!")

   districts <- rownames(map)

   neighbors <- list(length=length(regions))
   for(i in 1:length(regions)){
      ind <- which(districts == regions[i])  
      neighbors[[i]] <- districts[map[ind,] == -1]
   }
   
   names(neighbors) <- regions
   return(neighbors)
}
