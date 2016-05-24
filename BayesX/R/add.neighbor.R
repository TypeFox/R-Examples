add.neighbor <- function(map, region1, region2)
{
   if(! inherits(map,"gra"))
      stop("Argument 'map' is not an object of class 'gra'!")

   # names of districts
   districts <- rownames(map)

   # find index of the regions with changes
   ind1 <- which(districts == region1)
   ind2 <- which(districts == region2)
  
   # modify neighborhoopd structure
   map[ind1,ind2] <- map[ind2,ind1] <- -1
  
   # and adjust no. of neighbors
   map[ind1,ind1] <- -sum(map[ind1, -ind1])
   map[ind2,ind2] <- -sum(map[ind2, -ind2])

   class(map) <- "gra"
   return(map)
}
