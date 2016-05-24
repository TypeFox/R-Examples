
 #####################
 #### h.knn.funct ####
 #####################

 ## Description: Internal function.  setting a minimum bandwidth in order to
 ##              check that a candidate bandwidth h doesn't contains DB linear 
 ##              models with only one observation. Takes the distance that 
 ##              includes the 3 nearest neighbors for each individual row. 
 ## 
 ##        Inputs:  dist, k=3.
 ##        Outputs: h.knn.k
 ##

h.knn.funct <- function(dist,k=3){
   d.order <- apply(dist,1,order)
   h.knn.k <- sapply(1:dim(dist)[1],function(i){dist[i,d.order[k,i]]})
}