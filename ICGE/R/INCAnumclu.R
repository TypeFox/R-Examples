INCAnumclu <- function(d, K, method="pam",pert, L= NULL, noise=NULL){
############### Number of atypical individuals  #############################
# Input:
# d: distance matrix (n x n)
# K: partitions in 2:K clusters are considered
# method: method to perform the clustering. If method="partition" a collection
#         of M partitions must by given in pert.
# L: if NULL, there are no noise units. If L>=1, the units classified in clusters 
#     such that #C<=L are considered as noise units. If L="custom", the noise units
#     are custom-selected by the user and they must be indicated in argument "noise"
# noise: logic vector that indicates whether the corresponding individual should be 
#        considered as noise unit. Necessarily, argument L must be set as L="custom".
# Output:
# Totals : INCA index for k=2, ..., K
#############################################################################

d <- as.matrix(d)
n <- dim(d)[1]



METHODS <- c("pam","diana","fanny", "average", "single", "complete", "ward", "weighted","partition")
meth <- pmatch(method, METHODS)
if (is.na(meth)) 
    stop("invalid clustering method")
if (meth == -1) 
    stop("ambiguous clustering method")
                                                                
method <- METHODS[meth]



if (is.null(L)){                                                               
   aux <- INCAncinner(d, K, method, pert)
   out <- list(INCAindex= aux$INCAindex, method=method, noise=noise)
}else{
   if (is.character(L)){
      if (is.na(pmatch(L, "custom"))) 
          stop("invalid character name for argument L")
      if (is.null(noise))
          stop("with L=custom, noise argument must be specified")
      notnoise <- !noise
   }else{
       if(method=="partition")
          stop("with method=custom, noise units can not be established")
      
       p <- selectioncluster(d, K, method)
       f <- tabulate(p)
       quitarcluster <- (1:K)[f<=L]
       noise <- p %in% quitarcluster
       notnoise <- !noise
   }



  
  dp <- d[notnoise, notnoise]
  if(method=="partition"){
     pert <- as.matrix(pert)
     pertp <- pert[notnoise,]
  }
  aux1 <- INCAncinner(dp, K, method, pertp)
  aux2 <- INCAncinner(d, K, method, pert)
  out <- list(INCAindex=aux1$INCAindex, method=method, noise=list(indexnoise=aux2$INCAindex, 
              noiseunits=row.names(d)[noise]))
}                                                            
#out <- list(INCAindex=T, method=method)

class(out) <- "incanc"
return(out)

} #end of function
