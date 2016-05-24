selectioncluster <- function(d, K, method){
# To check whether there is some cluster with few units
pama <- function(d, K){
  p <- pam(d, K, diss=TRUE)$clustering
  return(p)
}  

dian <- function(d, K){
   aux <- diana(d, diss=TRUE)
   p <- cutree(aux, K)
   return(p)
}  


fannya <- function(d, K){
   p <- fanny(d, K, diss=TRUE)$clustering
   return(p)
} 


bestela <- function(d, K, method){
   aux <- agnes(d, diss=TRUE, method=method)
   p <- cutree(aux, K)
   return(p)
}  

 p <- switch(method,  pam=pama(d, K), diana=dian(d, K), fanny=fannya(d, K), bestela(d,K, method))

return(p)
}
