kPCA.kern <-
function(x, ...){
  objectPC <- getPrincipalComponents(x)
  object <- list(KPCs = objectPC$KPCs, Es = objectPC$Es, Vecs = objectPC$Vecs, K = x, theta = NULL, x  = NULL)
  class(object) <- c("kPCA.kern", "kPCA")   
  return(object)
}
