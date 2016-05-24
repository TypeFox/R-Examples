# Function to obtain the matrix of dyadic dominance indices -Dij- #

getDij <- function(X,names=NULL){
  if (nrow(X) != ncol(X)) 
    return("Error: Sociomatrix must be square");
  if ( is.na(X) || !is.numeric(X))
    return("Error: Sociomatrix must be numeric");

dyadc <- X + t(X);
Dij <- X/dyadc-(((X/dyadc)-0.5)/(dyadc+1))
Dij[is.nan(Dij)] <- 0.
if (is.null(names)) names <- paste('Ind.',1:nrow(X))
rownames(Dij) <- names
colnames(Dij) <- names
return(Dij)
}
