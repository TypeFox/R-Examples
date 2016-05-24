# Function to obtain the matrix of win proportions -Pij- #

getPij <- function(X,names=NULL){
  if (nrow(X) != ncol(X)) 
    return("Error: Sociomatrix must be square");
  if ( is.na(X) || !is.numeric(X))
    return("Error: Sociomatrix must be numeric");
dyadc <- X + t(X);
Pij <- X/dyadc
Pij[is.nan(Pij)] <- 0.
if (is.null(names)) names <- paste('Ind.',1:nrow(X))
rownames(Pij) <- names
colnames(Pij) <- names
return(Pij)
}