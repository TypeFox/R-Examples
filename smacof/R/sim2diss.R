# converts similarity matrix into dissimilarities

sim2diss <- function(similmat, method = "corr", to.dist = TRUE)
{
  # similmat... similarity matrix (not necessarily symmetric, nor quadratic)
  # method ... various methods allowed: "corr" for correlation matrices, "neglog", "counts" or an integer value
  # to.dist ... if TRUE, it creates an object of class "dist", if FALSE a matrix.

  if (method == "corr") dissmat <- sqrt(1-similmat)
  if (method == "neglog") dissmat <- -log(similmat)
  if (method == "counts") dissmat <- -log((similmat*t(similmat))/(outer(diag(similmat),diag(similmat))))
  if (is.numeric(method)) dissmat <- method - similmat

  if (to.dist) dissmat <- as.dist(dissmat)

  return(dissmat)
}
  
                       
                    
 
