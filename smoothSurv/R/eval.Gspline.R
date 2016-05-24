###########################################
#### AUTHOR:    Arnost Komarek         ####
####            23/07/2004             ####
####                                   ####
#### FILE:      eval.Gspline.R         ####
####                                   ####
#### FUNCTIONS: eval.Gspline.R         ####
###########################################
eval.Gspline <- function(Gspline, grid){
  if (!is.data.frame(Gspline)) stop("Gspline must be a data.frame")

  mus <- as.numeric(Gspline[["Knot"]])
  sig <- as.numeric(Gspline[["SD basis"]])
  ccoef <- as.numeric(Gspline[["c coef."]])

  if (is.null(mus)) stop("Vector of knots did not find in Gspline data.frame")
  if (is.null(sig)) stop("Vector of standard deviations did not find in Gspline data.frame")
  if (is.null(ccoef)) stop("Vector of weights did not find in Gspline data.frame")
  if (sum(is.na(mus))) stop("Incorrect knot vector supplied")
  if (sum(is.na(sig))) stop("Incorrect standard deviations vector supplied")
  if (sum(is.na(ccoef))) stop("Incorrect weights vector supplied")

  ## Function to evaluate a G-spline in one grid point
  dfitted <- function(u){
     normals <- dnorm(u, mean = mus, sd = sig)
     value <- t(ccoef) %*% normals
     return(value)
  }

  ## Compute it in a grid of values
  rooster <- matrix(grid, ncol = 1)  
  y.fitted <- apply(rooster, 1, "dfitted")  
  
  to.return <- data.frame(x = rooster, y = y.fitted)
  return(to.return)    
}  
