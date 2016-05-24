ALbw <-
function(type_kernel="n", vec_data)
#######################################################################
#               REQUIRED INPUTS FOR THE FUNCTION        		          #
#######################################################################
#   "type_kernel" kernel function:  "e" Epanechnikov,	"n" Normal, 
#   "vec_data" sample data to estimate distribution function
{
  n <- length(vec_data)
  orderr <- 0
  # pilot badwidth
  bp <- (n^(-0.3)) * sd(vec_data)
 
  # weight function
  ss <- quantile(vec_data, c(orderr, 1-orderr))
  w <- (ss[1]<=vec_data) & (vec_data<=ss[2])

  # plug-in estimator of the asymptotically optimal bandwidth
  aux <- outer(vec_data, vec_data, "-")/bp
  aux <- kernel_function(type_kernel, aux)
  diag(aux) <- 0
  D2_F <- sum(aux*w)/(bp*n*(n-1))
  V2 <- 2 * A1_k(type_kernel) * D2_F
  aux1 <- outer(vec_data, vec_data, "-")/bp
  aux1 <- derivative_kernel_function(type_kernel, aux1)
  aux2 <- apply(t(aux1)%*%(aux1*w), 1, sum)
  D3_F <- sum(aux2)/((n^3)*(bp^4))
  B3 <- 0.25 * (A2_k(type_kernel))^2 * D3_F

  ALbw <- (((0.25*V2)/B3)^(1/3)) * (n^(-1/3))
  
  return(ALbw)

}
