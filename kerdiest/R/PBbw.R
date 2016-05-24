PBbw <-
function(type_kernel="n", vec_data, num_stage=2)
#######################################################################
#               REQUIRED INPUTS FOR THE MAIN FUNCTION		              #
#######################################################################
#   "type_kernel" kernel function to estimate F:  "e" Epanechnikov,	"n" Normal, 
#                                              "b" Biweight, "t" Triweight    
#   "vec_data" data sample 
#   "num_stage" number of stages for estimating the optimal bandwidth. Usually num_stage=2 but it is possible to choose 2, 3 or 4 
{
  if(num_stage > 4 | num_stage < 2)
    {
      stop("The parameter num_stage must be 2, 3 or 4")
    }
  n <- length(vec_data)

  # dispersion parameter min{s,IQR/1.349}
  sigma <- min(sd(vec_data), IQR(vec_data)/1.349)

  # b-stage estimator of ho
  # Step 1
  psi_NR <- function(r)
  {
    result <- ((-1)^(r/2) * factorial(r))/((2*sigma)^(r + 1) * (factorial(r/2)) * (pi^(1/2)))
    return(result)
  }
  
  # Step 2
  # function Psi_r(g)
  psi <- function(r, g)
  {
    aux <- outer(vec_data,vec_data,"-")/g
    aux <- derivative_normal_kernel(r, aux)
    result <- sum(aux)*((g^(-r-1))*(n^(-2)))
    return(result)
  }

  # function g_r
  g_r <- function(r, psi_r2)
  {
		num <- 2*derivative_normal_kernel(r, 0)
		denom <- (-n) * mu2_k("n") * psi_r2
		result <- (num/denom)^(1/(r+3))
		return(result)
  }

  # plug-in estimator of the asymptotically optimal bandwidth
  psi_aux <- psi_NR((2*num_stage)+2)
  g_aux <- g_r(2*num_stage, psi_aux)

  for(i in 1:(num_stage-1))
  {
    psi_aux <- psi((2*(num_stage-i))+2, g_aux)
    g_aux <- g_r(2*(num_stage-i) , psi_aux) 
  }
  band_pilot <- g_aux
  psi2 <- psi(2, band_pilot) 
  ro <- ro_k(type_kernel)
  PBbw <- (ro/((-n)*((mu2_k(type_kernel))^2)*(psi2)))^(1/3)

  return(PBbw)

}
