#' @title Efficiency bounds 
#' 
#' @description
#' Finds upper A-efficiency bounds for regular block designs.
#' 
#' @details
#' Upper bounds for the A-efficiency factors of regular block designs 
#' (see Chapter 2.8 of John and Williams 1995). Non-trivial A-efficiency upper bounds
#' are calculated for regular block designs with equal block sizes 
#' and equal replication. All other designs return NA.    
#' 
#' @param n the total number of plots in the design.
#' 
#' @param v the total number of treatments in the design.
#' 
#' @param b the total number of blocks in the design.
#' 
#' @references
#' 
#' John, J. A. and Williams, E. R. (1995). Cyclic and Computer Generated Designs. Chapman and Hall, London.
#' 
#' @examples 
#' 
#' # 50 plots, 10 treatments and 10 blocks for a design with 5 replicates and blocks of size 5 
#' upper_bounds(n=50,v=10,b=10)
#'
#' @export
  upper_bounds=function(n,v,b) {
  if ( !isTRUE(all.equal(n%%v,0)) ||  !isTRUE(all.equal(n%%b,0)) || (v+b-1)>n ) return(NA) 
  r = n/v #replication
  if (isTRUE(all.equal(r%%b, 0))) return(1) 
  k = n/b #block size	
  # this bound is for non-binary designs - see John and Williams page 44
  if (k > v) {
    kp=k%%v
    U0=1 - kp*(v-kp)/(k*k*(v - 1)) 
    rp=b*kp/v
    phi=kp/k
    lambdap = rp*(kp - 1)/(v - 1)
    alphap = lambdap - floor(lambdap)
    s2=phi^4*v*(v-1)*alphap*(1-alphap)/((rp*kp)**2) # corrected second moment lower bound
    S=sqrt(s2/(v-1)/(v-2))
    U1 = U0 - (v - 2)*S*S/(U0 + (v - 3)*S)
    U2 = U0 - (1 - U0)*s2/((1 - U0)*(v - 1) - s2)
    bound=min(U0,U1,U2,na.rm = TRUE)	
  } else {
  # this bound is for binary designs 
  dual=v>b 
  if (dual) {
    temp = b
    b = v
    v = temp
    temp=r
    r = k
    k = temp
  }	
  bound =  v*(k - 1)/(k*(v - 1))
  lambda = r*(k - 1)/(v - 1)
  if  (!identical(lambda,floor(lambda))) {
    ebar=bound 
    alpha = lambda - floor(lambda)
    s2=v*(v-1)*alpha*(1-alpha)/(r*k*r*k) # corrected second moment lower bound
    S=sqrt(s2/(v-1)/(v-2))
    U1= ebar - (v - 2)*S*S/(ebar + (v - 3)*S)
    U2= ebar - (1 - ebar)*s2/((1 - ebar)*(v - 1) - s2)
    if ( alpha < v/(2*(v - 1)) ) 
      z = alpha*((v + 1)*alpha - 3) else
      z = (1 - alpha)*(v - (v + 1)*alpha)
    s3J = alpha*v*(v - 1)*z/((r*k)**3)
    U4= ebar - s2*s2/((v - 1)*(s3J+ ebar*s2))	
    if  (identical(0,floor(lambda)))
      s3P = alpha*v*(v - 1)*(  (v + 1)*alpha*alpha - 3*alpha - k + 2)/((r*k)**3) else s3P=s3J
    U5= ebar - s2*s2/((v - 1)*(s3P+ ebar*s2)) 
    bound=min(U1,U2,U4,U5,na.rm = TRUE)	
  }
  if (dual) 
    bound = (b - 1)/((b - v) + (v - 1)/bound)
  }
  return(round(bound,6))
  }	
