# adjfct.poly
#
# Polynomial adjustment terms
#
# distance  - perpendicular distance vector
# scaling     - scale parameter
# adj.order - vector of orders of polynomials to fit
# adj.parm  - vector of parameters (a_j)
adjfct.poly <- function(distance,scaling=1,adj.order,adj.parm=NULL,
                        adj.exp=FALSE){

  # Check the adjustment parameters
  if(is.null(adj.parm)){
    adj.parm <- as.vector(rep(1,length(adj.order)))
  }

  adj.order <- as.vector(adj.order)

  # Should have checked the order beforehand
  polysum <- 0

  for(i in 1:length(adj.order)){
    polysum <- polysum + (adj.parm[i]*(distance/scaling)^adj.order[i])
  }

  # if adj.exp return exp(polysum) to keep f(x)>0 
  if(adj.exp){
    return(exp(polysum))
  }else{
    return(polysum)
  }
}
