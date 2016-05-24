# adjfct.herm
#
# Hermite polynomial adjustment terms
#
# distance  - perpendicular distance vector
# scaling     - scale parameter
# adj.order - vector of orders of polynomials to fit
# adj.parm  - vector of parameters (a_j)
adjfct.herm <- function(distance,scaling=1,adj.order,adj.parm=NULL,
                        adj.exp=FALSE){

  # Check the adjustment parameters
  if(is.null(adj.parm)){
    adj.parm <- as.vector(rep(1,length(adj.order)))
  }

  adj.order <- as.vector(adj.order)

  hermsum <- 0

  for(i in 1:length(adj.order)){
    hermsum <- hermsum + (adj.parm[i]*hermite.poly((distance/scaling),adj.order[i]))
  }

  # if adj.exp return exp(hermsum) to keep f(x)>0 
  if(adj.exp){
    return(exp(hermsum))
  }else{
    return(hermsum)
  }
}
