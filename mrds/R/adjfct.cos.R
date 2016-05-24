# Cosine adjustment terms
#
# distance  - perpendicular distance vector
# scaling     - scale parameter
# adj.order - vector of orders of Cosine terms to fit
# adj.parm  - vector of parameters (a_j)
adjfct.cos <- function(distance,scaling=1,adj.order,adj.parm=NULL,
                       adj.exp=FALSE){

  # Check the adjustment parameters
  if(is.null(adj.parm)){
    adj.parm <- as.vector(rep(1,length(adj.order)))
  }

  adj.order <- as.vector(adj.order)

  cossum <- 0

  for(i in 1:length(adj.order)){
    cossum <- cossum + (adj.parm[i]*cos((adj.order[i]*pi*distance)/scaling))
  }

  # if adj.exp return exp(cossum) to keep f(x)>0  
  if(adj.exp){
    return(exp(cossum))
  }else{
    return(cossum)
  }
}
