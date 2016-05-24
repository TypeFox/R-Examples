dist.variables <-
function(data, associationFun=association, check.psd=TRUE){
# data: data.frame of original data 
# associationFun: function that calculates association measure for each pair of variables
# check.psd: check if resulting similarity matrix S is positive semi-definite?

  S <- similarity.variables(data, associationFun=associationFun, check.psd=check.psd)
  D <- as.dist(sqrt(1 - S))
  return(D)
}
