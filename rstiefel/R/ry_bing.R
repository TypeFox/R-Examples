ry_bing <-
function(y,l)
{
  #C function to perform a Gibbs update of 
  #y with invariant distribution 
  #p(y|l) \propto \exp{ sum(l*y^2) } 
  #with respect to the uniform measure on the sphere
  .C("ry_bing",y=as.double(y),l=as.double(l),n=as.integer(length(y)))$y
}
