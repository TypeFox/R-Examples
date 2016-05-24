# @title Sequence of combinatorial coefficients
# @description Compute the sequence of binomial coefficients \eqn{{m-1}\choose{r-1}}, for \eqn{r= 1, \dots, m}, 
# and then returns a vector of the same length as ordinal, whose i-th component is the corresponding binomial 
# coefficient \eqn{{m-1}\choose{r_i-1}}
# @aliases kkk
# @usage kkk(m, ordinal)
# @param m Number of ordinal categories
# @param Vector of ordinal responses
#' @keywords internal


kkk <-
function(m,ordinal){
  serie<-1:m
  vett<-choose(m-1,serie-1)
  return(vett[ordinal])
}
