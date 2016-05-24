#### derivative of the truncation region function as in Gross et al's write up
#### input:
####    - x = Point at which to evaluate derivative
####    - q, r = Coefficients of the function
gradient.truncation.region.function <-
function(x, q, r){
  q/2/sqrt(x) + r/2/sqrt(1+x)
}
