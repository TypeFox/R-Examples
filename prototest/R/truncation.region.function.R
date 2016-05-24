#### truncation region function as in Gross et al's write up
#### input:
####    - x = Point at which to evaluate function
####    - q, r, s = Coefficients of the function
truncation.region.function <-
function(x, q, r, s){
  q*sqrt(x) + r*sqrt(1 + x) + s
}
