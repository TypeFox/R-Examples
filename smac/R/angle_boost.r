angle_boost=function(x,y,weight,nlambda,lambda.min,lambda,standardize,epsilon)
{

if (is.null(lambda)) {z = boostway1(x, y, weight, nlambda, lambda.min, standardize, epsilon)}
if (!is.null(lambda)) {z = boostway2(x, y, weight, lambda, standardize, epsilon)}

return(z)
}
