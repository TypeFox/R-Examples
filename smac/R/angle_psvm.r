angle_psvm=function(x,y,weight,nlambda,lambda.min,lambda,standardize,epsilon)
{

if (is.null(lambda)) {z = psvmway1(x, y, weight, nlambda, lambda.min, standardize, epsilon)}
if (!is.null(lambda)) {z = psvmway2(x, y, weight, lambda, standardize, epsilon)}

return(z)
}
