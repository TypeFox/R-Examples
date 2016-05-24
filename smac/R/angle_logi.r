angle_logi=function(x,y,weight,nlambda,lambda.min,lambda,standardize,epsilon)
{

if (is.null(lambda)) {z = logiway1(x, y, weight, nlambda, lambda.min, standardize, epsilon)}
if (!is.null(lambda)) {z = logiway2(x, y, weight, lambda, standardize, epsilon)}

return(z)
}
