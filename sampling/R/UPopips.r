UPopips<-function(lambda, type=c("pareto","uniform","exponential"))
{
if(any(is.na(lambda))) stop("there are missing values in the lambda vector")
n=sum(lambda)
if(!(type %in% c("pareto","uniform","exponential"))) stop("the type argument is wrong")
omega=runif(n)
switch(type,pareto=order(omega*(1-lambda)/((1-omega)*lambda))[1:n],
uniform=order(omega/lambda)[1:n], exponential=order(log(1-omega)/log(1-lambda)))[1:n]
}




