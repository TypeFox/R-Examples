ci.src <- function(beta.k=NULL, SE.beta.k=NULL, N=NULL, K=NULL, R2.Y_X=NULL, R2.k_X.without.k=NULL, 
conf.level=.95, R2.Y_X.without.k=NULL, t.value=NULL, b.k=NULL, SE.b.k=NULL, s.Y=NULL, s.X=NULL, 
alpha.lower=NULL, alpha.upper=NULL, Suppress.Statement=FALSE, ...)
{

if(!is.null(b.k))
{
if(is.null(s.Y)) stop("Since you have specified the unstandardized regression coefficient, you must also specify the standard deviation of Y (so that the function can compute the standardized regression coefficient).")
if(is.null(s.X)) stop("Since you have specified the unstandardized regression coefficient, you must also specify the standard deviation of X (so that the function can compute the standardized regression coefficient).")

beta.k <- b.k*(s.X/s.Y)

if(!is.null(SE.b.k)) SE.beta.k <- SE.b.k*(s.X/s.Y)
}


if(beta.k > 1.1) warning("This function is only for standardized regression coefficients. Is your 'b.k' in standarized units (the observed value), although possible, seems quite large?", call.=FALSE)

result<-ci.reg.coef(b.j=beta.k, SE.b.j=SE.beta.k, s.Y=1, s.X=1, N=N, p=K, 
R2.Y_X=R2.Y_X, R2.j_X.without.j=R2.k_X.without.k, conf.level=conf.level, 
R2.Y_X.without.j=R2.Y_X.without.k, t.value=t.value, alpha.lower=alpha.lower, 
alpha.upper=alpha.upper, Noncentral=TRUE, Suppress.Statement=TRUE, ...)

if(Suppress.Statement!=TRUE)
print (paste(conf.level*100, "percent CI limits (with corresponding probability) for the kth population standard regression coefficient calculated using the (noncentral) t-distribution  with", N-K-1, "degrees of freedom follow"))


result.new <- list(Lower.Limit.for.beta.k=result$Lower.Limit.for.beta.j, Prob.Less.Lower=result$Prob.Less.Lower,
Upper.Limit.for.beta.k=result$Upper.Limit.for.beta.j, Prob.Greater.Upper=result$Prob.Greater.Upper) 

return(result.new)

}
