ci.rc <- function(b.k, SE.b.k=NULL, s.Y=NULL, s.X=NULL, N, K, R2.Y_X=NULL, R2.k_X.without.k=NULL, conf.level=.95, R2.Y_X.without.k=NULL, t.value=NULL, alpha.lower=NULL, alpha.upper=NULL, Noncentral=FALSE, Suppress.Statement=FALSE, ...)
{
result<-ci.reg.coef(b.j=b.k, SE.b.j=SE.b.k, s.Y=s.Y, s.X=s.X, N=N, p=K, R2.Y_X=R2.Y_X, R2.j_X.without.j=R2.k_X.without.k, conf.level=conf.level, R2.Y_X.without.j=R2.Y_X.without.k, t.value=t.value, alpha.lower=alpha.lower, alpha.upper=alpha.upper, Noncentral=Noncentral, Suppress.Statement=TRUE, ...)

if(Noncentral==FALSE)
{if(Suppress.Statement!=TRUE)
print (paste(conf.level*100, "percent CI limits (with corresponding probability) for the kth population regression coefficient calculated using the (central) t-distribution  with", N-K-1, "degrees of freedom follow"))
}

if(Noncentral==TRUE)
{if(Suppress.Statement!=TRUE)
print (paste(conf.level*100, "percent CI limits (with corresponding probability) for the kth population regression coefficient calculated using the (noncentral) t-distribution  with", N-K-1, "degrees of freedom follow"))
}

result.new <- list(Lower.Limit.for.b.k=result$Lower.Limit.for.beta.j, Prob.Less.Lower=result$Prob.Less.Lower,
Upper.Limit.for.b.k=result$Upper.Limit.for.beta.j, Prob.Greater.Upper=result$Prob.Greater.Upper) 

return(result.new)
}
