add_size.mean.two_stage <- function(L, alpha, presample)
{
n0=length(presample)
s0=sd(presample)
z=(L/2)^2/(qt(1-alpha/2,n0-1)^2)
n=max(n0+1,ceiling(s0^2/z))
list("Significance level"=alpha,
"Length of confidence interval"=L,
"Length of presample"=n0,
"Standard deviation of presample"=s0,
"Number of additional observations"=n-n0)
}

#add_size.mean.two_stage(0.9,0.05,rnorm(16,10,1.2))

