summary.bayes=function(object,coverage=.9,...)
{
x = as.numeric(names(object$prob))
p = object$prob
post.mean=sum(x*p)
post.sd=sqrt(sum((x-post.mean)^2*p))
names(p)=NULL
n = length(x)
sp = sort(p, index.return = TRUE)
ps = sp$x
i = sp$ix[seq(n, 1, -1)]
ps = p[i]
xs = x[i]
cp = cumsum(ps)
ii = 1:n
j = ii[cp >= coverage]
j = j[1]
eprob = cp[j]
set = sort(xs[1:j])
v = list(mean=post.mean,sd=post.sd,coverage = eprob, set = set)
return(v)
}
