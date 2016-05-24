mnormt.onesided=function(m0,normpar,data)
{
#
# mnormt.onesided Performs a test that a normal mean is <= certain value.
#    m0 = value to be tested
#    normpar = mean and standard deviation of normal prior on mu
#    data = (sample mean, sample size, known sampling standard deviation)

xbar=data[1]; n=data[2]; s=data[3]
prior.mean=normpar[1]
prior.sd=normpar[2]
prior.var=prior.sd^2

priorH=pnorm(m0,prior.mean,prior.sd)
priorA=1-priorH
prior.odds=priorH/priorA

post.precision=1/prior.var+n/s^2
post.var=1/post.precision
post.sd=sqrt(post.var)
post.mean=(xbar*n/s^2+prior.mean/prior.var)/post.precision
postH=pnorm(m0,post.mean,post.sd)
postA=1-postH
post.odds=postH/postA
BF=post.odds/prior.odds

return(list(BF=BF,prior.odds=prior.odds,post.odds=post.odds,postH=postH))
}
