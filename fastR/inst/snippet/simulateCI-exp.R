rate = 1/10
v = (1/rate)^2                                     # var of exponential
mu = 10                                            # mean of exponential
ci(rexp(20,rate),sd=sqrt(v))$conf.int              # an example CI
#
# 10,000 simulated samples of size 20
CIsim(n=20, samples=10000, rdist=rexp, args=list(rate=rate),
	estimand=mu, method=ci, method.args=list(sd=sqrt(v)))
#
# 10,000 simulated samples of size 5
CIsim(n=5, samples=10000, rdist=rexp, args=list(rate=rate),
	estimand=mu, method=ci, method.args=list(sd=sqrt(v)))
#
# 10,000 simulated samples of size 2
CIsim(n=2, samples=10000, rdist=rexp, args=list(rate=rate),
	estimand=mu, method=ci, method.args=list(sd=sqrt(v)))
