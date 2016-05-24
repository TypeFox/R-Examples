mu = 1/2; v = 1/12           # mean and variance
#
# 10,000 simulated samples of size 20
CIsim(n=20, samples=10000, rdist=runif, estimand=mu,
	method=ci, method.args=list(sd=sqrt(v)))
#
# 10,000 simulated samples of size 5
CIsim(n=5, samples=10000, rdist=runif, estimand=mu,
	method=ci, method.args=list(sd=sqrt(v)))
#
# 10,000 simulated samples of size 2
CIsim(n=2, samples=10000, rdist=runif, estimand=mu,
	method=ci, method.args=list(sd=sqrt(v)))
