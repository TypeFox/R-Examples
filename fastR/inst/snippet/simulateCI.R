# an example CI from a sample of size 20
ci(rnorm(20,500,100))$conf.int                     
# 10,000 simulated samples of size 20
CIsim(n=20, samples=10000, rdist=rnorm, args=list(mean=500,sd=100),
	estimand=500, method=ci, method.args=list(sd=100))
#
# an example CI from a sample of size 5
ci(rnorm(5,500,100))$conf.int
# 10,000 simulated samples of size 5
CIsim(n=5,samples=10000, rdist=rnorm, args=list(mean=500,sd=100),
	estimand=500, method=ci, method.args=list(sd=100))
#
# an example CI from a sample of size 2
ci(rnorm(2,500,100))$conf.int
# 10,000 simulated samples of size 2
CIsim(n=2,samples=10000, rdist=rnorm, args=list(mean=500,sd=100),
	estimand=500, method=ci, method.args=list(sd=100))
