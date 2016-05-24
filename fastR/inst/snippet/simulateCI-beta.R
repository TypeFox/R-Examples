mu <- 0.4 / (0.4 + 0.6); mu                          # mean for beta dist
v <- (0.4*0.6) / ((0.4 + 0.6)^2 * (0.4+0.6+1)); v    # var for beta dist
#
# 10,000 simulated samples of size 20
CIsim(n=20, samples=10000, rdist=rbeta, args=list(shape1=0.4,shape2=0.6),
	estimand=mu, method=ci, method.args=list(sd=sqrt(v)))
#
# 10,000 simulated samples of size 5
CIsim(n=5, samples=10000, rdist=rbeta, args=list(shape1=0.4,shape2=0.6),
	estimand=mu, method=ci, method.args=list(sd=sqrt(v)))
#
# 10,000 simulated samples of size 2
CIsim(n=2, samples=10000, rdist=rbeta, args=list(shape1=0.4,shape2=0.6),
	estimand=mu, method=ci, method.args=list(sd=sqrt(v)))
