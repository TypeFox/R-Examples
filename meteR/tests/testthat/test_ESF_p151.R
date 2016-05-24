context('ESF_p151')

s0=c(rep(4,8),rep(16,8),rep(64,8),rep(256,8))
n0.tmp=rep(c(2^2,2^4,2^8,2^12),8)
n0=s0*n0.tmp
e0.tmp=rep(c(rep(2^2,4),rep(2^10,4)),4)
e0=n0*e0.tmp

