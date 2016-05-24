numEff= 
function(x,m=as.integer(min(length(x),(100/sqrt(5000))*sqrt(length(x)))))
{
#
# P. Rossi
# revision history: 3/27/05
#
# purpose:
#  compute N-W std error and relative numerical efficiency
# 
# Arguments:
#  x is vector of draws
#  m is number of lags to truncate acf
#    def is such that m=100 if length(x)= 5000 and grows with sqrt(length)
#
# Output:
#  list with numerical std error and variance multiple (f)
#
wgt=as.vector(seq(m,1,-1))/(m+1)
z=acf(x,lag.max=m,plot=FALSE)
f=1+2*wgt%*%as.vector(z$acf[-1])
stderr=sqrt(var(x)*f/length(x))
list(stderr=stderr,f=f,m=m)
}

