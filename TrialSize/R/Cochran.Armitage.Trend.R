Cochran.Armitage.Trend <-
function(alpha, beta,pi,di,ni,delta){
xi=ni*pi
N=sum(ni)
qi=1-pi
p=sum(ni*pi)/N
q=1-p
ri=ni/ni[1]
dmean=sum(ni*di)/sum(ni)
A=sum(ri*pi*(di-dmean))
n0_star=(qnorm(1-alpha)*sqrt(p*q*sum(ri*(di-dmean)^2))+qnorm(1-beta)*sqrt(sum(pi*qi*ri*(di-dmean)^2)))^2/A^2
n=n0_star*(1+sqrt(1+2*delta/(A*n0_star)))^2/4
return(n)
}
