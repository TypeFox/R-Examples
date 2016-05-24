ISV.Equality <-
function(alpha, beta, sigma1, sigma2,m){
ratio=sigma1/sigma2
#ratio
for (i in 1:1000){
ratio.f=qf(p=(1-beta), i*(m-1),i*(m-1),lower.tail=FALSE)/qf(alpha/2,i*(m-1),i*(m-1),lower.tail=FALSE)
}
}
