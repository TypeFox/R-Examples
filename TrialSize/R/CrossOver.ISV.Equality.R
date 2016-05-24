CrossOver.ISV.Equality <-
function(alpha, beta, sigma1, sigma2,m){
ratio=sigma1/sigma2
ratio
for (i in 2:1000){

ratio.f=qf(p=(1-beta), (2*i-2)*(m-1),(2*i-2)*(m-1),lower.tail=FALSE)/qf(p=alpha/2,(2*i-2)*(m-1),(2*i-2)*(m-1),lower.tail=FALSE)
print(c(i,ratio.f))
}
}
