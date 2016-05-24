CrossOver.ISV.NIS <-
function(alpha, beta, sigma1, sigma2,m,margin){
ratio=sigma1/(margin^2*sigma2)
print(c("ratio",ratio))
for (i in 2:100){

ratio.f=qf(p=(1-beta), (2*i-2)*(m-1),(2*i-2)*(m-1),lower.tail=FALSE)/qf(alpha,(2*i-2)*(m-1),(2*i-2)*(m-1),lower.tail=FALSE)
print(c(i,ratio.f))
}
}
