ISV.NIS <-
function(alpha, beta, sigma1, sigma2,m,margin){
ratio=sigma1/(sigma2*margin^2)
#print(c("ratio",ratio))
for (i in 1:1000){
ratio.f=qf(p=(1-beta), i*(m-1),i*(m-1),lower.tail=FALSE)/qf(alpha,i*(m-1),i*(m-1),lower.tail=FALSE)
print(ratio.f)
}
}
