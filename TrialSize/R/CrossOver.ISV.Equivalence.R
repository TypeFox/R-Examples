CrossOver.ISV.Equivalence <-
function(alpha, beta, sigma1, sigma2,m,margin){
ratio=(margin^2*sigma1)/sigma2
print(c("ratio",ratio))
for (i in 2:100){

ratio.f=qf(p=(beta/2), (2*i-2)*(m-1),(2*i-2)*(m-1),lower.tail=FALSE)/qf((1-alpha),(2*i-2)*(m-1),(2*i-2)*(m-1),lower.tail=FALSE)
print(c(i,ratio.f))
}
}
