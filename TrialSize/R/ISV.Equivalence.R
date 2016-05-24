ISV.Equivalence <-
function(alpha, beta, sigma1, sigma2,m,margin){
ratio=margin^2*sigma1/sigma2
#ratio
for (i in 1:1000){

ratio.f=qf(beta/2, i*(m-1),i*(m-1),lower.tail=FALSE)/qf(1-alpha,i*(m-1),i*(m-1),lower.tail=FALSE)

}
}
