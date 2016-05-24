Propensity.Score.strata <-
function(alpha, beta,J,a,b,p1,phi){
p2=phi*p1/(phi*p1+1-p1)
delta=sum(a*b*(1-b)*(p1-p2))
sigma1=sum(a*b*(1-b)*((1-b)*p1*(1-p1)+b*p2*(1-p2)))
sigma0=sum(a*b*(1-b)*(b*p1+(1-b)*p2)*(b*(1-p1)+(1-b)*(1-p2)))
n=(sqrt(sigma0)*qnorm(1-alpha/2)+sqrt(sigma1)*qnorm(1-beta))^2/delta^2
}
