ANOVA.Repeat.Measure <-
function(alpha, beta,sigma,delta,m){
n=2*sigma^2*(qnorm(1-alpha/(m*2))+qnorm(1-beta))^2/delta^2
}
