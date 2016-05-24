pig.dist.glm <-
function(variable, lambda.start=1, delta.start=1, epsylon = 10^(-8), n=100){
lambda.old = c(); delta.old = c()
lambda.new = c(); delta.new = c()
g=glm(variable~1, family=poisson(log))
lambda.start=exp(g$coefficients[1])
lambda.old = lambda.start
delta.old = delta.start
diff = 0.1
n.iter=0
rules = laguerre.quadrature.rules(n)
rule = rules[[n]]

jpmf = function(variable, lambda, delta) {
mix = function(theta) {
y=(delta*exp(delta^2)*theta^(-3/2)*exp(-((delta^2)/2)*(1/theta + theta)))/(sqrt(2*pi))
for (i in 1:length(variable)) {
y=y*((lambda*theta)^(variable[i])*exp(-lambda*theta)/factorial(variable[i]))
return(y)
}
}
jpmf= laguerre.quadrature(mix, rule, lower = 0, upper = Inf, weighted = FALSE)
return(jpmf)
}
jpmf.loop = c()
jpmf.loop = jpmf(variable = variable, lambda = lambda.old, delta = delta.old)

while(diff>epsylon) { 
##########   E-step 
nominator = c(); denominator = c()
invGauss.density = function(theta) {
y=(delta.old*exp(delta.old^2)*theta^(-3/2)*exp(-((delta.old^2)/2)*(1/theta + theta)))/(sqrt(2*pi))
return(y)
}
for (i in 1:length(variable)) {
t.r.nominator = function(theta) {
y=(lambda.old*theta)^(variable[i]+1)*exp(-lambda.old*theta)*invGauss.density(theta)
return(y)
}
nominator[i] = laguerre.quadrature(t.r.nominator, rule, lower = 0, upper = Inf, weighted = FALSE)
t.r.denominator = function(theta) {
y=(lambda.old*theta)^(variable[i])*exp(-lambda.old*theta)*invGauss.density(theta)
return(y)
}
denominator[i] = laguerre.quadrature(t.r.denominator, rule, lower = 0, upper = Inf, weighted = FALSE)
}
t = nominator/denominator
##########   M-step
g=glm(variable~log(t),  family=poisson(log))
lambda.new=exp(g$coefficients[1])
loglikelihood.delta = function(delta) {
y=log(prod((delta*exp(delta^2)*t^(-3/2)*exp(-((delta^2)/2)*(1/t + t)))/(sqrt(2*pi))))
return(y)
}
delta.new = optimize(f = loglikelihood.delta, interval = c(0,3), maximum=TRUE)$maximum 
jpmf.loop = c(jpmf.loop, jpmf(variable = variable, lambda = lambda.new, delta = delta.new))
diff = abs(jpmf.loop[length(jpmf.loop)]-jpmf.loop[length(jpmf.loop)-1])
delta.old = delta.new 
lambda.old = lambda.new
n.iter = n.iter + 1
}
outlist = list(lambda=lambda.new, delta=delta.new, n.iter=n.iter, likelihood.values=jpmf.loop)
    class(outlist) = "pig.dist.glm"
    return(outlist)
}
