pln.dist.glm <-
function(variable, lambda.start=1, nu.start=1, epsylon = 10^(-8), n=100){
lambda.old = c(); nu.old = c()
lambda.new = c(); nu.new = c()
g=glm(variable~1, family=poisson(log))
lambda.start=exp(g$coefficients[1])
lambda.old = lambda.start
nu.old = nu.start
diff = 0.1
n.iter=0
rules = laguerre.quadrature.rules(n)
rule = rules[[n]]

jpmf = function(variable, lambda, nu) {
mix = function(theta) {
y=exp(-((log(theta)+(nu^2)/2)^2)/(2*nu^2))/(sqrt(2*pi)*nu*theta)
for (i in 1:length(variable)) {
y=y*((lambda*theta)^(variable[i])*exp(-lambda*theta)/factorial(variable[i]))
return(y)
}
}
jpmf= laguerre.quadrature(mix, rule, lower = 0, upper = Inf, weighted = FALSE)
return(jpmf)
}

jpmf.loop = c()
jpmf.loop = jpmf(variable = variable, lambda = lambda.old, nu = nu.old)

while(diff>epsylon) { 
##########   E-step 
nominator = c(); denominator = c()
lognorm.density = function(theta) {
y=exp(-((log(theta)+(nu.old^2)/2)^2)/(2*nu.old^2))/(sqrt(2*pi)*nu.old*theta)
return(y)
}
for (i in 1:length(variable)) {
t.r.nominator = function(theta) {
y=(lambda.old*theta)^(variable[i]+1)*exp(-lambda.old*theta)*lognorm.density(theta)
return(y)
}
nominator[i] = laguerre.quadrature(t.r.nominator, rule, lower = 0, upper = Inf, weighted = FALSE)
t.r.denominator = function(theta) {
y=(lambda.old*theta)^(variable[i])*exp(-lambda.old*theta)*lognorm.density(theta)
return(y)
}
denominator[i] = laguerre.quadrature(t.r.denominator, rule, lower = 0, upper = Inf, weighted = FALSE)
}
t = nominator/denominator
##########   M-step
g=glm(variable~log(t),  family=poisson(log))
lambda.new=exp(g$coefficients[1])
loglikelihood.nu = function(nu) {
y=log(prod(exp(-((log(t)+(nu^2)/2)^2)/(2*nu^2))/(sqrt(2*pi)*nu*t)))
return(y)
}
nu.new = optimize(f = loglikelihood.nu, interval = c(0,10), maximum=TRUE)$maximum
jpmf.loop = c(jpmf.loop, jpmf(variable = variable, lambda = lambda.new, nu = nu.new))
diff = abs(jpmf.loop[length(jpmf.loop)]-jpmf.loop[length(jpmf.loop)-1])
nu.old = nu.new 
lambda.old = lambda.new
n.iter = n.iter + 1
}
outlist = list(lambda=lambda.new, nu=nu.new, n.iter=n.iter, likelihood.values=jpmf.loop)
    class(outlist) = "pln.dist.glm"
    return(outlist)
}
