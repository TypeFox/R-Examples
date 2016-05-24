pg.dist.glm <-
function(variable, lambda.start=1, gamma.par.start=1, epsylon = 10^(-8), n=100){
lambda.old = c(); gamma.par.old = c()
lambda.new = c(); gamma.par.new = c()
g=glm(variable~1, family=poisson(log))
lambda.start=exp(g$coefficients[1])
lambda.old = lambda.start
gamma.par.old = gamma.par.start
diff = 0.1
n.iter=0
rules = laguerre.quadrature.rules(n)
rule = rules[[n]]
jpmf = function(variable, lambda, gamma.par) {
mix = function(theta) {
y=(gamma.par^(gamma.par)*theta^(gamma.par-1)*exp(-(gamma.par*theta)))/(gamma(gamma.par))
for (i in 1:length(variable)) {
y=y*((lambda*theta)^(variable[i])*exp(-lambda*theta)/factorial(variable[i]))
return(y)
}
}
jpmf = laguerre.quadrature(mix, rule, lower = 0, upper = Inf, weighted = FALSE)
return(jpmf)
}

jpmf.loop = c()
jpmf.loop = jpmf(variable = variable, lambda = lambda.old, gamma.par = gamma.par.old)

while(diff>epsylon) { 
##########   E-step 
nominator = c(); denominator = c()
Gamma.density = function(theta) {
y=(gamma.par.old^(gamma.par.old)*theta^(gamma.par.old-1)*exp(-(gamma.par.old*theta)))/(gamma(gamma.par.old))
return(y)
}
for (i in 1:length(variable)) {
t.r.nominator = function(theta) {
y=(lambda.old*theta)^(variable[i]+1)*exp(-lambda.old*theta)*Gamma.density(theta)
return(y)
}
nominator[i] = laguerre.quadrature(t.r.nominator, rule, lower = 0, upper = Inf, weighted = FALSE)
t.r.denominator = function(theta) {
y=(lambda.old*theta)^(variable[i])*exp(-lambda.old*theta)*Gamma.density(theta)
return(y)
}
denominator[i] = laguerre.quadrature(t.r.denominator, rule, lower = 0, upper = Inf, weighted = FALSE)
}
t = nominator/denominator
##########   M-step
g=glm(variable~log(t),  family=poisson(log))
lambda.new=exp(g$coefficients[1])
loglikelihood.gamma.par = function(gamma.par) {
y=log(prod((gamma.par^(gamma.par)*t^(gamma.par-1)*exp(-(gamma.par*t)))/(gamma(gamma.par))))
return(y)
}
gamma.par.new = optimize(f = loglikelihood.gamma.par, interval = c(0.1, 4), maximum=TRUE)$maximum 
jpmf.loop = c(jpmf.loop, jpmf(variable = variable, lambda = lambda.new, gamma.par = gamma.par.new))
diff = abs(jpmf.loop[length(jpmf.loop)]-jpmf.loop[length(jpmf.loop)-1])
gamma.par.old = gamma.par.new 
lambda.old = lambda.new
n.iter = n.iter + 1
}
outlist = list(lambda=lambda.new, gamma.par=gamma.par.new, n.iter=n.iter, likelihood.values=jpmf.loop)
    class(outlist) = "pg.dist.glm"
    return(outlist)
}
