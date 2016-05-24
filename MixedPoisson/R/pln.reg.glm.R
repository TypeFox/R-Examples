pln.reg.glm <-
function(variable, regressors, lambda.start=1, nu.start=1, epsylon = 10^(-8), n=100){
lambda.old = c(); nu.old = c()
lambda.new = c(); nu.new = c()
regr = as.data.frame(regressors)
if (ncol(regr)>0) {regressors.names=names(regr)} else regressors.names='1'
model.formula.X=" ~ "
for (k in 1:length(regressors.names)) {
if (k==1) {model.formula.X=paste(model.formula.X, regressors.names[k], sep="")} else {model.formula.X=paste(model.formula.X, "+", regressors.names[k], sep="")}
}
X=model.matrix(eval(parse(text=model.formula.X)), data=regr)
model.formula=paste("variable ", model.formula.X, sep="")
g=glm(eval(parse(text=model.formula)), family=poisson(log), data=cbind(variable, regr))
lambda.start=exp(X%*%g$coefficients) 
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
y=y*((lambda[i]*theta)^(variable[i])*exp(-lambda[i]*theta)/factorial(variable[i]))
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
y=(lambda.old[i]*theta)^(variable[i]+1)*exp(-lambda.old[i]*theta)*lognorm.density(theta)
return(y)
}
nominator[i] = laguerre.quadrature(t.r.nominator, rule, lower = 0, upper = Inf, weighted = FALSE)
t.r.denominator = function(theta) {
y=(lambda.old[i]*theta)^(variable[i])*exp(-lambda.old[i]*theta)*lognorm.density(theta)
return(y)
}
denominator[i] = laguerre.quadrature(t.r.denominator, rule, lower = 0, upper = Inf, weighted = FALSE)
}
t = nominator/denominator
##########   M-step
model.formula.Mstep=paste(model.formula, "+log(t)", sep="")
g=glm(eval(parse(text=model.formula.Mstep)), family=poisson(log), data=cbind(variable, regr, t))
lambda.new=exp(X%*%g$coefficients[-length(g$coefficients)])
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
outlist = list(lambda=lambda.new, nu=nu.new, regression.coefficients=g$coefficients[-length(g$coefficients)], n.iter=n.iter, likelihood.values=jpmf.loop)
    class(outlist) = "pln.reg.glm"
    return(outlist)
}
