pg.reg.glm <-
function(variable, regressors, lambda.start=1, gamma.par.start=1, epsylon = 10^(-8), n=100){
lambda.old = c(); gamma.par.old = c()
lambda.new = c(); gamma.par.new = c()
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
gamma.par.old = gamma.par.start
diff = 0.1
n.iter=0
rules = laguerre.quadrature.rules(n)
rule = rules[[n]]
jpmf = function(variable, lambda, gamma.par) {
mix = function(theta) {
y=(gamma.par^(gamma.par)*theta^(gamma.par-1)*exp(-(gamma.par*theta)))/(gamma(gamma.par))
for (i in 1:length(variable)) {
y=y*((lambda[i]*theta)^(variable[i])*exp(-lambda[i]*theta)/factorial(variable[i]))
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
y=(lambda.old[i]*theta)^(variable[i]+1)*exp(-lambda.old[i]*theta)*Gamma.density(theta)
return(y)
}
nominator[i] = laguerre.quadrature(t.r.nominator, rule, lower = 0, upper = Inf, weighted = FALSE)
t.r.denominator = function(theta) {
y=(lambda.old[i]*theta)^(variable[i])*exp(-lambda.old[i]*theta)*Gamma.density(theta)
return(y)
}
denominator[i] = laguerre.quadrature(t.r.denominator, rule, lower = 0, upper = Inf, weighted = FALSE)
}
t = nominator/denominator
##########   M-step
model.formula.Mstep=paste(model.formula, "+log(t)", sep="")
g=glm(eval(parse(text=model.formula.Mstep)), family=poisson(log), data=cbind(variable, regr, t))
lambda.new=exp(X%*%g$coefficients[-length(g$coefficients)])
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
outlist = list(lambda=lambda.new, gamma.par=gamma.par.new, regression.coefficients=g$coefficients[-length(g$coefficients)], n.iter=n.iter, likelihood.values=jpmf.loop)
    class(outlist) = "pg.reg.glm"
    return(outlist)
}
