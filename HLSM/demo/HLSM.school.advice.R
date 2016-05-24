#library(HLSM)
#data(schools-advice-data)

message ("Predictors: Edge variables")
print (ps.edge.vars.mat[[1]][1:10,1:10,])

message ("Outcome: Advice seeking")
print (ps.advice.mat[[1]][1:10,1:10])

#Random Regression Coefficients#
priors = NULL
tune = NULL
initialVals = NULL
niter = 10

random.fit = HLSM(X = ps.edge.vars.mat,Y = ps.advice.mat,
	initialVals = initialVals,priors = priors,
	tune = tune,tuneIn = FALSE,dd = 2,niter = niter,
	intervention = 0)


summary(random.fit)
names(random.fit)

fixed.fit = HLSMfixedEF(X = ps.edge.vars.mat,Y = ps.advice.mat,
	initialVals = initialVals,priors = priors,
	tune = tune,tuneIn = FALSE,dd = 2,niter = niter,
	intervention = 0)
summary(fixed.fit)
names(fixed.fit)





