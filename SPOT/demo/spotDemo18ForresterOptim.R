##
## use demo(... , ask=F) to run the demo all at once, without interruption
##
## Load package
require("SPOT")
## set up a list of settings for SPOT:
lower=c(-5,-5) #lower boundary of region of interest
upper=c(5,5) # upper boundary
config=list(
	alg.func=spotRastriginFunction, #select the target function to be optimized
	alg.roi=spotROI(lower,upper), #build region of interest
	seq.predictionModel.func="spotPredictForrester", #select the used surrogate model, here a kriging model
	seq.predictionOpt.func="spotModelOptim", #select a function to optimize on the surrogate model
	seq.predictionOpt.method="cmaes",	#select an optimization algorithm for the above function
	io.verbosity=3, #verbosity is set high, all output is given to user
	seq.design.oldBest.size=0,#no noise, so old design points will not be reevaluatied
	seq.design.new.size=2, #in each sequential step, two new points will be evaluated
	spot.ocba=FALSE,#target function is non noisy, so this has to be false
	seq.design.maxRepeats = 1, #each design point is evaluated once due to no noise
	init.design.repeats = 1, #each design point is evaluated once due to no noise
	auto.loop.nevals=25, #allowed budget of target function evaluations
	spot.seed=0, #seed for random number generator
	spot.fileMode=F) #write no file output
	
## start optimization with spot	
res1<-spot(spotConfig=config)	

## or use the optim-like interface (resembles the optim() function)
config2=list(
	seq.predictionOpt.func="spotModelOptim", #select a function to optimize on the surrogate model
	seq.predictionOpt.method="cmaes",	#select an optimization algorithm for the above function
	io.verbosity=3, #verbosity is set high, all output is given to user
	seq.design.oldBest.size=0,#no noise, so old design points will not be reevaluatied
	seq.design.new.size=2, #in each sequential step, two new points will be evaluated
	spot.ocba=FALSE,#target function is non noisy, so this has to be false
	seq.design.maxRepeats = 1, #each design point is evaluated once due to no noise
	init.design.repeats = 1, #each design point is evaluated once due to no noise
	maxit=25, #allowed budget of target function evaluations
	seed=0) #seed for random number generator	
res2<-spotOptim(par=c(NA,NA), #only used to define dimension of decision space
		upper=upper,
		lower=lower,
		method="spotPredictForrester",#select the used surrogate model, here a kriging model
		fn=spotRastriginFunction,  #select the target function to be optimized
		control=config2) #other settings
