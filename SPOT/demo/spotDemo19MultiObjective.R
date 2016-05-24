##
## use demo(... , ask=F) to run the demo all at once, without interruption
##
## Load package SPOT
require("SPOT")
## Load package MCO, which contains the target function (ZDT2)
require(mco)
## set up a list of settings for SPOT:
lower=c(0,0) #boundaries of decision space, in this case two dimensional
upper=c(1,1)
config=list(
	alg.func=zdt2, #target function, multi objective (2 objectives)
	alg.roi=spotROI(lower,upper), #constructor for ROI
	seq.predictionModel.func="spotPredictEarth", #multivariate adaptive regression splines, surrogate model
	io.verbosity=3,
	seq.design.oldBest.size=0,seq.design.new.size=10, 
	seq.predictionOpt.func="spotModelParetoOptim", #use pareto optimization function to optimize build surrogate models
	seq.predictionOpt.method="sms-emoa",  #use sms emoa for the above mentioned pareto optimization
	seq.predictionOpt.budget=2000, #the budget for sms emoa
	seq.predictionOpt.psize=20, #population size for sms emoa
	spot.ocba=FALSE,seq.design.maxRepeats = 1,init.design.repeats = 1, #some settings due to non noisy target
	auto.loop.nevals=40, #budget for SPOT, number of evaluations on zdt2
	spot.seed=0, #seed for random number generator
	spot.fileMode=F) #write no file output
	
res<-spot(spotConfig=config)
	
								   