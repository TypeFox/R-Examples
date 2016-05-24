##
## use demo(... , ask=F) to run the demo all at once, without interruption
##
## Load package
require("SPOT")
## set up a list of settings for SPOT:
lower=c(-2,-2,-2) #three dimensional region of interest
upper=c(2,2,2)
config=list(
	alg.func=function(x){sum(x^2)}, #target function is sphere
	alg.roi=spotROI(lower,upper),
	seq.forrester=list(budgetalgtheta=100),
	seq.predictionModel.func="spotPredictForrester",
	seq.predictionOpt.func="spotModelOptim",seq.predictionOpt.method="cmaes",	
	io.verbosity=3,report.interactive=TRUE,
	seq.design.oldBest.size=0,seq.design.new.size=2,
	spot.noise=0,spot.ocba=FALSE,seq.design.maxRepeats = 1,init.design.repeats = 1,
	auto.loop.nevals=30,spot.seed=0,
	spot.fileMode=F)#write no file output
	
## start optimization with spot	
res<-spot(spotConfig=config)	

## do an interactive 3d plot report
## chose variables and click eval in the GUI that pops up
## close the GUI to proceed
repres<-spot(spotConfig=append(list(report.func="spotReport3d"),res),spotTask="rep")

## do a non-interactive 3d plot report, specifiying which parameters to plot
repres<-spot(spotConfig=append(list(report.func="spotReport3d", report.main="This Is the Title!", report.interactive=F, report.aIndex=2, report.bIndex=3),res),spotTask="rep")

## do a non-interactive contour plot report, specifiying which parameters to plot
repres<-spot(spotConfig=append(list(report.func="spotReportContour", report.main="This Is the Title!", report.interactive=F, report.aIndex=2, report.bIndex=3),res),spotTask="rep")

