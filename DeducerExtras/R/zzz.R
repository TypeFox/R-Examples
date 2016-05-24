




.onLoad <- function(libname, pkgname){
	
	#if deducer gui is not running, do minimal load
	deducerLoaded <- try(.deducer == .jnull(),silent=TRUE)
	if(inherits(deducerLoaded,"try-error") || deducerLoaded)
		return(NULL)
	
	RFunction <<- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <<- J("org.rosuda.deducer.widgets.param.RFunctionDialog")
	
	Param<<- J("org.rosuda.deducer.widgets.param.Param")
	ParamAny <<- J("org.rosuda.deducer.widgets.param.ParamAny")
	ParamVariable <<- J("org.rosuda.deducer.widgets.param.ParamVariable")
	ParamMultipleVariables <<- J("org.rosuda.deducer.widgets.param.ParamMultipleVariables")
	ParamLogical <<- J("org.rosuda.deducer.widgets.param.ParamLogical")
	ParamCharacter<<- J("org.rosuda.deducer.widgets.param.ParamCharacter")
	ParamNumeric<<- J("org.rosuda.deducer.widgets.param.ParamNumeric")
	ParamRObject<<- J("org.rosuda.deducer.widgets.param.ParamRObject")
	ParamRFunctionResult<<- J("org.rosuda.deducer.widgets.param.ParamRFunctionResult")
	RFunctionList<<- J("org.rosuda.deducer.widgets.param.RFunctionList")
	RFunctionListDialog <<- J("org.rosuda.deducer.widgets.param.RFunctionListDialog")
	RFunctionDialog <<- J("org.rosuda.deducer.widgets.param.RFunctionDialog")
	
	
	.registerDialog("Distribution quantiles",
			function() .makeDistributionDialog("quantile"))
	
	.registerDialog("Distribution function values",
			function() .makeDistributionDialog("distribution"))
	
	.registerDialog("Cumulative distribution function",
			function() .makeDistributionDialog("CDF"))
	
	.registerDialog("Data summary",.makeSummaryDialog)
	
	.registerDialog("Paired test",.makePairedTestDialog)
	
	.registerDialog("Single proportion",.makeProportionDialog)
	
	.registerDialog("Single proportion: Exact",.makeExactProportionDialog)
	
	.registerDialog("k-sample proportion",.makeNProportionDialog)
	
	.registerDialog("k-sample variance test",.makeEqualVarianceDialog)
	
	.registerDialog("t-test power",.makeTTestPowerDialog)
	
	.registerDialog("k-means cluster",.makeKMeansDialog)
	
	.registerDialog("Apply k-means to data",.makeApplyKMeansDialog)
	
	.registerDialog("Hierarchical cluster",.makeHClustDialog)
	
	.registerDialog("Multi-dimensional scaling",.makeMDSDialog)
	
	.registerDialog("Ranking analysis",.makeRankingAnalysisDialog)
	
	#.registerDialog("3d Scatter plot",.make3dScatterPlotDialog)
	
	.registerDialog("Open Vignette",.makeVignetteDialog)
	
	gui.addSeperator <- function(){}
	if(.windowsGUI){
		winMenuAdd("Extras")
		gui.addMenuItem <- winMenuAddItem
	}else if(.jgr){
		DeducerMain$insertMenu(J("org.rosuda.JGR.JGR")$MAINRCONSOLE, "Extras",6L)
		#jgr.addMenu("Extras")
		gui.addMenuItem <- jgr.addMenuItem
		gui.addSeperator <- function () jgr.addMenuSeparator("Extras")
	}else
		gui.addMenuItem <- function(x,y,z){}
	
	deducer.addMenu("Extras")
	
	'%+%' <- function(x,y) paste(x,y,sep="")
	
	addMenuItem <- function(name){
		deducer.addMenuItem(name,,
				".getDialog('" %+% name %+% "')$run()","Extras")
		gui.addMenuItem("Extras",name,"deducer('"%+% name %+% "')")
	}
	
	addMenuItem("Distribution quantiles")
	addMenuItem("Distribution function values")
	addMenuItem("Cumulative distribution function")
	gui.addSeperator()
	deducer.addMenuItem("Load Data From Package",,
			"DeducerExtras:::.makePackageDataDialog()$run()","Extras")
	gui.addMenuItem("Extras","Load data from package","deducer('"%+% "Load Data From Package" %+% "')")
	addMenuItem("Data summary")
	gui.addSeperator()

	addMenuItem("Single proportion")
	addMenuItem("Single proportion: Exact")
	addMenuItem("k-sample proportion")
	gui.addSeperator()
	addMenuItem("k-sample variance test")
	gui.addSeperator()
	addMenuItem("t-test power")
	gui.addSeperator()	
	addMenuItem("k-means cluster")
	addMenuItem("Apply k-means to data")
	addMenuItem("Hierarchical cluster")
	gui.addSeperator()
	addMenuItem("Multi-dimensional scaling")
	gui.addSeperator()
	addMenuItem("Ranking analysis")
	#gui.addSeperator()
	#addMenuItem("3d Scatter plot")
	gui.addSeperator()
	addMenuItem("Open Vignette")
}




