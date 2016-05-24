# TODO: Add comment
# 
# Author: ianfellows
###############################################################################




########################################################################
#
#				paired test
#
########################################################################

.makePairedTestDialog <- function(){
	RFunction <- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")
	Param<- J("org.rosuda.deducer.widgets.param.Param")
	ParamAny <- J("org.rosuda.deducer.widgets.param.ParamAny")
	ParamVariable <- J("org.rosuda.deducer.widgets.param.ParamVariable")
	ParamMultipleVariables <- J("org.rosuda.deducer.widgets.param.ParamMultipleVariables")
	ParamLogical <- J("org.rosuda.deducer.widgets.param.ParamLogical")
	ParamCharacter<- J("org.rosuda.deducer.widgets.param.ParamCharacter")
	ParamNumeric<- J("org.rosuda.deducer.widgets.param.ParamNumeric")
	ParamRObject<- J("org.rosuda.deducer.widgets.param.ParamRObject")
	ParamRFunctionResult<- J("org.rosuda.deducer.widgets.param.ParamRFunctionResult")
	RFunctionList<- J("org.rosuda.deducer.widgets.param.RFunctionList")
	RFunctionListDialog <- J("org.rosuda.deducer.widgets.param.RFunctionListDialog")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")	
	
	
	
	
	#functions for dialog
	rf <- new(RFunction,"t.test")
	rfw <- new(RFunction,"wilcox.test")
	
	
	#make parameters and add them to the functions
	pv1 <- new(ParamVariable,"x")
	pv1$setTitle("First")
	rf$add(pv1)
	rfw$add(pv1)
	
	
	pv2 <- new(ParamVariable,"y")
	pv2$setTitle("Second")
	rf$add(pv2)
	rfw$add(pv2)
	
	
	p <- new(ParamCharacter,"alternative","two.sided")
	p$setTitle("Type")
	p$setOptions(c("two.sided","less","greater"))
	p$setViewType(p$VIEW_COMBO)
	rf$add(p)
	rfw$add(p$clone())
	
	p <- new(ParamNumeric,"conf.level",.95)
	p$setLowerBound(0)
	p$setUpperBound(1)
	rf$add(p)
	rfw$add(p$clone())
	
	p <- new(ParamLogical,"paired",TRUE)
	p$setDefaultValue(FALSE)
	p$setViewType(p$VIEW_HIDDEN)
	rf$add(p)
	rfw$add(p)
	
	#make function list and add functions
	prf <- new(RFunctionList,"Paired tests");
	prf$setViewType(p$VIEW_RFUNCTION_PANEL)
	prf$addRFunction("t-test",rf,FALSE)
	prf$addRFunction("Wilcoxon signed rank",rfw,FALSE)
	prf$setRequiresVariableSelector(TRUE)
	prf$setActiveFunctions("t-test")
	#parameters common to all functions should go at the top
	globals <- .jarray(list(pv1,pv2),"org.rosuda.deducer.widgets.param.Param")
	prf$setGlobalParams(globals)
	
	
	#make dialog and display
	rfd <- new(RFunctionListDialog, prf )
	rfd$setSize(500L,520L)
	rfd$setLocationRelativeTo(.jnull())
	rfd
}
