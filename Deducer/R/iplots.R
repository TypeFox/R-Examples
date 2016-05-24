# TODO: Add comment
# 
# Author: ianfellows
###############################################################################



########################################################################
#
#				iplots: Scatter
#
########################################################################

.makeIplotDialog <- function(){
	RFunction <- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")
	ParamVariable <- J("org.rosuda.deducer.widgets.param.ParamVariable")
	ParamCharacter<- J("org.rosuda.deducer.widgets.param.ParamCharacter")
	
	rf <- new(RFunction,"Interactive Scatter Plot")
	rf$setName("iplot")
	rf$setRequiresVariableSelector(TRUE)
	
	var1 <- new(ParamVariable,"x")
	rf$add(var1)
	
	var2 <- new(ParamVariable,"y")
	rf$add(var2)	
	
	pcx <- new(ParamCharacter,"xlab")
	pcx$setTitle("x label")
	rf$add(pcx)
	
	pcy <- new(ParamCharacter,"ylab")
	pcy$setTitle("y label")
	rf$add(pcy)	
	
	rfd <- new(RFunctionDialog, rf)
	rfd$setLocationRelativeTo(.jnull())
	rfd$setSize(460L,350L)
	rfd
}

########################################################################
#
#				iplots: Histogram
#
########################################################################

.makeIhistDialog <- function(){
	RFunction <- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")
	ParamVariable <- J("org.rosuda.deducer.widgets.param.ParamVariable")
	ParamNumeric<- J("org.rosuda.deducer.widgets.param.ParamNumeric")
	
	rf <- new(RFunction)
	rf$setName("ihist")
	rf$setTitle("Interactive Histogram")
	rf$setRequiresVariableSelector(TRUE)
	
	var1 <- new(ParamVariable,"var")
	rf$add(var1)
	
	
	pcx <- new(ParamNumeric,"binw")
	pcx$setTitle("Bin Width")
	pcx$setRequired(FALSE)
	pcx$setLowerBound(0)
	rf$add(pcx)
	
	rfd <- new(RFunctionDialog, rf)
	rfd$setLocationRelativeTo(.jnull())
	rfd$setSize(460L,300L)
	rfd
}

########################################################################
#
#				iplots: Bar
#
########################################################################

.makeIbarDialog <- function(){
	RFunction <- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")
	ParamVariable <- J("org.rosuda.deducer.widgets.param.ParamVariable")
	ParamLogical<- J("org.rosuda.deducer.widgets.param.ParamLogical")
	
	rf <- new(RFunction)
	rf$setName("ibar")
	rf$setTitle("Interactive Bar Chart")
	rf$setRequiresVariableSelector(TRUE)
	
	var1 <- new(ParamVariable,"var")
	rf$add(var1)
	
	
	pcx <- new(ParamLogical,"isSpine")
	pcx$setTitle("Spine Plot")
	rf$add(pcx)
	
	rfd <- new(RFunctionDialog, rf)
	rfd$setLocationRelativeTo(.jnull())
	rfd$setSize(460L,300L)
	rfd
}



########################################################################
#
#				iplots: Mosaic
#
########################################################################

.makeImosaicDialog <- function(){
	RFunction <- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")
	ParamMultipleVariables <- J("org.rosuda.deducer.widgets.param.ParamMultipleVariables")
	ParamCharacter<- J("org.rosuda.deducer.widgets.param.ParamCharacter")
	
	rf <- new(RFunction)
	rf$setName("imosaic")
	rf$setTitle("Interactive Mosaic Chart")
	rf$setRequiresVariableSelector(TRUE)
	
	var1 <- new(ParamMultipleVariables,"vars")
	rf$add(var1)
	
	
	pcx <- new(ParamCharacter,"type")
	pcx$setTitle("Type")
	pcx$setOptions(c("observed", "expected", "fluctuation","same.bin.size",
					"multiple.barchart"))
	pcx$setViewType(pcx$VIEW_COMBO)
	rf$add(pcx)
	
	rfd <- new(RFunctionDialog, rf)
	rfd$setLocationRelativeTo(.jnull())
	rfd$setSize(460L,400L)
	rfd
}


########################################################################
#
#				iplots: Box long
#
########################################################################

.makeIboxLongDialog <- function(){
	RFunction <- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")
	ParamVariable <- J("org.rosuda.deducer.widgets.param.ParamVariable")
	
	
	rf <- new(RFunction)
	rf$setName("ibox")
	rf$setTitle("Interactive Box Plot")
	rf$setRequiresVariableSelector(TRUE)
	
	var1 <- new(ParamVariable,"y")
	var1$setName(.jnull())
	rf$add(var1)
	
	var2 <- new(ParamVariable,"x")
	var2$setName(.jnull())
	rf$add(var2)	
	
	
	rfd <- new(RFunctionDialog, rf)
	rfd$setLocationRelativeTo(.jnull())
	rfd$setSize(460L,300L)
	rfd
}

########################################################################
#
#				iplots: Box wide
#
########################################################################

.makeIboxWideDialog <- function(){
	RFunction <- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")
	ParamMultipleVariables <- J("org.rosuda.deducer.widgets.param.ParamMultipleVariables")
	ParamCharacter<- J("org.rosuda.deducer.widgets.param.ParamCharacter")
	
	rf <- new(RFunction)
	rf$setName("ibox")
	rf$setTitle("Interactive Box Plot")
	rf$setRequiresVariableSelector(TRUE)
	
	var1 <- new(ParamMultipleVariables,"Variables")
	var1$setName(.jnull())
	rf$add(var1)
	
	rfd <- new(RFunctionDialog, rf)
	rfd$setLocationRelativeTo(.jnull())
	rfd$setSize(460L,400L)
	rfd
}




########################################################################
#
#				iplots: Par coord plot
#
########################################################################

.makeIpcpDialog <- function(){
	RFunction <- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")
	ParamMultipleVariables <- J("org.rosuda.deducer.widgets.param.ParamMultipleVariables")
	ParamCharacter<- J("org.rosuda.deducer.widgets.param.ParamCharacter")
	
	rf <- new(RFunction)
	rf$setName("ipcp")
	rf$setTitle("Interactive Parallel Coordinate Plot")
	rf$setRequiresVariableSelector(TRUE)
	
	var1 <- new(ParamMultipleVariables,"Variables")
	var1$setName(.jnull())
	rf$add(var1)
	
	rfd <- new(RFunctionDialog, rf)
	rfd$setLocationRelativeTo(.jnull())
	rfd$setSize(460L,400L)
	rfd
}
