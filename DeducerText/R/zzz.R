
.makeQuickWordcloud <- function(){
	RFunction <- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")	
	ParamVariable <- J("org.rosuda.deducer.widgets.param.ParamVariable")
	ParamLogical <- J("org.rosuda.deducer.widgets.param.ParamLogical")
	ParamNumeric	<- J("org.rosuda.deducer.widgets.param.ParamNumeric")
	func <- new(RFunction,"wordcloud")
	var1 <- new(ParamVariable,"words")
	func$add(var1)
	
	ord <- new(ParamLogical,"random.order",TRUE)
	ord$setTitle("Random order")
	ord$setValue(FALSE)
	func$add(ord)

	freq <- new(ParamNumeric,"min.freq",3)
	freq$setValue(0)
	freq$setTitle("Min count")
	freq$setLowerBound(0)
	func$add(freq)
	
	mw <- new(ParamNumeric,"max.words")
	mw$setValue(200)
	mw$setTitle("# of words")
	mw$setLowerBound(1)
	func$add(mw)	
	
	rfd <- new(RFunctionDialog, func)
	
	hb <- rfd$getHelpButton()
	hb$setVisible(TRUE)
	hb$setUrl("http://www.deducer.org/pmwiki/index.php?n=Main.TextQuickWordCount")
	
	rfd$setLocationRelativeTo(.jnull())
	rfd
}


.makeTextPlot <- function(){
	RFunction <- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")	
	ParamVariable <- J("org.rosuda.deducer.widgets.param.ParamVariable")
	ParamLogical <- J("org.rosuda.deducer.widgets.param.ParamLogical")
	ParamNumeric	<- J("org.rosuda.deducer.widgets.param.ParamNumeric")
	func <- new(RFunction,"textplot")
	var1 <- new(ParamVariable,"x")
	func$add(var1)
	var2 <- new(ParamVariable,"y")
	func$add(var2)
	var3 <- new(ParamVariable,"words")
	var3$setTitle("Text")
	func$add(var3)

	freq <- new(ParamNumeric,"cex",1)
	freq$setTitle("Size")
	freq$setLowerBound(0)
	func$add(freq)

	
	rfd <- new(RFunctionDialog, func)
	
	hb <- rfd$getHelpButton()
	hb$setVisible(TRUE)
	hb$setUrl("http://www.deducer.org/pmwiki/index.php?n=Main.TextTextPlot")
	
	rfd$setLocationRelativeTo(.jnull())
	rfd
}


.onLoad <- function(libname, pkgname) {
		if(!.jgr)
			return(TRUE)
        .jpackage("DeducerText");
        x <- .jnew(J("edu.cens.text.Text"));
        x <- try(x$initJGR(), silent=TRUE);
		
		.registerDialog("Quick Word Cloud",.makeQuickWordcloud)
		.registerDialog("Text Plot",.makeTextPlot)
		
		
		
        return(TRUE);
}
