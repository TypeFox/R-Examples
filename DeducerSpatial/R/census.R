# TODO: Add comment
# 
# Author: ianfellows
###############################################################################


loadCensusData <- function(state,level=c("county","tract","blkgrp","blk","cdp"),
		year=c("2010","2000"), verbose=TRUE, osmTransform=TRUE, envir = .GlobalEnv){
	level <- match.arg(level)
	year <- match.arg(year)
	pkg <- paste("UScensus",year,level,sep="")
	if(year=="2000" && level=="county")
		stop("County level data not available for the 2000 census")
	if(!suppressWarnings(require(pkg,character.only=TRUE))){
		if(year=="2010")
			library("UScensus2010")
		if(verbose)
			cat(pkg," not found. Installing (This may take some time)...\n")
		os <- if(length(grep("^darwin",R.version$os))>0)
			"osx"
		else if(.Platform$OS.type=="windows")
			"windows"
		else
			"linux"
		if(year=="2010")
			do.call(paste("install.",level,sep=""),list(x=os))
		else
			install.packages(pkg)
		if(!require(pkg,character.only=TRUE))
			stop("Package failed to install. Ensure that you are connected to the internet")
	}
	dataName <- paste(state,".",level,if(year=="2010") "10" else "",sep="")
	if(verbose)
		cat("Loading: ",dataName,"\n")
	data(list=dataName,envir=envir)
	if(osmTransform)
		eval(parse(text=paste(dataName,"<-spTransform(",dataName,",CRS=osm())")),envir=envir)
}




.makeCensusDialog <- function(){
	
	RFunction <- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")
	ParamCharacter <- J("org.rosuda.deducer.widgets.param.ParamCharacter")
	
	states <- c('Alabama', 'Alaska', 'Arizona', 'Arkansas', 
			'California', 'Colorado', 'Connecticut', 'Delaware', 
			'Florida', 'Georgia', 'Hawaii', 'Idaho', 'Illinois',
			'Indiana', 'Iowa', 'Kansas', 'Kentucky', 'Louisiana',
			'Maine', 'Maryland', 'Massachusetts', 'Michigan', 'Minnesota',
			'Mississippi', 'Missouri', 'Montana', 'Nebraska', 'Nevada', 
			'New Hampshire', 'New Jersey', 'New Mexico', 'New York', 
			'North Carolina', 'North Dakota', 'Ohio', 'Oklahoma', 'Oregon', 
			'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota', 
			'Tennessee', 'Texas', 'Utah', 'Vermont', 'Virginia', 'Washington', 
			'West Virginia', 'Wisconsin', 'Wyoming')
	
	rf <- new(RFunction,"loadCensusData")
	rf$setRequiresVariableSelector(FALSE)
	
	st <- new(ParamCharacter,"state")
	st$setValue("california")
	st$setTitle("State")
	st$setOptions(tolower(states))
	st$setViewType(st$VIEW_COMBO)
	rf$add(st)
	
	level <- new(ParamCharacter,"level","county")
	level$setTitle("Level")
	level$setOptions(c("county","tract","blkgrp","cdp"))
	level$setLabels(c("County","Tract","Block group","Census Designated Places"))
	level$setViewType(level$VIEW_COMBO)
	rf$add(level)
	
	year <- new(ParamCharacter,"year","2010")
	year$setTitle("Year")
	year$setOptions(c("2010","2000"))
	year$setViewType(year$VIEW_COMBO)
	rf$add(year)
	
	rfd <- new(RFunctionDialog, rf)
	rfd$setSize(280L,200L)
	rfd$setLocationRelativeTo(.jnull())
	#try(rfd$setLocationRelativeTo(J("org.rosuda.JGR.JGR")$MAINRCONSOLE))
	
	hb <- rfd$getHelpButton()
	hb$setVisible(TRUE)
	hb$setUrl("http://www.deducer.org/pmwiki/index.php?n=Main.SpatialCensus")
	
	rfd
}