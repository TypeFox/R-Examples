

GUI_WINDOW <- function(list.info=list()){
	
	##########################
	## PREAMBLE/INFORMATION ##
	##########################
	
	dialogtitle <- "This is the title of the window"

	usetabs <- TRUE
	tabnames <- c("Tab 1","Tab 2","Tab 3")
	helppage <- "plot" 
	
	# Do not change these lines
	if(usetabs){ntabs <- length(tabnames)} else {ntabs <- 1}
	new.frames <- .initialize.new.frames(ntabs)
	grid.config <- .initialize.grid.config(ntabs)
	grid.rows <- .initialize.grid.rows(ntabs)
	### end of "Do not change these lines"
	
	##################
	## GRID BUTTONS ##
	##################
	
	make.help.button <- TRUE
	make.setwd.button <- TRUE
	make.resetgws.button <- TRUE
	make.seed.button <- TRUE

	###########
	## TAB 1 ##
	###########
	Tab <- 1

	### 1. ADDING THE FRAMES ###
	
	# Add frames here
	

	### 2. CONFIGURING THE GRID ###
	grid.config <- .grid.matrix(Tab=Tab,c("frame1","frame2","frame3",NA),byrow=TRUE,nrow=2,ncol=2,grid.config=grid.config)


	### 3. COMBINING THE ROWS ###
	grid.rows <- .combine.rows(Tab=Tab,rows=c(1,2),title="A nice box: ",border=TRUE,grid.rows=grid.rows,grid.config=grid.config)
	
	#############
	### TAB 2 ###
	#############
	Tab <- 2

	# Repeat what you did for tab 1 for as many tabs as you like...
	
	
	###################################################################
	## USE ALL THE ARGUMENTS IN THE GENERAL GUI_TEMPLATE FUNCTION    ##
	###################################################################
	GUI_template(dialogtitle=dialogtitle,helppage=helppage,make.resetgws.button=make.resetgws.button,make.setwd.button=make.setwd.button,make.help.button=make.help.button,make.seed.button=make.seed.button,usetabs=usetabs,tabnames=tabnames,grid.config=grid.config,grid.rows=grid.rows,new.frames=new.frames)
	
}


