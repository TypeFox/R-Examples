## ----newwindow_script1,echo=TRUE,eval=FALSE------------------------------
#  GUI_WINDOW <- function(list.info=list()){
#  	
#  	##########################
#  	## PREAMBLE/INFORMATION ##
#  	##########################
#  	
#  	dialogtitle <- "This is the title of the window"
#  	
#  	usetabs <- TRUE
#  	tabnames <- c("Tab 1","Tab 2","Tab 3")
#  	helppage <- "plot"
#  	
#  	# Do not change these lines
#  	if(usetabs){ntabs <- length(tabnames)} else {ntabs <- 1}
#  	new.frames <- .initialize.new.frames(ntabs)
#  	grid.config <- .initialize.grid.config(ntabs)
#  	grid.rows <- .initialize.grid.rows(ntabs)
#  	### end of "Do not change these lines"
#  	
#  	##################
#  	## GRID BUTTONS ##
#  	##################
#  	
#  	make.help.button <- TRUE
#  	make.setwd.button <- TRUE
#  	make.resetgws.button <- TRUE
#  	make.seed.button <- TRUE
#  	
#  	# ... continuation of the script down below (these 2 parts are put here)	
#  
#  }	# Note: The curly bracket is placed here for syntax reasons.
#  	#       It should be placed after the call of GUI_template.

## ----newwindow_script2,eval=FALSE,echo=TRUE------------------------------
#  ###########
#  ## TAB 1 ##
#  ###########
#  Tab <- 1
#  
#  ### 1. ADDING THE FRAMES ###
#  
#  # Add frames here
#  
#  
#  ### 2. CONFIGURING THE GRID ###
#  grid.config <- .grid.matrix(Tab=Tab,c("frame1","frame2","frame3",NA),
#  		byrow=TRUE,nrow=2,ncol=2,grid.config=grid.config)
#  
#  
#  ### 3. COMBINING THE ROWS ###
#  grid.rows <- .combine.rows(Tab=Tab,rows=c(1,2),title="A nice box: ",
#  		border=TRUE,grid.rows=grid.rows,grid.config=grid.config)
#  

## ----newwindow_script3,echo=TRUE,eval=FALSE------------------------------
#  #############
#  ### TAB 2 ###
#  #############
#  Tab <- 2
#  
#  # Repeat the 3 steps of tab 1 for as many tabs as you like...
#  
#  
#  ###################################################################
#  ## USE ALL THE ARGUMENTS IN THE GENERAL GUI_TEMPLATE FUNCTION    ##
#  ###################################################################
#  GUI_template(dialogtitle=dialogtitle,helppage=helppage,make.resetgws.button=
#  		make.resetgws.button,make.setwd.button=make.setwd.button,
#  		make.help.button=make.help.button,make.seed.button=make.seed.button,
#  		usetabs=usetabs,tabnames=tabnames,grid.config=grid.config,grid.rows=grid.rows,
#  		new.frames=new.frames)
#  

## ----example_script1,echo=TRUE,eval=FALSE--------------------------------
#  plaid_WINDOW <- function(list.info=list()){
#  	
#      ##########################
#      ## PREAMBLE/INFORMATION ##
#      ##########################
#  	
#      dialogtitle <- "Plaid Biclustering"
#  		
#      usetabs <- TRUE
#      tabnames <- c("Biclustering","Plot & Diagnostics")
#  	
#      if(usetabs){ntabs <- length(tabnames)} else {ntabs <- 1}
#      new.frames <- .initialize.new.frames(ntabs)
#      grid.config <- .initialize.grid.config(ntabs)
#      grid.rows <- .initialize.grid.rows(ntabs)
#  		
#      helppage <- "BCPlaid"
#  		
#      ##################
#      ## GRID BUTTONS ##
#      ##################
#  	
#      make.help.button <- TRUE
#      make.setwd.button <- FALSE
#      make.resetgws.button <- FALSE
#      make.seed.button <- TRUE
#  	
#  	#... followed by tabs, frames,...
#  }

## ----example_script2,echo=TRUE,eval=FALSE--------------------------------
#  ### 2. CONFIGURING THE GRID ###
#  grid.config <- .grid.matrix(Tab=Tab,c("toclusterframe","modelframe","backgroundcheckframe"
#  		,NA,"backgroundentryframe1","backgroundentryframe2","plaidbutton",NA),
#  		byrow=TRUE,nrow=4,ncol=2,grid.config=grid.config)
#  
#  
#  ### 3. COMBINING THE ROWS ###
#  grid.rows <- .combine.rows(Tab=Tab,rows=c(1),title="Plaid Specifications",border=TRUE,
#  		grid.rows=grid.rows,grid.config=grid.config)
#  grid.rows <- .combine.rows(Tab=Tab,rows=c(2,3),title="Layer Specifications",border=TRUE,
#  		grid.rows=grid.rows,grid.config=grid.config)
#  
#  
#  #### Plot & Diagnostics Tab has been omitted ###
#  
#  
#  ###################################################################
#  ## USE ALL THE ARGUMENTS IN THE GENERAL GUI_TEMPLATE FUNCTION    ##
#  ###################################################################
#  GUI_template(dialogtitle=dialogtitle,helppage=helppage,make.resetgws.button=
#  	make.resetgws.button,make.setwd.button=make.setwd.button,make.help.button=
#  	make.help.button,make.seed.button=make.seed.button,usetabs=usetabs,tabnames=
#  	tabnames,grid.config=grid.config,grid.rows=grid.rows,new.frames=new.frames)
#  

## ----doitandprint_example,echo=TRUE,eval=FALSE---------------------------
#  CenterColumns <- function(type=c("mean","median")){
#  	command <- paste0("center.vector <- apply(",ActiveDataSet(),
#  			               ",MARGIN=2,FUN=",type,")")
#  	doItAndPrint(command)
#  	doItAndPrint("center.vector")
#  	
#  	doItAndPrint(paste0("plot(center.vector,main='Column Centers',xlab='Columns',
#                                 ylab='Value')"))
#  }	

