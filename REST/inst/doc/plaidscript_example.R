

plaid_WINDOW <- function(list.info=list()){
	
	##########################
	## PREAMBLE/INFORMATION ##
	##########################
	
	dialogtitle <- "Plaid Biclustering"
	
	
	usetabs <- TRUE
	
	tabnames <- c("Biclustering","Plot & Diagnostics")
	
	if(usetabs){ntabs <- length(tabnames)} else {ntabs <- 1}
	new.frames <- .initialize.new.frames(ntabs)
	grid.config <- .initialize.grid.config(ntabs)
	grid.rows <- .initialize.grid.rows(ntabs)
	
	
	helppage <- "BCPlaid"
	
	
	##################
	## GRID BUTTONS ##
	##################
	
	make.help.button <- TRUE
	make.setwd.button <- FALSE
	make.resetgws.button <- FALSE
	make.seed.button <- TRUE
	
	###########
	## TAB 1 ##
	###########
	Tab <- 1
	
	### 1. ADDING THE FRAMES ###
	
	####		RADIO BUTTONS FRAME 		 ####
	#                               			#
	
	type <- "radiobuttons"
	
	# Change variables accordingly:
	frame.name <- "toclusterframe"
	argument.names <- c("Rows","Columns","Rows & Columns")   
	arguments <- c("cluster")		
	argument.values <- c("r","c","b")
	argument.types <- "char"
	initial.values <- "b" 
	title <- "To Cluster"
	border <- FALSE
	
	# DO NOT CHANGE THIS LINE:
	new.frames <- .add.frame(Tab=Tab,type=type,frame.name=frame.name,argument.names=argument.names,arguments=arguments,argument.values=argument.values,initial.values=initial.values,title=title,border=border,new.frames=new.frames,argument.types=argument.types)	
	
	
	######		  ENTRY FIELDS FRAME 				#####
	#							    		 			#
	
	type <- "entryfields"
	
	# Change variables accordingly:
	frame.name <- "modelframe"  
	argument.names <- c("Model Formula") 
	argument.types <- c("num")
	arguments <- c("fit.model")
	initial.values <- c("y ~ m+a+b")
	title <- "Model"
	border <- FALSE
	entry.width <- c("10")  
	
	# Do not change this line:
	new.frames <- .add.frame(Tab=Tab,type=type,frame.name=frame.name,argument.names=argument.names,arguments=arguments,initial.values=initial.values,title=title,border=border,entry.width=entry.width,argument.types=argument.types  ,new.frames=new.frames)
	
	####		CHECK BOXES FRAME 			  ####
	#                               			 #
	
	type <- "checkboxes"
	
	# Change variables accordingly:
	frame.name <-  "backgroundcheckframe"
	argument.names <- c("Background Layer?") 
	arguments <- c("background") 
	initial.values <- c(1) 
	title <- ""
	border <- FALSE
	
	# DO NOT CHANGE THIS LINE:
	new.frames <- .add.frame(Tab=Tab,type=type,frame.name=frame.name,argument.names=argument.names,arguments=arguments,initial.values=initial.values,title=title,border=border,new.frames=new.frames)
	
	
	######		  ENTRY FIELDS FRAME 				#####	
	#							    		 			#
	
	type <- "entryfields"
	
	# Change variables accordingly:
	frame.name <- "backgroundentryframe1"  
	argument.names <- c("Shuffle","Back Fit","Max Layes") 
	argument.types <- c("num","num","num")
	arguments <- c("shuffle","back.fit","max.layers")
	initial.values <- c(3,0,20)
	title <- ""
	border <- FALSE
	entry.width <- c("3","3","3")  
	
	# Do not change this line:
	new.frames <- .add.frame(Tab=Tab,type=type,frame.name=frame.name,argument.names=argument.names,arguments=arguments,initial.values=initial.values,title=title,border=border,entry.width=entry.width,argument.types=argument.types  ,new.frames=new.frames)
	
	
	######		  ENTRY FIELDS FRAME 				#####
	#							    		 			#
	
	type <- "entryfields"
	
	# Change variables accordingly:
	frame.name <- "backgroundentryframe2"  
	argument.names <- c("Iteration Startup","Iteration Layer") 
	argument.types <- c("num","num")
	arguments <- c("iter.startup","iter.layer")
	initial.values <- c(5,10)
	title <- ""
	border <- FALSE
	entry.width <- c("3","3")  
	
	# Do not change this line:
	new.frames <- .add.frame(Tab=Tab,type=type,frame.name=frame.name,argument.names=argument.names,arguments=arguments,initial.values=initial.values,title=title,border=border,entry.width=entry.width,argument.types=argument.types  ,new.frames=new.frames)
	
	#### MANUAL BUTTONS FRAME ####

	type <- "buttons"
	
	# Change variables accordingly:
	frame.name <- "plaidbutton"  
	button.name <- "Plaid"  
	button.function <- "biclust" 
	button.data <- "x" 
	button.object <-  "PlaidResult" 
	button.width <- "12"
	button.data.transf <- "matrix" # only matrix available here !
	
	arg.frames <- c("toclusterframe","modelframe","backgroundcheckframe","backgroundentryframe1","backgroundentryframe2")
	
	save <- TRUE 
	show <- TRUE
	show.save <- TRUE
	button.otherarg <- ",method=BCPlaid()" 
	
	# Do not change this line: 
	new.frames <- .add.frame(Tab=Tab,frame.name=frame.name,
			type=type,button.name=button.name,button.width=button.width,
			button.data.transf=button.data.transf,
			button.function=button.function,button.data=button.data,
			button.object=button.object,button.otherarg=button.otherarg,
			arg.frames=arg.frames,save=save,show=show,new.frames=new.frames,show.save=show.save)
	
	### 2. CONFIGURING THE GRID ###
	grid.config <- .grid.matrix(Tab=Tab,c("toclusterframe","modelframe","backgroundcheckframe",NA,"backgroundentryframe1","backgroundentryframe2","plaidbutton",NA),byrow=TRUE,nrow=4,ncol=2,grid.config=grid.config)
	
	
	### 3. COMBINING THE ROWS ###
	grid.rows <- .combine.rows(Tab=Tab,rows=c(1),title="Plaid Specifications",border=TRUE,grid.rows=grid.rows,grid.config=grid.config)
	grid.rows <- .combine.rows(Tab=Tab,rows=c(2,3),title="Layer Specifications",border=TRUE,grid.rows=grid.rows,grid.config=grid.config)
	
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


