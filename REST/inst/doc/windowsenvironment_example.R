
main_WINDOW <- function(list.info=list()){
	
	##########################
	## PREAMBLE/INFORMATION ##
	##########################
	
	dialogtitle <- "Example"
		
	usetabs <- FALSE
	
	tabnames <- c("Biclustering","Plot & Diagnostics")
	
	if(usetabs){ntabs <- length(tabnames)} else {ntabs <- 1}
	new.frames <- .initialize.new.frames(ntabs)
	grid.config <- .initialize.grid.config(ntabs)
	grid.rows <- .initialize.grid.rows(ntabs)
		
	helppage <- ""
		
	##################
	## GRID BUTTONS ##
	##################
	
	make.help.button <- FALSE
	make.setwd.button <- FALSE
	make.resetgws.button <- FALSE
	make.seed.button <- FALSE
	
	###########
	## TAB 1 ##
	###########
	Tab <- 1
	
	### 1. ADDING THE FRAMES ###
	
	#### ENTRY FIELDS FRAME ####
	
	type <- "entryfields"
	
	# Change variables accordingly:
	frame.name <- "entryframe1"  
	argument.names <- c("Values?") 
	argument.types <- c("num") 
	arguments <- c("arg1") 
	initial.values <- c("")
	title <- ""
	border <- FALSE
	entry.width <- c("25")
	
	# Do not change this line:
	new.frames <- .add.frame(Tab=Tab,type=type
			,frame.name=frame.name,argument.names=argument.names
			,arguments=arguments,initial.values=initial.values
			,title=title,border=border,entry.width=entry.width
			,argument.types=argument.types  ,new.frames=new.frames)
	

	#### MANUAL BUTTONS FRAME ####
	
	type <- "buttons"
	
	# Change variables accordingly:
	frame.name <- "buttonframe1"  
	button.name <- "Choose"  
	button.function <- "choose_WINDOW" 
	button.data <- "" 
	button.object <-  "saveobject" 
	button.width <- "12"
	button.data.transf <- "matrix" 
	
	arg.frames <- c()
	
	save <- FALSE
	show.save <- FALSE
	show <- FALSE
	button.otherarg <- ""  # always start with a ,
	
	# Do not change this line: 
	new.frames <- .add.frame(Tab=Tab,frame.name=frame.name,
			type=type,button.name=button.name,button.width=button.width,
			button.data.transf=button.data.transf,
			button.function=button.function,button.data=button.data,
			button.object=button.object,button.otherarg=button.otherarg,
			arg.frames=arg.frames,save=save,show=show,show.save=show.save,
			new.frames=new.frames)
	
	### 2. CONFIGURING THE GRID ###
	grid.config <- .grid.matrix(Tab=Tab,c("entryframe1","buttonframe1"),byrow=TRUE,nrow=1,ncol=2,grid.config=grid.config)
		
	### 3. COMBINING THE ROWS ###
	
	###################################################################
	## USE ALL THE ARGUMENTS IN THE GENERAL GUI_TEMPLATE FUNCTION    ##
	###################################################################
	GUI_template(dialogtitle=dialogtitle,helppage=helppage,make.resetgws.button=make.resetgws.button,make.setwd.button=make.setwd.button,make.help.button=make.help.button,make.seed.button=make.seed.button,usetabs=usetabs,tabnames=tabnames,grid.config=grid.config,grid.rows=grid.rows,new.frames=new.frames)
	
}

choose_WINDOW <- function(){
	
	##########################
	## PREAMBLE/INFORMATION ##
	##########################
	
	dialogtitle <- "Choosing"
		
	usetabs <- FALSE
	
	tabnames <- c("Biclustering","Plot & Diagnostics")
	
	if(usetabs){ntabs <- length(tabnames)} else {ntabs <- 1}
	new.frames <- .initialize.new.frames(ntabs)
	grid.config <- .initialize.grid.config(ntabs)
	grid.rows <- .initialize.grid.rows(ntabs)
		
	helppage <- ""
	
	##################
	## GRID BUTTONS ##
	##################
	
	make.help.button <- FALSE
	make.setwd.button <- FALSE
	make.resetgws.button <- FALSE
	make.seed.button <- FALSE
	
	###########
	## TAB 1 ##
	###########
	Tab <- 1
	
	### 1. ADDING THE FRAMES ###
		
	#### LIST BOX FRAME ####
	
	type <- "listbox"
	
	# Change variables accordingly:
	frame.name <- "listboxframe1"
	arguments <- "values"		  # should only be 1
	argument.names <- c("Value 1","Value 2","Value 3")
	argument.values <- c("value1","value2","value3")
	argument.types <- "char"   		  # should be only 1
	initial.values <- c("value3")     # Can be 1 or multiple
	length <- 4  #no character , if not given, will take length of names
	select.multiple <- TRUE
	title <- "Possible Values:"
	border <- TRUE
	
	# DO NOT CHANGE THIS LINE:
	new.frames <- .add.frame(Tab=Tab,type=type,
			frame.name=frame.name,argument.names=argument.names,
			arguments=arguments,argument.values=argument.values,
			argument.types=argument.types, initial.values=initial.values,
			length=length,select.multiple=select.multiple,
			title=title,border=border,new.frames=new.frames)
	
	#### MANUAL BUTTONS FRAME ####
	
	type <- "buttons"
	
	# Change variables accordingly:
	frame.name <- "buttonframe1"  
	button.name <- "Ok"  
	button.function <- "setentry_example" 
	button.data <- "" 
	button.object <-  "saveobject" 
	button.width <- "12"
	button.data.transf <- "matrix" 
	
	arg.frames <- c("listboxframe1")
	
	save <- FALSE
	show.save <- FALSE
	show <- FALSE
	button.otherarg <- ""  # always start with a ,
	
	# Do not change this line: 
	new.frames <- .add.frame(Tab=Tab,frame.name=frame.name,
			type=type,button.name=button.name,button.width=button.width,
			button.data.transf=button.data.transf,
			button.function=button.function,button.data=button.data,
			button.object=button.object,button.otherarg=button.otherarg,
			arg.frames=arg.frames,save=save,show=show,show.save=show.save,
			new.frames=new.frames)
	
	
	### 2. CONFIGURING THE GRID ###
	grid.config <- .grid.matrix(Tab=Tab,c("listboxframe1","buttonframe1"),byrow=TRUE,nrow=2,ncol=1,grid.config=grid.config)
	
	
	### 3. COMBINING THE ROWS ###
	
	###################################################################
	## USE ALL THE ARGUMENTS IN THE GENERAL GUI_TEMPLATE FUNCTION    ##
	###################################################################
	GUI_template(dialogtitle=dialogtitle,helppage=helppage,make.resetgws.button=make.resetgws.button,make.setwd.button=make.setwd.button,make.help.button=make.help.button,make.seed.button=make.seed.button,usetabs=usetabs,tabnames=tabnames,grid.config=grid.config,grid.rows=grid.rows,new.frames=new.frames)
		
}

setentry_example <- function(values){
	
	new.value <- "c("
	for(i in 1:length(values)){
		new.value <- paste0(new.value,"'",values[i],"'")
		if(i!=length(values)){new.value <- paste0(new.value,",")}
	}
	new.value <- paste0(new.value,")")
	
	ChangeWindow("Example",tab=1,"entryframe1","arg1",new.value)
	CancelWindow("Choosing")
}

