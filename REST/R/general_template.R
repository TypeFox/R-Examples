

GUI_template <- function(dialogtitle="",helppage="",usetabs=FALSE,tabnames=c(),make.resetgws.button=FALSE,make.setwd.button=FALSE,make.help.button=FALSE,make.seed.button=FALSE,grid.config=grid.config,grid.rows=grid.rows,new.frames){
	
	
	if(dialogtitle==""){
		stop("No dialogtitle was defined!")
	}
	.updateEnvirObject(dialogtitle,ENVIR=environment())
	
	#########################################################################################################################################################
	## General preparation ##
	#########################
	
	
	
	# Making the frames #
	if(usetabs==TRUE){
		n.tabs <- length(tabnames)
		tabs <- c()
		for(t in 1:n.tabs){
			tabs <- c(tabs,paste("Tab",t,sep=""))
		}
		initializeDialog(title = gettextRcmdr(paste(dialogtitle,sep="")), use.tabs=TRUE,tabs=tabs)
		for(t in 1:n.tabs){
			tab.frame.command <- paste(tabs[t],"Frame <- tkframe(",tabs[t],")",sep="")	
			.eval.command(tab.frame.command)
		}
		
		
	} 
	else{
		n.tabs <- 1
		initializeDialog(title = gettextRcmdr(paste(dialogtitle,sep="")), use.tabs=FALSE)
		Tab1Frame <- tkframe(top)
	}
	
	
	
	
	
	
	
	
	##########################################################################################################################################################
	##########################################################################################################################################################
	
	
	##########################################################################################################################################################
	## THE CLUSTER FRAME & PLOTDIAG FRAME ##
	########################################
	
	
	for(Tab in 1:n.tabs){

	
	if(length(new.frames[[Tab]])!=0){
		
		##########################################################
		## Determining if there are rowframes, if so make them: ##
		##########################################################
		
		if(length(grid.rows[[Tab]])!=0){
			for(i in 1:length(grid.rows[[Tab]])){
				
				
				
				row.command <- paste("Tab",Tab,"Frame_row",i," <- .make.correct.frame(grid.rows[[",Tab,"]][[",i,"]]$title,grid.rows[[",Tab,"]][[",i,"]]$border,Tab",Tab,"Frame)",sep="")
				.eval.command(row.command)
					
				# Special case for a rowframe: NO border, but Title:
				if(grid.rows[[Tab]][[i]]$title!="" & grid.rows[[Tab]][[i]]$border==FALSE){
					temp.command <- paste("tkgrid(labelRcmdr(Tab",Tab,"Frame_row",i ,",fg=getRcmdr('title.color'),font='RcmdrTitleFont',text=gettextRcmdr(grid.rows[[",Tab,"]][[",i,"]]$title)),sticky='nw')"    ,sep="")
					.eval.command(temp.command)
				
					
					
				}
				
				
			}
			
		}
		
		
		
		##########################
		## Making of the frames ##
		##########################
		
		for(ii in 1:length(new.frames[[Tab]])){    # Had to use `ii` as radioButtons had an `i`for loop which was interfering with thisone
			
			current.frame <- new.frames[[Tab]][[ii]]
			frame.name <- current.frame$frame.name
			
			###### Determine if current.frame should be put in a rowframe or clusterframe ( + title,border options) #####
			
			if(length(grid.rows[[Tab]])!=0){			
				temp.names <- lapply(grid.rows[[Tab]],FUN=function(x){return(x$name.frames)})
				boolean <- sapply(temp.names,FUN=function(x){return( frame.name %in% x )})
				
				if(sum(boolean)==1){
					
					window <- paste("Tab",Tab,"Frame_row",which(boolean==TRUE),sep="")	
					
				}	
				else{
					
					window <- paste("Tab",Tab,"Frame",sep="")
					
				}			
			}else {
				
				window <- paste("Tab",Tab,"Frame",sep="")
				
			}
			
			
			##### Make the frame in which the current frame will be placed (for all kinds of types)
			frame.command <- paste("current.frame$frame <- .make.correct.frame(current.frame$title,current.frame$border,",window,")")
			.eval.command(frame.command)
			
		
			
			###### ENTRY FIELDS #########################################################################################################
			
			if(current.frame$type=="entryfields"){
				
#			# Make the entry frame (containing ALL entries)  (DEPENDING ON TITLE AND BORDER VALUES)
#			
				
				# Depending on other kind of boxes, this maybe can be put more upwards				
				
						
				number.entries <- length(current.frame$arguments)
				arguments <- current.frame$arguments
				argument.names <- current.frame$argument.names
				initial.values <- current.frame$initial.values
				
				# Making the element which will contain the separate entry frames
				
				current.frame$entry.frames <- list()
				current.frame$entry.vars <- list()
				current.frame$entry.fields <- list()
				
				# Make title when NOT using a border (special case of 'title = "a title"' & 'border=FALSE'
				
				if(current.frame$title!="" & current.frame$border==FALSE){
					tkgrid(labelRcmdr( current.frame$frame,fg=getRcmdr("title.color"),font="RcmdrTitleFont"   ,text=gettextRcmdr(current.frame$title)),sticky="nw")
					
				}
			
				
				# Make frames inside entry frame for each argument
				for(j in 1:number.entries){
					
					current.frame$entry.frames[[j]] <- tkframe(current.frame$frame)
					current.frame$entry.vars[[j]] <- tclVar(paste(initial.values[j]))
					current.frame$entry.fields[[j]] <- ttkentry(current.frame$entry.frames[[j]],width=current.frame$entry.width[j],textvariable=current.frame$entry.vars[[j]])
					
					
					
					tkgrid(labelRcmdr(current.frame$entry.frames[[j]],text=gettextRcmdr(paste(argument.names[j],": ",sep=""))),current.frame$entry.fields[[j]],sticky="nw")
					tkgrid(current.frame$entry.frames[[j]],sticky="ne")
					
					
				}
				
				new.frames[[Tab]][[ii]] <- current.frame
				
			}
			
			###### RADIO BUTTONS ##########################################################################################
			if(current.frame$type=="radiobuttons"){
				
				
				# Setting up the radio buttons
				current.frame$argument.values <- sapply(current.frame$argument.values,FUN=function(x){return(paste("BUTTONSTART",x,sep=""))})
				current.frame$initial.values <- paste("BUTTONSTART",current.frame$initial.values,sep="")
				
				radioButtons(current.frame$frame,name=frame.name,buttons=current.frame$argument.values,values=current.frame$argument.values,labels=gettextRcmdr(current.frame$argument.names),initialValue=current.frame$initial.values,title="")
				
				# Make title when NOT using a border (special case of 'title = "a title"' & 'border=FALSE')
				if(current.frame$title!="" & current.frame$border==FALSE){
					tkgrid(labelRcmdr( current.frame$frame,fg=getRcmdr("title.color"),font="RcmdrTitleFont"   ,text=gettextRcmdr(current.frame$title)),sticky="nw")
					
				}
				
				# Saving the tclVar variable of the radioframe
				eval(parse(text=paste("current.frame$radioVar <-",frame.name,"Variable",sep="" )))
				
				# Putting down the radio button frame (which is generated automatically by radioButtons)
				radio.command <- paste("tkgrid(",frame.name,"Frame,sticky='nw')",sep="")
				.eval.command(radio.command)
				
								
				new.frames[[Tab]][[ii]] <- current.frame
				
				
				
			}
			
			###### CHECK BOXES ##########################################################################################
			
			if(current.frame$type=="checkboxes"){
								
				# Setting up Check Boxes
				checkBoxes(current.frame$frame,frame=frame.name,boxes=paste(current.frame$arguments),initialValues=current.frame$initial.values,labels=sapply(current.frame$argument.names,FUN=gettextRcmdr))
				
				# Make title when NOT using a border (special case of 'title = "a title"' & 'border=FALSE')
				if(current.frame$title!="" & current.frame$border==FALSE){
					tkgrid(labelRcmdr( current.frame$frame,fg=getRcmdr("title.color"),font="RcmdrTitleFont"   ,text=gettextRcmdr(current.frame$title)),sticky="nw")
					
				}
				
				# Saving the tclVar variables of the checkbox frame
				current.frame$checkVar <- list()
				arguments <- current.frame$arguments
				
				for(j in 1:length(arguments)){
					temp.arg <- arguments[j]
					eval(parse(text=paste("current.frame$checkVar[[",j,"]] <- ",temp.arg,"Variable"  ,sep="")))
				}
				
				# Putting down the check button frame (which is generated automatically by checkBoxes)	
				
				check.command <- paste("tkgrid(",frame.name,",sticky='nw')",sep="")
				.eval.command(check.command)
				
				new.frames[[Tab]][[ii]] <- current.frame
				
				
			}
			###### VALUE SLIDERS ##########################################################################################
			
			if(current.frame$type=="valuesliders"){
								
				number.sliders <- length(current.frame$arguments)
				arguments <- current.frame$arguments
				argument.names <- current.frame$argument.names
				initial.values <- current.frame$initial.values
				
				# Making the element which will contain the separate slider frames
				
				current.frame$slider.frames <- list()
				current.frame$slider.vars <- list()
				current.frame$slider <- list()
				
				# Make title when NOT using a border (special case of 'title = "a title"' & 'border=FALSE')
				
				if(current.frame$title!="" & current.frame$border==FALSE){
					tkgrid(labelRcmdr( current.frame$frame,fg=getRcmdr("title.color"),font="RcmdrTitleFont"   ,text=gettextRcmdr(current.frame$title)),sticky="nw")
					
				}
				
				# Make frames inside valueslider frame for each slider
				
				for( j in 1:number.sliders){
					
					current.frame$slider.frames[[j]] <- tkframe(current.frame$frame)
					current.frame$slider.vars[[j]] <- tclVar(as.character(initial.values[j]))
					current.frame$slider[[j]] <- tkscale(current.frame$slider.frames[[j]],variable=current.frame$slider.vars[[j]],showvalue=TRUE,from=current.frame$from[j],to=current.frame$to[j],length=current.frame$length[j],resolution=current.frame$by[j],orient="horizontal")
					
					tkgrid(labelRcmdr(current.frame$slider.frames[[j]],text=gettextRcmdr(current.frame$argument.names[j])), current.frame$slider[[j]] ,sticky="sw")
					tkgrid(current.frame$slider.frames[[j]],sticky="nw")
					
					
				}
				
				new.frames[[Tab]][[ii]] <- current.frame
			}
			
			###### SPIN BOXES ##########################################################################################
			
			if(current.frame$type=="spinboxes"){
				
				number.spins <- length(current.frame$arguments)
				arguments <- current.frame$arguments
				argument.names <- current.frame$argument.names
				initial.values <- current.frame$initial.values
				
				# Making the element which will contain the separate spinbox frames
				
				current.frame$spin.frames <- list()
				current.frame$spin.vars <- list()
				current.frame$spin <- list()
				
				# Make title when NOT using a border (special case of 'title = "a title"' & 'border=FALSE')
				
				if(current.frame$title!="" & current.frame$border==FALSE){
					tkgrid(labelRcmdr( current.frame$frame,fg=getRcmdr("title.color"),font="RcmdrTitleFont"   ,text=gettextRcmdr(current.frame$title)),sticky="nw")
					
				}
				
				# Make frames inside spinboxes frame for each spinbox
				
				for(j in 1:number.spins){
					
					current.frame$spin.frames[[j]] <- tkframe(current.frame$frame)
					current.frame$spin.vars[[j]] <- tclVar(as.character(initial.values[j]))
					current.frame$spin[[j]] <- tkspinbox(current.frame$spin.frames[[j]],from=current.frame$from[j],to=current.frame$to[j],width=current.frame$entry.width,textvariable=current.frame$spin.vars[[j]] ,state="readonly",increment=current.frame$by[j])
					
					tkgrid(labelRcmdr(current.frame$spin.frames[[j]],text=gettextRcmdr(current.frame$argument.names[j])),current.frame$spin[[j]],sticky='nw')
					tkgrid(current.frame$spin.frames[[j]],sticky='nw')	
					
				}
				
				new.frames[[Tab]][[ii]] <- current.frame
				
			}
			
			######    LIST BOX    ########################################################################################
			
			if(current.frame$type=="listbox"){
				
				# Make title when NOT using a border (special case of 'title = "a title"' & 'border=FALSE')
				
				if(current.frame$title!="" & current.frame$border==FALSE){
					tkgrid(labelRcmdr( current.frame$frame,fg=getRcmdr("title.color"),font="RcmdrTitleFont"   ,text=gettextRcmdr(current.frame$title)),sticky="nw")
					
				}
				
				
				
				selectmode <- ifelse(current.frame$select.multiple,"multiple","single")
				height <- current.frame$length
								
				
				current.frame$listboxFrame <- tkframe(current.frame$frame)
				
				
				current.frame$listBox <- tklistbox(current.frame$listboxFrame,height=height,exportselection="FALSE",selectmode=selectmode,background="white")
				for(listname in current.frame$argument.names) tkinsert(current.frame$listBox,"end",listname)
				
				
				eval(parse(text=paste0("current.frame$listScroll <- ttkscrollbar(current.frame$listboxFrame,command=function(...){tkyview(new.frames[[",Tab,"]][[",ii,"]]$listBox, ...)})")))
				eval(parse(text=paste0("tkconfigure(current.frame$listBox, yscrollcommand=function(...){tkset(new.frames[[",Tab,"]][[",ii,"]]$listScroll, ...)})")))
				
				if(length(current.frame$argument.names)!=0){tkselection.set(current.frame$listBox,0)}
				
				tkgrid(current.frame$listBox,current.frame$listScroll)
				tkgrid.configure(current.frame$listScroll,sticky="ns")
				
				tkgrid(current.frame$listboxFrame,sticky="nw")
				
				if(current.frame$border==TRUE){tkgrid.configure(current.frame$listboxFrame,padx="5",pady="2")}
	
				new.frames[[Tab]][[ii]] <- current.frame
				
			}
			
			
			
			###### MANUAL BUTTONS ########################################################################################
			
			if(current.frame$type=="buttons"){

	
				button_result <- current.frame$button.object
				
				
				function.command <- paste(current.frame$button.function,"(",sep="")
				first.arg=TRUE
				
				if(current.frame$button.data!=""){
					
					input.data <- ActiveDataSet()		
					if(current.frame$button.data.transf=="matrix"){input.data <- paste("as.matrix(",input.data,")",sep="")}
					if(current.frame$button.data.transf=="ExprSet"){input.data <- paste0("as.ExprSet(",input.data,")")}
					
					
					function.command <- paste(function.command,current.frame$button.data,"=",input.data  ,sep="")
					first.arg=FALSE
				}
				
				
				
#				if(current.frame$button.data==""  & current.frame$button.otherarg==""){
#					function.command <- paste(function.command,"...",sep="")
#				}
				
				if(current.frame$button.otherarg!=""){
					# NOTE: AT THE MOMENT THE USER NEEDS TO DECIDE ITSELF IF A ',' IS NECESSARY, AUTOMATE?
					# NOTE2: You should always use a ','
					function.command <- paste(function.command,current.frame$button.otherarg,sep="")
					
				}
				
				arg.names <- .transform.vector2text(current.frame$arg.frames)

				
				save <- current.frame$save
				show <- current.frame$show
				show.save <- current.frame$show.save
				
				
				temp.command <- paste("function(){
								
								function.command <- .build.button.function(\"",function.command,"\",",arg.names,",\"",button_result,"\",new.frames,",save,",",Tab,")
								
								if(",show,"==TRUE){
								doItAndPrint(function.command)
								}
								
								if(",show,"!=TRUE){
								justDoIt(function.command)
								}
								if(",save,"==TRUE & ",show.save,"==TRUE){
								doItAndPrint('",button_result,"')
								}


								}",sep="")
				#if(!is.null(dev.list())){par(mfrow=c(1,1))}  # Code to reset the graphics window to c(1,1). Should be above after the 3rd if
				

				
				eval(parse(text=paste("button.command <- ",temp.command,sep="")))
				
				current.frame$button.command <- button.command
				
				current.frame$buttonRcmdr <- buttonRcmdr(current.frame$frame,command=current.frame$button.command,text=gettextRcmdr(current.frame$button.name),foreground="darkgreen",default="active",width=current.frame$button.width,borderwidth=3)
				
				tkgrid(current.frame$buttonRcmdr,sticky="s")
				
				new.frames[[Tab]][[ii]] <- current.frame
				
			}
		}
		
		########################################################################
		## Configuring of the frames (normal & combined rows) in CLUSTERFRAME ##
		########################################################################
		

		
		if(length(grid.rows[[Tab]])!=0){	special.rows <- lapply(grid.rows[[Tab]],FUN=function(x){return(x$rows)})
			special.rows.vector <- sort(unlist(special.rows))} else {special.rows.vector <- c()}
		row.index <- 1
		
		while(row.index <= dim(grid.config[[Tab]])[1]){
			
			## Placing 'rows' in grid####
			
			if((length(grid.rows[[Tab]])!=0) & (row.index %in% special.rows.vector) ){
				

				
				rowframe.index <- which(lapply(special.rows,FUN=function(x){return(row.index %in% x)}) ==TRUE)
				
				
				
				for(i in special.rows[[rowframe.index]]){
					grid.command <- paste("tkgrid(")
					for(column.index in 1:dim(grid.config[[Tab]])[2]){
						if(is.na(grid.config[[Tab]][i,column.index])){} else {grid.command <- paste(grid.command,"new.frames[[",Tab,"]][[",.find.frame(new.frames[[Tab]],grid.config[[Tab]][i,column.index]) ,"]]$frame,",sep="")}
						
					}
					grid.command <- paste(grid.command,"sticky='nw',padx=6,pady=6)")

					.eval.command(grid.command)
				}			
				
				grid.command <- paste("tkgrid(Tab",Tab,"Frame_row",rowframe.index,",sticky='nw')",sep="")

				
				.eval.command(grid.command)

				
				row.index <- i+1
				
			}
			
			
			###########################
			
			else{	

				## Normal grid matrix####
				grid.command <- "tkgrid("
				for(column.index in 1:dim(grid.config[[Tab]])[2]){
					if(is.na(grid.config[[Tab]][row.index,column.index])){} else {grid.command <- paste(grid.command,"new.frames[[",Tab,"]][[",.find.frame(new.frames[[Tab]],grid.config[[Tab]][row.index,column.index]) ,"]]$frame,",sep="")}
					
				}
				
				grid.command <- paste(grid.command,"sticky='nw',padx=6,pady=6)",sep="")

				.eval.command(grid.command)

				row.index <- row.index + 1
				#############################	
			}
		}
		

		
	}
	}
	
	
	
	
	

	############################################
	## Making ONOK, ONCANCEL, ONHELP,... funtions ##
	############################################
	
	onOK <- function(){}
	
	onCancel <- function() {
		if (GrabFocus()) 
			tkgrab.release(top)
		tkdestroy(top)
		tkfocus(CommanderWindow())
	}
	
	onHelp <- function() {
		tkgrab.release(window)
		print(help(helppage))
	}
	
	onSetwd <- function(){
		.Setwd()
	}
	
	onReset <- function(){
		doItAndPrint("rm(list=ls(,envir=.GlobalEnv))")
	}
	
	onSeed <- function(){
		doItAndPrint(paste0("set.seed(",tclvalue(seed_vars),")"))
	}
	
	##########################################################################################################################################################
	##########################################################################################################################################################
	
	
	###################################################################
	## MAKING THE BUTTONS WHICH WILL APPEAR BY STANDARD IN ALL FRAMES##
	###################################################################
	

		optional.buttons.vector <- c("helpButton","setButton","resetButton")
		chosen.buttons <- paste(optional.buttons.vector[c(make.help.button,make.setwd.button,make.resetgws.button)],collapse=",")

		if(usetabs){
			buttonsFrame <- tkframe(top)
		}
		else{
			buttonsFrame <- tkframe(Tab1Frame)
		}
	
		optional.buttons <- tkframe(buttonsFrame)
		exitButton <- buttonRcmdr(buttonsFrame, text = gettextRcmdr("Exit"), foreground = "red", width = "12", command = onCancel, borderwidth = 3)
		
		helpButton <- buttonRcmdr(optional.buttons, text=gettextRcmdr("Help"),foreground="red",width="12",command=onHelp,borderwidth=3)
		setButton <- buttonRcmdr(optional.buttons, text=gettextRcmdr("Set Working Directory"),foreground="red",width="20",command=onSetwd,borderwidth=3)
		resetButton <- buttonRcmdr(optional.buttons, text=gettextRcmdr("Reset Global Workspace"),foreground="red",width="22",command=onReset,borderwidth=3)
		
		
		if(sum(make.help.button,make.setwd.button,make.resetgws.button)>0){
			tk_command <- paste("tkgrid(",chosen.buttons,")"    ,sep="")
			.eval.command(tk_command)
		}
		
		tkgrid(optional.buttons,exitButton,pady=8)
		
		if(make.seed.button){
			seedboxFrame <- tkframe(buttonsFrame)
			seed_entry <- tkframe(seedboxFrame)
			seed_vars <- tclVar(paste0(round(runif(1,0,999))))
			seed_field <- ttkentry(seed_entry,width=4,textvariable=seed_vars)
			tkgrid(labelRcmdr(seed_entry,text=gettextRcmdr("Seed: ")),seed_field,sticky="nw")
			
			seedboxButton <- buttonRcmdr(seedboxFrame, text=gettextRcmdr("Set"),foreground="red",width="5",command=onSeed,borderwidth=3)
			
			
			tkgrid(seed_entry,seedboxButton,sticky="nw")
			tkgrid(seedboxFrame)
			
			tkgrid.configure(seedboxButton,padx="6")
			tkgrid.configure(seedboxFrame,sticky="w")
		}
		
		
		tkgrid.columnconfigure(buttonsFrame, 0, weight=1)
		tkgrid.columnconfigure(buttonsFrame, 1, weight=1)
		tkgrid.configure(optional.buttons, sticky="w")
		tkgrid.configure(exitButton, sticky="e")
		
		if(!usetabs){
			tkgrid(buttonsFrame,sticky="swe",pady="10")
		}


		
		
	##########################################################################################################################################################
	## Making the final tab grids ##
	################################
		
		for(Tab in 1:n.tabs){
			tk_command <- paste("tkgrid(Tab",Tab,"Frame,padx=5,pady=5,sticky='nw')",sep="")
			.eval.command(tk_command)
		}
		
		

	if(usetabs==TRUE){
		dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE,onOK=onOK,tabs=tabs,tab.names=tabnames,preventGrabFocus=TRUE)
		
	}
	else{
		dialogSuffix(use.tabs=FALSE,onOK=onOK,preventGrabFocus=TRUE)
	}
}




