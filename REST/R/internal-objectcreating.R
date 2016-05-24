

.add.frame <- function(Tab,type,frame.name,argument.names="",arguments="",initial.values=c(),title="",border=FALSE,entry.width="2",argument.values=c(),argument.types=c(),from=c(),to=c(),by=c(),length=c(),select.multiple=FALSE,button.name="",button.function="",button.data="",arg.frames=c(),button.otherarg="",button.object="",button.width="12",button.data.transf="matrix",save=TRUE,show.save=TRUE,show=TRUE ,new.frames=new.frames){
	
	
	# Entry Fields
	if(type=="entryfields"){
		
		new <-  list(type=type,frame.name=frame.name,argument.names=argument.names,arguments=arguments,initial.values=initial.values,title=title,border=border,entry.width=entry.width,argument.types=argument.types)
		
		new.frames[[Tab]][[length(new.frames[[Tab]])+1]] <- new
		new.frames <- .order.button.frames(new.frames,Tab)
		return(new.frames)

	   
	}
	
	# Radio buttons
	if(type=="radiobuttons"){
		
		new <-  list(type=type,frame.name=frame.name,argument.names=argument.names,arguments=arguments,initial.values=initial.values,title=title,border=border,argument.values=argument.values,argument.types=argument.types)
		
		new.frames[[Tab]][[length(new.frames[[Tab]])+1]] <- new
		new.frames <- .order.button.frames(new.frames,Tab)  # Make sure the manual button frames are the last ones in the list .
		return(new.frames)

	}
	
	# Check Boxes
	if(type=="checkboxes"){
		new <- list(type=type,frame.name=frame.name,argument.names=argument.names,arguments=arguments,initial.values=initial.values,title=title,border=border)
		
		new.frames[[Tab]][[length(new.frames[[Tab]])+1]] <- new
		new.frames <- .order.button.frames(new.frames,Tab)
		return(new.frames)

		
	}
	
	# Slider Values
	if(type=="valuesliders"){
		new <- list(type=type,title=title,border=border,frame.name=frame.name,argument.names=argument.names,arguments=arguments,initial.values=initial.values,from=from,to=to,by=by,length=length)
		
		new.frames[[Tab]][[length(new.frames[[Tab]])+1]] <- new
		new.frames <- .order.button.frames(new.frames,Tab)
		return(new.frames)

	}
	
	# Spin Boxes
	if(type=="spinboxes"){
		new <- list(type=type,frame.name=frame.name,argument.names=argument.names,arguments=arguments,initial.values=initial.values,from=from,to=to,by=by,entry.width=entry.width,title=title,border=border)
		
		new.frames[[Tab]][[length(new.frames[[Tab]])+1]] <- new
		new.frames <- .order.button.frames(new.frames,Tab)
		return(new.frames)

	}
	
	# List Boxes
	if(type=="listbox"){
		new <- list(type=type,frame.name=frame.name,argument.names=argument.names,arguments=arguments,argument.values=argument.values,argument.types=argument.types,initial.values=initial.values,length=length,select.multiple=select.multiple,title=title,border=border)
		
		new.frames[[Tab]][[length(new.frames[[Tab]])+1]] <- new
		new.frames <- .order.button.frames(new.frames,Tab)
		return(new.frames)

	}
	
	# Manual Buttons
	if(type=="buttons"){
		new <- list(frame.name=frame.name,type=type,button.name=button.name,button.function=button.function,button.data=button.data,button.otherarg=button.otherarg,arg.frames=arg.frames,button.object=button.object,title="",border=FALSE,save=save,show=show,show.save=show.save,button.width=button.width,button.data.transf=button.data.transf)
		new.frames[[Tab]][[length(new.frames[[Tab]])+1]] <- new
			
		new.frames <- .order.button.frames(new.frames,Tab)
		return(new.frames)
	}
	
	
}



.order.button.frames <- function(new.frames,Tab){
	
	boolean <- sapply(new.frames[[Tab]],FUN=function(x){x$type=="buttons"})
	new.frames[[Tab]] <- new.frames[[Tab]][order(boolean)]
	
	return(new.frames)
}


#.as.var <- function(x){ return(eval.parent( as.name(x)  ,n=1))}



.find.frame <- function(x,frame.name){
	temp.names <- (lapply(x,FUN=function(d){return(d$frame.name)}))
	
	find.boolean <- temp.names==frame.name
	
	if(sum(find.boolean)==1){
		return(which(find.boolean))
	}
	
	else if(sum(find.boolean)>1){ 
		stop(paste("'",frame.name,"' is used for multiple frames!",sep=""),call.=FALSE)
	}
	else{
		stop(paste("'",frame.name,"' is not recognised as a framename. Check your 'grid.config' matrix.",sep=""),call.=FALSE)

	}
	
	
}


.eval.command <- function(x){return(eval.parent(parse(text=x),n=1))}


.combine.rows <- function(Tab,rows,title,border,grid.rows,grid.config){
	
	all.grid.rows <- grid.rows
	all.grid.config <- grid.config
	
	
	eval(parse(text=paste("grid.rows <- all.grid.rows[[",Tab,"]]",sep="")))
	eval(parse(text=paste("grid.config <- all.grid.config[[",Tab,"]]",sep="")))
	
	# The names of the frames involved in this combined row are extracted. This information is needed in the template function.
	name.frames <- as.vector(grid.config[rows,])
	name.frames <- name.frames[!is.na(name.frames)]
	
	new <- list(rows=rows,title=title,border=border,name.frames=name.frames)
	grid.rows[[length(grid.rows)+1]] <- new
	
	# In order to keep the grid.rows object correct. Sort the list, based on the rows inside an element:
	# This will ensure making a correct grid, even if the rows were combined like first 3 & 4, then 1 & 2
	
	grid.rows <- grid.rows[order( unlist(lapply(grid.rows,FUN=function(x){return(min(x$rows))})) )]
	
	
	eval(parse(text=paste("all.grid.rows[[",Tab,"]] <- grid.rows",sep="")))
	
	return(all.grid.rows)
}


.make.correct.frame <- function(title,border,window){
	
	if(title!="" & border==TRUE){
		
		return(ttklabelframe(window,text=gettextRcmdr(title)))
		
	}
	else{
		
		if(border==TRUE){relief<-"groove"} else {relief <- "flat"}
		
		return(tkframe(window,relief=relief,borderwidth=2))
		
	}
	
}



.initialize.new.frames <- function(ntabs){
	new.frames <- list()
	
	for(i in 1:ntabs){
		new.frames[[i]] <- list()
	}
	return(new.frames)
}

.initialize.grid.config <- function(ntabs){
	grid.config <- list()
	for(i in 1:ntabs){
		grid.config[[i]] <- list()
	}

	return(grid.config)
	
}

.initialize.grid.rows <- function(ntabs){
	grid.rows <- list()
	for(i in 1:ntabs){
		grid.rows[[i]] <- list()
	}

	return(grid.rows)
	
}

.grid.matrix <- function(Tab,data,grid.config=grid.config,...){
	temp <- matrix(data=data,...)
	
	grid.config[[Tab]] <- temp
	return(grid.config)
	
}


.build.command.argument <- function(current.frame,command){
	
	if(current.frame$type=="entryfields"){
		
		number.entries <- length(current.frame$arguments)
		arguments <- current.frame$arguments
		
		for(j in 1:number.entries){
			
			if(current.frame$argument.types[j]=="num"){
				add.command <- if(tclvalue(current.frame$entry.vars[[j]])==""){""} else {paste(",",arguments[j],"=",tclvalue(current.frame$entry.vars[[j]]),sep="")}
			}
			if(current.frame$argument.types[j]=="char"){
				add.command <- if(tclvalue(current.frame$entry.vars[[j]])==""){""} else {paste(",",arguments[j],"='",tclvalue(current.frame$entry.vars[[j]]),"'",sep="")}
				
			}
			
			command <- paste(command,add.command,sep="")
			
		}
		return(command)
	}
	
	if(current.frame$type=="radiobuttons"){
		
				
		temp <- (tclvalue(current.frame$radioVar))
		
		if(grepl("BUTTONSTART",temp,fixed=TRUE)){
			temp <- gsub("BUTTONSTART","",temp,fixed=TRUE)
		}
		
		if(current.frame$argument.types=="char"){
			add.command <- paste( ",",current.frame$arguments,"='",temp,"'",sep=""   )
		}
		if(current.frame$argument.types=="num"){
			add.command <- paste( ",",current.frame$arguments,"=",temp,sep=""   )
		}	
			
		command <- paste(command,add.command,sep="")
		return(command)
	}
	
	if(current.frame$type=="checkboxes"){
		
		number.checks <- length(current.frame$arguments)
		arguments <- current.frame$arguments
		
		for(j in 1:number.checks){
			
#			temp.command <- paste("temp.var <- as.character(tclvalue(",arguments[[j]],"Variable))" ,sep="")
#			.eval.command(temp.command)
			temp.var <- as.character(tclvalue(current.frame$checkVar[[j]]))
			
			if(temp.var=="1"){check.var <- TRUE} else {check.var <- FALSE}
			
			add.command <- paste(",",arguments[j],"=",check.var,sep="")
			
			command <- paste(command,add.command,sep="")
			
		}
		return(command)
	}
	
	if(current.frame$type=="valuesliders"){
		number.sliders <- length(current.frame$arguments)
		arguments <- current.frame$arguments
		
		for(j in 1:number.sliders){
			
			add.command <- paste(",",arguments[j],"=",tclvalue(current.frame$slider.vars[[j]]),sep="")
			
			command <- paste(command,add.command,sep="")
			
		}
		return(command)
		
	}
	
	if(current.frame$type=="spinboxes"){
		number.spins <- length(current.frame$arguments)
		arguments <- current.frame$arguments
		
		for(j in 1:number.spins){
			
			add.command <- paste(",",arguments[j],"=",tclvalue(current.frame$spin.vars[[j]]),sep="")
			command <- paste(command,add.command,sep="")
			
		}
		return(command)
	}
	
	if(current.frame$type=="listbox"){
		
		sel <- as.integer(tkcurselection(current.frame$listBox))+1
		
		command <- paste0(command,",",current.frame$arguments,"=")
		
		
		if(current.frame$select.multiple==TRUE){
			add.command <- paste0("c(")
			if(current.frame$argument.types=="char"){chartype <- TRUE}else{chartype <- FALSE}
			
			for(i.sel in sel){
				if(chartype){
					add.command <- paste0(add.command,"'",current.frame$argument.values[i.sel],"'")
				}
				else{
					add.command <- paste0(add.command,current.frame$argument.values[i.sel])
				}
				
				if(!(i.sel==sel[length(sel)])){add.command <- paste0(add.command,",")}
			}
			add.command <- paste0(add.command,")")
			
		}
		else{
			if(current.frame$argument.types=="char"){
				add.command <- paste0("'",current.frame$argument.values[sel],"'")
			}
			else{
				add.command <- paste0(current.frame$argument.values[sel])
			}
		}
		
		command <- paste0(command,add.command)	
	}
	
	
}



.transform.vector2text <- function(x){
	if(length(x)==0){return("c()")}
	
	out <- "c("
	for(i.arg in 1:length(x)){
		out <- paste(out,"'",x[i.arg],"'",sep="")
		if(i.arg!=length(x)){out <- paste(out,",",sep="")}
	}
	out <- paste(out,")",sep="")
	return(out)
}



.build.button.function <- function(function.command,arg.names,button_result,new.frames,save,Tab){ 
	
	
	for(i.frame in arg.names){
	
		boolean <- sapply(new.frames[[Tab]],FUN=function(x){x$frame.name==i.frame})
		temp.index <- which(boolean==TRUE)
		if(sum(boolean) ==0){stop(paste("'",i.frame,"' is not defined in 'new.frames' object",sep=""))}
		if(sum(boolean)>1){stop(paste("'",i.frame,"' is defined multiple times in 'new.frames' object",sep=""))}
		if(sum(boolean)==1){
			current.arg.frame <- new.frames[[Tab]][[temp.index]]
		
			function.command <- .build.command.argument(current.arg.frame,function.command)
		}
	}

	function.command <- paste(function.command,")" ,sep="")
	function.command <- gsub("\\(,","\\(",function.command) # Fixing the case when no data or otherarg is used (function(,x=1))
	
	if(save==TRUE){function.command <- paste(button_result," <- ",function.command,sep="")}
	
	return(function.command)
}


.give.doublequote <- function(x){return(paste("\"",x,"\"",sep=""))}


.Setwd <- function (x=TRUE) 
{
	wd <- tclvalue(tkchooseDirectory(initialdir = getwd(), parent = CommanderWindow()))
	if (wd != "") 
		doItAndPrint(paste("setwd(\"", wd, "\")", sep = ""))
}





.rcmdr.warning <- function(x){
	warning.command <- paste("warning('",x,"',call.=FALSE)",sep="")
	justDoIt(warning.command)
}

SaveGUI <- function(object.names,init.name="result"){
	#object.names = vector of names of the objects which should be saved
		
	init.name <- paste0(init.name,".RData")
	
	filenameLoc <- tclvalue(tkgetSaveFile(initialfile=init.name,filetypes="{{RData Files} {.RData .rda}} {{All files} *}"))
	
	save.command <- paste0("save(list=",deparse(object.names),",file='",filenameLoc,"')")
	doItAndPrint(save.command)
}

LoadGUI <- function(){
	fileNameLoc <- tclvalue(tkgetOpenFile(filetypes="{{RData Files} {.RData .rda}} {{All files} *}")) 
	load.command <- paste0("load('",fileNameLoc,"')")
	doItAndPrint(load.command)
}


.updateEnvirObject <- function(name,ENVIR){
	
	
	if(is.null(.GetEnvREST("GUI_envir"))){
		GUI_envir <- list()
		.AssignEnvREST("GUI_envir",GUI_envir)
	}
	
	GUI_envir <- .GetEnvREST("GUI_envir")
	if(name %in% names(GUI_envir)){
		index.env <- which(name==names(GUI_envir))
		
		temp.env <- GUI_envir[[index.env]]
		rm(list=ls(temp.env),envir=temp.env)
	}
	else{
		index.env <- length(GUI_envir)+1
	}
	GUI_envir[[index.env]] <- ENVIR
	names(GUI_envir)[index.env] <- name
	
	.AssignEnvREST("GUI_envir",GUI_envir)

}


.EnvREST <- new.env()


.GetEnvREST <- function(x){
	if(!exists(x,envir=.EnvREST,inherits=FALSE)){
		return(NULL)
	}
	else{
		return(get(x=x,envir=.EnvREST,inherits=FALSE))
	}
	
}

.AssignEnvREST <- function(x,value){
	assign(x=x,value=value,envir=.EnvREST)
}


GetWindowsENVIR <- function(){
	return(.GetEnvREST("GUI_envir"))
}