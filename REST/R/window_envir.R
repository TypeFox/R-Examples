

ChangeWindow <- function(dialogtitle,tab=1,framename,argument,new.value){

	if(!(class(argument)=="character")){stop("argument should be a string",call.=FALSE)}
	if(!(class(framename)=="character")){stop("framename should be a string",call.=FALSE)}
	
	GUI_envir <- .GetEnvREST("GUI_envir")
	
	if(!(dialogtitle %in% names(GUI_envir))){stop("Dialog Title not found",call.=FALSE)}
	index.env <- which(dialogtitle == names(GUI_envir))
	new.frames <- get("new.frames",envir=GUI_envir[[index.env]])
	new.frames.tab <- new.frames[[tab]]
	
	vector.names <- sapply(new.frames.tab,FUN=function(x){x$frame.name})
	if(!(framename %in% vector.names)){stop("Frame Name not found",call.=FALSE)}
	index.frame <- which(framename==vector.names)
	new.frames.tab.frame <- new.frames.tab[[index.frame]]
	
	# checking new.value
	if(new.frames.tab.frame$type=="listbox"){
		if(!(class(new.value)=="data.frame")){stop("new.value should be a data.frame",call.=FALSE)}
		if(dim(new.value)[2]!=2){stop("new.value should have 2 columns",call.=FALSE)}
		for(i.c in 1:dim(new.value)[2]){
			new.value[,i.c] <- as.character(new.value[,i.c])
		}
	}
	else{
		if(!(class(new.value)=="character" )){stop("new.value should be a string")}
		
	}

	
	if(new.frames.tab.frame$type =="entryfields"){
		
		if(!(argument %in% new.frames.tab.frame$arguments)){stop("Argument not found",call.=FALSE)}
		index.arg <- which(argument==new.frames.tab.frame$arguments)
		
		eval(parse(text=paste0("tclvalue(new.frames[[",tab,"]][[",index.frame,"]]$entry.vars[[",index.arg,"]]) <- \"",new.value,"\"")),envir=GUI_envir[[index.env]])
	}
	
	if(new.frames.tab.frame$type=="checkboxes"){
		if(!(argument %in% new.frames.tab.frame$arguments)){stop("Argument not found",call.=FALSE)}
		if(!(new.value %in% c("0","1"))){stop("Argument should be either '1' or '0' for check boxes",call.=FALSE)}
		
		index.arg <- which(argument==new.frames.tab.frame$arguments)
		
		eval(parse(text=paste0("tclvalue(new.frames[[",tab,"]][[",index.frame,"]]$checkVar[[",index.arg,"]]) <- \"",new.value,"\"")),envir=GUI_envir[[index.env]])
		
	}
	
	if(new.frames.tab.frame$type=="radiobuttons"){
		if(!(argument %in% new.frames.tab.frame$arguments)){stop("Argument not found",call.=FALSE)}
		
		len.arg <- length(new.frames.tab.frame$argument.values)
		new.value <- as.numeric(new.value)
		if(!(new.value %in% 1:len.arg)){stop("new.value is not correct",call.=FALSE)}
		
		eval(parse(text=paste0("tclvalue(new.frames[[",tab,"]][[",index.frame,"]]$radioVar) <- new.frames[[",tab,"]][[",index.frame,"]]$argument.values[",new.value,"]")),envir=GUI_envir[[index.env]])
		
	}
	
	if(new.frames.tab.frame$type=="valuesliders"){
		if(!(argument %in% new.frames.tab.frame$arguments)){stop("Argument not found",call.=FALSE)}
		index.arg <- which(argument==new.frames.tab.frame$arguments)
		
		new.value <- as.numeric(new.value)
		if(!((new.value >= new.frames.tab.frame$from[index.arg]) & (new.value <= new.frames.tab.frame$to[index.arg]))){stop("new.value out of bounds",call.=FALSE)}
				
		eval(parse(text=paste0("tclvalue(new.frames[[",tab,"]][[",index.frame,"]]$slider.vars[[",index.arg,"]]) <- \"",new.value,"\"")),envir=GUI_envir[[index.env]])
		
	}
	
	if(new.frames.tab.frame$type=="spinboxes"){
		if(!(argument %in% new.frames.tab.frame$arguments)){stop("Argument not found",call.=FALSE)}
		index.arg <- which(argument==new.frames.tab.frame$arguments)
		
		new.value <- as.numeric(new.value)
		if(!((new.value >= new.frames.tab.frame$from[index.arg]) & (new.value <= new.frames.tab.frame$to[index.arg]))){stop("new.value out of bounds",call.=FALSE)}
		
		eval(parse(text=paste0("tclvalue(new.frames[[",tab,"]][[",index.frame,"]]$spin.vars[[",index.arg,"]]) <- \"",new.value,"\"")),envir=GUI_envir[[index.env]])
		
	}
	
	if(new.frames.tab.frame$type=="listbox"){
		if(!(argument %in% new.frames.tab.frame$arguments)){stop("Argument not found",call.=FALSE)}
		
		#new.value input should be a dataframe with column 1 the argument.names and column 2 the argument.values
		
		eval(parse(text=paste0("new.frames[[",tab,"]][[",index.frame,"]]$argument.names <- ",.transform.vector2text(new.value[,1]))),envir=GUI_envir[[index.env]])
		eval(parse(text=paste0("new.frames[[",tab,"]][[",index.frame,"]]$argument.values <- ",.transform.vector2text(new.value[,2]))),envir=GUI_envir[[index.env]])
	
		eval(parse(text=paste0("tkdelete(new.frames[[",tab,"]][[",index.frame,"]]$listBox,'0','end')")),envir=GUI_envir[[index.env]])
		eval(parse(text=paste0("for(name in ",.transform.vector2text(new.value[,1]),") tkinsert(new.frames[[",tab,"]][[",index.frame,"]]$listBox,'end',name)")),envir=GUI_envir[[index.env]])
	
		
	}
	
	
}

CancelWindow <- function(dialogtitle){
	GUI_envir <- .GetEnvREST("GUI_envir")
	
	if(!(dialogtitle %in% names(GUI_envir))){stop("Dialog Title not found",call.=FALSE)}
	index.env <- which(dialogtitle == names(GUI_envir))
	eval(parse(text="
	if (GrabFocus()){tkgrab.release(top)}\n
	tkdestroy(top)\n
	tkfocus(CommanderWindow())
	"),envir=GUI_envir[[index.env]])
}
