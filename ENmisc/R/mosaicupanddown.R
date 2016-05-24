### file: mosaicupanddown.R

raise <- function(w) tkraise(getToolkitWidget(w))


make_modal <- function(w) {
	if(!is(w@widget, "gWindowtcltk"))
		return()
	win <- w@widget@block
	tkwait.window(win)
}


loadedPackages <- function () {
    s <- grep("^package:", search(), value = TRUE)
    return(sub("^package:", "", s))
}

CommanderWindowP <- function(){
    if (!("Rcmdr" %in% loadedPackages())) return(FALSE)
	return(!is.null(CommanderWindow()))
}


getwinstate<-function(){
	ensurewinstate()
#	get("structablepermdialogwindowstate",envir=globalenv())
	getENmisc("structablepermdialogwindowstate")
}  

setwinstate <- function(winstate){
#	assign("structablepermdialogwindowstate",winstate,envir=globalenv())
	putENmisc("structablepermdialogwindowstate",winstate)
}

setwinstateval <- function(valnam,val){
	winstate<-getwinstate()
	winstate[[valnam]]<-val
	setwinstate(winstate)
}

getwinstateval <- function(valnam){
	getwinstate()[[valnam]]
}

ensurewinstate <- function(){
#   if (!exists("structablepermdialogwindowstate",envir=globalenv()))
   if (!exists("structablepermdialogwindowstate",envir=ENmiscEnv()))
   setwinstate(list())
} 

inv.perm <- function(perm){
    unlist(Map(function(i) which(i == perm), 1:length(perm)))
}

"%upto%" <- function(startind,endind){
	if(startind <= endind) startind:endind else NULL
}

upordownable <- function(i,splitVerticalVec,keepmarginorder=TRUE,checkdown=TRUE){
	compBoolVal <- ifelse(keepmarginorder,(!splitVerticalVec[i]),splitVerticalVec[i])
	baseSeq <- 	if (checkdown) ((i+1) %upto% length(splitVerticalVec)) else
								(1 %upto% i-1)
	candidates <- Filter(function(x) (splitVerticalVec[x]== compBoolVal),baseSeq)
		(length(candidates) > 0) 
} 


upable <- function(i,splitVerticalVec,keepmarginorder=TRUE){
	upordownable(i,splitVerticalVec,keepmarginorder,FALSE)
}

downable <- function(i,splitVerticalVec,keepmarginorder=TRUE){
	upordownable(i,splitVerticalVec,keepmarginorder,TRUE)
}


upperm <- function(i,splitVerticalVec,keepmarginorder=TRUE){
  if (!(upable(i,splitVerticalVec,keepmarginorder))) return(NULL)
  if (keepmarginorder) {
	newPos <- max(Filter(function(x) (splitVerticalVec[x]==!splitVerticalVec[i]),1 %upto% i-1))
	perm<-1:length(splitVerticalVec)
	perm[newPos:(i-1)]<-(newPos+1):i
	perm[i]<-newPos
  } else {
	newPos <- max(Filter(function(x) (splitVerticalVec[x]==splitVerticalVec[i]),1 %upto% i-1))
	perm<-1:length(splitVerticalVec)
	perm[newPos] <- i
	perm[i] <- newPos
  }
  perm  
}

downperm <- function(i,splitVerticalVec,keepmarginorder=TRUE){
  if (!(downable(i,splitVerticalVec,keepmarginorder))) return(NULL)
  if (keepmarginorder) {
	newPos <- min(Filter(function(x) (splitVerticalVec[x]==!splitVerticalVec[i]),
						(i+1) %upto% length(splitVerticalVec)))
	perm<-1:length(splitVerticalVec)
	perm[(i+1):newPos]<-i:(newPos-1)
	perm[i]<-newPos
  } else {
	newPos <- min(Filter(function(x) (splitVerticalVec[x]==splitVerticalVec[i]),
						(i+1) %upto% length(splitVerticalVec)))
	perm<-1:length(splitVerticalVec)
	perm[newPos] <- i
	perm[i] <- newPos     
  }
  perm  
}

makeFuncall <- function(funexpr,...){
	myargs<-list(...)
	funcall<-paste(funexpr,"(",sep="")
	nargs<-length(myargs)
	if (is.null(names(myargs))){
		namelist<-rep("",nargs)
	} else {
		namelist<-names(myargs)
	}
	
	for (i in 1:nargs){
	    if (nchar(namelist[i])>0)
			funcall<-paste(funcall,names(myargs)[i],"=",sep="")
		funcall<-paste(funcall,myargs[[i]],sep="")
		if (i < nargs)
			funcall<-paste(funcall,", ",sep="")	
	}
	funcall<-paste(funcall,")",sep="")
	funcall
}
  
vecArg <- function(vec){
    arg <- vec
	if (is.integer(vec)) arg <- as.numeric(vec)
	capture.output(dput(arg))
}   

find.perm <- function(vec1,vec2){
	unlist(Map(function(x)which(x == vec2),vec1))
}

margin.table.structable <- function(x,margin){
	mytable<-as.table(x)
	mytable<-margin.table(mytable,margin)
	mysplit<-attr(x,"split_vertical")[margin]
	structable(mytable,split_vertical=mysplit)
}

aperm.structable <- function(a,perm,...){
	mytable <- aperm(as.table(a),perm)
	if (is.numeric(perm)) {
		splitperm<-perm
	} else {
		splitperm<-find.perm(names(dimnames(mytable)),names(dimnames(a)))		
	}
	mysplit_vertical <- attr(a,"split_vertical")[splitperm]
	structable(mytable,split_vertical=mysplit_vertical)
}

#aperm.structable.for.mosaic <- function(mystructable,perm){
#	mytable <- aperm.table(as.table(mystructable),perm)
#	if (is.numeric(perm)) {
#		splitperm<-perm
#	} else {
#	splitperm<-find.perm(names(dimnames(mytable)),names(dimnames(mystructable)))		
#	}
#	mysplit_vertical <- attr(mystructable,"split_vertical")[splitperm]
#	structable(mytable,split_vertical=mysplit_vertical)
#}


structableExpression <- function(objExprString){
	if (exists(objExprString)) {
		myobj<-get(objExprString)
		if (inherits(myobj,"structable"))
			return(objExprString)
#		if (inherits(myobj,"table") | inherits(myobj,"ftable"))
		if (isMosaicable(myobj))
			return(makeFuncall("structable",objExprString))
		if (inherits(myobj,"data.frame"))
			return(makeFuncall("structable",makeFuncall("table",makeFuncall("extractFactors",objExprString))))
		return("")	
	} else {
		if (is.structable(eval(parse(text=objExprString)))) 
			return(objExprString)
		return("")
	}	
}


errorMessage <- function(msgText){
    if (CommanderWindowP()) {
		Message(msgText,type="error")
	} else {
		gmessage(msgText,type="error")
	}
}

warnMessage <- function(msgText){
    if (CommanderWindowP()) {
		Message(msgText,type="warning")
	} else {
		gmessage(msgText,type="warning")
	}
}


isValidMosaicPermArgument <- function(tableorname){
	if (!is.character(tableorname)) {
		if (!isMosaicable(tableorname)){
#			errorMessage("mosaicPermDialog called with invalid argument")
			return(FALSE)
		} else {
			return(TRUE)
		}
	}
	res <- try(eval(parse(text=tableorname)),silent=TRUE)
	if (inherits(res,"try-error")) {
#		errorMessage("mosaicPermDialog called with invalid argument")
		return(FALSE)
	} else {
		if (!isMosaicable(res)) {
#			errorMessage("mosaicPermDialog called with invalid argument")
			return(FALSE)
		} else {
			return(TRUE)		
		}
	}		
}

mosaicableObjectAndNamestring <- function(x){
	if (is.character(x)){
		obj <- eval(parse(text=x))
		structobj <- ifelse(inherits(obj,"structable"),obj,structable(obj))
		namestring <- x
		if (inherits(obj,"structable")) {
			structobj <- obj
		} else {
			structobj <- structable(obj)
		}	
		structstring <- ifelse(inherits(obj,"structable"),x,paste("structable(",x,")",sep=""))
	} else {
		if (inherits(x,"structable")) {
			structobj <- x
		} else {
			structobj <- structable(x)
		}	
		namestring<-deparse(substitute(x))
		structstring <- ifelse(inherits(x,"structable"),namestring,paste("structable(",namestring,")",sep=""))
	}
	return(list(structobject=structobj,namestring=namestring,structstring=structstring))
}


mosaicPermDialog <- function(tablename,allow.collapsing=TRUE,start.active=NULL,
							displayPermCommand=TRUE,extendedOptions=TRUE) {	

	require(gWidgets)
    require(gWidgetstcltk)
    require(Rcmdr) 	

	tablearg<-deparse(substitute(tablename))							
	mosaicPermDialogHelper(tablearg,callingwidget=NULL,allow.collapsing,start.active,
						   displayPermCommand,extendedOptions)				  
}

mosaicPermDialogHelper <- function(tablename,callingwidget,allow.collapsing=TRUE,start.active=NULL,
							displayPermCommand=TRUE,extendedOptions=TRUE){				  
	if (missing(tablename)){
		errorMessage("mosaicPermDialog called without a required argument")
		return(invisible())		
	}	
    if (!isValidMosaicPermArgument(tablename)) {
		errorMessage("mosaicPermDialog called with invalid argument")
		return(invisible())		
	}
	
	calledfromwidget <- FALSE
	
	if (!(missing(callingwidget) || is.null(callingwidget))) {
	    calledfromwidget <- TRUE
		visible(callingwidget) <- FALSE
	}
	
	ensurewinstate() 
	doAssoc <- FALSE
#	assign("myAssoc",doAssoc,envir=globalenv())
	putENmisc("myAssoc",doAssoc)

    resultType <- "command"
#	assign("myResultType",resultType,envir=globalenv())
	putENmisc("myResultType",resultType)

    if (is.character(tablename)) {
		tmpres <- mosaicableObjectAndNamestring(tablename)
	} else {
		tmpres <- mosaicableObjectAndNamestring(deparse(substitute(tablename)))	
	}
	
	mytable <- tmpres$structobject
	mytablename <- tmpres$namestring
	mystructablename <- tmpres$structstring
	
	setwinstateval("structable",mytable)
	setwinstateval("structablename",mytablename)
	setwinstateval("structableexpr",mystructablename)

    varnames<-names(dimnames(mytable))
   
	currvarnames <- names(dimnames(mytable))
	currsplits <- attr(mytable,"split_vertical")	
	currperm <- 1:(length(varnames))
	

	setwinstateval("varnames",currvarnames)
	setwinstateval("splits",currsplits)
	setwinstateval("splitsorig",currsplits)	
	setwinstateval("perm",currperm)
	
	if (is.null(start.active)) {
		activevars<-rep(T,length(varnames))
	} else {
		activevars<-start.active
	}
	
	setwinstateval("activevars",activevars)
	setwinstateval("plotcommand","")
	setwinstateval("tablecommand","")
	
	buttons<-list()
	 
	listindex <- function(i,j) 6*(i-1)+j 	

	updatedisplay<-function(updatemosaic=TRUE){
		winstate<-getwinstate()
	
		enabled(bcolrev) <- svalue(bcol)
		if (!svalue(bcol)) {
			svalue(bcolrev)<-FALSE
		}
		
		for (i in 1:length(winstate$varnames)) {
            if (winstate$splits[i]){
				svalue(buttons[[listindex(i+1,5)]])<- winstate$varnames[i]	
				svalue(buttons[[listindex(i+1,1)]])<- ""				
			} else {	
				svalue(buttons[[listindex(i+1,5)]])<- ""				
				svalue(buttons[[listindex(i+1,1)]])<- winstate$varnames[i]	
			}
		}
		activevarsCount <- length(Filter(function(x)x,winstate$activevars))
#		doAssoc <- get("myAssoc",envir=globalenv())
		doAssoc <- getENmisc("myAssoc")
		if (doAssoc) {
			allowDeactivate <- (activevarsCount > 2)
		} else {
			allowDeactivate <- (activevarsCount > 1)
		}

		for (i in 1:length(winstate$varnames)) {				
			enabled(buttons[[listindex(i+1,2)]]) <- 
				upable(i,winstate$splits,winstate$keepmarginorder) 
			if 	(upable(i,winstate$splits,winstate$keepmarginorder))
				showButton(i+1,2) else
				hideButton(i+1,2) 
				enabled(buttons[[listindex(i+1,4)]]) <- 
				downable(i,winstate$splits,winstate$keepmarginorder)
			if 	(downable(i,winstate$splits,winstate$keepmarginorder))
				showButton(i+1,4) else
				hideButton(i+1,4) 
			if (winstate$splits[i]) {
#					svalue(buttons[[listindex(i+1,3)]])<-"To_row_vars"
					svalue(buttons[[listindex(i+1,3)]])<-"backward"
			} else {
#					svalue(buttons[[listindex(i+1,3)]])<-"To_col_vars"
					svalue(buttons[[listindex(i+1,3)]])<-"forward"
			}
			if (allow.collapsing) {
				svalue(buttons[[listindex(i+1,6)]])<- winstate$activevars[i]	
				if (allowDeactivate) {
					enabled(buttons[[listindex(i+1,6)]]) <- TRUE
				} else {
					enabled(buttons[[listindex(i+1,6)]]) <- !winstate$activevars[i]				
				}
			}	
		}
		
		if (updatemosaic){
			if (!is.null(winstate$structable)){
				plotmosaic()
			}
		}	
	}	

	
	plotmosaic <- function(){
		winstate <- getwinstate()

		
		structableexpr <- winstate$structableexpr
		
        if (all(winstate$activevars)) {
			command <- structableexpr
			if (!all(winstate$perm==(1:length(winstate$perm)))) {	
				command <- makeFuncall("aperm",
					command,
					vecArg(winstate$perm))
			}
			if (any(winstate$splitsorig[winstate$perm] != winstate$splits))
				command <- makeFuncall("structable",
					command,
					split_vertical=vecArg(winstate$splits))										
		} else {	
			varselect<- winstate$activevars[inv.perm(winstate$perm)]				
			reducedperm<-rank(winstate$perm[winstate$activevars])
			command <- makeFuncall("margin.table.structable",
				structableexpr,
				vecArg(which(varselect)))
			if (!all(reducedperm==(1:length(reducedperm)))) {	
				command <- makeFuncall("aperm",
					command,
					vecArg(reducedperm))
			}
			if (any((winstate$splitsorig[winstate$activevars])[reducedperm] != winstate$splits[winstate$activevars]))
			
			command <- makeFuncall("structable",
				command,
				split_vertical=vecArg(winstate$splits[winstate$activevars]))
		}
		
		
#		doAssoc <- get("myAssoc",envir=globalenv())
		doAssoc <- getENmisc("myAssoc")
#		resultType <- get("myResultType",envir=globalenv())
		resultType <- getENmisc("myResultType")
		graphCommand <- ifelse(doAssoc,"assoc","mosaic")
        tablecommand <- command
		
		if (svalue(bcol) & !doAssoc) {
			tmpTab <- eval(parse(text=command))
			colornum <- length(last.element(attr(tmpTab,"dnames")))
			lastInd <- length(attr(tmpTab,"dnames"))
			palettename <- getOption("ENmisc")$mosaicpalette
			if (is.null(palettename)) palettename <- "RdYlGn"
			reversepalette <- getOption("ENmisc")$mosaicpalettereverse
			if (is.null(reversepalette)) reversepalette <- FALSE
			
			if (svalue(bcolrev)) reversepalette <- !reversepalette
			
			colorString <- paste("highlighting_fill=",
								 ifelse((colornum < 3) || reversepalette ,"brewer.pal.ext","brewer.pal"),
								 "(",as.character(colornum),",'",palettename,"'",
								 ifelse(reversepalette,",reverse=TRUE",""),")",sep="")
			hiString <- paste("highlighting",as.character(lastInd),sep="=")			

			command<-makeFuncall(graphCommand,command,hiString,colorString)
			

		} else {
			command<-makeFuncall(graphCommand,command) 	
		}
		if (displayPermCommand)
			svalue(permtextfield)<-command
			
		
		setwinstateval("plotcommand",command)
		setwinstateval("tablecommand",tablecommand)

		eval(parse(text=command))
	}	
    
    permhandler<-function(butnum,dodown){z<-butnum;zz<-dodown;
		function(h,...){
			winstate<-getwinstate()
			currvarnames<-getwinstateval("varnames")
			currsplits<-getwinstateval("splits")
			currperm<-getwinstateval("perm")
			curractivevars<-getwinstateval("activevars")
			keepmarginorder<-getwinstateval("keepmarginorder")
			if (dodown) {
				perm<-downperm(butnum,currsplits,keepmarginorder) 
			} else {
				perm<-upperm(butnum,currsplits,keepmarginorder) 
			}	
			currsplits<-currsplits[perm]
			currvarnames<-currvarnames[perm]
			curractivevars<-curractivevars[perm]
			currperm<-currperm[perm]

			setwinstateval("varnames",currvarnames)
			setwinstateval("splits",currsplits)
			setwinstateval("perm",currperm)
			setwinstateval("activevars",curractivevars)
			
			updatedisplay()		
		}
   }
	
	downhandler <- function(x)permhandler(x,TRUE)
	uphandler <-   function(x)permhandler(x,FALSE)

	rchandler<-function(butnum){z<-butnum;
		function(h,...){
			currsplits<-getwinstateval("splits")
			currsplits[butnum]<-!currsplits[butnum]
			setwinstateval("splits",currsplits)
			updatedisplay()
		}
	}			

	checkbuttonhandler <- function(h,...) {
		ensurewinstate()
		winstate<-getwinstate()
     	activevars <- unlist(Map(function(i)svalue(buttons[[listindex(i+1,6)]]),1:length(varnames)))
		winstate$activevars <- activevars
		setwinstate(winstate)
		updatedisplay()
	}

	keepMarginHandler <- function(h,...){
		setwinstateval("keepmarginorder",svalue(h$obj,index=TRUE)==1)
		updatedisplay()				
	}
	  
	hideButton <- function(i,j){
		delete(tbl[i,j],buttons[[listindex(i,j)]])
	}
	 
	showButton <- function(i,j){
		add(tbl[i,j],buttons[[listindex(i,j)]])
	}
	
	setLastPlotCommand <- function(){
	}
	
	
	colorbuttonhandler <- function(h,...){
		updatedisplay()
	}
	
	plottypebuttonhandler <- function(h,...){
		if (svalue(h$obj) == "assoc plot") {
#			assign("myAssoc",TRUE,envir=globalenv())
			putENmisc("myAssoc",TRUE)
			svalue(bcol)<-FALSE
			enabled(bcol)<-FALSE
		} else {
#			assign("myAssoc",FALSE,envir=globalenv())
			putENmisc("myAssoc",FALSE)
			enabled(bcol)<-TRUE
		}
		updatedisplay()
	}
	
	OKbuttonhandler <- function(h,...){
		if (CommanderWindowP()){ 
#			resultType <- get("myResultType",envir=globalenv())
			resultType <- getENmisc("myResultType")
			if (resultType == "structable"){
				doItAndPrint(getwinstateval("tablecommand"))
			} else {	
				doItAndPrint(getwinstateval("plotcommand"))
			}	
		}	
	}
	
	resulttypebuttonhandler <- function(h,...){
		if (svalue(h$obj) == "return plot command") {
			resultType <- 'command'
		} else {
			resultType <- 'structable'
		}
#	   assign("myResultType",resultType,envir=globalenv())
	   putENmisc("myResultType",resultType)
	}
	
	
	if (doAssoc) {
		mytitle <- "Reorganize structable for assoc plot"
	} else {
		mytitle <- "Reorganize structable for mosaic plot"
	}
	win <- gbasicdialog(title=mytitle,handler=OKbuttonhandler)	
	tgroup <- ggroup(container=win,expand=TRUE,fill="x")
	addSpring(tgroup)
	tbl<-glayout(container=tgroup,spacing=0)
	addSpring(tgroup)
	tbl[1,1] <- (g <- ggroup(container=tbl,spacing=1))
	buttons[[listindex(1,1)]] <- (gl <- glabel("Row vars",container=g))
	tbl[1,5] <- (g <- ggroup(container=tbl))
	buttons[[listindex(1,5)]] <- (gl <- glabel("Col vars",container=g))
	if (allow.collapsing){
		tbl[1,6] <- (g <- ggroup(container=tbl))
		buttons[[listindex(1,6)]] <- (gl <- glabel("Active vars",container=g))
	}
	
	for (i in 1:length(varnames)){
		tbl[i+1,1] <- (gg <- ggroup(container=tbl,spacing=1))
		buttons[[listindex(i+1,1)]]<- (gl <- glabel(varnames[i],container=gg))
		tbl[i+1,5] <- (gg <- ggroup(container=tbl,spacing=1))
		buttons[[listindex(i+1,5)]]<- (gl <- glabel(varnames[i],container=gg))
		tbl[i+1,2] <- (gg <- ggroup(container=tbl,spacing=1))
		buttons[[listindex(i+1,2)]]<-gbutton("up",container=gg,
			handler=uphandler(i),compound="image")
		tbl[i+1,3] <- (gg <- ggroup(container=tbl,spacing=1))
#		buttons[[listindex(i+1,3)]]<-gbutton(ifelse(currsplits[i],"To_col_vars","To_row_vars"),
		buttons[[listindex(i+1,3)]]<-gbutton(ifelse(currsplits[i],"backward","forward"),
			container=gg,
			handler=rchandler(i),compound="image")
		tbl[i+1,4] <- (gg <- ggroup(container=tbl,spacing=1))
		buttons[[listindex(i+1,4)]]<-gbutton("down",container=gg,
			handler=downhandler(i),compound="image")
		if (allow.collapsing){
			tbl[i+1,6] <- (gg <- ggroup(container=tbl,spacing=1))
			addSpring(gg)
			buttons[[listindex(i+1,6)]] <- (gcb <- gcheckbox(checked=start.active[i],
				handler=checkbuttonhandler,container=gg))
			addSpring(gg)
		}
	}
#	assign("buttons",buttons,envir=globalenv())
	putENmisc("buttons",buttons)
#	assign("tbl",tbl,envir=globalenv())
	putENmisc("tbl",tbl)
	for (i in 1:length(varnames)){
		if (currsplits[i]){
			svalue(buttons[[listindex(i+1,5)]])<-""
		} else {
			svalue(buttons[[listindex(i+1,1)]])<-""
		}
	}
		
	gg<-ggroup(container=win,expand=TRUE,fill="x")
	addSpring(gg)
	gf <- gframe("Arrow button action",container=gg)
	grb<-gradio(c("Keep var order within margins","Reorder vars within margins"),container=gf)	
	addSpring(gg)
	addHandlerClicked(grb, handler=keepMarginHandler)
	addSpring(gg)
	gg<-ggroup(container=win,horizontal=FALSE)
	addSpace(gg,45)
	bcol<-gcheckbox("Colorize last variable",container=gg,handler=colorbuttonhandler)
	bcolrev<-gcheckbox("Reverse color scheme",container=gg,handler=colorbuttonhandler)

	
	if (extendedOptions){
		gg<-ggroup(container=win,horizontal=TRUE,expand=TRUE,fill="x")
        addSpring(gg)
		gf<-gframe(text='Graphics type',container=gg)
		addSpring(gg)
		grg <- gradio(c("mosaic plot","assoc plot"),container=gf)
	    addHandlerClicked(grg, handler=plottypebuttonhandler)
		addSpring(gg)
		gf<-gframe(text='Returned object',container=gg)
		grret <- gradio(c("return plot command","return structable object"),container=gf,handler=resulttypebuttonhandler)
		addSpring(gg)		
	}
	
	if (displayPermCommand) {
		gg<-ggroup(container=win,expand=TRUE,fill="x",horizontal=FALSE)
		gseparator(horizontal=TRUE,container=gg,expand=TRUE,fill="x")
		gg<-ggroup(container=win,horizontal=FALSE)
		gg<-ggroup(container=win,horizontal=FALSE)
		addSpace(gg,15)
		glabel("Command to create mosaic plot (can be copied):",container=gg)
		gg<-ggroup(container=win,horizontal=FALSE)
		addSpace(gg,10)
		permtextfield<-gedit(container=gg)
		size(permtextfield)<-c(300)
		editable(permtextfield)<-FALSE
	}	
	setwinstateval("keepmarginorder",svalue(grb,index=TRUE)==1)
 	updatedisplay()
#	raise(win)
#	focus(win,set=TRUE)
	visible(win,set=TRUE)
	plotcommand<-getwinstateval("plotcommand")
	tablecommand<-getwinstateval("tablecommand")
#    cat(paste("Final call:",plotcommand,tablecommand,sep="\n"))	

#		print(getwinstateval("tablecommand"))
	
#	assign(".lastMosaicOrAssocPlotCommandVar",plotcommand,envir=globalenv())
	putENmisc(".lastMosaicOrAssocPlotCommandVar",plotcommand)
	
#    rm(myAssoc,envir=globalenv())

    rm(list="myAssoc",envir=ENmiscEnv())

#	if (exists("structablepermdialogwindowstate",envir=globalenv()))	
#		rm("structablepermdialogwindowstate",envir=globalenv())

	if (exists("structablepermdialogwindowstate",envir=ENmiscEnv()))	
		rm(list="structablepermdialogwindowstate",envir=ENmiscEnv())

#		rm(buttons,tbl,envir=globalenv())	
		rm(buttons,tbl,envir=ENmiscEnv())	
 #   resultType <- get("myResultType",envir=globalenv())
 #   rm(myResultType,envir=globalenv())
   resultType <- getENmisc("myResultType")
   rm(list="myResultType",envir=ENmiscEnv())

	
#	print(resultType)
	
  	
    if (resultType == "structable") {
	    myResult <- eval(parse(text=tablecommand))
	} else {
		myResult <- plotcommand
	}
	if (calledfromwidget) {
#		assign("gwidgetscallresult",myResult,envir=globalenv())
		putENmisc("gwidgetscallresult",myResult)
		dispose(callingwidget)
	}	
	return(myResult)
}

##############################################################################################

hasMoreFactors <- function(mydataframe){
	factorind<-unlist(Map(function(i)is.factor(mydataframe[,i]),1:length(mydataframe)))
	length(factorind)>1
}

extractFactors <- function(mydataframe){
	factorind<-unlist(Map(function(i)is.factor(mydataframe[,i]),1:length(mydataframe)))
	mydataframe[,factorind,drop=FALSE]
}


isMosaicable <- function(x){
    mosaicableTypes <- c("table","ftable","structable")
	any(sapply(mosaicableTypes,function(z) inherits(x,z)))
}

existMosaicAbles <- function(){
	length(mosaicAbles())>0
}

mosaicAbles <- function(){
#	mytables<-Filter(function(x)any(class(get(x)) %in% c("table","ftable","structable")),ls(envir=globalenv()))
	mytables<-Filter(function(x)isMosaicable(get(x)),ls(envir=globalenv()))
	myfactordfs <- Filter(function(x)hasMoreFactors(get(x)),Filter(function(x)inherits(get(x),"data.frame"),ls(envir=globalenv())))
	c(mytables,myfactordfs)
}


mosaicSelectDialog <- function(){
    require(gWidgets)
    require(gWidgetstcltk)
    require(Rcmdr) 	
	dotable <- function(h,...){
		tablename<-svalue(h$obj)
	}
    if (!existMosaicAbles()) return(invisible(NULL))
	tlist <- mosaicAbles()
	mytitle <- paste("Select table for plot")

	mywin <- gwindow(title=mytitle,visible=FALSE,handler=function(h,...){
#		if (exists("structablepermdialogwindowstate",envir=globalenv())) 
#			rm(structablepermdialogwindowstate,envir=globalenv())
		if (exists("structablepermdialogwindowstate",envir=ENmiscEnv())) 
			rm(list="structablepermdialogwindowstate",envir=ENmiscEnv())
	})


	gg<-ggroup(horizontal=FALSE,container=mywin)
	dummy<-glabel("Select table",container=gg)
	cb<-gcombobox(tlist,selected=ifelse(length(tlist)==1,1,0),container=gg,handler=dotable)

	doOKbutton <- function(h,...){
		tablename <- svalue(cb)
		
		if (is.na(tablename)) {
			warnMessage("No table has been selected for plotting")
			dispose(mywin)
		}		
		if(!is.na(tablename)) mosaicPermDialogHelper(tablename,callingwidget=mywin)	
	}

	addSpace(gg,0,horizontal=FALSE)
	okbut <- gbutton('Create plot',container=gg,handler=doOKbutton)
	focus(cb)<-TRUE
	visible(mywin)<-TRUE
	dummy <- raise(mywin)
	make_modal(mywin)
#	if (exists("gwidgetscallresult",envir=globalenv())) {
	if (exists("gwidgetscallresult",envir=ENmiscEnv())) {
		res <- getENmisc("gwidgetscallresult")
#		rm(gwidgetscallresult,envir=globalenv())
		rm(list="gwidgetscallresult",envir=ENmiscEnv())
		return(res)
	}
}



last.element <- function(l){
   l[[length(l)]]
}

setMosaicPalette <- function(palettename,reverse=FALSE){
	require(RColorBrewer)
	if (!missing(palettename)) {
		if (!palettename %in% rownames(brewer.pal.info)) {
			warning(paste(palettename,"is not a valid name of a ColorBrewer pallette"))
			return()
		}
		enopt <- getOption("ENmisc")
		enopt$mosaicpalette <- palettename
		options(ENmisc=enopt)
	}
	enopt <- getOption("ENmisc")
	enopt$mosaicpalettereverse <- reverse
	options(ENmisc=enopt)
}

brewer.pal.ext <- function(n,name,reverse=FALSE){
   nup <- max(3,n)
   pal <- brewer.pal(nup,name)
   if (n == 1) pal <- pal[2]
   if (n == 2) pal <- pal[c(1,3)]
   if (reverse) pal  <- pal[n:1]
   return(pal)
}

