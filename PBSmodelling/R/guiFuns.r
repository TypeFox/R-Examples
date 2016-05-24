##################################################
#                 PBS Modelling                  #
# -----------------------------------------------#
# This file aims to include functions specific   #
# to createWin and other GUI functions           #
#                                                #
# Authors:                                       #
#  Jon T. Schnute <schnutej-dfo@shaw.ca>,        #
#  Alex Couture-Beil <alex@mofo.ca>, and         #
#  Rowan Haigh <rowan.haigh@dfo-mpo.gc.ca>       #
#                                                #
# Hidden object code line numbers:               #
#   1151 - .create functions                     #
#   3312 - .map functions                        #
#   3459 - .convert functions                    #
#   3775 - .extract functions                    #
#   4387 - helper functions                      #
##################################################


#createWin------------------------------2012-12-17
# Create a GUI window from a given file,
# or GUI description list.
#-------------------------------------------ACB/RH
createWin <- function( fname, astext=FALSE, env=NULL )  #parent.frame() ) #globalenv() )
{
	#must be called here for examples in rd to pass check
	.initPBSoptions()
	if(is.null(env)) env = parent.frame()

	#parse window description into a valid widget tree
	if (is.character(fname)) {
		guiDesc <- parseWinFile(fname, astext=astext)
	}
	else if (is.list(fname)) {
		guiDesc <- .validateWindowDescList(fname)
		# Disabled the grid check in `.validateWindowDescWidgets` because the structure of a 
		# windows description list must have changed since this little-used function was written. (RH 12/17)
		# Note that the list result of `parseWinFile` in line 27 never goes through this validation algorithm.
		# This call only occurs if a list from `parseWinFile` (or hand-made) is passed directly into `createWin`.
		# Note: A list from `parseWinFile` never appears to create `.widgets` under `.widgets`. ***** SHOULD BE CHECKED *****
	}
	else {
		cat("ERROR, supplied argument is wrong type\n")
		return()
	}

        #uncomment the following code to test colours
        ## if (!exists ("counter")) {
        ##         counter <<- 0
        ## }
        ## counter <<- counter + 1
        ## colourset <- c("#D4D0C8", "#FF0000", "#00FF00")
        ## guiDesc[[1]]$windowname <- paste(guiDesc[[1]]$windowname, counter, sep="")
        ## guiDesc[[1]]$winBackground <- colourset[counter %% 3 + 1]
        #end of colour test code
        
	if (is.null(guiDesc)) {
		return()
	}

	#iterate over all possible windows
	for(i in 1:length(guiDesc)) {
		winName <- guiDesc[[i]]$windowname
		
		if (is.null(winName))
			stop("No window name given.")

		#destroy any existing windows with the same name
		#tt <- .PBSmod[[winName]]$tkwindow
		tget(.PBSmod)
		tt <- .PBSmod[[winName]]$tkwindow
		if (!is.null(tt))
			tkdestroy(tt)

		#clear the storage for this window
		.map.init(guiDesc[[i]]$windowname)  # alters .PBSmod

		#store windowname as most recent active window
		#eval(parse(text=".PBSmod$.activeWin <<- guiDesc[[i]]$windowname"))
		#eval(parse(text = ".PBSmodEnv$.PBSmod$.activeWin <<- guiDesc[[i]]$windowname")) 
		tget(.PBSmod)
		.PBSmod$.activeWin <- guiDesc[[i]]$windowname
		tput(.PBSmod)

		#Had to re-write tktoplevel to allow binding to the destroy
		#action without breaking the cleanup process
		mytktoplevel <- function(widget, parent = .TkRoot, ...)
		{
			command <- widget$onclose
			w <- tkwidget(parent, "toplevel", ...)
			ID <- .Tk.ID(w)
			tkbind(w, "<Destroy>", function() {

				#call the user function (if supplied)
				if (is.null(command)) {}
				else if (command=="") {}
				else if (exists(command,mode="function")) {
					.do.gui.call(command, list())
				}
				else {
					cat(paste("Warning: cannot find function '", command, "'.\n", sep=""))
				}

				#finish up with tcltk clean up (from real tktoplevel func)
				if (exists(ID, envir = parent$env, inherits = FALSE)) 
					rm(list = ID, envir = parent$env)
				tkbind(w, "<Destroy>", "")
				
				if( !is.null( widget[["remove"]] ) && widget$remove == TRUE ){
					#eval(parse(text=".PBSmod[[ widget$windowname ]] <<- NULL"))
					#eval(parse(text=".PBSmodEnv$.PBSmod[[ widget$windowname ]] <<- NULL"))
					tget(.PBSmod)
					.PBSmod[[ widget$windowname ]] <- NULL
					tput(.PBSmod)
				}
			})
			w
		}

		#set Tcl/Tk color palette with window bg and fg colour;
		#tk_setPalette will derive a number of other colours,
		#e.g., highlight, based on these colours
		#
		#note that calls to tk_setPalette normally change
		#colours on ALL existing widgets (if they use the
		#default colours), so after creating a widget, we'll
		#change its colour very slightly to ensure that the
		#colour isn't the default
		#
		#for the background, we'll change the colour when we
		#call setPalette so that we can use the proper value
		#when setting the actual background

		bgCol <- .getSimilarColour (guiDesc[[i]]$winBackground)
		fgCol <- .getSimilarColour (guiDesc[[i]]$winForeground)
		tcl("tk_setPalette", 
		    "background", bgCol,
		    "activebackground", bgCol,
		    "foreground", fgCol,
		    "activeforeground", fgCol,
		    "selectColor", "white" #inner colour of checkboxes
		    )
                
		#create TK window (blank canvas)
		tt <- mytktoplevel( guiDesc[[i]] )
                # if the bg color is the same as setPalette, future calls to setPalette
                # will change the background color
                tkconfigure (tt, bg=guiDesc[[i]]$winBackground)

		#store the TK handle (so we can destroy it at a later time via closeWin)
		#eval(parse(text=".PBSmod[[winName]]$tkwindow <<- tt"))
		#eval(parse(text=".PBSmodEnv$.PBSmod[[winName]]$tkwindow <<- tt"))
		tget(.PBSmod)
		.PBSmod[[winName]]$tkwindow <- tt

		#store environment to look for functions under
		#eval(parse(text=".PBSmod[[winName]]$env <<- env"))
		#eval(parse(text=".PBSmodEnv$.PBSmod[[winName]]$env <<- env"))
		.PBSmod[[winName]]$env <- env
		tput(.PBSmod)

		#set window title
		tkwm.title(tt,guiDesc[[i]]$title)

		#remove any old history data
		if (exists("PBS.history",envir=.PBSmodEnv)) {
			tget(PBS.history)
			j<-grep(paste("^", winName, "\\.", sep=""), names(PBS.history))
			for(n in names(PBS.history)[j]) {
				#eval(parse(text="PBS.history[[n]] <<- NULL"))
				PBS.history[[n]] <- NULL
				}
			tput(PBS.history)
		}

		#create menus
		if (length(guiDesc[[i]]$.menus) > 0) {
			#create the top space where drop down menus are attached
			topMenu <- tkmenu(tt)
			tkconfigure(tt,menu=topMenu)
                        .adjustAllColours (topMenu)

			#define function to create menu widgets
			.createMetaWidget.menu <- function(topMenu, widget)
			{
				label <- widget$label
				if (widget$nitems < 1)
					stop("menu nitems must have atleast one menuitem.")

				subMenu <- tkmenu(topMenu, tearoff=FALSE)
                                .adjustAllColours (subMenu)

				for(i in 1:widget$nitems) {
					if (widget$.widgets[[i]]$type=="menu") {
						.createMetaWidget.menu(subMenu, widget$.widgets[[i]])
						#stop("submenus need work.")
					}
					else if (widget$.widgets[[i]]$type=="menuitem")
					{
						#if eval isn't used, all callbacks call the same callback -> the last one which is defined
						#(probably since widget is redefined, but the value is only ever checked when clicked on, which is after all widgets are created)
						eval(parse(text=paste('command=function(...) .extractData("',widget$.widgets[[i]][["function"]],'", "',widget$.widgets[[i]][["action"]],'", "',winName,'")', sep="")))
						
							argList <- list( subMenu, "command", label=widget$.widgets[[i]]$label, command=command )
							if( any( widget$.widgets[[i]]$font != "" ) )
								argList$font <- .createTkFont( widget$.widgets[[i]]$font )
							if( widget$.widgets[[i]]$fg != "" )
								argList$foreground <- widget$.widgets[[i]]$fg
							if( widget$.widgets[[i]]$bg != "" )
								argList$background <- widget$.widgets[[i]]$bg
							.do.gui.call( "tkadd", argList )
					}
					else {
						stop(paste("widget type", widget$.widgets[[i]]$type, "found when expecting a menu or menuitem widget"))
					}

				}
				argList <- list( topMenu, "cascade", label=label, menu=subMenu )
				if( any( widget$font != "" ) )
					argList$font <- .createTkFont( widget$font )
				if( widget$fg != "" )
					argList$foreground <- widget$fg
				if( widget$bg != "" )
					argList$background <- widget$bg				
				.do.gui.call( "tkadd", argList )
			}
			#.createMetaWidget.menu function finished


			#create menus
			for(menu_i in 1:length(guiDesc[[i]]$.menus)) {
				.createMetaWidget.menu(topMenu, guiDesc[[i]]$.menus[[menu_i]])
			}
		}

		#create the widgets
		widgetList <- guiDesc[[i]]$.widgets
		grid_i <- 1 #used for placing widgets into the grid
		while( length(widgetList) > 0 ) {
			tmp <- .createWidget( tt, widgetList, guiDesc[[i]]$windowname )
			sticky = ""
			if( is.null( widgetList[[ 1 ]][[ "sticky" ]] ) == FALSE )
				sticky = widgetList[[ 1 ]][[ "sticky" ]]

			if( is.null( widgetList[[ 1 ]][[ "padx" ]] ) )
				padx = 0
			else
				padx = widgetList[[ 1 ]][[ "padx" ]]

			if( is.null( widgetList[[ 1 ]][[ "pady" ]] ) )
				pady = 0
			else
				pady = widgetList[[ 1 ]][[ "pady" ]]
			
			vert <- guiDesc[[i]]$vertical
			tkgrid( 
					tmp[[ "widget" ]], 
					sticky = sticky, 
					row=ifelse(vert, grid_i, 0), 
					column=ifelse(vert, 0, grid_i ),
				 	pady = pady, 
				 	padx = padx )
			widgetList <- tmp[[ "widgetList" ]]
			grid_i <- grid_i + 1
		}

		#finish setup to any history widgets which have default imports
		if (exists("PBS.history",envir=.PBSmodEnv)) {
			tget(PBS.history)
			j<-grep(paste("^", winName, "\\.", sep=""), names(PBS.history))
			for(n in names(PBS.history)[j]) {
				if (length(PBS.history[[n]]) > 1)
					jumpHistory(n, 1)
			}
		}
	}
	return(invisible(NULL))
}
#----------------------------------------createWin


#closeWin-------------------------------2012-12-06
# Closes a window.
# Arguments:
#   name - window name to close
#----------------------------------------------ACB
closeWin <- function(name)
{
	tget(.PBSmod)
	if (missing(name))
		name <- names(.PBSmod)
	name <- grep("^[^\\.]", name, value=TRUE)
	for(n in name) {
		if(!is.null(.PBSmod[[ n ]])) {
			tt <- .PBSmod[[n]]$tkwindow
			tkdestroy(tt)
		}
	}
}
#-----------------------------------------closeWin


#focusWin-------------------------------2012-12-06
# brings focus to a window (doesn't work from R console)
# args:  winName   - window to focus
#        winVal - if T, make this the active window too
#----------------------------------------------ACB
focusWin <- function(winName, winVal=TRUE)
{
	tget(.PBSmod)
	if( is.null( .PBSmod[[ winName ]] ) )
		stop(paste("supplied winName \"", winName, "\" is not a valid window", sep=""))
	tkfocus(.PBSmod[[winName]]$tkwindow)
	if (winVal) {
		#packList(".activeWin",".PBSmod",winName) #.PBSmod$.activeWin <<- winName
		tget(.PBSmod)
		.PBSmod$.activeWin <- winName
		tput(.PBSmod)
	}
}
#-----------------------------------------focusWin


#getWinVal------------------------------2012-12-06
# All variables starting with "PBS." will not be returned by default
# since they should really be hidden by the user in most cases.
# Arguments:
#  v        - values to get
#  scope    - "L" for local, "P" for .PBSmodEnv, "G" for global, "" for return list only
#  asvector - if T return a vector, if F, return list
#  winName  - specify a specific window if more than one are in use
#-------------------------------------------ACB/RH
getWinVal <- function(v=NULL, scope="", asvector=FALSE, winName="")
{
	if (!exists(".PBSmod",envir=.PBSmodEnv)) {
		stop(".PBSmod was not found")
	}
	tget(.PBSmod)
	if (winName=="") {
		winName <- .PBSmod$.activeWin
		if( is.null( winName ) )
			return( list() )
	}

	if( is.null( .PBSmod[[ winName ]] ) )
		stop(paste("supplied window \"",winName,"\" name not found", sep=""))

	#extract all variables regardless if asked for by user
	vars <- .extractVar(winName)

	#get list of all vars (if user didnt supply any)
	if (is.null(v)) {
		v <- names(vars)
		if (is.null(v))
			return(list()) #no widgets with values found
		v <- v[substr(v,1,4)!="PBS."]
		if (!length(v))
			return(list()) #no widgets with values found
	}

	if (asvector)
		vals <- vector()
	else
		vals <- list()

	#iterate over all var names
	for(key in v) {
		if (asvector)
			vals[key] <- vars[[key]]
		else
			vals[[key]] <- vars[[key]]

		if (scope=="L")
			assign(key,vars[[key]],pos=parent.frame(1))
		else if (scope=="P")
			assign(key, vars[[key]], envir = .PBSmodEnv)
		else if (scope=="G")
			eval(parse(text="assign(key, vars[[key]], envir = .GlobalEnv)"))
	}
	return(vals)
}
#----------------------------------------getWinVal


#setWinVal------------------------------2012-12-04
# Updates a widget with a new value
# Arguments:
#  vars       - named list or vector specifying new values
#  winName - which window to update if multiple are active
#-------------------------------------------ACB/RH
setWinVal <- function(vars, winName="")
{
	tget(.PBSmod)
	if (winName=="")
		winName <- .PBSmod$.activeWin
	if( is.null( .PBSmod[[ winName ]] ) )
		stop(paste("unable to find .PBSmod$", winName))

	if (!length(vars))
		return(vars)

	name <- names(vars)
	for(i in 1:length(vars)) {

		if (is.list(vars))
			.setWinValHelper(name[i], vars[[i]], winName)
		else if (is.vector(vars))
			.setWinValHelper(name[i], vars[i], winName)
	}
}
#----------------------------------------setWinVal


#clearWinVal----------------------------2012-12-05
#   removes any global variables that have a name
#   which corresponds to a name in the window desc file
#-------------------------------------------ACB/RH
clearWinVal <- function() 
{
	objs <- names(getWinVal())
	globs <- ls(all.names=TRUE,pos=".PBSmodEnv") #.GlobalEnv")
	rmlist <- intersect(objs,globs)
	rm(list=rmlist,pos=".PBSmodEnv") #.GlobalEnv")
	invisible(rmlist)
}
#--------------------------------------clearWinVal


#chooseWinVal---------------------------2008-09-05
# Allows user to choose a string value from choices and write 
# chosen string into specified variable of specified window.
# Arguments:
#    choice  - vector of strings to choose from
#    varname - variable name to which choice is assigned in the target GUI.
#    winname - window name for getChoice
#-----------------------------------------------RH
chooseWinVal <- function(choice,varname,winname="window") {
	setPBSoptions("setChoice",NULL);
	setPBSoptions("setChoice",
		paste(";\nsetWinVal(list(",varname,"=chosen),winName=\"",winname,"\");}",sep="",collapse=""));
	getChoice(choice=choice,question="Select from:",horizontal=FALSE,radio=TRUE,qcolor="red3",gui=TRUE,quiet=TRUE);
	setPBSoptions("setChoice",NULL); }
#-------------------------------------chooseWinVal


#getChoice------------------------------2008-09-05
# Prompts user for an input from choices displayed in a GUI.
# The default getChoice() yields TRUE or FALSE.
# Answer is stored in .PBSmod$options$getChoice (or whatever winname is supplied).
# Arguments:
#   choice     - vector of strings to choose from
#   question   - question or prompting statement
#   winname    - window name for getChoice (default="getChoice")
#   horizontal - if T, display the choices horizontally, else vertically 
#   radio      - if T, display the choices as radio buttons, else buttons
#   qcolor     - colour for question
#   gui        - if T, functional when called from a GUI, else functional from command lines
#   quiet      - if T, don't print choice on command line.
# Examples:
#   getChoice("What do you want?",c("Everything","Nothing","Lunch","Money","Fame"),qcolor="red",gui=F)
#   getChoice("Who`s your daddy?",c("Stephen Harper","Homer Simpson","Jon Schnute"),horiz=F,radio=T,gui=F)
#-----------------------------------------------RH
getChoice <- function(choice=c("Yes","No"),question="Make a choice: ",winname="getChoice",
                      horizontal=TRUE, radio=FALSE,qcolor="blue",gui=FALSE,quiet=FALSE) {

	#Construct the hidden choice function
	fn1 <- paste(".makeChoice <- function(){
		act <- getWinAct(winName=\"",winname,"\")[1];
		if (act==\"Yes\") answer <- TRUE
		else if (act==\"No\") answer <- FALSE
		else answer <- act;
		setPBSoptions(\"",winname,"\",answer);
		closeWin(\"",winname,"\") }",sep="",collapse="");
	eval(parse(text=fn1))
	assign(".makeChoice",.makeChoice,envir=.PBSmodEnv)

	#Construct an onClose function
	fn2 <- paste(
		".closeChoice <- function() {\n",
		"chosen <- getPBSoptions(\"",winname,"\");\n",
		"active <- getPBSoptions(\"activeWin\");\n",
		ifelse(quiet,"","print(chosen);\n"),
		"if (is.null(chosen)) setPBSoptions(\"",winname,"\",\"abort\")\n",
		"if (gui && !is.null(active)) focusWin(winName=active);\n",
		"invisible(chosen)}",sep="",collapse="");
	setChoice <- getPBSoptions("setChoice");
	if (!is.null(setChoice))
		fn2 <- sub(";\\\ninvisible\\(chosen\\)}",setChoice,fn2) # only used by chooseWinVal
	eval(parse(text=fn2));
	assign(".closeChoice",.closeChoice,envir=.PBSmodEnv)

	#Construct the Window Description file
	n <- length(choice); ni <- 0;
	nrow <- ifelse(horizontal,1,n); ncol <- ifelse(horizontal,n,1);
	btype <- ifelse(radio,"radio","button");
	btext <- paste("blist <- c(\"window name=\\\"",winname,"\\\" title=Choice",sep="",collapse="");
	btext <- paste(btext," onClose=.win.closeChoice\",",sep="",collapse="")
	#if(!gui) btext <- paste(btext,"\", ",sep="",collapse="");
	qtext <- paste("\"label text=\\\"",question,"\\\" font=\\\"bold 10\\\" fg=\\\"",
		qcolor,"\\\" sticky=W\",",sep="",collapse="");
	btext <- paste(btext,qtext,"\"grid ",nrow," ",ncol," sticky=W\",",sep="",collapse="")

	for (i in choice) {
		if (radio) {
			ni <- ni + 1;
			btext <- paste(btext,
				paste("\"radio text=\\\"",i,"\\\" name=myC sticky=W value=",ni,
				" function=.win.makeChoice action=\\\"",i,"\\\"\",",sep=""),sep="",collapse="") }
		else {
			btext <- paste(btext,paste("\"button text=\\\"",i,"\\\" action=\\\"",i,
				"\\\" function=.win.makeChoice sticky=W\",",sep=""),sep="",collapse="") }
	}
	btext <- paste(btext,"\"\")",sep="",collapse="");
	eval(parse(text=btext));
	#if (exists(".PBSmod")) {
	#	setPBSoptions(winname,NULL); setPBSoptions("activeWin",.PBSmod$.activeWin) }
	if (exists(".PBSmod",envir=.PBSmodEnv)) {
		setPBSoptions(winname,NULL)
		setPBSoptions("activeWin",.PBSmodEnv$.PBSmod$.activeWin) } # changes .PBSMod

	#Create the Window Description file
	createWin(blist,astext=TRUE)
	if (radio) setWinVal(list(myC=0),winName=winname)
	answer <- NULL
	if (!gui) {
		while(is.null(answer)) {answer <- getPBSoptions(winname) } } 
	invisible(answer)
}
#----------------------------------------getChoice


#getWinAct------------------------------2012-12-20
getWinAct <- function(winName)
{
	if (!exists(".PBSmod",envir=.PBSmodEnv)) {
		stop(".PBSmod was not found")
	}
	tget(.PBSmod)
	if (missing(winName))
		winName <- .PBSmod$.activeWin
	return(.PBSmod[[winName]]$action)
}
#----------------------------------------getWinAct


#setWinAct------------------------------2012-12-20
setWinAct <- function(winName, action)
{
	if (is.null(action))
		return()
	tget(.PBSmod)
	if (length(.PBSmod[[winName]]$actions) >= .maxActionSize)
		#eval(parse(text=".PBSmod[[winName]]$actions <<- .PBSmod[[winName]]$actions[1:(.maxActionSize-1)]"))
		.PBSmod[[winName]]$actions <- .PBSmod[[winName]]$actions[1:(.maxActionSize-1)]
	#eval(parse(text=".PBSmod[[winName]]$actions <<- c(action, .PBSmod[[winName]]$actions)"))
	.PBSmod[[winName]]$actions <- c(action, .PBSmod[[winName]]$actions)
	tput(.PBSmod)
	invisible()
}
#----------------------------------------setWinAct


#getWinFun------------------------------2012-12-20
getWinFun <- function(winName)
{
	if (!exists(".PBSmod",envir=.PBSmodEnv)) {
		stop(".PBSmod was not found")
	}
	tget(.PBSmod)
	if (missing(winName))
		winName <- .PBSmod$.activeWin
	return(.PBSmod[[winName]]$functions)
}
#----------------------------------------getWinFun


#parseWinFile---------------------------2012-12-19
#   parse window description file into a list
# Arguments:
#   fname - filename or vector of strings
#   astext - if F, treat fname as a filename
#            if T, treat it as the contents of a file
#----------------------------------------------ACB
parseWinFile <- function(fname, astext=FALSE)
{
	if (astext) {
		#treat "\n" in astext mode as normal whitespace
		srcfile <- orgfile <- sub("\n", " ", fname)
		fname <- "read as text"
	}
	else {
		#read in file
		if (fname=="")
			stop("No filename given")
		srcfile <- orgfile <- scan(fname, what=character(), sep="\n", quiet=TRUE, blank.lines.skip=FALSE)
	}
	data <- list()
	j <- 0
	halt <- FALSE
	extendLine <- FALSE #used for extending a single line into lines with \
	extendLineNumber <- 0 #where a new widget starts - used for error messages
	str <- ""

	if (!length(srcfile)) {
		stop("Input file is empty\n")
	}

	#if comments were striped out earlier, we would lose the line count.
	for(i in 1:length(srcfile)) {
		if (!any(grep("^[[:space:]]*(#.*)?$", srcfile[i]))) {

			srcfile[i] <- .stripComments(srcfile[i])

			#append last string onto new string if applicable
			if (extendLine == TRUE)
				str <- paste(str, srcfile[i], sep=" ")
			else {
				str <- srcfile[i]
				extendLineNumber <- i
			}

			#determine if this string is extended by a \ at the end.
			tmp <- sub('\\\\$', '', str)
			if (tmp==str) #no sub took place
				extendLine = FALSE
			else
				extendLine = TRUE
			str <- tmp

			#parse the line once it is complete (no \)
			if (extendLine == FALSE) {
				tmp <- .getParamFromStr(str, fname, extendLineNumber, i, orgfile)
				if (is.null(tmp)) {
					halt <- TRUE
				}
				else if(halt==FALSE) {
					j <- j + 1
					data[[j]]<-tmp
				}
			}
		}
	}
	if (halt==TRUE) {
		stop("Errors were found in the GUI description file. Unable to continue\n")
	}

	#by this point all widgets from the text file have been converted into
	#an appropriate list of widgets, we will need to setup the nested grids
	#this will result in a more recursive tree-like list.

	#create a blank window if data is empty
	if (!length(data)) {
		data[[1]] <- .getParamFromStr("window") 
	}

	#we must make sure the first element is a window, if not we will insert one
	#as the head of the list
	if (data[[1]]$type != "window") {
		data <- c(1, data)
		data[[1]] <- .getParamFromStr("window") #pull in all defaults from defs.R
	}

	#data[[1]] is now guarenteed to be a window type

	#start parsing the read data - this mostly setups grid data
	parsedData<-list()
	j <- 0; #widget index
	i <- 0; #window index
	k <- 0; #menu index
	while(length(data)) {
		#pull out any options
		if (data[[1]]$type=="window") {
			i <- i + 1 #increment new window index
			j <- 0 #reset widget index
			k <- 0 #reset menu index
			parsedData[[i]] <- list()
			parsedData[[i]]$title <- data[[1]]$title
			parsedData[[i]]$windowname <- data[[1]]$name
			parsedData[[i]]$vertical <- data[[1]]$vertical
			parsedData[[i]]$onclose <- data[[1]]$onclose
			parsedData[[i]]$remove <- data[[1]]$remove
			parsedData[[i]]$winBackground <- data[[1]]$bg
			parsedData[[i]]$winForeground <- data[[1]]$fg
			parsedData[[i]]$.widgets <- list() #holds all widgets
			parsedData[[i]]$.menus <- list() #holds all menu widgets

			data <- data[-1]
		}
		else {
			#look for menu widgets
			if (data[[1]]$type=="menu") {
				k <- k + 1 #increment menu index

				#save menu widget
				parsedData[[i]]$.menus[[k]] <- data[[1]]

				#pull out n menuitem
				tmp <- .parsemenu(data[-1], data[[1]]$nitems)
				parsedData[[i]]$.menus[[k]]$.widgets <- tmp$menuData

				#parse remaining widgets
				data <- tmp$unparsedData
			}

			#look for regular widgets
			else {
				j <- j + 1 #incrememnt widget index

				#save widget
				parsedData[[i]]$.widgets[[j]] <- data[[1]]
				data <- data[-1] #remove widget from to parse list

#				#associate child widgets if grid
#				if (data[[1]]$type=="grid") {
#					tmp <- .parsegrid(data[-1], data[[1]]$nrow, data[[1]]$ncol)
#					parsedData[[i]]$.widgets[[j]]$.widgets <- tmp$gridData
#					#.parsegrid returns all left over widgets
#					data <- tmp$unparsedData
#				}
#				else if (data[[1]]$type=="notebook") {
#					tmp <- .parsegrid(data[-1], data[[1]]$nrow, data[[1]]$ncol)
#					parsedData[[i]]$.widgets[[j]]$.widgets <- tmp$gridData
#					#.parsegrid returns all left over widgets
#					data <- tmp$unparsedData
#				}
#				else {
#					data <- data[-1] #remove widget from to parse list
#				}
			}
		}
	}
	return(parsedData)
}
#-------------------------------------parseWinFile


#updateGUI------------------------------2012-12-06
# Update the active GUI with local values 
#-------------------------------------------ARK/RH
updateGUI <- function(scope="L") {
	# Translate the scope argument into a target environment
	if (!is.environment(scope) && scope=="L")      tenv=parent.frame(n=1)
	else if (!is.environment(scope) && scope=="P") tenv=.PBSmodEnv
	else if (!is.environment(scope) && scope=="G") tenv=.GlobalEnv
	else tenv=scope
	if (!is.environment(tenv)) stop("'scope' must be 'L', 'G', or a vaild R environment")
	parentList = ls( name=tenv )

	#if (!exists(".PBSmod",envir=.GlobalEnv)) return (invisible("'.PBSmod' does not exist"))
	if (!exists(".PBSmod",envir=.PBSmodEnv)) return (invisible("'.PBSmod' does not exist"))
	tget(.PBSmod)
	win = .PBSmod$.activeWin                 # Get the current active window name
	if (is.null(.PBSmod[[win]])) return (invisible("No active window"))
	guiList=.extractVar(win)  # GUI information from .PBSmod[[win]]

	# Check for parent environment variables that match the GUI list.
	isMatch = is.element( parentList,names(guiList) )
	if (any(isMatch)) {
		parentList = parentList[isMatch]
		# Now evaluate the variables into a list.
		nVals = length( parentList )
		vals  = as.list( 1:nVals )
		names( vals ) = parentList
		for (i in parentList) 
			vals[[i]]=get(i,envir=tenv)
		setWinVal( vals ) }
	invisible(isMatch)
}
#----------------------------------------updateGUI


#doAction-------------------------------2011-11-08
# Executes the action created by a widget.
#-----------------------------------------------RH
doAction=function(act){
	tget(.PBSmod)
	if (missing(act)) {
		if(is.null(.PBSmod$.activeWin)) return()
		act=getWinAct()[1] }
	if(is.null(act) || act=="") return()
	
	#get win's environment
	winName <- .PBSmod$.activeWin
	if( !is.null( winName ) )
		envir <- .PBSmod[[ winName ]]$env
	else
		#envir <- globalenv() #maybe parent.frame() is better
		envir <- .PBSmodEnv #maybe parent.frame() is better

	# Translation symbols used in Window Description File to create R-code:
	expr=gsub("`","\"",act)                # convert back-tick to double-quote
	expr=gsub("(_\\.)","\\\\\\\\.",expr)   # convert underscore period to four backslahes and one period
	eval(parse(text=expr),envir=envir)
	invisible(expr) }
#-----------------------------------------doAction


#setWidgetColor-------------------------2012-12-20
# TODO might want to rename this (?)
# Colours can be set (fg, bg) for and droplist, 
# check widgets, (entryfg, entrybg) for entry widget
# these are passed as `...`
#----------------------------------------------ACB
#setWidgetColor <- function( name, radioValue, winName = .PBSmod$.activeWin, ... )
setWidgetColor <- function( name, radioValue, winName = .PBSmodEnv$.PBSmod$.activeWin, ... )
{
	configure.entry <- function( ptr, entryfg, entrybg, noeditfg, noeditbg )
	{
		if( !missing( entryfg ) )
			tkconfigure( ptr, fg = entryfg )
		if( !missing( entrybg ) )
			tkconfigure( ptr, bg = entrybg )
		if( !missing( noeditfg ) )
			tkconfigure( ptr, disabledforeground=noeditfg )
		if( !missing( noeditbg ) ) {
			tkconfigure( ptr, disabledbackground=noeditbg )
			tkconfigure( ptr, readonlybackground=noeditbg )
		}
	}

	configure.droplist <- function( ptr, fg, bg )
	{
		if( !missing( fg ) )
			tkconfigure( ptr, foreground = fg, selectforeground = fg )
		if( !missing( bg ) )
			tkconfigure( ptr, selectbackground = bg, insertbackground = bg, highlightbackground = bg, entrybg = bg )
	}

	configure.check <- function( ptr, fg, bg, disablefg, entryfg, entrybg )
	{
		if( !missing( entryfg ) )
			fg <- entryfg
		if( !missing( entrybg ) )
			bg <- entrybg
		if( !missing( fg ) )
			tkconfigure( ptr, fg = fg )
		if( !missing( bg ) )
			tkconfigure( ptr, bg = bg )
		if( !missing( disablefg ) )
			tkconfigure( ptr, disabledforeground = disablefg )
	}

	configure.slide <- function( ptr, fg, bg )
	{
		if( !missing( fg ) )
			tkconfigure( ptr, fg = fg )
		if( !missing( bg ) )
			tkconfigure( ptr, bg = bg )
	}
	
	configure.label <- function( ptr, fg, bg )
	{
		if( !missing( fg ) )
			tkconfigure( ptr, fg = fg )
		if( !missing( bg ) )
			tkconfigure( ptr, bg = bg )
	}
	
	configure.button <- function( ptr, fg, bg, disablefg )
	{
		if( !missing( fg ) )
			tkconfigure( ptr, fg = fg )
		if( !missing( bg ) )
			tkconfigure( ptr, bg = bg )
		if( !missing( disablefg ) )
			tkconfigure( ptr, disabledforeground = disablefg )
	}

	configure.progressbar <- function( ptr, fg, bg )
	{
		if( !missing( fg ) )
			tcl( ptr, "configure", fg = fg )
		if( !missing( bg ) )
			tcl( ptr, "configure", troughcolor = bg )
	}

	configure.text <- function( ptr, fg, bg )
	{
		if( !missing( fg ) )
			tkconfigure( ptr, fg = fg )
		if( !missing( bg ) )
			tkconfigure( ptr, bg = bg )
	}

	configure.spinbox <- function( ptr, entryfg, entrybg )
	{
		if( !missing( entryfg ) )
			tkconfigure( ptr, foreground = entryfg, selectforeground = entryfg, entryfg = entryfg )
		if( !missing( entrybg ) )
			tkconfigure( ptr, bg = entrybg, insertbackground = entrybg, selectbackground = entrybg, entrybg = entrybg )
	}

	configure.radio <- function( name, radioValue, winName, fg, bg )
	{
		widget_list <- .map.get( winName, name )$tclwidgetlist
		if( missing( radioValue ) )
			radioValue <- names( widget_list ) #use all values (change state for every matching var name)
		else
			radioValue <- as.character( radioValue )
		for( key in radioValue ) {
			ptr = widget_list[[ key ]]
			if( !missing( fg ) )
				tkconfigure( ptr, fg = fg )
			if( !missing( bg ) )
				tkconfigure( ptr, bg = bg )
		}
	}

	#### function starts here ####

	tget(.PBSmod)
	#get window
	winwidget <- .PBSmod[[ winName ]]
	if( is.null( winwidget ) ) 
		stop( paste( "unable to find window:", winName ) )
	
	#get widget
	widget <- winwidget$widgets[[ name ]]
	if( is.null( widget ) )
		stop( paste( "unable to find widget: ", name ) )

	#HACK to change initial values (as set in window desc file)
	myargs = list( ... )
	if( !is.null( myargs[[ "noeditfg" ]] ) ) {
		widget$noeditfg = myargs[[ "noeditfg" ]]
		#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ name ]] <<- widget"))
		tget(.PBSmod)
		.PBSmod[[ winName ]]$widgets[[ name ]] <- widget
		tput(.PBSmod)
	}
	if( !is.null( myargs[[ "fg" ]] ) ) {
		widget$fg = myargs[[ "fg" ]]
		#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ name ]] <<- widget"))
		tget(.PBSmod)
		.PBSmod[[ winName ]]$widgets[[ name ]] <- widget
		tput(.PBSmod)
	}
	if( !is.null( myargs[[ "entryfg" ]] ) ) {
		widget$entryfg = myargs[[ "entryfg" ]]
		#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ name ]] <<- widget"))
		tget(.PBSmod)
		.PBSmod[[ winName ]]$widgets[[ name ]] <- widget
		tput(.PBSmod)
	}

	#get tcl ptr to tk widget
	tget(.PBSmod)
	widget_ptr <- .PBSmod[[ winName ]]$widgetPtrs[[ name ]]$tclwidget

	#special case for radio widgets
	if( widget$type == "radio" ) {
		configure.radio( name = name, radioValue = radioValue, winName = winName, ... )
		return(invisible())
	}

	#special case for matrix
	if (widget$type=="matrix") {
		if (length(widget$names)==1) {
			for(i in 1:widget$nrow)
				for(j in 1:widget$ncol) {
					argList <- list(paste(name,"[",i,",",j,"]",sep=""), winName = winName, ... ) 
					if( widget$mode == "logical" ) {
						argList[[ "entrybg" ]] <- NULL
						argList[[ "noeditbg" ]] <- NULL
						if( is.null( argList[[ "entryfg" ]] ) == FALSE ) {
							argList[[ "fg" ]] <- argList[[ "entryfg" ]]
							argList[[ "entryfg" ]] <- NULL
						}
						if( is.null( argList[[ "noeditfg" ]] ) == FALSE ) {
							argList[[ "disablefg" ]] <- argList[[ "noeditfg" ]]
							argList[[ "noeditfg" ]] <- NULL
						}
					}
					.do.gui.call( setWidgetColor, argList )
				}
			return(NULL)
		}
	}

	#special case for data
	if (widget$type=="data") {
		if (length(widget$names)==1) {
			for(i in 1:widget$nrow)
				for(j in 1:widget$ncol) {
					argList <- list(paste(name,"[",i,",",j,"]d",sep=""), winName = winName, ... ) 
					if( widget$modes[ j ] == "logical" ) {
						argList[[ "entrybg" ]] <- NULL
						argList[[ "noeditbg" ]] <- NULL
						if( is.null( argList[[ "entryfg" ]] ) == FALSE ) {
							argList[[ "fg" ]] <- argList[[ "entryfg" ]]
							argList[[ "entryfg" ]] <- NULL
						}
						if( is.null( argList[[ "noeditfg" ]] ) == FALSE ) {
							argList[[ "disablefg" ]] <- argList[[ "noeditfg" ]]
							argList[[ "noeditfg" ]] <- NULL
						}
					}
					.do.gui.call( setWidgetColor, argList )
				}
			return(NULL)
		}
	}

	#special case for vector
	if (widget$type=="vector") {
		if (length(widget$names)==1) {
			for(i in 1:widget$length) {
				argList <- list(paste(name,"[",i,"]",sep=""), winName = winName, ... ) 
				if( widget$mode == "logical" ) {
					argList[[ "entrybg" ]] <- NULL
					argList[[ "noeditbg" ]] <- NULL
					if( is.null( argList[[ "entryfg" ]] ) == FALSE ) {
						argList[[ "fg" ]] <- argList[[ "entryfg" ]]
						argList[[ "entryfg" ]] <- NULL
					}
					if( is.null( argList[[ "noeditfg" ]] ) == FALSE ) {
						argList[[ "disablefg" ]] <- argList[[ "noeditfg" ]]
						argList[[ "noeditfg" ]] <- NULL
					}
				}
				.do.gui.call( setWidgetColor, argList )
			}
			return(NULL)
		}
	}

	#scrolling object
	if (widget$type=="object") {
		return( setWidgetColor(paste("[superobject]", name,sep=""), winName = winName, ... ) )
	}



	#call specific config method based on widget type
	func <- paste( "configure.", widget$type, sep="" )
	if( exists( func ) == FALSE )
		stop( paste( "not supported for widget type:", widget$type ) )

	.do.gui.call( func, list( ptr=widget_ptr, ... ) )
	return(invisible())
}
#-----------------------------------setWidgetColor


#setWidgetState-------------------------2012-12-20
setWidgetState <- function( varname, state, radiovalue, winname, warn = TRUE )
{
	if (!exists(".PBSmod",envir=.PBSmodEnv))
		stop(".PBSmod was not found")
	tget(.PBSmod)
	if( missing( winname ) )
		winname <- .PBSmod$.activeWin
	if( is.null( .PBSmod[[ winname ]] ) )
		stop(paste("supplied window \"",winname,"\" name not found"))
	
	if( any( state == c( "disabled", "normal", "readonly", "active" ) ) == FALSE ) 
		stop( "state must be disabled, normal, readonly (for entry), or active( for radio)" )

	x  <- .map.get(winname, varname)
	tget(.PBSmod)
	wid <- .PBSmod[[winname]]$widgets[[varname]]
	if( is.null( wid ) ) stop(paste("supplied widget \"",varname,"\" name not found", sep=""))

	if( any( wid$type == c( "notebook" ) ) )
		stop( paste( wid$type, "widget is not supported" ) )

	#change readonly -> disabled for widgets which dont support readonly
	if( any( wid$type == c( "check", "radio", "droplist", "spinbox", "table", "text" ) ) && state == "readonly" ) {
		state <- "disabled"
		if( warn )
			warning( paste( "setting readonly to disabled (readonly not supported by ", wid$type, " widgets)", sep="" ) )
	}
	
	#special case since radio has several widgets with the same name
	if( wid$type == "radio" ) {
		widget_list <- .map.get( winname, varname )$tclwidgetlist
		if( missing( radiovalue ) )
			radiovalue <- names( widget_list ) #use all values (change state for every matching var name)
		else
			radiovalue <- as.character( radiovalue )
		for( key in radiovalue )
			tkconfigure( widget_list[[ key ]], state=state )
		return(invisible(NULL))
	}

	#if tclwidget is known - set it directly here
	if( !is.null( x[["tclwidget"]] ) ) {
		tkconfigure( x$tclwidget, state=state )
		#get correct foreground
		if( !is.null( wid[[ "noeditfg" ]] ) && !is.null( wid[[ "noeditbg" ]] ) ) {
			fg <- ifelse( any( wid$type == c( "entry", "spinbox" ) ), wid$entryfg, wid$fg )
			fg <- ifelse( state == "normal", fg, wid$noeditfg )
			#reset widget colors
			tkconfigure( x$tclwidget, disabledforeground=fg )
			tkconfigure( x$tclwidget, foreground=fg ) #no readonlyforeground
		}
		return(invisible(NULL))
	}
	
	#deal with more complicated objects
	if (is.null(wid)) {
		stop(paste('unable to set "', varname, '": not found.', sep=""))
	}

	#special case for matrix
	if (wid$type=="matrix") {
		if (length(wid$names)==1) {
			for(i in 1:wid$nrow)
				for(j in 1:wid$ncol)
					setWidgetState(paste(varname,"[",i,",",j,"]",sep=""), state, winname, warn=FALSE)
			return(NULL)
		}
	}

	#special case for data
	if (wid$type=="data") {
		if (length(wid$names)==1) {
			for(i in 1:wid$nrow)
				for(j in 1:wid$ncol)
					setWidgetState(paste(varname,"[",i,",",j,"]d",sep=""), state, winname, warn=FALSE)
			return(NULL)
		}
	}

	#special case for vector
	if (wid$type=="vector") {
		if (length(wid$names)==1) {
			for(i in 1:wid$length)
				setWidgetState(paste(varname,"[",i,"]",sep=""), state, winname, warn=FALSE)
			return(NULL)
		}
	}

	#scrolling object
	if (wid$type=="object") {
		return( setWidgetState(paste("[superobject]", varname,sep=""), state, winname, warn=FALSE) )
	}

	stop(paste("unable to update \"", varname, "\" - unable to handle type:", wid$type, sep=""))
}
#-----------------------------------setWidgetState


#=================================================
#               HIDDEN FUNCTIONS
#=================================================



#===== Create Functions ==========================
#----- (.createWidget, etc.) ---------------------


#.createWidget--------------------------2012-12-19
# Generic function to create most widgets, which
# calls appropriate .createWidget.xxxWidgetTypexxx() func
# Arguments:
#  tk      - frame to attach widget to
#  widget  - widget list
#  winName - active window name
#----------------------------------------------ACB
.createWidget <- function(tk, widgetList, winName)
{
	tget(.PBSmod)
	widget <- widgetList[[ 1 ]]
	#save functions
	if (!is.null(widget[["function"]])) {
		if (widget[["function"]]!="") {
			if (!any(.PBSmod[[winName]]$functions==widget[["function"]]))
				#eval(parse(text=".PBSmod[[winName]]$functions <<- c(.PBSmod[[winName]]$functions, widget[[\"function\"]])"))
				.PBSmod[[winName]]$functions <- c(.PBSmod[[winName]]$functions, widget[["function"]])
		}
	}
	#save widget information by name parameter (and not widget$name)
	#widget name can sometimes be "foo[1,2,3]" or some such combo.
	#where as type="vector" name="foo" is never seen in the regular map
	if (!is.null(widget$name)) {
		if (length(widget$name)==1 && widget$name != "" ) {
			if (is.null(.PBSmod[[winName]]$widgets[[ widget$name ]])) {
				#eval(parse(text=".PBSmod[[winName]]$widgets[[widget$name]] <<- widget"))
				.PBSmod[[winName]]$widgets[[widget$name]] <- widget
			} else {
				#duplicate widget name found -> likely a radio with many options but only one var name
				#save additional widget information under .duplicate
				if( is.null( .PBSmod[[winName]]$widgets[[widget$name]]$.duplicate ) )
					#eval(parse(text=".PBSmod[[winName]]$widgets[[widget$name]]$.duplicate <<- list()"))
					.PBSmod[[winName]]$widgets[[widget$name]]$.duplicate <- list()
				i <- length( .PBSmod[[winName]]$widgets[[widget$name]]$.duplicate ) + 1
				#eval(parse(text=".PBSmod[[winName]]$widgets[[widget$name]]$.duplicate[[ i ]] <<- widget"))
				.PBSmod[[winName]]$widgets[[widget$name]]$.duplicate[[ i ]] <- widget
			}
		}
	}
	tput(.PBSmod)

	#look for a function called .createWidget.WIDGETTYPE
	#all of these functions have the same parameters: (tk, widget, winName)
	func <- paste(".createWidget.", widget$type, sep="")
	if (exists(func,mode="function")) {
		ret <- .do.gui.call(func, list(tk, widgetList, winName))
		if( is.list( ret ) == FALSE )
			stop( paste( func, "didn't return a list( widget = <tkpointer>, widgetList = <uncreated remaining widgets> )" ) )
		return( ret )
	}
	stop(paste("Don't know how to create '", widget$type, "' widget\n", sep=""))
	return()
}
#------------------------------------.createWidget

#.createWidget.button-------------------2012-12-20
.createWidget.button <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	param <- list(parent=tk, text=widget$text)
	if ( any( widget$font != "" ) )
		param$font=.createTkFont(widget$font)
	if (!is.null(widget[["fg"]]) && widget$fg!="")
		param$foreground=widget$fg
	if (!is.null(widget[["bg"]]) && widget$bg!="")
		param$background=widget$bg
	if (is.numeric(widget$width) && widget$width > 0)
		param$width=widget$width
	if (widget[["function"]]!="")
		param$command=function(...) { .extractData(widget[["function"]], widget$action, winName) }
	if (!is.null(widget[["disablefg"]]) && widget$disablefg!="")
		param$disabledforeground = widget$disablefg

	button <- .do.gui.call("tkbutton", param)

	if( !is.null( widget[[ "name" ]] ) ) 
		.map.add(winName, widget$name, tclwidget=button)
	return(list( widget = button, widgetList = widgetList[ -1 ] ) )
}
#-----------------------------.createWidget.button


#.createWidget.check--------------------2012-12-19
.createWidget.check <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	widget$mode = "logical"

	if (widget$checked==TRUE)
		val <- 1
	else
		val <- 0

	argList <- list(parent=tk, text=widget$text)
	if (!is.null(widget[["fg"]]) && widget$fg!="")
		argList$foreground=widget$fg
	if (!is.null(widget[["bg"]]) && widget$bg!="")
		argList$background=widget$bg
	if (!is.null(widget[["font"]]) && any( widget$font!=""))
		argList$font <- .createTkFont(widget$font)
	if (!is.null(widget[["disablefg"]]) && widget$disablefg!="")
		argList$disabledforeground = widget$disablefg

	argList$variable <- .map.add(winName, widget$name, tclvar=tclVar(val))$tclvar
	argList$command=function(...) { .extractData(widget[["function"]], widget$action, winName)}

	tkWidget<-.do.gui.call("tkcheckbutton", argList)
	.map.set( winName, widget$name, tclwidget=tkWidget )
	if( widget$edit == FALSE )
		tkconfigure( tkWidget, state="disabled" )

	return(list( widget = tkWidget, widgetList = widgetList[ -1 ] ) )
}
#------------------------------.createWidget.check


#.createWidget.data---------------------2012-12-19
.createWidget.data <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]

	nrow <- widget$nrow
	ncol <- widget$ncol

	names <- widget$names
	modes <- widget$modes

	rowlabels <- widget$rowlabels
	rownames <- widget$rownames
	collabels <- widget$collabels
	colnames <- widget$colnames

	if (all(widget$values==""))
		values <- ""
	else {
		values <- widget$values
		#dim(values) <- c(nrow, ncol)
	}

	wid <- list(type="grid", bg=widget$bg, borderwidth=widget$borderwidth) #new grid widget to create
	wid$.widgets <- list() #elements in the grid
	wid$byrow <- TRUE

	nNames <- length(names)
	nModes <- length(modes)
	nRowlabels <- length(rowlabels)
	nRowNames <- length(rownames)
	nValues <- length(values)
	nCollabels <- length(collabels)
	nColNames <- length(colnames)

	#count names
	if (nNames!=1 && nNames!=(ncol*nrow))
		.stopWidget(paste('names argument must contain 1 or',ncol*nrow,'names seperated by whitespace.'), widget$.debug, winName)

	#count modes
	if (nModes!=1 && nModes!=ncol)
		.stopWidget(paste('modes argument must contain 1 or',ncol,'modes seperated by whitespace.'), widget$.debug, winName)
    
	#count rowlabels
	if (nRowlabels!=1 && nRowlabels!=0 && nRowlabels!=nrow)
		.stopWidget(paste('rowlabels should contain 1 or',nrow,'labels.'), widget$.debug, winName)

	#count rownames
	if (nRowNames!=1 && nRowNames!=nrow)
		.stopWidget(paste('rownames argument should contain 1 or',nrow,'labels.'), widget$.debug, winName)

	#count collabels
	if (nCollabels!=1 && nCollabels!=0 && nCollabels!=ncol)
		.stopWidget(paste('collabels argument should contain 1 or',ncol,'labels.'), widget$.debug, winName)

	#count colnames
	if (nColNames!=1 && nColNames!=ncol)
		.stopWidget(paste('colnames argument should contain 1 or',ncol,'labels.'), widget$.debug, winName)

	if( is.null( rowlabels ) )
			colLabelOffset <- 0
		else
			colLabelOffset <- 1

	#single column labels should be displayed as the title
	if( is.null( collabels ) ) {
		#nothing
	} else if (nCollabels==1 && ncol>1) {
		wid$toptitle<-collabels[1]
		wid$topfont<-widget$font
		wid$toptitle.offset<-1 #to help center the label
		#have counting labels above each column
		wid$.widgets[[1]] <- list()
		wid$.widgets[[1]][[1]] <- list(type='label', text="", bg=widget$bg, fg=widget$fg)
		for(j in 1:ncol) {
			wid$.widgets[[1]][[j+colLabelOffset]] <- list(type='label', text=j, font=widget$font, bg=widget$bg, fg=widget$fg)
		}
	} else {
		wid$.widgets[[1]] <- list()
		wid$.widgets[[1]][[1]] <- list(type='label', text="", bg=widget$bg, fg=widget$fg)
		for(j in 1:ncol) {
			wid$.widgets[[1]][[j+colLabelOffset]] <- list(type='label', text=collabels[j], font=widget$font, bg=widget$bg, fg=widget$fg)
		}
	}

	#row title
	if( is.null( rowlabels ) ) {
		wid$sidetitle.offset<-1 #to help center the label
	} else if (nRowlabels==1 && nrow>1) {
		wid$sidetitle<-rowlabels[1]
		wid$sidefont<-widget$font
		wid$sidetitle.offset<-1 #to help center the label
	}

	for(i in 1:nrow) {
		rowCount <- i #the first row of inputs should be 1 (even if there are labels ontop)
		i <- i + 1 #first row has labels
		wid$.widgets[[i]] <- list()

		for(j in 1:(ncol+1)) {
			#first col is for labels
			if (j==1) {
				if (!is.null(rowlabels)) {
					if (nRowlabels==1 && nrow>1) {
						text <- as.character(rowCount)
					}
					else
						text <- rowlabels[rowCount]

					tmp_i <- i
					#define a row label (per each row)
					row_number = tmp_i - 1 #for renaming row labels
					label_name <- paste( widget$name, "[rowlabel][", row_number, "]", sep="" )

					if (is.null( collabels ))
						tmp_i <- tmp_i - 1

					if( nNames != 1 ) label_name <- ""
					wid$.widgets[[tmp_i]][[j]] <- list(type='label', text=text, name=label_name, mode="character", font=widget$font, bg=widget$bg, fg=widget$fg )
					if( !is.null( widget[[ ".rowlabelwidth" ]] ) )
						wid$.widgets[[tmp_i]][[j]]$width <- widget$.rowlabelwidth
				}
			}
			else {
				if (nNames==1) #single name given
					name <- paste(names,'[',rowCount,',',j-1,']d',sep="")
				else #many names given
					if (widget$byrow)
						name <- names[j-1+(ncol*(i-2))]
					else
						name <- names[i-1+(nrow*(j-2))]
				if (nValues==1) {
					value <- values
				}
				else {
					if (widget$byrow)
						value <- values[j-1+(ncol*(i-2))]
					else
						value <- values[i-1+(nrow*(j-2))]
				}
				if (nModes==1)
					mode <- modes[1]
				else
					mode <- modes[j-1] #columns are offset by one

				#get width
				width <- rep( widget$width, length=ncol+1)[ j - 1 ]

				#tweak offset if labels are disabled
				# ****must be un-tweaked after creating the entry or check widget
				if (is.null( rowlabels ))
					j <- j - 1
				if (is.null( collabels ))
					i <- i - 1

				if (mode=="logical") {
					#display a checkbox
					if (is.na(as.logical(value)))
						checked=FALSE
					else if (as.logical(value)==TRUE)
						checked=TRUE
					else
						checked=FALSE
					wid$.widgets[[i]][[j]] <- list(
						  type="check",
						  mode="logical",
						  name=name,
						  text="",
						  "function"=widget[["function"]],
						  action=widget$action,
						  checked=checked,
						  bg=widget$bg,
						  fg=widget$entryfg,
						  noeditbg=widget$noeditbg,
						  noeditfg=widget$noeditfg,
						  edit=widget$edit
					)
				}
				else {
					#display a entry box
					wid$.widgets[[i]][[j]] <- list(
						  type='entry', 
						  name=name,
						  "function"=widget[["function"]],
						  .up_func=widget[[".up_func"]],
						  .down_func=widget[[".down_func"]],
						  .pageup_func=widget[[".pageup_func"]],
						  .pagedown_func=widget[[".pagedown_func"]],
						  action=widget$action,
						  enter=widget$enter,
						  value=value,
						  width=width,
						  mode=mode,
						  entryfont=widget$entryfont,
						  entrybg=widget$entrybg,
						  entryfg=widget$entryfg,
						  noeditbg=widget$noeditbg,
						  noeditfg=widget$noeditfg,
						  edit=widget$edit
					)
				}

				# ***untweak offset for special case of no labels
				if (is.null( rowlabels ))
					j <- j + 1
				if (is.null( collabels ))
					i <- i + 1
			}
		}
	}

	#look out for a trailing list (only happens if rowlabels=NULL)
	i <- length(wid$.widgets)
	if (!length(wid$.widgets[[i]]))
		wid$.widgets[[i]] <- NULL

	#remove titles if applicable
	if (is.null(rowlabels)) {
		wid$sidetitle <- ""
		wid$toptitle.offset <- NULL
	}
	if (is.null(collabels)) {
		wid$toptitle <- ""
		wid$sidetitle.offset <- NULL
	}

	widgets <- .convertOldGridToNewGrid( wid )
	tmp <- .createWidget.grid(tk, widgets, winName)
	return(list( widget = tmp$widget, widgetList = widgetList[ -1 ] ) )
}
#-------------------------------.createWidget.data


#.createWidget.droplist-----------------2012-12-20
.createWidget.droplist <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	#assert( choices != NULL xor values != NULL )
	if( is.null( widget[[ "choices" ]] ) && is.null( widget[[ "values" ]] ) )
		.stopWidget( paste( "either choices or values must be specified", sep="" ), widget$.debug, winName )

	if( !is.null( widget[[ "choices" ]] ) && !is.null( widget[[ "values" ]] ) )
		.stopWidget( paste( "only one of choices or values can be specified", sep="" ), widget$.debug, winName )

	if( !is.null( widget[[ "choices" ]] ) )
		values <- .getValueForWidgetSetup( widget$choices, widget, winName )
	else
		values <- widget$values
	labels <- widget[[ "labels" ]]
	if( is.null( labels ) ) labels <- values
	if ( length( labels ) != length( values ) )
		.stopWidget( paste( "values and labels must be the same size (or NULL)", sep="" ), widget$.debug, winName )
	#create real tk widget below

	#tk splits on " " if only one string is given - but the \ escape isn't working great.
	if( length( labels ) == 1 ) 
		new_labels <- gsub( " ", "\\\\ ", labels ) #should replace with "\\\\ ", but it still doesn't work well
	else
		new_labels <- labels

	argList <- list(parent=tk, type="ComboBox", editable=widget$add,values=new_labels)
	if (!is.null(widget[["fg"]]) && widget$fg!="") {
		#see http://tcltk.free.fr/Bwidget/ComboBox.html for possible options
		argList$foreground=widget$fg
		#argList$entryfg=widget$fg
		#argList$selectforeground=widget$fg
	}
	if (!is.null(widget[["bg"]]) && widget$bg!="") {
		#argList$background=widget$bg #this affects the colour of the drop down arrow
		#argList$selectbackground=widget$bg #covers what's selected - but stays after the item is selected
		#argList$insertbackground=widget$bg #color of insert cursor - leave as default
		argList$highlightbackground=widget$bg
		argList$entrybg=widget$bg
	}
	if (!is.null(widget[["font"]]) && any(widget$font!=""))
		argList$font <- .createTkFont(widget$font)
	argList$textvariable<-.map.add(winName, widget$name, tclvar=tclVar(labels[ widget$selected ]))$tclvar
	argList$width<-widget$width


	#callback
	argList$modifycmd = function(...) { .extractData(widget[["function"]], widget$action, winName)}
	drop_widget <- .do.gui.call( "tkwidget", argList )
	if( length( labels ) == 1 )
		tclvalue( argList$textvariable ) <- labels


	#save widget - so we can use tcl( drop_widget, "getvalue" ) at a later time
	.map.set(winName, widget$name, tclwidget=drop_widget )
	.map.set(winName, widget$name, droplist_values=values )

	.map.set( winName, paste( widget$name, ".values", sep="" ), droplist_widget=drop_widget )
	.map.set( winName, paste( widget$name, ".id", sep="" ), droplist_widget=FALSE )
	#eval(parse(text=".PBSmod[[winName]]$widgets[[ paste( widget$name, \".values\", sep=\"\" ) ]]$labels <<- values"))
	tget(.PBSmod)
	.PBSmod[[winName]]$widgets[[ paste( widget$name, ".values", sep="" ) ]]$labels <- values
	tput(.PBSmod)

	if( widget$edit == FALSE )
		tkconfigure( drop_widget, state="disabled" )


	enter <- !is.null(widget$enter)
	if (enter) enter <- widget$enter
	if (enter == FALSE) {
		#not sure what to bind to get a key by key update, validatecommand looks like an interesting hack, but return value isn't working
		.stopWidget( paste( "enter=F is not implemented", sep="" ), widget$.debug, winName )
	}
	#command will return after user hits enter
	tkconfigure( drop_widget, command=function(...) { .extractData(widget[["function"]], widget$action, winName) } )


	return(list( widget = drop_widget, widgetList = widgetList[ -1 ] ) )
}
#---------------------------.createWidget.droplist


#.createWidget.entry--------------------2012-12-19
.createWidget.entry <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	if (!is.null(widget[["label"]]) && widget$label!="") {
		#if label is set, then create a 2x1 grid
		label <- widget$label
		widget$label <- "" #blank it out, inf loop if not.
		widgets <- list( 
			list(type="grid", nrow=1, ncol=2, font="", byrow=TRUE, borderwidth=1, relief="flat", padx=0, pady=0, fg=widget$fg, bg=widget$bg ),
			list(type="label", text=label, padx=0, pady=0, font=widget$font, fg=widget$fg, bg=widget$bg),
			widget
		)
		tmp <- .createWidget.grid(tk, widgets, winName)
		return(list( widget = tmp$widget, widgetList = widgetList[ -1 ] ) )
	}
	#create real tk widget below
	argList <- list(parent=tk)
	if (!is.null(widget[["entryfg"]]) && widget$entryfg!="")
		argList$foreground=widget$entryfg
	if (!is.null(widget[["entrybg"]]) && widget$entrybg!="")
		argList$background=widget$entrybg
	if (!is.null(widget[["entryfont"]]) && any(widget$entryfont!=""))
		argList$font <- .createTkFont(widget$entryfont)
	argList$textvariable <- .map.add(winName, widget$name, tclvar=tclVar(widget$value))$tclvar
	argList$width <- widget$width

#if (widget$name=="something") browser()

	tkWidget<-.do.gui.call("tkentry", argList)
	.map.set( winName, widget$name, tclwidget=tkWidget )

	if( widget$edit == FALSE ) {
		tkconfigure( tkWidget, state="readonly" )
		if (!is.null(widget[["noeditfg"]]) && widget$noeditfg!="") {
			tkconfigure( tkWidget, disabledforeground=widget$noeditfg )
			tkconfigure( tkWidget, foreground=widget$noeditfg ) #no readonlyforeground
		}
	}
	if (!is.null(widget[["noeditbg"]]) && widget$noeditbg!="") {
		tkconfigure( tkWidget, disabledbackground=widget$noeditbg )
		tkconfigure( tkWidget, readonlybackground=widget$noeditbg )
	}

	#if (!is.null(widget[["noeditfg"]]) && widget$noeditfg!="") {
	#	tkconfigure( tkWidget, disabledforeground="red" )
	#	tkconfigure( tkWidget, foreground="red" )
	#}
		
	if( !is.null( widget[["password"]] ) && widget$password == TRUE )
		tkconfigure( tkWidget, show="*")

	enter <- !is.null(widget[["enter"]])
	if (enter)
		enter <- widget$enter
	if (enter) {
		#dont update it (unless an return was pressed) as it can slow it down a lot
		tkbind(tkWidget,"<KeyPress-Return>",function(...) { .extractData(widget[["function"]], widget$action, winName)})
	}
	else
		tkbind(tkWidget,"<KeyRelease>",function(...) { .extractData(widget[["function"]], widget$action, winName)})
	if( !is.null( widget[[".up_func"]] ) ) {
		tkbind(tkWidget,"<Up>",function(...) { .do.gui.call( widget[[".up_func"]], list( selected_widget_name = widget$name, ... ) ) } )
		tkbind(tkWidget,"<Prior>",function(...) { .do.gui.call( widget[[".pageup_func"]], list( selected_widget_name = widget$name, ... ) ) } )
	}
	if( !is.null( widget[[".down_func"]] ) ) {
		tkbind(tkWidget,"<Down>",function(...) { .do.gui.call( widget[[".down_func"]], list( selected_widget_name = widget$name, ... ) ) } )
		tkbind(tkWidget,"<Next>",function(...) { .do.gui.call( widget[[".pagedown_func"]], list( selected_widget_name = widget$name, ... ) ) } )
	}
	return(list( widget = tkWidget, widgetList = widgetList[ -1 ] ) )
}
#------------------------------.createWidget.entry


#.createWidget.grid---------------------2012-12-19
.createWidget.grid <- function(tk, widgetList, winName)
{
	#these "defaults" only apply to the first layer grid
	#because it is added as padding, and not parsed.
	#all other defaults are set in paramOrder list

	widget <- widgetList[[ 1 ]]

	if (is.null(widget[["borderwidth"]]))
		widget$borderwidth <- 5

	if (is.null(widget[["relief"]]))
		widget$relief <- "flat"

	argList <- list(parent=tk, borderwidth=widget$borderwidth,relief=widget$relief)
	if (!is.null(widget[["bg"]]) && widget$bg!="")
		argList$background=widget$bg
	if (!is.null(widget[["font"]]) && any(widget$font!=""))
		argList$font <- .createTkFont(widget$font)

	tkWidget<-.do.gui.call("tkframe", argList)
        
	#call buildgrid to attach all children widgets to grid
	tmp <- .buildgrid(tkWidget, widget, winName, widgetList[ - 1 ])

	return( tmp )
}
#-------------------------------.createWidget.grid


#.createWidget.image--------------------2012-12-19
.createWidget.image <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	argList <- list(parent=tk)

	#one must be given, the other must be null
	if( is.null( widget[["file"]] ) == is.null( widget[["varname"]] ) )
    	.stopWidget("a value must be specified for one of either `file' or `varname'", widget$.debug, winName)
	
	if( is.null( widget[["varname"]] ) == FALSE ) {
		#get image from variable
		tget(.PBSmod)
		if( !exists( widget$varname, envir=.PBSmod[[ winName ]]$env ) ) {
			msg <- paste( "unable to load image from variable \"", widget[["varname"]], "\" - variable not found", sep="" )
			#R is crashing here.... wtf? .stopWidget( msg, widget$.debug, winName ) #maybe its crashing when the window is closed by tkclose
			stop( msg )
		}
		image <- get( widget$varname, envir = .PBSmod[[ winName ]]$env )
	} else {
		image = widget$file
	}

	#check file exists
	if( file.exists( image ) == FALSE ) {
		msg <- paste( "unable to open image \"", image, "\" - file does not exist", sep="" )
		stop( msg )
	}

	image.id = paste( "PBSmodelling::", image, sep="" )

	x = try(tcl("image","create","photo", image.id, file=image ),silent=TRUE)
	if( inherits( x, "try-error" ) ) {
    	.stopWidget(paste("unable to open file", widget$file, " - only gif is supported"), widget$.debug, winName)
	}
	if( !is.null( widget[["subsample"]] ) && widget$subsample > 1 ) {
		#subsample the image (not a true resize, but the closest we'll get without using any fancy img libs
		image.sized = paste( image.id, "::subsample", sep="" )
		sized_img <- tcl("image","create","photo", image.sized )
		tcl( sized_img, "copy", x, subsample=widget$subsample )
		argList$image = image.sized #use subsampled image
	} else {
		argList$image = image.id #use original image
	}

	tkWidget<-.do.gui.call("tklabel", argList)
	return( list( widget = tkWidget, widgetList = widgetList[ -1 ] ) )
}
#------------------------------.createWidget.image


#.createWidget.include------------------2012-12-19
.createWidget.include <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	if( !is.null( widget[[ "file" ]] ) && !is.null( widget[[ "name" ]] ) )
    		.stopWidget( "both file and name can not be set at the same time", widget$.debug, winName)
	if( is.null( widget[[ "file" ]] ) && is.null( widget[[ "name" ]] ) )
    		.stopWidget( "either file or name must be supplied", widget$.debug, winName)

	file <- widget[[ "file" ]]
	if( is.null( file ) ) {
		if( exists( widget[[ "name" ]] ) == FALSE )
			.stopWidget( paste( "unable to load included filename from variable \"", widget[[ "name" ]], "\" - variable not found", sep="" ), widget$.debug, winName )
		file <- get( widget[[ "name" ]] )
	}

	if( file.exists( file ) == FALSE )
		.stopWidget( paste( "unable to load included filename \"", file, "\" - file does not exist", sep="" ), widget$.debug, winName )

	gui_desc <- parseWinFile( file )
	if( length( gui_desc ) == 0 ) {
                frm <- tkframe( tk )
                .adjustAllColours( frm )
		return( frm )
        }
	if( length( gui_desc ) > 1 ) warning( "Multiple windows found in the window description file - only the first will be included (and displayed)" )
	grid <- .packWidgetsIntoGrid( gui_desc[[ 1 ]]$.widgets, gui_desc[[ 1 ]]$vertical )
	return( list( widget = .createWidget( tk, grid, winName )$widget, widgetList = widgetList[ -1 ] )  )
}
#----------------------------.createWidget.include


#.createWidget.label--------------------2012-12-19
.createWidget.label <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	argList <- list(parent=tk)
	if( !is.null(widget[["name"]]) && widget$name != "" ) {
		argList$text<-tclvalue( .map.add(winName, widget$name, tclvar=tclVar(widget$text))$tclvar )
	} else {
		argList$text = widget$text
	}
	if (!is.null(widget[["fg"]]) && widget$fg!="")
		argList$foreground=widget$fg
	if (!is.null(widget[["bg"]]) && widget$bg!="")
		argList$background=widget$bg
	if (!is.null(widget[["font"]]) && any(widget$font!=""))
		argList$font <- .createTkFont(widget$font)
	if (!is.null(widget[["anchor"]]) && widget$anchor!="")
		argList$anchor <- tolower( widget$anchor )
	if (!is.null(widget[["justify"]]) && widget$justify!="")
		argList$justify <- widget$justify
	if (!is.null(widget[["wraplength"]]) && widget$wraplength > 0)
		argList$wraplength <- widget$wraplength 
	if (!is.null(widget[["width"]]) && widget$width > 0)
		argList$width <- widget$width 

	tkWidget<-.do.gui.call("tklabel", argList)
	if( !is.null(widget[["name"]]) && widget$name != "" ) {
		tkconfigure( tkWidget,textvariable = .map.get(winName, widget$name )$tclvar )
		#eval(parse(text=".PBSmod[[ winName ]]$widgetPtrs[[ widget$name ]]$tclwidget <<- tkWidget"))
		tget(.PBSmod)
		.PBSmod[[ winName ]]$widgetPtrs[[ widget$name ]]$tclwidget <- tkWidget
		tput(.PBSmod)
	}
	return( list( widget = tkWidget, widgetList = widgetList[ -1 ] ) )
}
#------------------------------.createWidget.label


#.createWidget.matrix-------------------2012-12-19
.createWidget.matrix <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	nrow <- widget$nrow
	ncol <- widget$ncol

	names <- widget$names
	#TODO - check all names are valid

	rowlabels <- widget$rowlabels
	rownames <- widget$rownames
	collabels <- widget$collabels
	colnames <- widget$colnames

	if (is.null(rownames)) rownames <- ""
	if (is.null(colnames)) colnames <- ""

	if (all(widget$values==""))
		values <- ""
	else {
		values <- widget$values
	}

	wid <- list(type="grid", bg=widget$bg, fg=widget$fg, borderwidth=widget$borderwidth) #new grid widget to create
	wid$.widgets <- list() #elements in the grid
	wid$byrow <- TRUE

	nNames <- length(names)
	nRowlabels <- length(rowlabels)

	nRowNames <- length(rownames)
	nValues <- length(values)
	nCollabels <- length(collabels)
	nColNames <- length(colnames)

	#count names
	if (nNames!=1 && nNames!=(ncol*nrow))
    	.stopWidget(paste('"names" argument must contain 1 or', ncol*nrow, 'names seperated by whitespace.'), widget$.debug, winName)
  
	#count rowlabels
	if (nRowlabels!=1 && nRowlabels!=0 && nRowlabels!=nrow)
		.stopWidget(paste('"rowlabels" argument should contain 1 or', nrow, 'labels.'), widget$.debug, winName)

	#count collabels
	if (nCollabels!=1 && nCollabels!=0 && nCollabels!=ncol)
		.stopWidget(paste('"collabels" argument should contain 1 or',ncol,'labels.'), widget$.debug, winName)

	#count rownames
	if (nRowNames!=1 && nRowNames!=nrow)
		.stopWidget(paste('"rownames" argument should contain 1 or',nrow,'labels.'), widget$.debug, winName)

	#count colnames
	if (nColNames!=1 && nColNames!=ncol)
		.stopWidget(paste('"colnames" argument should contain 1 or',ncol,'labels.'), widget$.debug, winName)

	if (is.null(rowlabels))
			colLabelOffset <- 0
		else
			colLabelOffset <- 1

	#single labels should be displayed as the title
	if( is.null( collabels ) ) {
		#nothing
	} else if (nCollabels==1 && ncol>1) {
		wid$toptitle<-collabels[1]
		wid$topfont<-widget$font
		wid$toptitle.offset<-1 #to help center the label
		#have counting labels above each column
		wid$.widgets[[1]] <- list()
		wid$.widgets[[1]][[1]] <- list(type='label', text="", bg=widget$bg, fg=widget$fg)
		for(j in 1:ncol) {
			wid$.widgets[[1]][[j+colLabelOffset]] <- list(type='label', text=j, font=widget$font, bg=widget$bg, fg=widget$fg)
		}
	}
	else {
		wid$.widgets[[1]] <- list()
		wid$.widgets[[1]][[1]] <- list(type='label', text="", bg=widget$bg, fg=widget$fg)
		for(j in 1:ncol) {
			wid$.widgets[[1]][[j+colLabelOffset]] <- list(type='label', text=collabels[j], font=widget$font, bg=widget$bg, fg=widget$fg)
		}
	}

	#row title
	if( is.null( rowlabels ) ) {
		wid$sidetitle.offset<-1 #to help center the label
	} else if (nRowlabels==1 && nrow>1) {
		wid$sidetitle<-rowlabels[1]
		wid$sidefont<-widget$font
		wid$sidetitle.offset<-1 #to help center the label
	}

	for(i in 1:nrow) {
		rowCount <- i #the first row of inputs should be 1 (even if there are labels ontop)
		i <- i + 1 #first row has labels
		wid$.widgets[[i]] <- list()
		for(j in 1:(ncol+1)) {
			#first row is for labels
			if (j==1) {
				if (!is.null(rowlabels)) {
					if (nRowlabels==1 && nrow>1) {
						text <- as.character(rowCount)
					}
					else
						text <- rowlabels[rowCount]

					tmp_i <- i
					if( is.null( collabels ) )
						tmp_i <- tmp_i - 1
					wid$.widgets[[tmp_i]][[j]] <- list(type='label', text=text, font=widget$font, bg=widget$bg, fg=widget$fg)
				}
			}
			else {
				if (nNames==1) #single name given
					name <- paste(names,'[',rowCount,',',j-1,']',sep="")
				else #many names given
					if (widget$byrow)
						name <- names[j-1+(ncol*(i-2))]
					else
						name <- names[i-1+(nrow*(j-2))]
				if (nValues==1) {
					value <- values
				}
				else {
					if (widget$byrow)
						value <- values[j-1+(ncol*(i-2))]
					else
						value <- values[i-1+(nrow*(j-2))]
				}

				#tweak offset if labels are disabled
				# ****must be un-tweaked after creating the entry or check widget
				if( is.null(rowlabels) )
					j <- j - 1
				if( is.null(collabels) )
					i <- i - 1

				if (widget$mode=="logical") {
					#display a checkbox
					if (is.na(as.logical(value)))
						checked=FALSE
					else if (as.logical(value)==TRUE)
						checked=TRUE
					else
						checked=FALSE
					wid$.widgets[[i]][[j]] <- list(
						  type='check',
						  mode="logical",
						  name=name,
						  text="",
						  "function"=widget[["function"]],
						  action=widget$action,
						  checked=checked,
						  bg=widget$bg,
						  fg=widget$entryfg,
						  noeditbg=widget$noeditbg,
						  noeditfg=widget$noeditfg,
						  edit=widget$edit
					)
				}
				else {
					#display a entry box
					wid$.widgets[[i]][[j]] <- list(
						  type='entry', 
						  name=name,
						  "function"=widget[["function"]],
						  action=widget$action,
						  enter=widget$enter,
						  value=value,
						  width=widget$width,
						  mode=widget$mode,
						  entryfont=widget$entryfont,
						  entrybg=widget$entrybg,
						  entryfg=widget$entryfg,
						  noeditbg=widget$noeditbg,
						  noeditfg=widget$noeditfg,
						  edit=widget$edit
					)

				}
				# ***untweak offset for special case of no labels
				if( is.null(rowlabels) )
					j <- j + 1
				if( is.null(collabels) )
					i <- i + 1
			}
		}
	}

	#look out for a trailing list (only happens if rowlabels=NULL)
	i <- length(wid$.widgets)
	if (!length(wid$.widgets[[i]]))
		wid$.widgets[[i]] <- NULL

	if (is.null(rowlabels)) {
		wid$sidetitle <- ""
		wid$toptitle.offset <- NULL
	}
	if (is.null(collabels)) {
		wid$toptitle <- ""
		wid$sidetitle.offset <- NULL
	}

	widgets <- .convertOldGridToNewGrid( wid )
	tmp <- .createWidget.grid(tk, widgets, winName)
	return(list( widget = tmp$widget, widgetList = widgetList[ -1 ] ) )
}
#-----------------------------.createWidget.matrix


#.createWidget.notebook-----------------2012-12-19
.createWidget.notebook <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]

	#create the notebook
	argList <- list( tk, "NoteBook" )

	if (!is.null(widget[["fg"]]) && widget$fg!="")
		argList$foreground=widget$fg
	if (!is.null(widget[["bg"]]) && widget$bg!="")
		argList$bg=widget$bg
	if (!is.null(widget[["arcradius"]]) && widget$arcradius!="")
		argList$arcradius=widget$arcradius
	if (!is.null(widget[["tabbevelsize"]]) && widget$tabbevelsize!="")
		argList$tabbevelsize=widget$tabbevelsize
	if (!is.null(widget[["homogeneous"]]) && widget$homogeneous!="")
		argList$homogeneous=widget$homogeneous
	if (!is.null(widget[["tabpos"]]) && widget$tabpos!="")
		argList$side=widget$tabpos
	if (!is.null(widget[["font"]]) && any(widget$font!=""))
		argList$font <- .createTkFont(widget$font)

	notebook <- .do.gui.call( "tkwidget", argList )
	if( !is.null( widget[["name"]] ) )
	.map.set( winName, widget$name, tclwidget=notebook ) # changes .PBSmod

	#create tabs
	childWidgets <- widgetList[ -1 ]
	tab_i <- 1
	for( tab in widget$tabs ) {

		.makeRaiseCmd <- function( tab )
		{
			tab <- tab #If tab isn't assigned here, the value isn't correctly copied - WTF - stupid R environments?
			#save most recently raised tab
			return( 
			function() { 
				if( !is.null( widget[["name"]] ) ) {
					#eval(parse(text=".PBSmod[[winName]]$widgets[[ paste( widget$name, \".raised\", sep=\"\" ) ]] <<- tab"))
					tget(.PBSmod)
					.PBSmod[[winName]]$widgets[[ paste( widget$name, ".raised", sep="" ) ]] <- tab
					tput(.PBSmod)
				}
				#callback
				.extractData(widget[["function"]], widget$action, winName)
			} )
		}

		tab.button <- tclvalue(tkinsert(notebook, "end", tab_i, "-text", tab, "-raisecmd", .makeRaiseCmd( tab_i ) ))
		tab.win <- .Tk.newwin(tab.button)
		tab.frame <- tkframe(tab.win)
                .adjustAllColours(tab.frame)

		#embed the next widget from the list as the *only* object in the frame (users should use a grid for multi items)
		embedded <- .createWidget( tab.frame, childWidgets, winName )
		tkpack( embedded[[ "widget" ]], padx = childWidgets[[1]]$padx, pady = childWidgets[[1]]$pady )
		tkpack( tab.frame )

		childWidgets <- embedded[[ "widgetList" ]]
		tab_i <- tab_i + 1
	}

	if( widget$selected > length( widget$tabs ) )
		.stopWidget( "given selected index is greater than the number of tabs", widget$.debug, winName )


	#raise selected tab
	tcl( notebook, "raise", widget$selected )

	#set height
	if( any( c( widget$width, widget$height ) != 0 ) && any( c( widget$width, widget$height ) == 0 ) )
		.stopWidget( "both width and height must be non-zero (to manually set the size) or both zero to automatically set the size", widget$.debug, winName)
	tcl( notebook, "configure", width = widget$width, height = widget$height )
	
	return( list( widget = notebook, widgetList = childWidgets ) )
}
#---------------------------.createWidget.notebook


#.createWidget.null---------------------2012-12-19
.createWidget.null <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	argList <- list( parent = tk, text="" )
	if( !is.null(widget[["bg"]]) && widget$bg != "" )
		argList$bg <- widget$bg
	tkWidget <- .do.gui.call( "tklabel", argList )

	return(list( widget = tkWidget, widgetList = widgetList[ -1 ] ) )
}
#-------------------------------.createWidget.null


#.createWidget.object-------------------2012-12-20
.createWidget.object <- function(tk, widgetList, winName, userObject = NULL )
{
	widget <- widgetList[[ 1 ]]
	if( is.null( userObject ) ) {
		tmp <- .check.object.exists( tk, widget, winName )
		if( !is.null( tmp ) )
			return( tmp )
	}
	
	if( !is.null( widget[["rowshow"]] ) && widget$rowshow > 0 )
		return( .createWidget.object.scrolling( tk, widgetList, winName ) )

	if( is.null( userObject ) ) {
		tget(.PBSmod)
		userObject <- get( widget$name, envir = .PBSmod[[ winName ]]$env )
	}

	#matrix
	if (is.matrix(userObject)) {
		wid <- list(type="matrix",
		            nrow=dim(userObject)[1],
		            ncol=dim(userObject)[2],
		            names=widget$name,
		            rowlabels=rownames(userObject),
		            collabels=colnames(userObject),
		            rownames=rownames(userObject),
		            colnames=colnames(userObject),
		            values=userObject,
		            byrow=FALSE,
		            font=widget$font,
		            fg=widget$fg,
		            bg=widget$bg,
		            entryfont=widget$entryfont,
		            entryfg=widget$entryfg,
		            entrybg=widget$entrybg,
		            noeditfg=widget$noeditfg,
		            noeditbg=widget$noeditbg,
		            "function"=widget[["function"]],
		            enter=widget$enter,
		            action=widget$action,
		            width=widget$width,
		            mode=mode(userObject),
		            sticky=widget$sticky,
		            borderwidth=widget$borderwidth,
		            padx=widget$padx,
		            pady=widget$pady,
		            edit=widget$edit
		            );
		if( widget$rowlabels == FALSE ) wid$rowlabels <- NULL
		if( widget$collabels == FALSE ) wid$collabels <- NULL
		#eval(parse(text=".PBSmod[[winName]]$widgets[[widget$name]] <<- wid"))
		tget(.PBSmod)
		.PBSmod[[winName]]$widgets[[widget$name]] <- wid
		tput(.PBSmod)
		tmp <- .createWidget(tk, list(wid), winName)
		return( list( widget = tmp$widget, widgetList = widgetList[ -1 ] ) )
	}

	#data.frame
	if (is.data.frame(userObject)) {
		dataModes <- c()
		dataValues <- c()
		for(i in 1:length(userObject)) {
			dataModes <- c(dataModes, mode(as.vector(userObject[[i]])))
			for(v in as.vector(userObject[[i]]))
				dataValues <- c(dataValues, v)
		}
		
		wid <- list(type="data",
		            nrow=dim(userObject)[1],
		            ncol=dim(userObject)[2],
		            names=widget$name,
		            rowlabels=rownames(userObject),
		            collabels=colnames(userObject),
		            rownames=rownames(userObject),
		            colnames=colnames(userObject),
		            values=dataValues,
		            byrow=FALSE,
		            font=widget$font,
		            fg=widget$fg,
		            bg=widget$bg,
		            entryfont=widget$entryfont,
		            entryfg=widget$entryfg,
		            entrybg=widget$entrybg,
		            noeditfg=widget$noeditfg,
		            noeditbg=widget$noeditbg,
		            "function"=widget[["function"]],
		            ".up_func"=widget[[".up_func"]],
		            ".down_func"=widget[[".down_func"]],
		            ".pageup_func"=widget[[".pageup_func"]],
		            ".pagedown_func"=widget[[".pagedown_func"]],
		            enter=widget$enter,
		            action=widget$action,
		            width=widget$width,
		            modes=dataModes,
		            sticky=widget$sticky,
		            borderwidth=widget$borderwidth,
		            padx=widget$padx,
		            pady=widget$pady,
		            edit=widget$edit,
		            .rowlabelwidth=widget$.rowlabelwidth
		            );
		if( widget$rowlabels == FALSE ) wid$rowlabels <- NULL
		if( widget$collabels == FALSE ) wid$collabels <- NULL
		#eval(parse(text=".PBSmod[[winName]]$widgets[[widget$name]] <<- wid"))
		tget(.PBSmod)
		.PBSmod[[winName]]$widgets[[widget$name]] <- wid
		tput(.PBSmod)
		tmp <- .createWidget(tk, list(wid), winName)
		return( list( widget = tmp$widget, widgetList = widgetList[ -1 ] ) )
	}

	#vector
	if (is.vector(userObject)) {
		wid <- list(type="vector",
		            names=widget$name,
		            length=length(userObject),
		            labels=names(userObject),
		            vecnames=names(userObject),
		            values=userObject,
		            font=widget$font,
		            fg=widget$fg,
		            bg=widget$bg,
		            entryfont=widget$entryfont,
		            entryfg=widget$entryfg,
		            entrybg=widget$entrybg,
		            noeditfg=widget$noeditfg,
		            noeditbg=widget$noeditbg,
		            vertical=widget$vertical,
		            "function"=widget[["function"]],
		            enter=widget$enter,
		            action=widget$action,
		            width=widget$width,
		            mode=mode(userObject),
		            sticky=widget$sticky,
		            borderwidth=widget$borderwidth,
		            padx=widget$padx,
		            pady=widget$pady,
		            edit=widget$edit
		            );
		if( widget$collabels == FALSE ) wid$labels <- NULL
		#eval(parse(text=".PBSmod[[winName]]$widgets[[widget$name]] <<- wid"))
		tget(.PBSmod)
		.PBSmod[[winName]]$widgets[[widget$name]] <- wid
		tput(.PBSmod)
		tmp <- .createWidget(tk, list(wid), winName)
		return( list( widget = tmp$widget, widgetList = widgetList[ -1 ] ) )
	}

	stop(paste("Error: variable \"", widget$name, "\" is an incompatible mode.", sep=""))
}
#-----------------------------.createWidget.object


#.createWidget.object.scrolling---------2012-12-19
.createWidget.object.scrolling <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	#check for existence
	tmp <- .check.object.exists( tk, widget, winName )
	if( !is.null( tmp ) )
		return( list( widget = tmp$widget, widgetList = widgetList[ -1 ] ) )

	widget_name <- widget$name	
	userObject <- .getValueForWidgetSetup( widget$name, widget, winName )

	if( is.matrix( userObject ) ) {
		userObject <- as.data.frame( userObject )
		#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ widget_name ]]$class <<- \"matrix\""))
		tget(.PBSmod)
		.PBSmod[[ winName ]]$widgets[[ widget_name ]]$class <- "matrix"
		tput(.PBSmod)
	} else {
		#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ widget_name ]]$class <<- \"\""))
		tget(.PBSmod)
		.PBSmod[[ winName ]]$widgets[[ widget_name ]]$class <- ""
		tput(.PBSmod)
	}
	if( !is.data.frame( userObject ) ) {
		stop( "superobjects only support data.frames" )
	}

	#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ widget_name ]]$display_top <<- 1"))
	tget(.PBSmod)
	.PBSmod[[ winName ]]$widgets[[ widget_name ]]$display_top <- 1
	tput(.PBSmod)
	rows_to_display <- widget$rowshow #num of rows visible
	enable_scrolling <- TRUE
	if( rows_to_display <= 0 || rows_to_display >= nrow( userObject ) ) {
		rows_to_display = nrow( userObject )
		enable_scrolling <- FALSE
	}
	#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ widget_name ]]$rows_to_display <<- rows_to_display"))
	tget(.PBSmod)
	.PBSmod[[ winName ]]$widgets[[ widget_name ]]$rows_to_display <- rows_to_display
	tput(.PBSmod)
	ncols <- ncol( userObject )
	nrows <- nrow( userObject )

	widget$.rowlabelwidth <- max( nchar( rownames( userObject ) ) )

	new_widget_name <- paste( "[superobject]", widget$name, sep="" ) 
	cols <- 1:ncol( userObject )
	sub_object_value <- userObject[1:rows_to_display, cols, drop = FALSE ]
	
	widget$name <- new_widget_name
	#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ widget_name ]]$.data <<- userObject"))
	tget(.PBSmod)
	.PBSmod[[ winName ]]$widgets[[ widget_name ]]$.data <- userObject
	tput(.PBSmod)
	rm( userObject )

	scroll_callback <- function( ... )
	{
		.superobject.saveValues( winName, widget_name )
		
		tget(.PBSmod)
		display_top <- .PBSmod[[ winName ]]$widgets[[ widget_name ]]$display_top
		userObject <- .PBSmod[[ winName ]]$widgets[[ widget_name ]]$.data
		
		args = list( ... )
		if( args[[1]] == "scroll" ) {
			display_top <- display_top + as.integer( args[[2]] )
		} else {
			#warning( "scroll not handled" )
			display_top = ceiling( as.numeric( args[[2]] ) * ( nrows - rows_to_display ) ) + 1
		}
		#enforce range
		if( display_top < 1 ) 
			display_top <- 1
		if( display_top + rows_to_display - 1 > nrows ) 
			display_top <- nrows - rows_to_display + 1
		#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ widget_name ]]$display_top <<- display_top"))
		#eval(parse(text=".PBSmodEnv$.PBSmod[[ winName ]]$widgets[[ widget_name ]]$display_top <<- display_top"))
		tget(.PBSmod)
		.PBSmod[[ winName ]]$widgets[[ widget_name ]]$display_top <- display_top
		tput(.PBSmod)
		
		#update row labels
		for( i in 1:rows_to_display ) {
			var_name = paste( new_widget_name, "[rowlabel][", i ,"]", sep="" )
			tmp_ptr <- .map.get( winName, var_name )
			tclvalue( tmp_ptr$tclvar ) <- rownames( userObject )[ display_top + i - 1 ]
		}

		#update row values
		for( i in 1:rows_to_display ) {
			for( j in 1:ncols ) {
				var_name = paste( new_widget_name, "[", i, ",", j ,"]d", sep="" )
				tmp_ptr <- .map.get( winName, var_name )
				tclvalue( tmp_ptr$tclvar ) <- userObject[ display_top + i - 1, j ]			
			}
		}
		
		#if using a scrollbar, use the following code to reposition scrollbar
		#update position of scroll bar
		#beg <- (display_top-1) / ( nrows - rows_to_display )
		#end <- beg
		#tkset( scroll, beg , end )
		tkconfigure( button_up, state=ifelse( display_top == 1, "disabled", "normal" ) )
		tkconfigure( button_pageup, state=ifelse( display_top == 1, "disabled", "normal" ) )

		tkconfigure( button_down, state=ifelse( display_top + rows_to_display - 1 >= nrows, "disabled", "normal" ) )
		tkconfigure( button_pagedown, state=ifelse( display_top + rows_to_display - 1 >= nrows, "disabled", "normal" ) )
	}

	#selected_widget_name is of the form "[superobject]somewidgetname[i,j]d"
	#move the tk focus up/down by y_offset (negative values up, postive down)
	set_widget_row_focus <- function( selected_widget_name, y_offset )
	{
		row <- as.integer( sub( "^\\[superobject\\].+\\[([0-9]+),[0-9]+\\]d$", "\\1", selected_widget_name ) )
		col <- as.integer( sub( "^\\[superobject\\].+\\[[0-9]+,([0-9]+)\\]d$", "\\1", selected_widget_name ) )

		#get new row to move to
		new_row <- row + y_offset
		if( new_row < 1 ) new_row <- 1
		if( new_row > rows_to_display ) new_row <- rows_to_display
		
		focus_to <- paste( "[superobject]", widget_name, "[", new_row, ",", col, "]d", sep="" )
		#tkfocus( .PBSmod[[ winName ]]$widgetPtrs[[ focus_to ]]$tclwidget )
		tget(.PBSmod)
		tkfocus( .PBSmod[[ winName ]]$widgetPtrs[[ focus_to ]]$tclwidget )
	}

	frame <- tkframe( tk )
        .adjustAllColours(frame)
	rowshow <- ceiling( widget$rowshow / 2 )
	widget$rowshow <- 0 #now we are just creating a regular object, if this was > 0, then we would get inf recursion
	widget$.up_func <- function( selected_widget_name, ...) { 
		#if( .PBSmod[[ winName ]]$widgets[[ widget_name ]]$display_top == 1 ) {
		tget(.PBSmod)
		if( .PBSmod[[ winName ]]$widgets[[ widget_name ]]$display_top == 1 ) {
			#no more hidden rows to scroll, change focus
			set_widget_row_focus( selected_widget_name, -1 )
		}
		scroll_callback( "scroll", "-1", "units" )
	}
	widget$.down_func <- function( selected_widget_name,...) { 
		#display_top <- .PBSmod[[ winName ]]$widgets[[ widget_name ]]$display_top
		tget(.PBSmod)
		display_top <- .PBSmod[[ winName ]]$widgets[[ widget_name ]]$display_top
		if( display_top + rows_to_display - 1 >= nrows ) {
			#no more hidden rows to scroll, change focus
			set_widget_row_focus( selected_widget_name, 1 )
		}
		scroll_callback( "scroll", "1", "units" )
	}
	widget$.pageup_func <- function(...) { scroll_callback( "scroll", as.character( -rowshow ), "units" ) }
	widget$.pagedown_func <- function(...) { scroll_callback( "scroll", as.character( rowshow ), "units" ) }
	obj_tk <- .createWidget.object( frame, list(widget), winName, userObject = sub_object_value )$widget
	
	#TODO make scroll grow to size of object
	#scroll <- tkscrollbar( frame, repeatinterval=5, command=function(...)scroll_callback(...))

	if( enable_scrolling == TRUE ) {
		#create left side button grid
		#scroll <- tkframe( frame )

		#two frames - one for the top two buttson, one for the bottom two
		f1 <- tkframe( frame )
		f2 <- tkframe( frame )
                # to deal with square artifacts above the buttons
                .adjustAllColours(f1)
                .adjustAllColours(f2)
                
		#keys for pageup, up, down, pagedown keys
		keys <- c( "<<", "<", ">", ">>" )
		font <- .createTkFont( "7" )

		button_pageup <- tkbutton( parent = f1, text = keys[1], font = font, width=1, command=function(...) { scroll_callback( "scroll", as.character( -rowshow ), "units" ) } )
		button_up <- tkbutton( parent = f1, text = keys[2], font = font, width=1, command=function(...) { scroll_callback( "scroll", "-1", "units" ) } )
		button_down <- tkbutton( parent = f2, text = keys[3], font = font, width=1, command=function(...) { scroll_callback( "scroll", "1", "units" ) } )
		button_pagedown <- tkbutton( parent = f2, text = keys[4], font = font, width=1, command=function(...) { scroll_callback( "scroll", as.character( rowshow ), "units" ) } )
                .adjustAllColours(button_pageup)
                .adjustAllColours(button_up)
                .adjustAllColours(button_down)
                .adjustAllColours(button_pagedown)
                
		#attach buttons to correct frame
		if( widget[["collabels"]] != FALSE ) {
                        lab <- tklabel( parent = f1, text = "" )
                        .adjustAllColours(lab)
			tkgrid( lab ) #align with rows
                }
		tkgrid( button_pageup ); tkgrid( button_up );
		tkgrid( button_down ); tkgrid( button_pagedown );

		tkgrid( obj_tk, column = 0, row = 0, rowspan = 2 )
		tkgrid( f1, column = 1, row = 0, sticky = "N" )
		tkgrid( f2, column = 1, row = 1, sticky = "S" )

		#disable up arrows (since we start on row 1)
		tkconfigure( button_up, state="disabled" )
		tkconfigure( button_pageup, state="disabled" )
	} else {
		tkgrid( obj_tk )
	}
	
	return(list( widget = frame, widgetList = widgetList[ -1 ] ) )
}
#-------------------.createWidget.object.scrolling


#.createWidget.progressBar--------------2012-12-19
.createWidget.progressbar <- function(tk, widgetList, winName)
{
	#example using ttk instead of bwidget
	#widget <- widgetList[[ 1 ]]
	#variable <- .map.add(winName, widget$name, tclvar=tclVar(widget$value))$tclvar
	#tkwidget <- tkwidget(tk, "ttk::progressbar", variable = variable )
	#return( list( widget = tkwidget, widgetList = widgetList[ -1 ] ) )

	widget <- widgetList[[ 1 ]]
	#usually, I would use tkwidget(.....), but we need to pass `type' to the widget creation (unfortunately, type is used by tkwidget too)
	win <- .Tk.subwin( tk )
	argList <- list( "ProgressBar", win )

	#used to change progress value
	if( widget$style == "incremental" )
		widget$value = widget$value / 2 #there's a but in the bwidgets/tk that doubles this value
	argList$variable <- .map.add(winName, widget$name, tclvar=tclVar( widget$value ))$tclvar

	if (!is.null(widget[["fg"]]) && widget$fg!="")
		argList$foreground=widget$fg
	if (!is.null(widget[["bg"]]) && widget$bg!="")
		argList$troughcolor=widget$bg
	if (!is.null(widget[["maximum"]]) && widget$maximum > 0)
		argList$maximum <- widget$maximum
	if (!is.null(widget[["height"]]) && widget$height > 0)
		argList$height <- widget$height 
	if (!is.null(widget[["width"]]) && widget$width > 0)
		argList$width <- widget$width 
	if (!is.null(widget[["borderwidth"]]))
		argList$borderwidth <- widget$borderwidth 
	if (!is.null(widget[["relief"]]) && widget$relief != "")
		argList$relief <- widget$relief 
	if (!is.null(widget[["vertical"]]) && widget$vertical == TRUE)
		argList$orient <- "vertical"
	if (!is.null(widget[["style"]]) && widget$style!="")
		argList$type <- widget$style

	tmp <- .do.gui.call( "tcl", argList )
	.map.set( winName, widget$name, tclwidget=tmp )

	return( list( widget = win, widgetList = widgetList[ -1 ] ) )
}
#------------------------.createWidget.progressBar


#.createWidget.radio--------------------2012-12-20
.createWidget.radio <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	argList <- list(parent=tk, text=widget$text, value=widget$value)
	if (!is.null(widget[["fg"]]) && widget$fg!="")
		argList$foreground=widget$fg
	if (!is.null(widget[["bg"]]) && widget$bg!="")
		argList$background=widget$bg
	if (!is.null(widget[["font"]]) && any(widget$font!=""))
		argList$font <- .createTkFont(widget$font)
	argList$variable<-.map.add(winName, widget$name, tclvar=tclVar(widget$value))$tclvar
	if (!is.null(widget[["selected"]]) && widget$selected==TRUE)
		tclvalue(argList$variable) <- widget$value
	argList$command=function(...) { .extractData(widget[["function"]], widget$action, winName)}

	tkWidget<-.do.gui.call("tkradiobutton", argList)
	
	#save widget pointer - radio can have many widgets for ONE varname, so store in a list
	widget_list <- .map.get( winName, widget$name )$tclwidgetlist
	if( is.null( widget_list ) )
		widget_list <- list()
	widget_list[[ as.character( widget$value ) ]] <- tkWidget
	.map.set( winName, widget$name, tclwidgetlist=widget_list )
	
	if( widget$edit == FALSE )
		tkconfigure( tkWidget, state="disabled" )

	return(list( widget = tkWidget, widgetList = widgetList[ -1 ] ) )
}
#------------------------------.createWidget.radio


#.createWidget.slide--------------------2012-12-20
.createWidget.slide <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	if (is.null(widget[["value"]]))
		widget$value <- widget$to
	argList <- list(parent=tk, from=widget$from, to=widget$to, orient=widget$orientation, showvalue=widget$showvalue)
	if (!is.null(widget[["fg"]]) && widget$fg!="")
		argList$foreground<-widget$fg
	if (!is.null(widget[["bg"]]) && widget$bg!="")
		argList$background<-widget$bg
	if (!is.null(widget[["font"]]) && any(widget$font!=""))
		argList$font <- .createTkFont(widget$font)
	argList$variable<-.map.add(winName, widget$name, tclvar=tclVar(widget$value))$tclvar
	argList$command<-function(...) { .extractData(widget[["function"]], widget$action, winName)}
	tkWidget<-.do.gui.call("tkscale", argList)
	.map.set( winName, widget$name, tclwidget=tkWidget )
	return(list( widget = tkWidget, widgetList = widgetList[ -1 ] ) )
}
#------------------------------.createWidget.slide


#.createWidget.slideplus----------------2012-12-20
.createWidget.slideplus <- function(tk, widgetList, winName)
{
        floatRegExp <- "^-?(([0-9]+\\.?[0-9]*)|([0-9]*\\.[0-9]+))$"
	widget <- widgetList[[ 1 ]]

	#initial widget$value defaults to <from> argument
	if (is.na(widget$value))
		widget$value <- widget$from

	#to remember last valid number; used to reset values upon an
	#invalid value
	lastMinVal <- ""
	lastCurVal <- ""
	lastMaxVal <- ""

	#update the slider
	updateSlideBounds <- function(slider, slideVar, curVar, minVar, maxVar, widget, winName)
	{
                #convert Tcl variables to values
		minVal <- tclvalue(minVar)
		curVal <- tclvalue(curVar)
		maxVal <- tclvalue(maxVar)

		#change min
		if (any(grep(floatRegExp,minVal))) {
                        # new value = valid: update
			tkconfigure(slider,from=as.numeric(minVal)/widget$by);
		} else {
                        # reset to last valid
                        tclvalue(minVar)<-lastMinVal 
		}

		#change max
		if (any(grep(floatRegExp,maxVal))) {
			tkconfigure(slider,to=as.numeric(maxVal)/widget$by);
		} else {
                        tclvalue(maxVar)<-lastMaxVal
		}

		#change current
		if (any(grep(floatRegExp,curVal))) {
			tclvalue(slideVar)<-round(as.numeric(curVal)/widget$by)
		} else {
			tclvalue(curVar)<-lastCurVal
		}
	}

        # with the saving below, we'll save values on key down and possibly
        # restore them on key release (if key was invalid); with such a short
        # lifetime, saving the values in .PBSmodEnv should be safe
	saveSlideBounds <- function(slider, curVar, minVar, maxVar)
	{
                # convert Tcl variables to values
		minVal <- tclvalue(minVar)
		curVal <- tclvalue(curVar)
		maxVal <- tclvalue(maxVar)

		if (any(grep(floatRegExp,minVal)))
			assign("lastMinVal",minVal,envir=.PBSmodEnv)
		if (any(grep(floatRegExp,curVal)))
			assign("lastCurVal",curVal,envir=.PBSmodEnv)
		if (any(grep(floatRegExp,maxVal)))
			assign("lastMaxVal",maxVal,envir=.PBSmodEnv)
	}

	convertCurVal <- function(widget, slideVar, curVar)
	{
		tclvalue(curVar) <- as.numeric(tclvalue(slideVar)) * widget$by
	}

        # create the widget using the specified font and fg/bg colours
	getColourfulWidget <- function( parent, type, widget, ... )
	{
		argList <- list( parent = parent, ... )
		if( !is.null( widget[["fg"]] ) && widget$fg != "" )
			argList$fg <- widget$fg
		if( !is.null( widget[["bg"]] ) && widget$bg != "" )
			argList$bg <- widget$bg
		if ( !is.null( widget[["font"]] ) && any( widget$font != "" ) )
			argList$font=.createTkFont(widget$font)
		return( .do.gui.call( type, argList ) )
	}

	#calculate fractional values
	from <- widget$from / widget$by
	to <- widget$to / widget$by

	if (is.null(widget[["value"]])) {
		value <- to
		widget$value <- widget$to
	} else {
		value <- widget$value / widget$by
	}

	slideVar<-.map.add(winName, paste(".", widget$name, ".slide", sep=""), tclvar=tclVar(value))$tclvar
	minVar<-.map.add(winName, paste(widget$name, ".min", sep=""), tclvar=tclVar(widget$from))$tclvar
	curVar<-.map.add(winName, widget$name,                        tclvar=tclVar(widget$value))$tclvar
	maxVar<-.map.add(winName, paste(widget$name, ".max", sep=""), tclvar=tclVar(widget$to))$tclvar

	#hold the widgets in this frame
	tkWidget <- getColourfulWidget( tk, "tkframe", list( bg = widget$bg ) )

        #create slider
	slider <- getColourfulWidget(tkWidget, "tkscale", widget,
                                     from=from, to=to, orient="horizontal", showvalue=FALSE, variable=slideVar,
                                     command=function(...) {
                                             # update the current value entry box
                                             convertCurVal(widget, slideVar, curVar)
                                             if (!widget$enter)
                                                     # call the user function automatically if enter = FALSE
                                                     .extractData(widget[["function"]], widget$action, winName)
                                     })

	#insert slider
	tkgrid(slider, columnspan=5, row=1, column=1)

	#create entries
	minWid <- getColourfulWidget( tkWidget, "tkentry", list( fg=widget$entryfg, bg=widget$entrybg, font=widget$entryfont ), textvariable=minVar, width=5 )
	curWid <- getColourfulWidget( tkWidget, "tkentry", list( fg=widget$entryfg, bg=widget$entrybg, font=widget$entryfont ), textvariable=curVar, width=5 )
	maxWid <- getColourfulWidget( tkWidget, "tkentry", list( fg=widget$entryfg, bg=widget$entrybg, font=widget$entryfont ), textvariable=maxVar, width=5 )

	if (widget$enter) {
                #capture value after return/enter is released; evaluate expressions like...
                #  tkbind(minWid,"<KeyRelease-KP_Enter>",function() updateSlideBounds(slider, slideVar, curVar, minVar, maxVar, widget, winName))
                for (wid in c("minWid", "curWid", "maxWid")) {
                        for (key in c("<KeyRelease-Return>", "<KeyRelease-KP_Enter>")) {
                                expr <- paste('tkbind(', wid, ', "', key, '", ',
                                              'function() { updateSlideBounds(slider, slideVar, curVar, minVar, maxVar, widget, winName); .extractData(widget[["function"]], widget$action, winName) })',
                                              sep="")
                                eval(parse(text=expr))
                        }
                }
	} else {
                #capture value before press (in case new value isn't
                #valid); then capture value after key is received
                for (wid in c("minWid", "curWid", "maxWid")) {
                        expr <- paste('tkbind(', wid, ', "<KeyPress>",function() { saveSlideBounds(slider, curVar, minVar, maxVar) }); ',
                                      'tkbind(', wid, ', "<KeyRelease>",function() { updateSlideBounds(slider, slideVar, curVar, minVar, maxVar, widget, winName); .extractData(widget[["function"]], widget$action, winName) })',
                                      sep ="")
                        eval(parse(text=expr))
                }
        }

	#bind functions for setWinVal() changes
	.map.set(winName, paste(widget$name, ".min", sep=""), onChange=function() {
                updateSlideBounds(slider, slideVar, curVar, minVar, maxVar, widget, winName)
                if (!widget$enter) .extractData(widget[["function"]], widget$action, winName) })
	.map.set(winName, widget$name, onChange=function() {
                updateSlideBounds(slider, slideVar, curVar, minVar, maxVar, widget, winName)
                if (!widget$enter) .extractData(widget[["function"]], widget$action, winName) })
	.map.set(winName, paste(widget$name, ".max", sep=""), onChange=function() {
                updateSlideBounds(slider, slideVar, curVar, minVar, maxVar, widget, winName)
                if (!widget$enter) .extractData(widget[["function"]], widget$action, winName) })

	#place widgets in grid
	tkgrid(getColourfulWidget( tkWidget, "tklabel", widget, text="Min->"), row=2, column=1)
	tkgrid(minWid, row=2, column=2)
	tkgrid(curWid, row=2, column=3)
	tkgrid(maxWid, row=2, column=4)
	tkgrid(getColourfulWidget( tkWidget, "tklabel", widget, text="<-Max"), row=2, column=5)

	return(list( widget = tkWidget, widgetList = widgetList[ -1 ] ) )
}
#--------------------------.createWidget.slideplus


#.createWidget.spinbox------------------2012-12-20
.createWidget.spinbox <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	if (!is.null(widget[["label"]]))
	if (widget$label!="") {
		#if label is set, then create a 2x1 grid
		label <- widget$label
		widget$label <- "" #blank it out, inf loop if not.
		widgets <- list(
			list(type="grid", nrow=1, ncol=2, font="", byrow=TRUE, borderwidth=1, relief="flat", padx=0, pady=0, fg=widget$fg, bg=widget$bg ),
			list(type="label", text=label, padx=0, pady=0, font=widget$font, fg=widget$fg, bg=widget$bg),
			widget
			)
		tmp <- .createWidget.grid(tk, widgets, winName)
		return(list( widget = tmp$widget, widgetList = widgetList[ -1 ] ) )
	}
	
	#create real tk widget spinbox below
	argList <- list(parent=tk, type="SpinBox", editable=TRUE, range=c( widget$from, widget$to, widget$by ) )
	if (!is.null(widget[["entryfg"]]) && widget$entryfg!="") {
		argList$foreground=widget$entryfg
		argList$selectforeground=widget$entryfg
		argList$entryfg=widget$entryfg
	}
	if (!is.null(widget[["entrybg"]]) && widget$entrybg!="") {
		argList$background=widget$entrybg
		argList$insertbackground=widget$entrybg
		argList$selectbackground=widget$entrybg
		argList$entrybg=widget$entrybg
		
	}
	if (!is.null(widget$entryfont) && any(widget$entryfont!=""))
		argList$font <- .createTkFont(widget$entryfont)
		
	if( is.na( widget$value ) )
		widget$value = widget$from
		

	argList$textvariable<-.map.add(winName, widget$name, tclvar=tclVar(widget$value))$tclvar
	argList$width<-widget$width
	argList$validate = "all"
	
	tkWidget<-.do.gui.call("tkwidget", argList)
	.map.set( winName, widget$name, tclwidget=tkWidget )

	#bug in spinbox -> up/down arrows still modify value in disabled mode
	if( widget$edit == FALSE )
		tkconfigure( tkWidget, state="disabled" )

	enter <- !is.null(widget[["enter"]])
	if (enter)
		enter <- widget$enter
        
	if (enter) {
		# don't update it unless Return/Enter was pressed (updating can slow it down a lot)
                keylist <- c("<KeyRelease-Return>", "<KeyRelease-KP_Enter>")
	}
	else {
                # update on every key release
                keylist <- c("<KeyRelease>")
                # update when otherwise modified (clicking the up/down arrows)
                tkconfigure( tkWidget, modifycmd=function(...) {
                        .extractData(widget[["function"]], widget$action, winName)
                })
        }
        # NOTE: tkbind does not work with SpinBox (from BWidget); instead, we need to use
        # Tcl directly with 'pathName bind ?arg...?' to set bindings for this widget
        for (key in keylist) {
                cmd <- paste(tkWidget$ID, 'bind ', key,
                             .Tcl.args(function(...) {
                                     .extractData(widget[["function"]],
                                                  widget$action, winName)}),
                             sep=" ")
                .Tcl(cmd)
        }
	return(list( widget = tkWidget, widgetList = widgetList[ -1 ] ) )
}
#----------------------------.createWidget.spinbox

#.createWidget.table--------------------2012-12-19
.createWidget.table <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]

	tmp <- .check.object.exists( tk, widget, winName )
	if( !is.null( tmp ) )
		return( list( widget = tmp$widget, widgetList = widgetList[ -1 ] ) )
	
	userObject <- .getValueForWidgetSetup( widget$name, widget, winName )
	if( !is.matrix( userObject ) && !is.data.frame( userObject ) )
		stop( "table only supports matrix" )
		
	#convert any factors to characters
	for( i in 1:ncol(userObject) ) {
		userObject[[ i ]] <- as.vector( userObject[[ i ]] ) #will convert factors into characters
	}

	widget_name <- widget$name

	tget(.PBSmod)
	#to help us getWinVal the correct size/mode/names
	#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ widget_name ]]$.dim <<- dim( userObject )"))
	.PBSmod[[ winName ]]$widgets[[ widget_name ]]$.dim <- dim( userObject )
	modes <- c()
	for( i in 1:ncol( userObject ) ) modes <- c( modes, mode( userObject[[ i ]] ) )
	#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ widget_name ]]$.modes <<- modes"))
	.PBSmod[[ winName ]]$widgets[[ widget_name ]]$.modes <- modes
	#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ widget_name ]]$.dimnames <<- dimnames( userObject )"))
	.PBSmod[[ winName ]]$widgets[[ widget_name ]]$.dimnames <- dimnames( userObject )
	#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ widget_name ]]$.class <<- ifelse( is.matrix( userObject ), \"matrix\", \"data.frame\" )"))
	.PBSmod[[ winName ]]$widgets[[ widget_name ]]$.class <- ifelse( is.matrix( userObject ), "matrix", "data.frame" )
	tput(.PBSmod)
	nrows <- nrow( userObject )
	ncols <- ncol( userObject )
	table_nrows <- nrows
	table_ncols <- ncols

	frame <- tkframe( tk )
        .adjustAllColours(frame)
	
	#create tcl storage for matrix
	tcl_array <- .map.add( winName, widget$name, tclarray=tclArray() )$tclarray
	.map.set( winName, widget$name, widgetname = widget$name )
	
	#array also holds labels - if they are wanted, then put labels in [[0,*]] and [[*,0]]
	#and start elements at [[1,*]] and [[*,1]]. If labels aren't required, start data at [[0,*]] and [[*,0]]
	row_label_offset <- ifelse( is.null( widget[[ "collabels" ]] ), 1, 0 )[1]
	col_label_offset <- ifelse( is.null( widget[[ "rowlabels" ]] ), 1, 0 )[1]
	show_rowtitle <- ifelse( is.null( widget[[ "rowlabels" ]] ), 0, 1 )[1]
	show_coltitle <- ifelse( is.null( widget[[ "collabels" ]] ), 0, 1 )[1]
	
	#fill in titles
	if( !is.null( widget$rowlabels ) ) {
		for (i in (1:nrows)) {
			if( all( widget$rowlabels == "" ) )
				tcl_array[[i-row_label_offset,0]] <- rownames( userObject )[i]
			else
				tcl_array[[i-row_label_offset,0]] <- widget$rowlabels[i]
		}
		#reserve left col for row titles
		table_ncols <- table_ncols + 1
	}
	if( widget$collabels[1]!="NULL" ) {
		for (i in (1:ncols)) {
			if( all( widget$collabels == "" ) )
				tcl_array[[0,i-col_label_offset]] <- colnames( userObject )[i]
			else
				tcl_array[[0,i-col_label_offset]] <- widget$collabels[i]
		}
		#reserve top row for column titles
		table_nrows <- table_nrows + 1
	}
	
	#fill in data
	for (i in (1:nrows))
		for (j in (1:ncols)) {
			tcl_array[[i-row_label_offset,j-col_label_offset]] <- userObject[i,j] #tcl arrays start at 0
		}
	
	#example from: http://bioinf.wehi.edu.au/~wettenhall/RTclTkExamples/tktable.html
	argList <- list( parent = frame, type = "table", rows=table_nrows,cols=table_ncols,titlerows=show_coltitle,titlecols=show_rowtitle,
	             height=widget$rowshow,width=10,multiline=0,
	             xscrollcommand=function(...) tkset(xscr,...),yscrollcommand=function(...) tkset(yscr,...))

	#if( length( widget$width ) > 1 )
	#	argList$width <- 

	if (!is.null(widget[["font"]]) && any(widget$font!=""))
		argList$font <- .createTkFont(widget$font)

	

	table1 <- .do.gui.call( "tkwidget", argList )
	.map.set( winName, widget$name, tclwidget=table1 )
	             
	#bug with tktable for "moveto" scrollbar scrolling - doesn't hit last element - must push down on keyboard, or click down arrow
	xscr <- tkscrollbar( frame,orient="horizontal", command=function(...)tkxview(table1,...))
	yscr <- tkscrollbar( frame,command=function(...)tkyview(table1,...))
        .adjustAllColours(xscr)
        .adjustAllColours(yscr)

	tkgrid(table1,yscr)
	tkgrid.configure(yscr,sticky="nsw")
	tkgrid(xscr,sticky="new")
	tkconfigure(table1,variable=tcl_array,background=widget$bg,foreground=widget$fg,selectmode="extended")
	
	tmp_width <- rep( widget$width, times = table_ncols )
	for( i in 1:table_ncols )
		tcl( table1, "width", i-1, tmp_width[ i ] )
	
	#callback func
	tkconfigure(table1, browsecmd=function(...) { .extractData(widget[["function"]], widget$action, winName)} )
	
	
	
	if( widget$edit == FALSE )
		tkconfigure( table1, state="disabled" )
	
	
	return(list( widget = frame, widgetList = widgetList[ -1 ] ) )
}
#------------------------------.createWidget.table


#.createWidget.text---------------------2012-12-20
.createWidget.text <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	tk <- tkframe(tk)
        .adjustAllColours(tk)

	param <- list(
	              parent=tk, 
	              bg=widget$bg, 
	              fg=widget$fg, 
	              height=widget$height, 
	              width=widget$width,
	              relief=widget$relief,
	              yscrollcommand=function(...)tkset(scrollBar,...)
	              )
	if ( any( widget$font != "" ) )
		param$font=.createTkFont(widget$font)

	scrollBar <- tkscrollbar(tk, repeatinterval=5, command=function(...)tkyview(txtBox,...))
        .adjustAllColours(scrollBar)
	txtBox <- .do.gui.call("tktext", param)

	.map.add(winName, widget$name, tclwidget=txtBox)
	tkinsert(txtBox,"end",widget$value)

	if (widget$edit==FALSE)
		tkconfigure(txtBox, state="disabled")

	if (widget$scrollbar == TRUE) {
		tkgrid(txtBox,scrollBar)
		tkgrid.configure(scrollBar,sticky="ns")
	} else {
		tkgrid(txtBox) }

	return(list( widget = tk, widgetList = widgetList[ -1 ] ) )
}
#-------------------------------.createWidget.text


#.createWidget.vector-------------------2012-12-19
.createWidget.vector <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	names <- widget$names
	labels <- widget$labels

	if (all(widget$values==""))
		values <- ""
	else
		values <- widget$values

	if( is.null( widget[["vecnames"]] ) ) widget$vecnames <- ""

	n <- widget$length
	nNames <- length(names)
	nLabels <- length(labels)
	nVecnames <- length(widget$vecnames)

	if (n==0) {
		if (nNames == 1 && nLabels == 1) {
			n<-1
		}
		else if (nNames != 1)
			n<-nNames
		else
			.stopWidget('missing "length" argument', widget$.debug, winName)
	}

	#count names
	if (nNames!=1 && nNames!=n)
		.stopWidget(paste("names argument must contain 1 or",n,"names seperated by spaces.\nreceived:", widget$names), widget$.debug, winName)

	#count vecnames
	if (nVecnames!=n && widget$vecnames!="")
		.stopWidget(paste('vecnames argument should contain',n,'vector names.'), widget$.debug, winName)

	#count labels
	if( !is.null( labels ) ) {
		if( labels[1] != "" || nLabels > 1 )
			labels <- rep( labels, length = n )
		else if( nNames > 1 || n == 1 )
			labels <- rep( names, length = n )
		else
			labels <- 1:n
	}

	#build grid
	wid <- list(type="grid", bg=widget$bg, fg=widget$fg, borderwidth=widget$borderwidth ) #new grid widget to create
	if (widget$vertical) {
		wid$nrow <- n
		wid$ncol <- ifelse( is.null( labels ), 1, 2 )
		wid$toptitle.offset=1
		wid$toptitle=""
		wid$byrow = TRUE
	}
	else {
		wid$nrow <- ifelse( is.null( labels ), 1, 2 )
		wid$ncol <- n
		wid$byrow = FALSE
	}

	nValues <- length(values)
	if (nValues!=1 && nValues!=n)
		.stopWidget(paste('values argument should contain 1 or',n,'values seperated by whitespace.'), widget$.debug, winName)

	#create list for grid, and children
	widgets <- list( wid )
	for(i in 1:n) {
		#create label
		if( is.null( labels ) ) {
			entryIndex <- 1
		} else {
			entryIndex <- 2
			text <- labels[i]
			widgets[[ length(widgets)+1]] <- list(type='label', text=text, font=widget$font, fg=widget$fg, bg=widget$bg)
		}

		#create entry
		if (nNames==1)
			name <- paste(names, '[', i, ']',sep="")
		else
			name <- names[i]

		if (nValues==1)
			value <- values[1]
		else
			value <- values[i]

		if (all(widget$vecnames==""))
			vname <- NULL
		else
			vname <- widget$vecnames[i]

		if (widget$mode=="logical") {
			#display a checkbox
			if (is.na(as.logical(value)))
				checked=FALSE
			else if (as.logical(value)==TRUE)
				checked=TRUE
			else
				checked=FALSE
			widgets[[ length(widgets)+1 ]] <- list(
				  type='check',
				  mode="logical",
				  name=name,
				  text="",
				  "function"=widget[["function"]],
				  action=widget$action,
				  checked=checked,
				  bg=widget$bg,
				  fg=widget$entryfg,
				  noeditbg=widget$noeditbg,
				  noeditfg=widget$noeditfg,
				  edit=widget$edit,
				  .name=vname
			)
		}
		else {
			#display a entry box
			widgets[[ length(widgets)+1 ]] <- list(
				  type='entry', 
				  name=name,
				  "function"=widget[["function"]],
				  action=widget$action,
				  enter=widget$enter,
				  value=value,
				  width=widget$width,
				  mode=widget$mode,
				  entryfont=widget$entryfont,
				  entrybg=widget$entrybg,
				  entryfg=widget$entryfg,
				  noeditbg=widget$noeditbg,
				  noeditfg=widget$noeditfg,
				  edit=widget$edit,
				  .name=vname
			)

		}

	}
	tmp <- .createWidget.grid(tk, widgets, winName)
	return(list( widget = tmp$widget, widgetList = widgetList[ -1 ] ) )

#	wid$ncol <- NULL
#	wid$nrow <- NULL
#	widgets <- .convertOldGridToNewGrid( wid )
#	tmp <- .createWidget.grid(tk, widgets, winName)
#	return(list( widget = tmp$widget, widgetList = widgetList[ -1 ] ) )
}
#-----------------------------.createWidget.vector


#.createTkFont--------------------------2012-12-20
#   creates a usable TK font from a given string
# Arguments:
#   fontStr - string describing a font and colour
#----------------------------------------------ACB
.createTkFont <- function( fontstr )
{
	#print( fontStr )
	#fontstr <- .convertParamStrToVector(casefold(fontStr))
	#print( fontstr )

	#default options
	fontparam<-list()

	for(i in 1:length(fontstr)) {
		if (fontstr[i]=="bold")
			fontparam$weight="bold"
		else if (fontstr[i]=="italic")
			fontparam$slant="italic"
		else if (fontstr[i]=="underline")
			fontparam$underline=TRUE
		else if (fontstr[i]=="overstrike")
			fontparam$overstrike=TRUE
		else if (fontstr[i]=="times")
			fontparam$family="Times"
		else if (fontstr[i]=="courier")
			fontparam$family="Courier"
		else if (fontstr[i]=="helvetica")
			fontparam$family="Helvetica"
		else if (any(grep("^[0-9]+$", fontstr[i])))
			fontparam$size=fontstr[i]
		else if (fontstr[i] != "") {
			fontparam$family = fontstr[i]
			#cat(paste("warning: font familly \"", fontstr[i], "\" is not guarenteed to work with TK on all platforms\n", sep=""))
		}
	}
	return(.do.gui.call(tkfont.create, fontparam))
}
#------------------------------------.createTkFont


#.buildgrid-----------------------------2012-12-20
#   used to create a grid on a window
# Arguments:
#   tk      - parent tk frame to attach widget to
#   grid    - widget list describing the grid
#   winName - active window name
#----------------------------------------------ACB
.buildgrid <- function(tk, grid, winName, childWidgets)
{
	toptitle <- grid$toptitle
	sidetitle <- grid$sidetitle

	if (is.null(toptitle))
		toptitle <- ""

	if (is.null(sidetitle))
		sidetitle <- ""

	if (is.null(grid[["ncol"]])) {
		stop( "ncol not given to grid" )
		grid$ncol=length(grid$.widgets[[1]])
	}

	if (is.null(grid[["nrow"]])) {
		stop( "nrow not given to grid" )
		grid$nrow=length(grid$.widgets)
	}

	#offset the title (useful for centering titles over a certain part)
	#like over the 3 columns of a matrix, but not row labels
	if (is.null(grid[["toptitle.offset"]]))
		grid$toptitle.offset<-0
	if (is.null(grid[["sidetitle.offset"]]))
		grid$sidetitle.offset<-0

	#set byrow
	if (is.null(grid[["byrow"]])) {
		grid$byrow=TRUE
	}

	#set font options
	if( is.null(grid[["topfont"]]) )
		topfont = ""
	else
		topfont = grid$topfont
	if( is.null(grid[["sidefont"]]) )
		sidefont = ""
	else
		sidefont = grid$sidefont
	
	#set fg
	if( is.null(grid[["topfg"]]) ) grid$topfg <- grid$fg
	if( is.null(grid[["sidefg"]]) ) grid$sidefg <- grid$fg
	#set bg
	if( is.null(grid[["topbg"]]) ) grid$topbg <- grid$bg
	if( is.null(grid[["sidebg"]]) ) grid$sidebg <- grid$bg

	#display title (if set)
	if (toptitle!="") {
		colspan=as.integer(grid$ncol)-grid$toptitle.offset
		argList <- list(parent=tk, text=toptitle)
		if (!is.null(grid[["topfg"]]) && grid$topfg!="") {
			argList$foreground <- grid$topfg
		}
		if (!is.null(grid[["topbg"]]) && grid$topbg!="")
			argList$background <- grid$topbg
		if (any( topfont!="" ) )
			argList$font <- .createTkFont(topfont)
		mytklabel<-.do.gui.call("tklabel", argList)
		tkgrid(mytklabel, columnspan=colspan, row=0, column=1+grid$toptitle.offset)
	}
	
	#display column title (if set)
	if (sidetitle!="") {
		rowspan=as.integer(grid$nrow)-grid$sidetitle.offset
		argList <- list(parent=tk, text=sidetitle)
		if (!is.null(grid[["sidefg"]]) && grid$sidefg!="")
			argList$foreground=grid$sidefg
		if (!is.null(grid[["sidebg"]]) && grid$sidebg!="")
			argList$background=grid$sidebg
		if (any(topfont!=""))
			argList$font <- .createTkFont(sidefont)
		mytklabel<-.do.gui.call("tklabel", argList)
		tkgrid(mytklabel, rowspan=rowspan, row=1+grid$sidetitle.offset, column=0)
		showsidetitle<-TRUE
	}

	#loop over all children widget of the grid.
	#these are stored as grid$.widgets[[row_id]][[col_id]]
#	for(i in 1:grid$nrow) {
#		for(j in 1:grid$ncol) {
	for( i in 0:(grid$nrow * grid$ncol - 1) ) {
		#create Widget (here's the recursive widget creation call)
		if( length( childWidgets ) == 0 )
			.stopWidget( "not enough children widgets to build grid", grid$.debug, winName )
		widget_def <- childWidgets[[ 1 ]]
		tmp <- .createWidget(tk, childWidgets, winName)
		widget <- tmp[[ "widget" ]]
		childWidgets <- tmp[[ "widgetList" ]]

		#set row and column position
		if (grid$byrow==TRUE) {
			row = floor( i / grid$ncol ) + 1
			column = i %% grid$ncol + 1
		} else {
			column = floor( i / grid$nrow ) + 1
			row = i %% grid$nrow + 1
		}

		#y padding
		if (is.null(widget_def$pady))
			pady <- 0
		else
			pady <- widget_def$pady

		#x padding
		if (is.null(widget_def$padx))
			padx <- 0
		else
			padx <- widget_def$padx

		#build begining of argument list for tkgrid() function
		argList <- list(widget, row=row, column=column, padx=padx, pady=pady)

		#append sticky flag argument if set
		if (is.character(widget_def$sticky)) {
			argList$sticky <- widget_def$sticky
		}
		.do.gui.call("tkgrid", argList)
	}
	return( list( widget = tk, widgetList = childWidgets ) )
}
#---------------------------------------.buildgrid


#===== Map Functions ============================
#----- (.map, etc.) ------------------------------


#.map.add-------------------------------2012-12-19
#   save a new value for a given key.
#   if a previous exists ignore the new value, and return previous value
# Arguments:
#   winName - map to extract values from
#   key     - name of item to extract (i.e. widget name)
#   ...     - named items to save (in a list)
#----------------------------------------------ACB
.map.add <- function(winName, key, ...)
{
	#if( is.null( .PBSmod[[ winName ]] ) )
	tget(.PBSmod)
	if( is.null(.PBSmod[[ winName ]]) )
		.map.init(winName)
	if (!is.character(key)) {
		stop("map error - key must be a string")
	}
	if (key=="") {
		stop("map error - key must be at least 1 character long")
	}
	tget(.PBSmod)
	if (!is.null(.PBSmod[[winName]]$widgetPtrs[[ key ]]))
		return(.PBSmod[[winName]]$widgetPtrs[[key]])
	#eval(parse(text=".PBSmod[[winName]]$widgetPtrs[[key]] <<- list(...)"))

	### The following caused days of confusion as to why values were not being displayed in the GUI:
	#.PBSmod[[winName]]$widgetPtrs[[key]] <- list(...)   # This was the orginal final line in the function.
	#tput(.PBSmod)                                       # This additional line disrupted the return of the final line.
	#return(list(...))                                   # Therefore need an explicit return.

	#Alex's suggestion (2012-12-07)
	newlist <- list(...)
	.PBSmod[[winName]]$widgetPtrs[[key]] <- newlist
	tput(.PBSmod)
	return(newlist)

	#.PBSmod[[winName]]$widgetPtrs[[key]] <- list(...)   # This was the orginal final line in the function.
	#tput(.PBSmod)                                       # This additional line disrupted the return of the final line.
	#return(list(...))                                   # Therefore need an explicit return.
}
#-----------------------------------------.map.add


#.map.get-------------------------------2012-12-19
#   Returns a value associated with a key
# Arguments:
#   winName - map to extract values from
#   key     - name of item to extract (i.e. widget name)
#----------------------------------------------ACB
.map.get <- function(winName, key)
{
	#return(.PBSmod[[winName]]$widgetPtrs[[key]])
	tget(.PBSmod)
	return(.PBSmod[[winName]]$widgetPtrs[[key]])
}
#-----------------------------------------.map.get


#.map.getAll----------------------------2012-12-19
#   Returns all visible items of a map of a certain window
# Arguments:
#   winName - map to extract values from
#----------------------------------------------ACB
.map.getAll <- function(winName)
{
	#return(.PBSmod[[winName]]$widgetPtrs)
	tget(.PBSmod)
	return(.PBSmod[[winName]]$widgetPtrs)
}
#--------------------------------------.map.getAll


#.map.init------------------------------2012-12-19
#   initialize the datastructure that holds the map(s)
# Arguments:
#   winName - name of map to initialize
#----------------------------------------------ACB
.map.init <- function(winName) {
	tget(.PBSmod)
	.PBSmod[[winName]] <- list()
	#to hold tclvar pointers
	.PBSmod[[winName]]$widgetPtrs <- list()
	#to hold widget definition lists (i.e. from win desc file)
	.PBSmod[[winName]]$widgets <- list()
	tput(.PBSmod)
	#packList(winName,".PBSmod",list()) #.PBSmod[[winName]] <<- list()
	#to hold tclvar pointers
	#target=paste(".PBSmod[[\"",winName,"\"]]",sep="")
	#eval(parse(text=paste("packList(\"widgetPtrs\",",target,",list())",sep=""))) #.PBSmod[[winName]]$widgetPtrs <<- list()
	#to hold widget definition lists (i.e. from win desc file)
	#eval(parse(text=paste("packList(\"widgets\",",target,",list())",sep=""))) #.PBSmod[[winName]]$widgets <<- list()
	invisible()
}
#----------------------------------------.map.init


#.map.set-------------------------------2012-12-19
#   save a new value for a given key, even if it involves
#   overwriting a previously stored value
# Arguments:
#   winName - map to extract values from
#   key     - name of item to extract (i.e. widget name)
#   ...     - named items to save (in a list)
#----------------------------------------------ACB
.map.set <- function(winName, key, ...)
{
	#if( is.null( .PBSmod[[ winName ]] ) )
	tget(.PBSmod)
	if( is.null( .PBSmod[[ winName ]] ) )
		.map.init(winName)

	if (!is.character(key))
		stop("map error - key must be a string")
	if (key=="")
		stop("map error - key must be at least 1 character long")

	tget(.PBSmod)
	if (!is.list(.PBSmod[[winName]]$widgetPtrs[[key]]))
		.PBSmod[[winName]]$widgetPtrs[[key]] <- list()
		#eval(parse(text=".PBSmod[[winName]]$widgetPtrs[[key]] <<- list()"))

	#set additional keys
	tmp <- list(...)
	tmpNames <- names(tmp)
	if (length(tmp)>0) {
		for (i in 1:length(tmp)) {
			if (is.null(tmpNames[i]))
				#eval(parse(text=".PBSmod[[winName]]$widgetPtrs[[key]][[i]] <<- tmp[[i]]"))
				.PBSmod[[winName]]$widgetPtrs[[key]][[i]] <- tmp[[i]]
			else if (tmpNames[i]=="")
				#eval(parse(text=".PBSmod[[winName]]$widgetPtrs[[key]][[i]] <<- tmp[[i]]"))
				.PBSmod[[winName]]$widgetPtrs[[key]][[i]] <- tmp[[i]]
			else
				#eval(parse(text=".PBSmod[[winName]]$widgetPtrs[[key]][[tmpNames[i]]] <<- tmp[[i]]"))
				.PBSmod[[winName]]$widgetPtrs[[key]][[tmpNames[i]]] <- tmp[[i]]
		}
	}
	tput(.PBSmod)
	return(.PBSmod[[winName]]$widgetPtrs[[key]])
}
#-----------------------------------------.map.set


#===== Convert Functions =============================
#----- (.convert, etc.) ------------------------------


#.convertOldGridToNewGrid---------------2012-12-20
#This is use to convert old style (where grid has a grid$.widgets[[i]][[j]] data structure)
#to the new format which is list( grid, child, child, ..., child )
#this should *never* be used by new code -- hopefuly it will be removed one day
#----------------------------------------------ACB
.convertOldGridToNewGrid <- function( grid )
{
	widgets <- list( grid )
	widgets[[ 1 ]]$.widgets <- NULL #remove old widgets

	if( is.null( widgets[[1]][[ "nrow" ]] ) )
		widgets[[1]][[ "nrow" ]] = length( grid$.widgets )
	if( is.null( widgets[[1]][[ "ncol" ]] ) )
		widgets[[1]][[ "ncol" ]] = length( grid$.widgets[[ 1 ]] )

	for( row in grid$.widgets ) {
		if( length( row ) != widgets[[1]][["ncol"]] ) stop( "missmatching ncol vs data in .widgets during old to new conversion" )
		for( item in row ) {
			if( item$type == "grid" ) {
				tmp <- .convertOldGridToNewGrid( item )
				widgets <- c( widgets, tmp )
			} else {
				widgets[[ length( widgets ) + 1 ]] <- item
			}
		}
	}
	return( widgets )
}
#-------------------------.convertOldGridToNewGrid


#.convertParamStrToVector---------------2012-12-20
#  Function to convert a string, x, into a vector of 
#  elements seperated by whitespace.
#  Whitespace can be interupted as values if it is enclosed by quotes.
#  special characters (newline, tab, \, ', ") must be escaped with \
# Arguments:
#   x     - string
#   fname - filename string for warning messages
#   line  - line number for warning messages
# Output:  
#   vector of values
#----------------------------------------------ACB
.convertParamStrToVector <- function(x, fname="", line=0)
{
	#TODO - this function needs optimization for readlist preformance
	#profilling shows a slow spot to be paste() calls
	escape<-0
	quoted<-0	#0=none, 1=", 2='
	word<-""
	j <- 0 #counter for words array
	equal <- 0 #counter for equal char
	quotefound <- 0 #used to capture ""
	words <- NULL
	#todo: triple check it's not needed - 
	# I dont think so since stripcomments is used prior
	#x<-.trimWhiteSpace(x)

	myNewEnv <- new.env()
	environment(myNewEnv) <- asNamespace("PBSmodelling")
	words <- .Call("strToVector", 
	               x, 
	               myNewEnv,
	               fname,
	               line,
	               PACKAGE="PBSmodelling")
	if (is.null(words))
		stop("Errors were found in the P-formated list. Unable to continue\n", call.=FALSE)
	return(words)
}
#-------------------------.convertParamStrToVector


#.convertParamStrToList-----------------2012-12-20
#  Converts a given string of values seperated by spaces into a list
#  while preserving space and escaped quotes within quotes 
#  (kindof - the value must still be stripped depending if its a single string, or vector of strings)
#----------------------------------------------ACB
.convertParamStrToList <- function(x, fname="", line.start=0, line.end=0, sourcefile=list())
{
	escape<-0
	quoted<-0	#0=none, 1=", 2='
	word<-""
	j <- 0 #counter for words array
	equal <- 0 #counter for equal char
	quotefound <- 0 #used to capture ""
	words <- list()
	#don't think its needed since stripcomments should do this too
	#x<-.trimWhiteSpace(x)

	myNewEnv <- new.env()
	environment(myNewEnv) <- asNamespace("PBSmodelling")
	words <- .Call("strToList", 
	               x, 
	               myNewEnv,
	               fname,
	               line.start,
	               PACKAGE="PBSmodelling")
	if (is.null(words))
		stop("Errors were found in the GUI description file. Unable to continue\n", call.=FALSE)
	return(words)
}
#---------------------------.convertParamStrToList


#.convertMatrixListToMatrix-------------2012-12-20
#  Converts a list into an N-dim array
# Arguments:
#   mList = z[[1]][[1]]...[[1]]=x
#           z[[1]][[1]]...[[2]]=x
#           ...
#           z[[1]][[1]]...[[1]]=x
#           ...
#           z[[i]][[j]]...[[k]]=x
#
# output an N-dim array
#----------------------------------------------ACB
.convertMatrixListToMatrix <- function(mList)
{
	size <- .getMatrixListSize(mList)
	arr <- array(dim=size)
	arr <- .setMatrixElement(mList, arr)
	return(arr)
}
#-----------------------.convertMatrixListToMatrix


#.convertMatrixListToDataFrame----------2012-12-20
#  Similar to toArray but to data.frame
# Arguments:
#   mList - see .convertMatrixListToMatrix:
#----------------------------------------------ACB
.convertMatrixListToDataFrame <- function(mList, colName="Y", rowNames="X")
{
	size <- .getMatrixListSize(mList)
	arr <- array(dim=size)
	arr <- .setMatrixElement(mList, arr)

	x<-list()
	for(i in 1:size[2]) {
		x[[i]]<-list()
	}

	for(i in 1:length(mList)) {
		for(j in 1:length(mList[[i]])) {
			x[[j]][[i]] <- mList[[i]][[j]]
		}
	}

	if (length(rowNames)==0) {
		rowNames=NULL
	}
	else {
		if (length(rowNames)==1) {
			if (rowNames=="") {
				rowNames=NULL
			}
			rowNames <- paste(rowNames, 1:size[1], sep="")
		}
		else if (length(rowNames)!=size[1])
			stop(paste("rowNames should be NULL, or a vector of size 1 or", size[1], ".\nGot rowNames=", rowNames, sep=""))
	}

	#create a data.frame
	argList <- list(row.names=rowNames)
	for(i in 1:size[2]) { #foreach column
		name <- paste("X", i, sep="")
		argList[[name]] <- unlist(x[[i]])
	}
	return(do.call("data.frame", argList))
}
#--------------------.convertMatrixListToDataFrame


#.setMatrixElement----------------------2012-12-20
#  Helper function used by .convertMatrixListToMatrix
#  to assign values from the matrix list into the array
#----------------------------------------------ACB
.setMatrixElement <- function(m, a, ind=NULL)
{
	if (is.null(m))
		return(a)
	if (!is.list(m)) {
		eval(parse(text=paste("a[", paste(ind, collapse=','), "] <- m", sep="")))
		return(a)
	}

	for(i in 1:length(m)) {
		a<-.setMatrixElement(m[[i]], a, c(ind,i))
	}
	return(a)
}
#--------------------------------.setMatrixElement


#.getMatrixListSize---------------------2012-12-20
#  Helper function used by .convertMatrixListToMatrix
#  to determine the minumum required size of the array
#  needed to create to convert the list into an array
#----------------------------------------------ACB
.getMatrixListSize <- function(m, d=NULL, big=0)
{
	if (!is.list(m)) {
		return(pmax(d, big))
	}

	for(i in 1:length(m)) {
		big <- .getMatrixListSize(m[[i]], c(d,i), big)
	}
	return(big)
}
#-------------------------------.getMatrixListSize


#.matrixHelp----------------------------2012-12-20
#  Used to help .extractVar deal with N-dim maticies
#  firstly it is converted into a "matrix list"
#  once the matrix list is completed (and size known)
#  it should be converted into a true array
#----------------------------------------------ACB
.matrixHelp <- function(matrixList, ind, value)
{
	if (length(ind)>1) {
		if (length(matrixList)<ind[1])
			matrixList[[ind[1]]]<-list()
		else if(!is.list(matrixList[[ind[1]]]))
			matrixList[[ind[1]]]<-list()

		matrixList[[ind[1]]] <- .matrixHelp(matrixList[[ind[1]]], ind[-1], value)
		return(matrixList)
	}
	else if(length(ind)==1) {
		matrixList[[ind[1]]]<-value
		return(matrixList)
	}
	else {
		stop(".matrixHelp() was called with no indices.")
	}
}
#--------------------------------------.matrixHelp


#.autoConvertMode-----------------------2012-12-20
#  Cconverts x into a numeric mode, if it looks like a valid number
# Arguments:
#   x - string to convert
#----------------------------------------------ACB
.autoConvertMode <- function(x)
{
	#nice regular expression to see if it could be logical
	if (length(grep(.regex.logical, x))==length(x)) {
		x <- as.logical(x)
	}
	#ugly regular expression to see if it looks like a numeric
	else if (length(grep(.regex.numeric, x))==length(x) 
	         && all(x!="-")) {

		x <- as.numeric(x)
	}
	#uglier regular expression to see if its complex
	else if (length(grep(.regex.complex, x))==length(x)
	         && all(x!="-")) {

		x <- as.complex(x)
	}
	return(x)
}
#---------------------------------.autoConvertMode


#.convertMode---------------------------2012-12-20
# This funcion is deprecated (see code in supportFuns).
# Converts a variable into a mode without showing any warnings.
# Arguments:
#   x    - variable to convert
#   mode - mode to convert to
#----------------------------------------------ACB
.convertMode <- function(x, mode)
{
	#cat( "converting: <<", x, ">> to", mode,"\n")
	#TODO - fix regexs
	#they dont slow down, but will mess up with NAs

	if (mode=="logical") {
		# "1" will be TRUE
		x[x=="1"] <- TRUE
		x[x=="0"] <- FALSE
		# anything else that fails the regex will be NA
		tmp <- grep(.regex.logical, x)
		x[-tmp] <- NA
		#anything else is now valid
		x <- as.logical(x)
	}
	else if (mode=="numeric" || mode=="integer") {
		tmp <- grep(.regex.numeric, x)
		x[-tmp] <- NA
		x[x=="-"] <- NA
		x <- as.numeric(x)
	}
	else if (mode=="complex") {
		tmp <- grep(.regex.complex, x)
		x[-tmp] <- NA
		x[x=="-"] <- NA
		x <- as.complex(x)
	}
	#print("convert mode ended at"); print(date());
	#otherwise it is character and nothing needs converting
	return(x)
}
#-------------------------------------.convertMode


#===== Extract Functions =========================
#----- (.extract, etc.) --------------------------


#.extractData---------------------------2012-12-20
#  Called by TK on button presses (or binds, onchanges, slides, ...)
# Arguments:
#   command - user command to call (ie function argument of widget)
#   action  - action value
#   winName - active window name
#----------------------------------------------ACB
.extractData <- function(command, action, winName)
{
	tget(.PBSmod)
	.PBSmod$.activeWin <- winName
	tput(.PBSmod)
	#packList(".activeWin",".PBSmod",winName) # default tenv=.PBSmodEnv #.PBSmod$.activeWin <<- winName
	setWinAct(winName, action)
	if (is.null(command))
		return()
	if (command=="")
		return()
	#if (exists(command,mode="function", envir = .PBSmod[[ winName ]]$env))
	#	.do.gui.call(command, list(), envir = .PBSmod[[ winName ]]$env )
	tget(.PBSmod)
	if (exists(command,mode="function", envir = .PBSmod[[ winName ]]$env))
		.do.gui.call(command, list(), envir = .PBSmod[[ winName ]]$env )
	else
		cat(paste("Warning: cannot find function '", command, "'.\n", sep=""))
}
#-------------------------------------.extractData


#.extractFuns---------------------------2012-12-20
#  Get a list of called functions.
# Arguments:
#   data - widget lists
#----------------------------------------------ACB
.extractFuns <- function(data)
{
	retData <- c()
	for(i in 1:length(data)) {
		if (!is.null(data[[i]][["widget"]][["function"]]))
			retData <- c(retData, data[[i]][["widget"]][["function"]])
	}
	return(retData)
}
#-------------------------------------.extractFuns


#.extractVar----------------------------2012-12-20
#  Extracts values from the tclvar ptrs
# Arguments:
#   winName - name of target window to extract data from
#----------------------------------------------ACB
.extractVar <- function(winName)
{
	# list of regular expressions of keys which should NOT be returned to the user
	keys_to_skip <- c( "\\[rowlabel\\]\\[[0-9]+\\]$" )

	#data is a list containing sub-lists in the form:
	#list(type="tcl", tclvar="tcl_var_ptr", mode="numeric")
	data <- .map.getAll(winName)
	#tget(.PBSmod) perhaps best not to use static copy; call object directly with .PBSmodEnv$.PBSmod

	values <- list()
	keys <- names(data)
	#if (length(data)<1) return(NULL)

	superobjects_to_process <- list() #superobjects (which scroll) which are stored elsewhere
	tables_to_process <- list() #tktable matrix objects - stored in tclarray

	#extract values from tcl into an R list whose index corresponds to the data list
	if( length(data) > 0 )
	for(i in 1:length(data)) {
		skip <- FALSE
		for( ignore_pattern in keys_to_skip ) {
			if( any( grep( ignore_pattern, keys[i] ) ) ) {
				skip <- TRUE
				break
			}
		}
		if( skip == TRUE ) next
		if( any( grep( "^\\[superobject\\]",  keys[i] ) ) ) {
			#skip superobjects
			#get superobject name
			
			tmp <- gsub( "]d$", "]", keys[i] )
			tmp <- gsub( "\\[[a-z0-9,]*\\]", "", tmp )
			superobjects_to_process[[ tmp ]] = TRUE
			next
		}
		#wid <- .PBSmod[[winName]]$widgets[[keys[i]]]
		tget(.PBSmod)
		wid <- .PBSmod[[winName]]$widgets[[keys[i]]]
		if( is.null( wid ) || is.null( wid[["type"]] ) )
			next
		if( wid$type == "button" )
			next #no data to extract
		if (!is.null(data[[i]][[ "tclvar" ]])) {
			values[[i]] <- tclvalue(data[[i]]$tclvar)
		}
		else if (!is.null(data[[i]][["tclarray"]])) {
			#special case for table matrix
			tables_to_process[[ data[[i]]$widgetname ]] = TRUE
			next
		}
		else if( !is.null( wid[["type"]] ) && wid$type == "text" ) {
			#special case for text widgets
			values[[i]] <- tclvalue(tkget(data[[i]]$tclwidget,"0.0","end"))
			wid$mode <- "character"
		}
		else if( !is.null( data[[i]][[ "droplist_widget" ]] ) ) {
			#nothing to extract for this one - only here for setting vars
			next
		} else if( !is.null( data[[i]][[ "tclwidget" ]] ) ) {
			#this isn't extracted, but rather set from a bind function whenever a tab is raised (see createwidget.notebook - and a few lines down in this func)
			next
		} else {
			stop(paste("unknown type:", data[[i]]))
		}

		#convert data to propper type
		if (is.null(wid[["mode"]]))
			mode <- "numeric"
		else
			mode <- wid[["mode"]]

		values[[i]] <- .convertMode(values[[i]], mode)
	}

	retData <- list()
	if( length( values ) )
	for(i in 1:length(values)) {
		#look for any vectors (arrays, matrices)
		#vector names end in [1,4,2...]
		if (any(grep("^[^\\[]+\\[([0-9,]+)\\]$", keys[i]))) {
			#extract the indicies (ind) and name of vector
			ind<-gsub("^[^\\[]+\\[([0-9,]+)\\]$", "\\1", keys[i])
			ind<-as.numeric(unlist(strsplit(ind, ",")))
			name <- gsub("\\[[0-9,]+\\]", "", keys[i])
			if (length(ind)>1) {
				#process multiple idicies (matrix, or array)

				#values from matricies are stored into a list 
				#and then converted into a matrix at a later stage
				if (!exists("matrixTmp"))
					matrixTmp <- list()

				#create a list for the new matrix
				if( is.null( matrixTmp[[ name ]] ) ) {
					matrixTmp[[name]] <- list()
				}

				#call matrixhelper to build a list, and then save the new changes
				matrixTmp[[name]] <- .matrixHelp(matrixTmp[[name]], ind, values[[i]])
			}
			else {
				#single index found (vector)
				if ( is.null( retData[[ name ]] ) )
					retData[[name]] <- NA
				retData[[name]][ind] <- values[[i]]
			}
		}
		#any var ending with indicies and a "d" EX: var[3,5]d is an element of a data.frame
		else if (any(grep("^[^\\[]+\\[([0-9,]+)\\]d$", keys[i]))) {
			ind<-gsub("^[^\\[]+\\[([0-9,]+)\\]d$", "\\1", keys[i])
			ind<-as.numeric(unlist(strsplit(ind, ",")))
			name <- gsub("\\[[0-9,]+\\]d", "", keys[i])
			if (length(ind)>1) {
				#store into a list just like we do with a matrix
				if (!exists("dataframeTmp"))
					dataframeTmp <- list()

				#create a list for the new matrix
				if ( is.null( dataframeTmp[[ name ]] ) ) {
					dataframeTmp[[name]] <- list()
				}

				#call matrixhelper to build a list, and then save the new changes
				dataframeTmp[[name]] <- .matrixHelp(dataframeTmp[[name]], ind, values[[i]])
			}
			else {
				#single index
				retData[[name]][ind] <- values[[i]]
			}
		}
		else {
			#no index (ie var is of standard type: [a-z0-9_]+)
			retData[[keys[i]]] <- values[[i]]
			#if (!is.null(.PBSmod[[ winName ]]$widgets[[ keys[i] ]][[ ".name" ]]))
			#	names(retData[[keys[i]]]) <- .PBSmod[[winName]]$widgets[[keys[i]]]$.name
			tget(.PBSmod)
			if (!is.null(.PBSmod[[ winName ]]$widgets[[ keys[i] ]][[ ".name" ]]))
				names(retData[[keys[i]]]) <- .PBSmod[[winName]]$widgets[[keys[i]]]$.name
		}
	}

	#convert all collected matrix lists into real n-dim arrays.
	if (exists("matrixTmp")) {
		keys <- names(matrixTmp)
		for(i in 1:length(matrixTmp)) {
			#colnames <- .PBSmod[[winName]]$widgets[[keys[i]]]$colnames
			tget(.PBSmod)
			colnames <- .PBSmod[[winName]]$widgets[[keys[i]]]$colnames
			if (is.null(colnames)) 
				colnames <- ""

			#rownames <- .PBSmod[[winName]]$widgets[[keys[i]]]$rownames
			rownames <- .PBSmod[[winName]]$widgets[[keys[i]]]$rownames
			if (is.null(rownames)) 
				rownames <- ""

			retData[[keys[i]]] <- .convertMatrixListToMatrix(matrixTmp[[i]])
			#can't use dimnames incase of 3 dimension or higher arrays
			tmpNames <- .PBSdimnameHelper(rownames, colnames, dim(retData[[keys[i]]]))
			#base::rownames <- does not work, use the following work-around
			retData[[keys[i]]] <- base::"rownames<-"(retData[[keys[i]]], tmpNames[[1]])
			retData[[keys[i]]] <- base::"colnames<-"(retData[[keys[i]]], tmpNames[[2]])
		}
	}

	#convert dataframe lists into dataframes
	if (exists("dataframeTmp")) {
		keys <- names(dataframeTmp)
		for(i in 1:length(dataframeTmp)) {
			#colnames <- .PBSmod[[winName]]$widgets[[keys[i]]]$colnames
			#rownames <- .PBSmod[[winName]]$widgets[[keys[i]]]$rownames
			tget(.PBSmod)
			colnames <- .PBSmod[[winName]]$widgets[[keys[i]]]$colnames
			rownames <- .PBSmod[[winName]]$widgets[[keys[i]]]$rownames

			retData[[keys[i]]] <- .convertMatrixListToDataFrame(dataframeTmp[[i]], colnames, rownames)
			#can't use dimnames incase of 3 dimension or higher arrays
			tmpNames <- .PBSdimnameHelper(rownames, colnames, dim(retData[[keys[i]]]))
			#base::rownames <- does not work, use the following work-around
			retData[[keys[i]]] <- base::"rownames<-"(retData[[keys[i]]], tmpNames[[1]])
			retData[[keys[i]]] <- base::"colnames<-"(retData[[keys[i]]], tmpNames[[2]])
		}
	}

	# look for special widgets (which don't use a tclvar)
	# assign vecnames to any vectors
	# droplist widget - get position of selected item (and possibly the complete set of possible choices)
	# for(wid in .PBSmod[[winName]]$widgets) {
	tget(.PBSmod)
	for(wid in .PBSmod[[winName]]$widgets) {
		if( !is.list( wid ) ) next
		if( is.null( wid[["type"]] ) ) next
		if (wid$type=="vector") {
			if (!is.null(retData[[ wid$names ]])) {
				if (any(wid$vecnames!=""))
					names(retData[[wid$names]]) <- wid$vecnames
			}
		}
		else if (wid$type=="droplist") {
			wid_name <- paste( wid$name, ".id", sep="" )
			
			#get widget ptr (not widget's variable ptr)
			tk_widget <- data[[ wid$name ]]$tclwidget
			
			#get selected index (not value)
			selected_i <- as.integer( tcl( tk_widget, "getvalue" ) ) + 1
			retData[[wid_name]] <- selected_i

			#extract values directly from widget (won't work if labels are used)
			#values <- tclvalue( tcl( tk_widget, "cget", "-values" ) )
			#values <- .tclArrayToVector( values )

			#get stored values (useful if labels were applied)
			wid_name <- paste( wid$name, ".values", sep="" )
			#values <- .PBSmod[[winName]]$widgets[[ wid_name ]]$labels
			values <- .PBSmod[[winName]]$widgets[[ wid_name ]]$labels
			retData[[wid_name]] <- values 

			#overwrite label values with real values (except when the user input their own choice)
			if( selected_i > 0 )
				retData[[ wid$name ]] <- values[ selected_i ]
		} else if( wid$type == "notebook" ) {
			if( !is.null( wid[["name"]] ) )
				#retData[[ wid$name ]] <- .PBSmod[[winName]]$widgets[[ paste( wid$name, ".raised", sep="" ) ]]
				retData[[ wid$name ]] <- .PBSmod[[winName]]$widgets[[ paste( wid$name, ".raised", sep="" ) ]]
		}
	}
	
	#superobjects to process
	for( k in names( superobjects_to_process ) ) {
		.superobject.saveValues( winName, k ) #save currently visible data
		#retData[[ k ]] <- .PBSmod[[ winName ]]$widgets[[ k ]]$.data
		tget(.PBSmod)
		retData[[ k ]] <- .PBSmod[[ winName ]]$widgets[[ k ]]$.data
		#if( .PBSmod[[ winName ]]$widgets[[ k ]]$class == "matrix" )
		if( .PBSmod[[ winName ]]$widgets[[ k ]]$class == "matrix" )
			retData[[ k ]] <- as.matrix( retData[[ k ]] )
	}
	
	#tables to process
	
	for( k in names( tables_to_process ) ) {
		retData[[ k ]] <- .table.getvalue( winName, k )
	}

	#convert factors to characters
	if( length( retData ) > 0 )
	for( i in 1:length( retData ) ) {
		if( is.factor( retData[[ i ]] ) )
			retData[[ i ]] <- as.character( retData[[ i ]] )
		else if( is.data.frame( retData[[ i ]] ) ) {
			for( j in 1:length( retData[[ i ]] ) ) {
				if( is.factor( retData[[ i ]][,j] ) )
					retData[[ i ]][,j] <- as.character( retData[[ i ]][,j] )
			}
		}
	}
	return(retData)
}
#--------------------------------------.extractVar


#.getParamFromStr-----------------------2012-12-20
#  Returns a list with all parameters extracted from a list.
# Arguments:
#   inputStr   - string from win desc file describing a widget
#   fname      - filename to display with error messages
#   line.start - line number where widget is first found
#   line.end   - line number of last line of widget (ie extended line)
#   sourcefile - 
#   paramOrder
#----------------------------------------------ACB
.getParamFromStr <- function(inputStr, fname="", line.start=0, line.end=0, 
                             sourcefile=list(), paramOrder=.widgetDefs)
{
	#now passed in function - to enable overriding
	# a "constant" defines how the parameters should look.

	namedArguments <- 0 # determines if we can accept unnamed arguments
	paramData <- list() #extracted params from a line to return
	typeDefined <- 0    #this must be set before returning
	if (inputStr=="") {
		.catError(paste("input line is empty", sep=""), fname, line.start)
		return(NULL) }
	#split input string into seperate arguments
	s<-.convertParamStrToList(inputStr, fname, line.start, line.end, sourcefile)
	for(j in 1:length(s)) {
		value<-s[[j]]$value
		quoted <- s[[j]]$quoted == "Y"
		#argument is named, unnamed arguments are no longer valid
		if (!is.null(s[[j]][["key"]])) {
			namedArguments <- 1 
			key<-casefold(s[[j]]$key)
			if (typeDefined==0) {
				if (key=="type") {
					typeDefined <- 1;
					#case of type is ignored
					paramData[[key]] <- value <- casefold(value, upper=FALSE) 
					#fetch argument Ordering
					argOrder <- paramOrder[[value]]
				}
			} else {
				#special case where unquoted NULL is converted to a real NULL, but can't work for ALL values because of NULL widget type.
				if( !quoted && value == "NULL" ) {
					paramData[ key ] <- list( NULL )
				} else {
					paramData[[ key ]] <- value
				}
			}
		}
		#argument is not named (no key was given)
		else if(namedArguments==0) {
			if (j==1) {
				#first argument must be "type"
				widgetType <- paramData$type <- casefold(value)
				#fetch argument Ordering
				argOrder <- paramOrder[[widgetType]]
				if (is.null(argOrder)) {
					#given widget type is not valid
					.catError(paste("unknown widget type '", paramData$type,"'", sep=""), fname, line.start, line.end, sourcefile)
					return(NULL)
				}
				typeDefined <- 1;
			}
			else if(j > length(argOrder)) {
				.catError(paste("more arguments given than supported by widget '",paramData$type,"'", sep=""), fname, line.start, line.end, sourcefile)
				return(NULL)
			}
			else {
				#determine the name of the second argument
				argName <- argOrder[[j]]$param
				#special case where unquoted NULL is converted to a real NULL, but can't work for ALL values because of NULL widget type.
				if( !quoted && value == "NULL" ) {
					paramData[ argName ] <- list( NULL )
				} else {
					paramData[[ argName ]] <- value
				}
			}
		}
		#error - unnamed arg given after named arg
		else {
			.catError(paste("unnamed argument with value \"",value,"\" given after a named argument.", sep=""), fname, line.start, line.end, sourcefile)
			return(NULL)
		}
	}
	#test if a type has been defined
	if (typeDefined==0) {
		.catError(paste("no widget type given", sep=""), fname, line.start, line.end, sourcefile)
		return(NULL) }
	#check that widget type is valid
	if( is.null( paramOrder[[ paramData[["type"]] ]] ) ) {
		.catError(paste("unknown widget type '", paramData$type,"'", sep=""), fname, line.start, line.end, sourcefile)
		return(NULL) }
	#test if all given arguments are valid arguments of the given type
	#and then convert from character to another type if applicable
	errorFound <- 0
	givenNames <- names(paramData)
	for(i in 1:length(givenNames)) {
		pos <- .searchCollection(argOrder, givenNames[i])
		if (pos==-1) {
			#argument name not valid
			.catError(paste('argument \'', givenNames[i], '\' is not supported by widget \'', paramData$type, '\'', sep=""), fname, line.start, line.end, sourcefile)
			errorFound <- 1
		}
		else if (pos == -2) {
			.catError(paste('argument \'', givenNames[i], '\' of widget \'', paramData$type, '\' matches multiple formal arguments.', sep=""), fname, line.start, line.end, sourcefile)
			errorFound <- 1
		}
		else {
			#sometimes, only a bit of the paramater name is given
			#if so, we would like to know the full name
			fullParamName <- argOrder[[pos]]$param
			if (fullParamName != givenNames[i]) {
				names(paramData)[i] <- fullParamName
			}
			#check supplied argument data matches grep pattern (if defined)
			if (!is.null(argOrder[[pos]]$grep) && !all( is.null( paramData[[i]] ) ) ) {
				if (!any(grep(argOrder[[pos]]$grep, paramData[[i]]))) {
					#supplied data is not formatted correctly
					.catError(paste('argument \'', givenNames[i], '\': value \'', paramData[[i]], '\' is not accepted. It should match', argOrder[[pos]]$grep , sep=""), 	fname, line.start, line.end, sourcefile)
					errorFound <- 1
				}
			}
			#check argument for user-specified NULL options (represented by an NA)
			if( all( is.null( paramData[[i]] ) ) ) {
				if ( is.null( argOrder[[pos]][["allow_null"]] ) || argOrder[[pos]]$allow_null == FALSE ) {
					#NULL is not valid for this option
					.catError(paste('argument \'', givenNames[i], '\': NULL is not accepted. To specify the word use "NULL"', sep=""), fname, line.start, line.end, sourcefile)
					errorFound <- 1
				}
			}
			#some strings - if character need to be stripped of slashes, or sub-divided
			if (errorFound == 0 && !is.null(argOrder[[pos]][["class"]]) && !is.null(paramData[[i]]) ) {
				if (argOrder[[pos]]$class=="character") {
					#convert value to either a single string (.stripSlashes)
					tmp <- .stripSlashes(paramData[[i]], fname, line.start, line.end, sourcefile)
					if (is.null(tmp)) {
						#most likely an unescaped quote was found - .catError called from .stripSlashes, not here
						errorFound <- 1
					}
					else
						paramData[[i]] <- tmp
				}
				else if (argOrder[[pos]]$class=="characterVector") {
					#convert value to a vector of strings
					tmp <- .stripSlashesVec(paramData[[i]], fname, line.start, line.end, sourcefile)
					if (is.null(tmp)) {
						errorFound <- 1
					}
					else
						paramData[[i]] <- tmp
				}
				else if (argOrder[[pos]]$class=="integerVector") {
					#convert value to a vector of strings
					tmp <- .stripSlashesVec(paramData[[i]], fname, line.start, line.end, sourcefile)
					if (is.null(tmp)) {
						errorFound <- 1
					}
					else
						paramData[[i]] <- as.integer( tmp )
				}
				else if ((argOrder[[pos]]$class=="numeric" || argOrder[[pos]]$class=="integer") && any(grep("[a-z]", paramData[[i]], ignore.case=TRUE))) {
					paramData[[i]] <- 0
				}
				else if (argOrder[[pos]]$class=="logical" && !any(grep("^(F|T|FALSE|TRUE|false|true)$", paramData[[i]]))) {
					paramData[[i]] <- FALSE
				}
				else {
					paramData[[i]] <- as(paramData[[i]], argOrder[[pos]]$class)
				}
			}
		}
	}
	if (errorFound != 0)
		return(NULL)
	#check that all required arguments have been supplied, and fill in any missing values with defaults
	for(i in 1:length(argOrder)) {
		if (argOrder[[i]]$required) {
			if( is.null( paramData[[ argOrder[[i]][["param"]] ]] ) ) {
				.catError(paste('required argument \'', argOrder[[i]]$param, '\' is missing from widget \'', paramData$type, '\'', sep=""), fname, line.start, line.end, sourcefile)
				errorFound <- 1
			}
		}
		else if (!is.null(argOrder[[i]][["default"]])) {
			#if argument name (from widgetdefs) != any names set in paramData, then no value was given by the user.
			#note: can't use is.null since a user can set a value to equal null (which will set a arg name)
			if( all( argOrder[[ i ]][[ "param" ]] != names( paramData ) ) ) {
				#set it to default value
				paramData[[argOrder[[i]]$param]] <- argOrder[[i]]$default
			}
		}
	}
	if (errorFound != 0)
		return(NULL)

	#convert default values from character class to specified class ($mode in widgetDefs.r) - but why is this only specific to radio?
	if ( !is.null( paramData[[ "value" ]] ) && !is.null( paramData[[ "mode" ]] ) && paramData$type=="radio" ) {
		paramData$value <- as(paramData$value, paramData$mode)
	}
	#store debug information to be used if there are any errors while building the GUI
	sourceCode <- c()
	if (line.start<=line.end && !missing(sourcefile)) {
		for(i in line.start:line.end) {
			sourceCode <- c(sourceCode, sourcefile[[i]])
		}
	}
	
	paramData$.debug <- list(sourceCode=sourceCode, fname=fname, line.start=line.start, line.end=line.end)
	return(paramData)
}
#---------------------------------.getParamFromStr


#.table.getvalue------------------------2012-12-20
.table.getvalue <- function( winName, widgetName )
{
	tget(.PBSmod)
	widget <- .PBSmod[[ winName ]]$widgets[[ widgetName ]]
	tcl_array <- .map.get( winName, widgetName)$tclarray
	nrows <- widget$.dim[1]
	ncols <- widget$.dim[2]
	
	mat <- list() #data.frame( nrow = nrows, ncol = ncols, dimnames = widget$.dimnames )
	
	row_label_offset <- ifelse( is.null( widget[["collabels"]] ), 1, 0 )
	col_label_offset <- ifelse( is.null( widget[["rowlabels"]] ), 1, 0 )
	
	#extract data
	for (j in (1:ncols)) {
		mat[[ j ]] <- vector( widget$.modes[j], length=nrows )
		for (i in (1:nrows)) {
			tmp <- tcl_array[[i-row_label_offset,j-col_label_offset]]
			mat[[ j ]][ i ] <- as( tclvalue( tmp ), widget$.modes[j] )
		}
	}
	
	tmp <- data.frame( mat )
	dimnames( tmp ) <- widget$.dimnames
	if( widget$.class == "matrix" )
		tmp <- as.matrix( tmp )
	return( tmp )
}
#----------------------------------.table.getvalue


#.table.setvalue------------------------2012-12-20
.table.setvalue <- function( winName, widgetName, value )
{
	tget(.PBSmod)
	widget <- .PBSmod[[ winName ]]$widgets[[ widgetName ]]
	tcl_array <- .map.get( winName, widgetName)$tclarray
	nrows <- widget$.dim[1]
	ncols <- widget$.dim[2]
	
	if( any( dim( value ) != widget$.dim ) ) stop( "value for table is not the correct dimension" )
	
	row_label_offset <- ifelse( is.null( widget[["collabels"]] ), 1, 0 )
	col_label_offset <- ifelse( is.null( widget[["rowlabels"]] ), 1, 0 )
	
	#extract data
	for (i in (1:nrows))
		for (j in (1:ncols)) {
			tcl_array[[i-row_label_offset,j-col_label_offset]] <- value[ i, j ]
		}
	return( value )
}
#----------------------------------.table.setvalue


#.getValueForWidgetSetup----------------2012-12-20
# Returns the value(s) to use for initializing a widget
# -> from widget$value (if set)
# -> or global var matching widget$name when value is NULL
# widget: widget list passed to .creatWidget.XXX
# winName: name of window being created
#----------------------------------------------ACB
.getValueForWidgetSetup <- function( varname, widget, winName )
{
	tget(.PBSmod)
	if( !exists( varname, envir = .PBSmod[[ winName ]]$env ) )
		.stopWidget( paste( "unable to find variable \"", varname, "\" in global memory - this search happend since value=NULL", sep="" ), widget$.debug, winName )
	var <- get( varname, envir = .PBSmod[[ winName ]]$env )
	if( is.factor( var ) )
		var <- as.character( var )
	else if( is.data.frame( var ) ) {
		for( i in 1:length( var ) ) {
			if( is.factor( var[,i] ) )
				var[,i] <- as.character( var[,i] )
		}
	}
	return( var )
}
#--------------------------.getValueForWidgetSetup


#===== Helper Functions ==========================
#----- (.parse/.strip/.search/.validate,etc.) ----

#.getSimilarColour-------------------------2013-02-06
# obtains a similar colour (slightly brighter/darker) than the argument
# args:  col - the colour to adjust
#----------------------------------------------NMB
.getSimilarColour <- function (col) {
        # convert the current colour to an RGB vector
        col <- col2rgb(col)[, 1]

        # slightly lighten the colour unless it's already white, in which case
        # we'll slightly darken it
        col <- rgb(ifelse(col[1] == 255, col[1] - 1, col[1] + 1),
                   ifelse(col[2] == 255, col[2] - 1, col[2] + 1),
                   ifelse(col[3] == 255, col[3] - 1, col[3] + 1),
                   maxColorValue = 255)

        return (col)
}
#-----------------------------------.getSimilarColour

#.adjustAllColours----------------------------2013-02-06
# for the passed Tk widget, attempts to slightly modify all applicable colours
# args:  tObj - Tk object
#----------------------------------------------NMB
.adjustAllColours <- function(tObj)
{
        # obtain the configuration for the Tk widget
        config <- as.character(tkconfigure(tObj))

        # property list obtained from http://www.tcl.tk/man/tcl/TkCmd/palette.htm
        for (attr in c("activebackground", "activeforeground", "background",
                       "disabledforeground", "foreground", "highlightbackground",
                       "highlightcolor", "insertbackground", "selectbackground",
                       "selectcolor", "selectforeground", "troughcolor")) {
                # if the attribute exists in the configuration...
                if (length(grep(paste("^-", attr, " ", sep=""), config,
                                ignore.case = TRUE)) > 0) {
                        # ... let's set it to a similar colour
                        expr <- paste(
                                "tkconfigure(tObj, ", attr,
                                '=.getSimilarColour(as.character(tkconfigure(tObj,',
                                ' paste("-", attr, sep="")))[5]))', sep="")
                        eval(parse(text=expr))
                }
        }
}
#--------------------------------------.adjustAllColours

#.do.gui.call------------------------------2013-02-07
# extends do.call, which is used to create most Tk widgets, to
# immediately adjust a new widget's colours (to get around how
# tk_setPalette changes the colour of existing widgets)
#-------------------------------------------NMB
.do.gui.call <- function (what, args, quote = FALSE, envir = parent.frame()) 
{
        rval <- do.call (what, args, quote, envir)

        # if it's likely that we've just created a widget, try setting
        # the colours for it
        if (is.character(what)
            && is.element(what, c("tcl", "tkbutton", "tkcheckbutton", "tkentry",
                                  "tkframe", "tklabel", "tkradiobutton",
                                  "tkscale", "tktext", "tkwidget"))) {
                .adjustAllColours(rval)
        }

        return (rval)
}

#.catError------------------------------2012-12-20
#  Used to display parsing errors.
# Arguments:
#   err        - error string to display
#   fname      - file name where error was found
#   line.start - starting line of widget with error
#   line.end   - end line of widget with error
#   sourcefile - source code of the file in question
#   errorType  - type of error to display
#----------------------------------------------ACB
.catError <- function(err, fname, line.start, line.end, sourcefile=list(), errorType="GUI parse error")
{
	err <- paste(errorType, " (", fname, ":", line.start, ") : ", err, "\n", sep="")
	cat(err)
	if (length(sourcefile)>0) {
		for(i in line.start:line.end) {
			cat(paste(i, ": ", sourcefile[[i]], "\n", sep=""))
		}
		cat("\n")
	}
}
#----------------------------------------.catError


#.catError2-----------------------------2012-12-20
.catError2 <- function(err, fname, line.start)
{
	err <- paste("Parse error", " (", fname, ":", line.start, ") : ", err, "\n", sep="")
	cat(err)
}
#---------------------------------------.catError2


#.check.object.exists-------------------2012-12-20
#  Call to test for existence of dynamically loaded 
#  object for "object", "table". Returns NULL on 
#  no errors; on errors, a tk display error is 
#  created which can be embeded in the window
#----------------------------------------------ACB
.check.object.exists <- function( tk, widget, winName )
{
	.dispError <- function(errorTxt)
	{
		wid <- list(type="label", 
		            text=errorTxt,
		            bg="white",
		            fg="red",
		            font="bold"
		            )
		return(.createWidget(tk, list(wid), winName))
	}

	#if (!exists(widget$name, envir = .PBSmod[[ winName ]]$env )) {
	tget(.PBSmod)
	if (!exists(widget$name, envir = .PBSmod[[ winName ]]$env )) {
		return(.dispError(paste("Error: variable \"", widget$name, "\" could not be found.", sep="")))
	}
	return( NULL )
}
#-----------------------------.check.object.exists


#.inCollection--------------------------2012-12-20
#   returns true if needle occurs in haystack
# Input: 
#   haystack - a vector to search
#   needle   - a single element to search for
#----------------------------------------------ACB
.inCollection <- function(haystack, needle)
{
	if (is.null(haystack)) {
		return(FALSE)
	}
	if (is.vector(haystack)) {
		return(any(haystack==needle))
	}
	stop("only vectors are supported")
	return(FALSE)
}
#------------------------------------.inCollection


#.packWidgetsIntoGrid-------------------2012-12-20
#  Pack all widgets into a grid with ncol=1 nrow=<number of widgets>
.packWidgetsIntoGrid <- function( widgets, vertical = TRUE ) {
	gridWidget = list(
	             type="grid", 
	             font="",
	             borderwidth=0,
	             relief="flat",
	             padx=0,
	             pady=0,
	             .widgets=list(),
	             nrow = length( widgets ),
	             ncol = 1,
	             byrow = vertical
	)
	return( c( list( gridWidget ), widgets ) )
#	#add all widgets to this grid, each in a new row [[j]] and 1st column [[1]]
#	for(j in 1:length(widgets)) {
#		gridWidget$.widgets[[j]] <- list()
#		gridWidget$.widgets[[j]][[1]] <- widgets[[j]]
#	}
#	return( gridWidget )
}
#-----------------------------.packWidgetsIntoGrid


#.parsegrid-----------------------------2012-12-20
#  Returns two items in a list:
#   - $gridData which is a list of lists representing columns
#   - $unparsedData - which is left over from the grid and 
#     still needs parsing
# Arguments:
#   data - list of widget lists
#   nRow - num of grid rows
#   nCol - num of grid columns
#----------------------------------------------ACB
.parsegrid <- function(data, nRow, nCol)
{
	parsedData=list()
	rows=list()
	cols=list()
	row <- 0;
	col <- 0;
	while(length(data)) {
		#add item into column
		col <- col + 1
		cols[[col]] <- data[[1]]

		#add a nested grid type
		if (data[[1]]$type=="grid") {
			tmp <- .parsegrid(data[-1], data[[1]]$nrow, data[[1]]$ncol)
			cols[[col]]$.widgets <- tmp$gridData
			data <- tmp$unparsedData
		}
		else {
			data <- data[-1]
		}

		#check for a filled row
		if (col == nCol) {
			row <- row + 1
			col <- 0
			rows[[row]] <- cols

			#any more rows left?
			if (row == nRow) {
				#return two parts of the data
				tmp <- list()
				tmp$gridData <- rows
				tmp$unparsedData <- data #left over data
				return(tmp)
			}
		}
	}
	stop("Grid did not have enough child objects.")
}
#---------------------------------------.parsegrid


#.parsemenu-----------------------------2012-12-20
# func: - very similar to .parsegrid but for menus
#   set up a menu with children menus or menuitems
# Arguments:
#  data   - list of widgets to be used as child of menu
#  nItems - how many children to select for the menu
#----------------------------------------------ACB
.parsemenu <- function(data, nItems)
{
	menuitems <- list()
	itemCount <- 0
	while(length(data)) {
		#increment count
		itemCount <- itemCount + 1


		if (data[[1]]$type!="menuitem" && data[[1]]$type!="menu")
			stop("non menu, or menuitem widget found, when expecting one. Check your menu nitems count.")


		#add/associate widget with menu
		menuitems[[itemCount]] <- data[[1]]

		#add a nested menu type
		if (data[[1]]$type=="menu") {
			tmp <- .parsemenu(data[-1], data[[1]]$nitems)
			menuitems[[itemCount]]$.widgets <- tmp$menuData
			data <- tmp$unparsedData
		}
		else {
			data <- data[-1] #remove widget
		}

		#return menu sub items
		if (itemCount == nItems) {
			tmp <- list()
			tmp$menuData <- menuitems
			tmp$unparsedData <- data #left over data
			return(tmp)
		}
	}
	stop("menu did not have enough child menuitems. Check your menu nitems count.")
}
#---------------------------------------.parsemenu


#.PBSdimnameHelper----------------------2012-12-20
#  Adds dimnames to stuff (matrix, data.frame)
# Arguments:
#   rownames - vector of size 1, or dim[1] nameing the rows.
#              if only one name is given, a number (1..dim[1]) will be appended to the name
#   colnames - vector of size 1 or dim[2] naming columns
#              if only one name is given, then (1..dim[2]) is appended
#   dim      - vector of size 2, dim[1] is nRows, dim[2] is nCols
#----------------------------------------------ACB
.PBSdimnameHelper <- function(rownames, colnames, dim)
{
	if (is.null(rownames)) rownames <- ""
	if (is.null(colnames)) colnames <- ""

	nRows <- dim[1]
	if( length(rownames)>1 || nRows == 1 )
		rName <- rownames
	else if (length(rownames)==0)
		rName <- NULL
	else if (rownames=="")
		rName <- NULL
	else {
		rName <- paste(rownames, 1:nRows, sep="")
	}

	nCols <- dim[2]
	if (length(colnames)>1 || nCols == 1)
		cName <- colnames
	else if (length(colnames)==0)
		cName <- NULL
	else if (colnames=="")
		cName <- NULL
	else {
		cName <- paste(colnames, 1:nCols, sep="")
	}
	return(list(rName, cName))
}
#--------------------------------.PBSdimnameHelper


#.searchCollection----------------------2012-12-20
#  Searches a haystack for a needle, or a similar longer needle.
# Arguments:
#   haystack - list to search
#   needle = scaler to search for
# Output: 
#   position of needle in list (non-negative)
#   -1 if none are found
#   -2 if two similar needles are found.
#      ex) -2 for "nee" is similar to "need", and "needle"
#----------------------------------------------ACB
.searchCollection <- function(haystack, needle)
{
	similar <- -1
	for(i in 1:length(haystack)) {
		if (haystack[[i]]$param == needle) {
			return(i)
		}
		else if (any(grep(paste("^", needle, sep=""), haystack[[i]]$param))) {
			#this is used to find any similar matches
			#if more than two are similar, then it is impossible
			#to determine which needle we are after
			if (similar == -1)
				similar <- i #this is the first similar needle
			else
				similar <- -2 #two similar needles were found
		}
	}
	return(similar)
}
#--------------------------------.searchCollection


#.setWinValHelper-----------------------2012-12-20
.setWinValHelper <- function(varname, value, winName)
{
	x  <- .map.get(winName, varname)
	tget(.PBSmod)
	wid<- .PBSmod[[winName]]$widgets[[varname]]

	#if tclvar is known, we can set it directly.
	if (!is.null(x[["tclvar"]])) {

		#value should only be length 1
		if (length(value)!=1)
			stop(paste('unable to set "', varname, '": value given should be of length 1. given length: ',length(value), sep=""))

		if (is.na(value))
			value <- ""
		else if (!is.logical(value))
			value <- as.character(value)

		tclvalue(x$tclvar) <- value

		#some widgets must update other widgets(or functions) when they change
		if (!is.null(x[["onChange"]])) {
			if (is.function(x$onChange)) {
				.do.gui.call(x$onChange, list())
			}
			if (is.character(x$onChange)) #function names are accepted as strings
				if (exists(x$onChange,mode="function"))
					.do.gui.call(x$onChange, list())
		}
		return(value)
	}

	#special case for table arrays
	else if( !is.null(x[["tclarray"]]) ) {
		return( .table.setvalue( winName, wid$name, value ) )
	}

	#otherwise if tclwidget is known (only text widget)
	else if (!is.null(x[["tclwidget"]])) {
		#special case for text boxes
		if (wid$type=="text") {
			if (wid$edit==FALSE)
				tkconfigure(x$tclwidget, state="normal")

			if (length(value)>1)
				value = paste(value,collapse="\n")

			tkdelete(x$tclwidget, "0.0", "end") #clear text widget
			tkinsert(x$tclwidget,"0.0",value) #update

			if (wid$edit==FALSE)
				tkconfigure(x$tclwidget, state="disabled")
			return(value)
		} else if( wid$type == "notebook" ) {
			tcl( x$tclwidget, "raise", value )
			return(value)
		}
		stop(paste("unhandled widget type", x$tclwidget))
	}

	else if( !is.null( x[[ "droplist_widget" ]] ) ) {
		if( is.logical( x[[ "droplist_widget" ]] ) ) #hack to only set values, and not .id
			return( value )
		tkconfigure( x[[ "droplist_widget" ]], values = value )
		#eval(parse(text=".PBSmod[[winName]]$widgets[[varname]]$labels <<- value")) #there's no way to specify different labels via setWinVal, so assume the same
		tget(.PBSmod)
		.PBSmod[[winName]]$widgets[[varname]]$labels <- value #there's no way to specify different labels via setWinVal, so assume the same
		tput(.PBSmod)
		return( value )
	}
	
	#catch any special "high level" widgets that do not have
	#tclvar or tclwidget. If however no wid is defined, we are doomed
	if (is.null(wid)) {
		stop(paste('unable to set "', varname, '": not found.', sep=""))
	}

	#special case for matrix
	if (wid$type=="matrix") {
		if (length(wid$names)==1) {
			if (!is.matrix(value))
				stop(paste('unable to set "', varname, '": supplied value is not a matrix.', sep=""))
			for(i in 1:nrow(value))
				for(j in 1:ncol(value))
					.setWinValHelper(paste(varname,"[",i,",",j,"]",sep=""), value[i,j], winName)
			return(value)
		}
	}

	#special case for data
	if (wid$type=="data") {
		if (length(wid$names)==1) {
			#todo: if not a data.frame
			#	stop(paste('unable to set "', varname, '": supplied value is not a dataframe.', sep=""))
			for(i in 1:nrow(value))
				for(j in 1:ncol(value))
					.setWinValHelper(paste(varname,"[",i,",",j,"]d",sep=""), value[i,j], winName)
			return(value)
		}
	}

	#special case for vector
	if (wid$type=="vector") {
		if (length(wid$names)==1) {
			if (length(value)!=wid$length)
				stop(paste('unable to set "', varname, '": supplied vector should have length ', wid$length, sep=""))
			for(i in 1:length(value))
				.setWinValHelper(paste(varname,"[",i,"]",sep=""), value[i], winName)
			return(value)
		}
	}
	
	#special case for superobject
	if( wid$type == "superobject" || wid$type == "object" ) {
		#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ wid$name ]]$.data <<- value"))
		tget(.PBSmod)
		.PBSmod[[ winName ]]$widgets[[ wid$name ]]$.data <- value
		tput(.PBSmod)
		.superobject.redraw( winName, wid$name )
		return( value )
	}

	print( wid )
	stop(paste("unable to update \"", varname, "\" - no widget found.", sep=""))
}
#---------------------------------.setWinValHelper


#.stopWidget----------------------------2012-12-20
#  Fatal error during window creation (not parse)
# Arguments:
#   err       - error string to display
#   wid.debug - list of widget code (created in parsing process
#   winName   - active window name
#----------------------------------------------ACB
.stopWidget <- function(err, wid.debug, winName)
{
	err <- paste("\nGUI parse error (", wid.debug$fname, ":", 
	             wid.debug$line.start, ") : ", err, "\n\n", sep="")

	if (length(wid.debug$sourceCode)>0) {
		j <- 0;
		for(i in wid.debug$line.start:wid.debug$line.end) {
			j <- j + 1
			err <- paste(err, i, ": ", wid.debug$sourceCode[j], "\n", sep="")
		}
	}
	tget(.PBSmod)
	tt <- .PBSmod[[winName]]$tkwindow
	tkdestroy(tt)
	stop(err, call.=FALSE)
}
#--------------------------------------.stopWidget


#.stripComments-------------------------2012-12-20
#  removes any trailing comments from a line, 
#  but ignores #'s in quoted strings
# Arguments:
#  x - a string with or without comments
# Output:
#  string without comments
# Example: 
#   x='type="label" text="I am #1" #comment'
#   returns 'type="label" text="I am #1"'
#----------------------------------------------ACB
.stripComments <- function(x)
{
	if (length(x)>1) {
		retVal <- c()
		for(i in 1:length(x)) {
			retVal[i] <- .stripComments(x[i])
		}
		return(retVal)
	}
	return(.Call("stripComments", as.character(x), PACKAGE="PBSmodelling"))
}
#-----------------------------------.stripComments


#.stripSlashes--------------------------2012-12-20
#  Removes slashes from a string.
# Arguments:
#   TODO
#----------------------------------------------ACB
.stripSlashes <- function(x, fname="", line.start=0, line.end=0, sourcefile=list())
{
	word<-""
	escape<-0
	for(i in 1:nchar(x)) {
		ch<-substr(x,i,i)

		#escaped char is expected
		if (escape!=0) {
			if (ch=="n")
				ch <- "\n"
			else if (ch=="t")
				ch <- "\t"
			else if (ch=="r")
				ch <- "\r"
			word <- paste(word, ch, sep="")
			escape <- 0
		}
		#next char will be escapped
		else if (ch=="\\") {
			escape <- 1
		}
		#shouldnt find any singlequotes - if we did it should be a vector of strings
		else if (ch=="'" || ch=="\"") {
			.catError("unexpected singlequote found.", fname, line.start, line.end, sourcefile)
			return(NULL)
		}
		#any other character
		else {
			word <- paste(word, ch, sep="")
		}
	}
	return(word)
}
#------------------------------------.stripSlashes


#.stripSlashesVec-----------------------2012-12-20
#  Given a string x, x is split into a vector of
#  words, which were seperated by spaces however,
#  if single quotes are used, space is perserved
#  x="a b 'c d'" converts into "a" "b" "c d"
#----------------------------------------------ACB
.stripSlashesVec <- function(x, fname="", line.start=0, line.end=0, sourcefile=list())
{
	word<-""
	words=c()
	escape<-0
	quoted<-0
	quoteFound<-0
	j <- 0
	for(i in 1:nchar(x)) {
		ch<-substr(x,i,i)

		#escaped char is expected
		if (escape!=0) {
			if (ch=="n")
				ch <- "\n"
			else if (ch=="t")
				ch <- "\t"
			else if (ch=="r")
				ch <- "\r"

			word <- paste(word, ch, sep="")
			escape <- 0
		}

		#next char will be escapped
		else if (ch=="\\") {
			escape <- 1
		}
		#shouldnt find any doublequotes anywhere
		else if (ch=="\"") {
			.catError("unexpected doublequote found.", fname, line.start, line.end, sourcefile)
			return(NULL)
		}
		else if (ch=="'") {
			if (quoted==0) {
				quoted<-1
				quoteFound<-1
			}
			else {
				quoted <- 0
			}
		}
		#space found
		else if (ch==" " || ch=="\t") {
			if (quoted==0) {
				if (word!="" || quoteFound==1) {
					##save key and value
					j <- j + 1
					words[j] <- word

					#reset variables for next loop
					word <- ""
					quoteFound <- 0
				}
			}
			else {
				word <- paste(word, ch, sep="")
			}
		}
		#any other character
		else {
			word <- paste(word, ch, sep="")
		}
	}
	#look for last word
	if (quoted==0) {
		if (word!="" || quoteFound==1) {
			##save key and value
			j <- j + 1
			words[j] <- word

			#reset variables for next loop
			word <- ""
			quoteFound <- 0
		}
	}
	else {
		.catError("unterminated quote found.", fname, line.start, line.end, sourcefile)
		return(NULL)
	}
	if (is.null(words))
		words[1]<-""
	return(words)
}
#---------------------------------.stripSlashesVec


#.superobject.redraw--------------------2012-12-20
.superobject.redraw <- function( winName, widget_name )
{
	tget(.PBSmod)
	userObject <- .PBSmod[[ winName ]]$widgets[[ widget_name ]]$.data
	rows_to_display <- .PBSmod[[ winName ]]$widgets[[ widget_name ]]$rows_to_display
	display_top <- .PBSmod[[ winName ]]$widgets[[ widget_name ]]$display_top
	ncols <- ncol( userObject )
	new_widget_name <- paste( "[superobject]", widget_name, sep="" ) 
	
	#reload screen data with userObject values
	for( i in 1:rows_to_display ) {
		for( j in 1:ncols ) {
			var_name = paste( new_widget_name, "[", i, ",", j ,"]d", sep="" )
			tmp_ptr <- .map.get( winName, var_name )
			tclvalue( tmp_ptr$tclvar ) <- userObject[ display_top + i - 1, j ]
		}
	}
}
#------------------------------.superobject.redraw


#.superobject.saveValues----------------2012-12-20
#  Must be outside of .createwidget.superobject 
#  since getWinVal must be able to call this.
#----------------------------------------------ACB
.superobject.saveValues <- function( winName, widget_name )
{
	tget(.PBSmod)
	modes <- .PBSmod[[ winName ]]$widgets[[ paste( "[superobject]", widget_name, sep="" ) ]]$modes
	
	userObject <- .PBSmod[[ winName ]]$widgets[[ widget_name ]]$.data
	rows_to_display <- .PBSmod[[ winName ]]$widgets[[ widget_name ]]$rows_to_display
	display_top <- .PBSmod[[ winName ]]$widgets[[ widget_name ]]$display_top
	ncols <- ncol( userObject )
	new_widget_name <- paste( "[superobject]", widget_name, sep="" ) 
	
	#Save data viewable on screen (which user might have changed)
	for( i in 1:rows_to_display ) {
		for( j in 1:ncols ) {
			var_name = paste( new_widget_name, "[", i, ",", j ,"]d", sep="" )
			tmp_ptr <- .map.get( winName, var_name )
			userObject[ display_top + i - 1, j ] <- .convertMode( tclvalue( tmp_ptr$tclvar ), modes[ j ] )
		}
	}
		
	#save back to global memory (so getWinVal can access it)
	#eval(parse(text=".PBSmod[[ winName ]]$widgets[[ widget_name ]]$.data <<- userObject"))
	tget(.PBSmod)
	.PBSmod[[ winName ]]$widgets[[ widget_name ]]$.data <- userObject
	tput(.PBSmod)
	invisible()
}
#--------------------------.superobject.saveValues


#.trimWhiteSpace------------------------2012-12-20
#  remove leading and trailing whitespace
# Arguments:
#  x - string to trim
# Example:
#  "   foo bar " becomes "foo bar"
#----------------------------------------------ACB
.trimWhiteSpace <- function(x)
{
	return(sub("[[:space:]]+$", "", sub("^[[:space:]]+", "", x)))
}
#----------------------------------.trimWhiteSpace


#.validateWindowDescList----------------2012-12-20
#   determines if the list represents a valid PBS Modelling description list
#   if any required fields are missing, it will halt via stop()
#   if any fields are omitted which have default values defined in the
#   .widgetDefs list, then those fields and values will be set
# Arguments:
#   x - list to validate
#----------------------------------------------ACB
.validateWindowDescList <- function(x)
{
	if (!is.list(x))
		stop("no list was given")
	if (length(x)==0)
		stop("No windows were given")

	paramOrder <- .widgetDefs

	for (i in 1:length(x)) {
		#validate each window

		#check for a window title
		if (is.null(x[[i]][["title"]]))
			x[[i]]$title <- paramOrder$window$title

		#check for widgets
		if (!is.list(x[[i]]$.widgets))
			stop("The widget list is missing")

		x[[i]]$.widgets <- .validateWindowDescWidgets(x[[i]]$.widgets)

		#TODO - check .menu too?

		if (is.null(x[[i]][["winBackground"]]))
			x[[i]]$winBackground = "#D4D0C8"

		if (is.null(x[[i]][["winForeground"]]))
			x[[i]]$winForeground = "#000000"
	}

	return(x)
}
#--------------------------.validateWindowDescList


#.validateWindowDescWidgets-------------2012-12-20
#  Used by .validateWindowDescList to validate each widget
# Arguments:
#   x - widget list to validate
# Note: this function is similar to .getParamFromStr but is
#       only designed for lists, and is not as robust.
#       ex: -no error messages for filename/line number
#           -no support for expanding shortened param names
#           -no type conversion
#----------------------------------------------ACB
#***** THIS FUNCTION APPEARS TO BE BUGGY ***** RH
#  It tests for $.widgets under grids, which doesn't seem to
#  occur any more. (I could be wrong.)
#-----------------------------------------------RH
.validateWindowDescWidgets <- function(x)
{
	paramOrder <- .widgetDefs
	for(i in 1:length(x)) {
		type <- x[[i]]$type
		if (is.null(paramOrder[[type]]))
			stop(paste("unknown widget type found:", type))
		# RH (2012-12-17)-----------------
		# I've disabled this grid check because the structure of a
		# windows description list must have changed since this 
		# little-used function was written.
		# --------------------------------
		#check children widgets of grid
		if (type=="oldgrid") {
			if (!is.list(x[[i]]$.widgets))
				stop("grid needs a .widgets list")
			for(j in 1:length(x[[i]]$.widgets))
				x[[i]]$.widgets[[j]] <- .validateWindowDescWidgets(x[[i]]$.widgets[[j]])
		}
		#look for all options, if any are missing assign the default value
		#unless they are absolutely required
		args <- paramOrder[[type]]
		for(j in 1:length(args)) {
			if (is.null(x[[i]][[args[[j]][["param"]] ]] )) {
				#a paramater is missing from the list.

				#is the paramater required?
				if (args[[j]]$required)
					stop(paste("missing argument", args[[j]]$param, "from widget", type))

				#is there a default value?
				if (!is.null(args[[j]][["default"]]))
					x[[i]][[args[[j]]$param]] <- args[[j]]$default
			}
			else {
				#the argument was found. let's check that its the right type
				#and matches the grep

				#check grep if applicable
				if (!is.null(args[[j]][["grep"]])) {
					#check grep from .widgetDefs with supplied value from list
					if (!any(grep(args[[j]]$grep, x[[i]][[args[[j]]$param]])))
						stop(paste("given value \"", x[[i]][[args[[j]]$param]], 
						"\" does not match grep:", args[[j]]$grep, sep=""))
				}
			}
		}
	}
	return(x)
}
#-----------------------.validateWindowDescWidgets


#===== THE END ===================================

