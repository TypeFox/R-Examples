##################################################
#                 PBS Modelling                  #
# -----------------------------------------------#
# This file aims to include functions specific   #
# to the history widget and other GUI functions  #
#                                                #
# Authors:                                       #
#  Jon T. Schnute <schnutej-dfo@shaw.ca>,        #
#  Alex Couture-Beil <alex@mofo.ca>, and         #
#  Rowan Haigh <rowan.haigh@dfo-mpo.gc.ca>       #
#                                                #
# The history widget creates a bunch of standard #
# PBS widgets, which call the appropriate func-  #
# tions as if a power-user created history.      #
##################################################

#=================================================
#              HISTORY FUNCTIONS
#=================================================

#addHistory-----------------------------2012-12-19
# save history
# Arguments:
#  hisname - history instance name if multiple are active
#----------------------------------------------ACB
addHistory <- function(hisname="")
{
	if (hisname=="") 
		hisname <- getWinAct()[1]

	tget(PBS.history)
	if (!is.list(PBS.history)) 
		stop("History not intialized - see initHistory")
	if (!is.list(PBS.history[[hisname]])) 
		stop(paste("History \"", hisname,"\" not intialized - see initHistory", sep=""))

	x <- PBS.history[[hisname]]                 # old history
	itemLen <- length(x)-1                      # don't count header
	index <- PBS.history[[hisname]][[1]]$index  # make it a real index
	insertMode <- getWinVal(PBS.history[[hisname]][[1]]$modename)[[PBS.history[[hisname]][[1]]$modename]]
	
	if (is.null(insertMode) || index==0) {
		insertMode <- "a"
	}
	if (insertMode=="a") {
		#insert to the right of current index
		#eval(parse(text="PBS.history[[hisname]][[index+2]] <<- getWinVal()"))
		PBS.history[[hisname]][[index+2]] <- getWinVal()
		if (index < itemLen) {
			for(i in (index+2):(itemLen+1)) {
				#eval(parse(text="PBS.history[[hisname]][[i+1]] <<- x[[i]]"))
				PBS.history[[hisname]][[i+1]] <- x[[i]]
			}
		}
		#point index to inserted pos
		#eval(parse(text="PBS.history[[hisname]][[1]]$index <<- index+1"))
		PBS.history[[hisname]][[1]]$index <- index + 1
		tput(PBS.history)
	}
	else if (insertMode=="b") {
		#insert to the left of current index
		#eval(parse(text="PBS.history[[hisname]][[index+1]] <<- getWinVal()"))
		PBS.history[[hisname]][[index+1]] <- getWinVal()
		for(i in (index+1):(itemLen+1)) {
			#eval(parse(text="PBS.history[[hisname]][[i+1]] <<- x[[i]]"))
			PBS.history[[hisname]][[i+1]] <- x[[i]]
		}
		tput(PBS.history)
	}
	else if (insertMode=="o") {
		#overwrite the current index
		#eval(parse(text="PBS.history[[hisname]][[index+1]] <<- getWinVal()"))
		PBS.history[[hisname]][[index+1]] <- getWinVal()
		tput(PBS.history)
	}
	else {
		stop(paste("unknown insert mode:", insertMode))
	}
	.updateHistory(hisname)
	.updateHistoryButtons( hisname )
}
#---------------------------------------addHistory


#backHistory----------------------------2012-12-20
#  Move back in history.
# Arguments:
#   hisname   - history instance name if multiple are active
#----------------------------------------------ACB
backHistory <- function(hisname="")
{
	if (hisname=="") 
		hisname <- getWinAct()[1]
	win <- strsplit(hisname, "\\.")[[1]][1]

	tget(PBS.history)
	if (!is.list(PBS.history)) 
		stop("History not intialized - see initHistory function help")
	if (!is.list(PBS.history[[hisname]])) 
		stop(paste("History \"", hisname,"\" not intialized - see initHistory", sep=""))

	i <- PBS.history[[hisname]][[1]]$index
	if (i < 2) {
		#cat("history widget: warning, current position is already at front of history list.\n")
		return()
	}
	#PBS.history[[hisname]][[1]]$index <<- i <- i-1
	i <- i-1
	#eval(parse(text="PBS.history[[hisname]][[1]]$index <<- i"))
	PBS.history[[hisname]][[1]]$index <- i
	tput(PBS.history)
	.updateHistoryButtons( hisname )                      # does not modify PBS.history
	setWinVal(PBS.history[[hisname]][[i+1]], winName=win) #i is always one lower
	.updateHistory(hisname)                               # does not modify PBS.history
	if (!is.null(PBS.history[[hisname]][[1]]$func))
		do.call(PBS.history[[hisname]][[1]]$func, list(),envir=.PBSmodEnv$.PBSmod[[.PBSmodEnv$.PBSmod$.activeWin]]$env)
}
#--------------------------------------backHistory


#clearHistory---------------------------2012-12-20
#  Remove all history elements from
# Arguments:
#   hisname - history instance name if multiple are active
#----------------------------------------------ACB
clearHistory <- function(hisname="")
{
	if (hisname=="") 
		hisname <- getWinAct()[1]

	tget(PBS.history)
	if (!is.list(PBS.history)) 
		stop("History not intialized - see initHistory")
	if (!is.list(PBS.history[[hisname]])) 
		stop(paste("History \"", hisname,"\" not intialized - see initHistory", sep=""))

	tmp <- PBS.history[[hisname]][[1]]
	tmp$index = 0
	#eval(parse(text="PBS.history[[hisname]] <<- list(0)"))
	PBS.history[[hisname]] <- list(0)
	#eval(parse(text="PBS.history[[hisname]][[1]] <<- tmp"))
	PBS.history[[hisname]][[1]] <- tmp
	tput(PBS.history)

#	len <- length(PBS.history[[hisname]])
#	if (len > 1) {
#		for(i in 2:len) {
			#PBS.history[[hisname]][[i]] <<- NULL #something weird is happening here
#			rmHistory(hisname)
#		}
#	}
	.updateHistory(hisname)
	.updateHistoryButtons( hisname )
}
#-------------------------------------clearHistory


#exportHistory--------------------------2012-12-20
#  Save PBS history to a file
# Arguments:
#   hisname - history instance name if multiple are active
#   fname   - initial filename to save under
#----------------------------------------------ACB
exportHistory <- function(hisname="", fname="")
{
	if (hisname=="") hisname <- getWinAct()[1]

	tget(PBS.history)
	if (!is.list(PBS.history[[hisname]]))
		stop("unable to export history. Incorect history name given.")

	if (fname=="")
		fname <- selectFile( initialfile=paste(hisname,".History.r", sep=""), mode="save" )
	if (fname=="")
		stop("no filename given.")

	x = PBS.history[[hisname]]
	x[[1]] <- NULL # remove history widget info
	writeList(x, fname)
}
#------------------------------------exportHistory


#firstHistory---------------------------2012-12-20
#  Move to first history slide.
# Arguments:
#   hisname   - history instance name if multiple are active
#----------------------------------------------ACB
firstHistory <- function(hisname="")
{
	if (hisname=="")
		hisname <- getWinAct()[1]
	tget(PBS.history)
	if(length(PBS.history[[hisname]])>1)
		jumpHistory(hisname, 1)
	.updateHistoryButtons( hisname )
}
#-------------------------------------firstHistory


#forwHistory----------------------------2012-12-20
#  Move forward in history.
# Arguments:
#   hisname   - history instance name if multiple are active
#----------------------------------------------ACB
forwHistory <- function(hisname="")
{
	if (hisname=="") 
		hisname <- getWinAct()[1]
	win <- strsplit(hisname, "\\.")[[1]][1]

	tget(PBS.history)
	if (!is.list(PBS.history)) 
		stop("History not intialized - see initHistory")
	if (!is.list(PBS.history[[hisname]])) 
		stop(paste("History \"", hisname,"\" not intialized - see initHistory", sep=""))

	i <- PBS.history[[hisname]][[1]]$index
	if (i >= (length(PBS.history[[hisname]])-1)) {
		#cat("history widget: warning, current position is already at end of history list.\n")
		return()
	}
	#PBS.history[[hisname]][[1]]$index <<- i <- i+1
	i <- i+1
	#eval(parse(text="PBS.history[[hisname]][[1]]$index <<- i"))
	PBS.history[[hisname]][[1]]$index <- i
	tput(PBS.history)
	.updateHistoryButtons( hisname )
	setWinVal(PBS.history[[hisname]][[i+1]], winName=win) #i is always one lower
	.updateHistory(hisname)
	if (!is.null(PBS.history[[hisname]][[1]]$func))
		do.call(PBS.history[[hisname]][[1]]$func, list(),envir=.PBSmodEnv$.PBSmod[[.PBSmodEnv$.PBSmod$.activeWin]]$env)
}
#--------------------------------------forwHistory


#importHistory--------------------------2012-12-20
#  Import PBS history from a file.
# Arguments:
#   hisname - history instance name if multiple are active
#   fname   - initial filename to open from
#----------------------------------------------ACB
importHistory <- function(hisname="", fname="", updateHis=TRUE)
{
	if (hisname=="") hisname <- getWinAct()[1]
	win <- strsplit(hisname, "\\.")[[1]][1]
	
	tget(PBS.history)
	if (!is.list(PBS.history[[hisname]]))
		stop("unable to import history. Incorect history name given.")

	if (fname=="")
		fname <- selectFile( mode="open" )
	if ( is.null( fname ) || fname=="" )
		stop("no filename given.")

	newHist <- readList(fname)
	insertMode <- getWinVal(PBS.history[[hisname]][[1]]$modename)[[PBS.history[[hisname]][[1]]$modename]]

	a <- PBS.history[[hisname]]
	#eval(parse(text="PBS.history[[hisname]] <<- list()"))
	PBS.history[[hisname]] <- list()
	index <- max(0, min(a$index, length(a)-1))
	if (insertMode!="b" || index==0)
		index <- index + 1
	i <- 1

	repeat {
		if( !length(a) && !length(newHist) )
			break
		if( i > index && length(newHist) ) {
			#eval(parse(text="PBS.history[[hisname]][[i]] <<- newHist[[1]]"))
			PBS.history[[hisname]][[i]] <- newHist[[1]]
			newHist[[1]] <- NULL
		} else {
			#eval(parse(text="PBS.history[[hisname]][[i]] <<- a[[1]]"))
			PBS.history[[hisname]][[i]] <- a[[1]]
			a[[1]] <- NULL
		}
		i <- i + 1
	}
	#eval(parse(text="PBS.history[[hisname]][[1]]$index <<- 1"))
	PBS.history[[hisname]][[1]]$index <- 1
	tput(PBS.history)

	#update with new history settings
	if (updateHis)
		jumpHistory(hisname, index)

	return(invisible(PBS.history[[hisname]]))
}
#------------------------------------importHistory


#initHistory----------------------------2012-12-20
#  Setup the History "list".
# Arguments:
#   hisname   - history instance name if multiple are active
#   indexname - customized index widget name
#   sizename  - customized size widget name
#   overwrite - retain old history?
#----------------------------------------------ACB
initHistory <- function(hisname, indexname=NULL, sizename=NULL, buttonnames=NULL, modename=NULL, func=NULL, overwrite=TRUE)
{

	if (!exists("PBS.history", envir = .PBSmodEnv)) #.GlobalEnv))
		PBS.history <- list()
	else
		PBS.history <- get("PBS.history", envir = .PBSmodEnv) #.GlobalEnv)

	if (is.null(func) || func=="")
		func <- NULL
	if (!is.null(func)) {
		if (!exists(func,mode="function",envir=.PBSmodEnv$.PBSmod[[.PBSmodEnv$.PBSmod$.activeWin]]$env)) {  #look for function in the original environment
			cat(paste("Warning: cannot find function '", func, "'.\n", sep=""))
			func <- NULL
		}
	}

	if (!is.list(PBS.history))
		assign("PBS.history", list(), envir = .PBSmodEnv) #.GlobalEnv)

	if (!is.list(PBS.history[[hisname]]) || overwrite) {
		PBS.history[[hisname]] <- list(0)
		PBS.history[[hisname]][[1]] <- list(index=0) #the first element is the index, all other elements are history items
	}
	#save names of entry boxes
	PBS.history[[hisname]][[1]]$indexname <- indexname
	PBS.history[[hisname]][[1]]$buttonnames <- buttonnames
	PBS.history[[hisname]][[1]]$sizename <- sizename
	PBS.history[[hisname]][[1]]$modename <- modename
	PBS.history[[hisname]][[1]]$func <- func
	assign("PBS.history", PBS.history, envir = .PBSmodEnv) #.GlobalEnv)
}
#--------------------------------------initHistory


#jumpHistory----------------------------2012-12-20
#  Need history name
#   and what index to jump to - or what entry to pull it out of
# Arguments:
#   hisname   - history instance name if multiple are active
#----------------------------------------------ACB
jumpHistory <- function(hisname="", index="")
{
	if (hisname=="") 
		hisname <- getWinAct()[1]
	win <- strsplit(hisname, "\\.")[[1]][1]

	tget(PBS.history)
	if (!is.list(PBS.history)) 
		stop("History not intialized - see initHistory")
	if (!is.list(PBS.history[[hisname]])) 
		stop(paste("History \"", hisname,"\" not intialized - see initHistory", sep=""))

	if (is.numeric(index))
		i <- index
	else if (index=="")
		i <- as.numeric(getWinVal(PBS.history[[hisname]][[1]]$indexname))
	else
		i <- as.numeric(getWinVal(index))

	if (i > length(PBS.history[[hisname]])-1 || i <= 0) {
		cat("Error: history index is out of bounds.\n")
		return()
	}
	#eval(parse(text="PBS.history[[hisname]][[1]]$index <<- i")) #update index
	PBS.history[[hisname]][[1]]$index <- i # update index
	tput(PBS.history)
	setWinVal(PBS.history[[hisname]][[i+1]], winName=win)       #i is always one lower
	.updateHistory(hisname)
	if (!is.null(PBS.history[[hisname]][[1]]$func))
		do.call(PBS.history[[hisname]][[1]]$func, list(),envir=.PBSmodEnv$.PBSmod[[.PBSmodEnv$.PBSmod$.activeWin]]$env)
	.updateHistoryButtons( hisname )
}
#--------------------------------------jumpHistory


#lastHistory----------------------------2012-12-20
#  Move to last history slide.
# Arguments:
#   hisname   - history instance name if multiple are active
#----------------------------------------------ACB
lastHistory <- function(hisname="")
{
	if (hisname=="")
		hisname <- getWinAct()[1]
	tget(PBS.history)
	if(length(PBS.history[[hisname]])-1 > 0)
		jumpHistory(hisname, length(PBS.history[[hisname]])-1)
	.updateHistoryButtons( hisname )
}
#--------------------------------------lastHistory


#rmHistory------------------------------2012-12-19
#  If index is numeric - delete history in that spot
#  else delete the history where the current index points to 
#  (and not the value of the current index box - as a user might not have pushed enter)
# Arguments:
#   hisname   - history instance name if multiple are active
#----------------------------------------------ACB
rmHistory <- function(hisname="", index="")
{
	if (hisname=="") 
		hisname <- getWinAct()[1]
	win <- strsplit(hisname, "\\.")[[1]][1]

	tget(PBS.history)
	if (!is.list(PBS.history)) 
		stop("History not intialized - see initHistory")
	if (!is.list(PBS.history[[hisname]])) 
		stop(paste("History \"", hisname,"\" not intialized - see initHistory", sep=""))

	if (is.numeric(index))
		i <- index
	else
		i <- PBS.history[[hisname]][[1]]$index

	if (length(PBS.history[[hisname]]) == 1) {
		cat("History list is already empty.\n")
		return()
	}

	#eval(parse(text="PBS.history[[hisname]] <<- PBS.history[[hisname]][-(i+1)]"))
	PBS.history[[hisname]] <- PBS.history[[hisname]][-(i+1)]

	#change index if it was the last element
	if (i > length(PBS.history[[hisname]])-1)
		#eval(parse(text="PBS.history[[hisname]][[1]]$index <<- length(PBS.history[[hisname]])-1")) #set index to size
		PBS.history[[hisname]][[1]]$index <- length(PBS.history[[hisname]])-1 #set index to size
	tput(PBS.history)
	#change values to current index
	i <- PBS.history[[hisname]][[1]]$index
	if (i > 0) {
		setWinVal(PBS.history[[hisname]][[i+1]], winName=win) #i is always one lower
		if (!is.null(PBS.history[[hisname]][[1]]$func))
			do.call(PBS.history[[hisname]][[1]]$func, list(),envir=.PBSmodEnv$.PBSmod[[.PBSmodEnv$.PBSmod$.activeWin]]$env)
	}
	.updateHistory(hisname)
	.updateHistoryButtons( hisname )
}
#----------------------------------------rmHistory


#sortHistory----------------------------2012-12-19
#  Sort history
# Arguments:
#  hisname - history instance name if multiple are active
#----------------------------------------------ACB
sortHistory <- function(file="",outfile=file,hisname="")
{
	if (file!="") {
		return(.sortHelperFile(file, outfile))
	}
	if (hisname!="") {
		return(.sortHelperActive(hisname))
	}
	currHist <- NULL
	if (exists("PBS.history",envir=.PBSmodEnv)) {
		tget(PBS.history)
		currHist <- names(PBS.history)
	}
	if (!is.null(currHist)) {
		radios <- list(list(type="label", text="Select an active window history to sort", sticky="w", padx=12))
		i <- 2
		for(h in currHist) {
			len <- length(PBS.history[[h]])-1
			if (len==1)
				items <- "(1 item)"
			else
				items <- paste("(", len, " items)", sep="")
			radios[[i]] <- list(type="radio",
			                    name="hisname",
			                    value=h,
			                    text=paste(h, items),
			                    mode="character",
			                    sticky="w",
			                    padx=12)
			i <- i+1
		}
		radios[[i]] <- list(type="button", "function"=".sortHelper", action="active", text="Sort Active History", sticky="w", padx=12)
	} 
	else {
		radios <- list(list(type="label", text="No active history widgets could be found.\nTry creating a window with a history widget.", sticky="w", padx=12))
	}
	win <- list(title = "Sort History",
	            windowname = "sort.History",
	            .widgets = c(
	                list(
	                  list(type="label", text="Active History                                            ", font="bold underline", fg="red3", sticky="w")
	                ),
	                radios,
	                list(
	                  list(type="null"),
	                  list(type="label", text="Saved History                                            ", font="bold underline", fg="red3", sticky="w"),
	                  list(type="label", text="Open a saved history file to sort and select a file \nto save to. (which can be the same file)", sticky="w", padx=12),
                         list(type = "grid", .widgets = list(
                           list(
                             list(type="entry", name="openfile", sticky="e", padx=12, width=30, mode="character"),
                             list(type="button", "function"=".updateFile", action="open", text="Open From", sticky="w", width=10)
                           ),
                           list(
                             list(type="entry", name="savefile", sticky="e", padx=12, width=30, mode="character"),
                             list(type="button", "function"=".updateFile", action="save", text="Save To", sticky="w", width=10)
                           ))),
                         list(type="button", "function"=".sortHelper", action="file", text="Sort Saved History", sticky="w", padx=12)
                       )
	            ))
	createWin(list(win))
}
#--------------------------------------sortHistory


#=================================================
#               HIDDEN FUNCTIONS
#=================================================


#.createWidget.history------------------2012-12-20
.createWidget.history <- function(tk, widgetList, winName)
{
	widget <- widgetList[[ 1 ]]
	#FIXME TODO figure out a better naming convension for what vars should be returned, or hidden by getWinVal
	#currently, anything with "PBS." is NOT returned. maybe we want PBS.hidden.XXXXX
	indexname=paste("PBS.history.", widget$name, ".index", sep="") #widget name that stores/displays the index number
	sizename=paste("PBS.history.", widget$name, ".size", sep="") #widget name that displays the size of history
	modename=paste("PBS.history.", widget$name, ".mode", sep="") #widget name that displays the size of history
	textname=paste( "PBShistory.", widget$name, ".caption", sep="") #widget name that displays text

	button_names <- list( 
		first = paste("PBS.history.", widget$name, ".button.first", sep=""),
		back = paste("PBS.history.", widget$name, ".button.back", sep=""),
		"next" = paste("PBS.history.", widget$name, ".button.next", sep=""),
		last = paste("PBS.history.", widget$name, ".button.last", sep="") )

	widget$name <- paste(winName, widget$name, sep=".")

	#initialize a list to be used once the window is created
	initHistory(widget$name, indexname=indexname, sizename=sizename, 
		buttonnames=button_names, modename=modename, func=widget[["function"]])

	historyGrid <- 
	list(type="grid", nrow=2, ncol=1, font="", fg=widget$fg, bg=widget$bg, byrow=TRUE, borderwidth=1, relief="sunken", padx=widget$padx, pady=widget$pady, .widgets=
		list(
			list(
				list(type="grid", nrow=3, ncol=4, font="", fg=widget$fg, bg=widget$bg, byrow=TRUE, borderwidth=1, relief="flat", padx=0, pady=0, sticky="we", .widgets=
					list(
						list(
							list(type="button", text="<<", name=button_names[["first"]], font="", fg=widget$fg, bg=widget$bg, width=5, "function"="firstHistory", action=widget$name, sticky="", padx=0, pady=0),
							list(type="button", text="<",  name=button_names[["back"]], font="", fg=widget$fg, bg=widget$bg, width=5, "function"="backHistory", action=widget$name, sticky="", padx=0, pady=0),
							list(type="button", text=">",  name=button_names[["next"]], font="", fg=widget$fg, bg=widget$bg, width=5, "function"="forwHistory", action=widget$name, sticky="", padx=0, pady=0),
							list(type="button", text=">>", name=button_names[["last"]], font="", fg=widget$fg, bg=widget$bg, width=5, "function"="lastHistory", action=widget$name, sticky="", padx=0, pady=0)
						),
						list(
							#list(type="label", text="Index", font="", fg=widget$fg, bg=widget$bg, sticky="", padx=0, pady=0),
							list(type="button", text="Sort", font="", fg=widget$fg, bg=widget$bg, width=5, "function"="doAction", action=paste("sortHistory(hisname=\"",widget$name,"\")",sep=""), sticky="", padx=0, pady=0),
							list(type="entry", name=indexname, value="0", width=5, label="", font="", entryfg=widget$entryfg, entrybg=widget$entrybg, "function"="jumpHistory", action=widget$name, enter=TRUE, mode="numeric", padx=0, pady=0, entrybg="white", edit=T),
							list(type="entry", name=sizename, value="0", width=5, label="", font="", action="", enter=TRUE, mode="numeric", padx=0, pady=0, edit=F),
							list(type="button", text="Empty", font="", fg=widget$fg, bg=widget$bg, width=5, "function"="clearHistory", action=widget$name, sticky="", padx=0, pady=0)
						),
						list(
							list(type="button", text="Insert", font="", fg=widget$fg, bg=widget$bg, width=5, "function"="addHistory", action=widget$name, sticky="", padx=0, pady=0),
							list(type="button", text="Delete", font="", fg=widget$fg, bg=widget$bg, width=5, "function"="rmHistory", action=widget$name, sticky="", padx=0, pady=0),
							list(type="button", text="Import", font="", fg=widget$fg, bg=widget$bg, width=5, "function"="importHistory", action=widget$name, sticky="", padx=0, pady=0),
							list(type="button", text="Export", font="", fg=widget$fg, bg=widget$bg, width=5, "function"="exportHistory", action=widget$name, sticky="", padx=0, pady=0)
						)
					)
				)
			),
			list(
				list(type="grid", nrow=1, ncol=3, font="", fg=widget$fg, bg=widget$bg, byrow=TRUE, borderwidth=1, relief="raised", padx=0, pady=0, sticky="we", .widgets=
					list(
						list(
							list(type="radio", name=modename, value="b", text="before", font="7", fg=widget$fg, bg=widget$bg, "function"="", action=widget$name, mode="character", sticky="w", padx=0, pady=0, edit=T),
							list(type="radio", name=modename, value="a", selected=TRUE, text="after", font="7", fg=widget$fg, bg=widget$bg, "function"="", action=widget$name, mode="character", sticky="w", padx=0, pady=0, edit=T),
							list(type="radio", name=modename, value="o", text="ovr", font="7", fg=widget$fg, bg=widget$bg, "function"="", action=widget$name, mode="character", sticky="w", padx=0, pady=0, edit=T)
						)
					)
				)
			)
		)
	)
	
	if( !is.null( widget[[ "text" ]] ) ) {
		widget$text <- tolower( widget$text )
		text_wid <- list(type = "text", name = textname, height = 6, width = 18, 
                       edit = TRUE, scrollbar = TRUE, fg = "black", bg = "white", 
                       mode = "character", font = "", value = "", borderwidth = 0, 
                       relief = "sunken", sticky = "", padx = 0, pady = 0 )
		if( widget[[ "textsize" ]] > 0 ) {
			if( widget$text == "n" || widget$text == "s" )
				text_wid$height = widget$textsize
			else
				text_wid$width = widget$textsize
		}
		if( widget$text == "n" || widget$text == "s" ) {
			historyGrid$nrow <- historyGrid$nrow + 1
			if( widget$text == "s" ) {
				historyGrid$.widgets[[3]] <- list( text_wid )
			} else { #North
				historyGrid$.widgets[[ 3 ]] <- historyGrid$.widgets[[ 2 ]]
				historyGrid$.widgets[[ 2 ]] <- historyGrid$.widgets[[ 1 ]]
				historyGrid$.widgets[[ 1 ]] <- list( text_wid )
			}
		} else { # E or W
			historyGrid$borderwidth <- 0
			historyGrid$padx <- 0
			historyGrid$pady <- 0
			newGrid <- list(type="grid", nrow=1, ncol=2, font="", fg=widget$fg, bg=widget$bg, byrow=TRUE, borderwidth=1, relief="sunken", padx=widget$padx, pady=widget$pady )
			
			if( widget$text == "w" )
				newGrid$.widgets = list( list( text_wid, historyGrid ) )
			else
				newGrid$.widgets = list( list( historyGrid, text_wid ) )

			historyGrid <- newGrid
		}
	}
	widgets <- .convertOldGridToNewGrid( historyGrid )
	tmp <- .createWidget.grid(tk, widgets, winName)
	if (widget$import!="")
		importHistory(widget$name, widget$import, FALSE)
	.updateHistoryButtons( widget$name )
	return(list( widget = tmp$widget, widgetList = widgetList[ -1 ] ) )
}
#----------------------------.createWidget.history


#.updateHistory-------------------------2012-12-19
# update widget values
# Arguments:
#  hisname - history instance name if multiple are active
#----------------------------------------------ACB
.updateHistory <- function(hisname)
{
	tget(PBS.history)
	indexname <- PBS.history[[hisname]][[1]]$indexname
	sizename  <- PBS.history[[hisname]][[1]]$sizename
	x<-list()

	if (!is.null(indexname))
		x[[indexname]] <- PBS.history[[hisname]][[1]]$index

	if (!is.null(sizename))
		x[[sizename]] <- length(PBS.history[[hisname]])-1

	win <- strsplit(hisname, "\\.")[[1]][1]

	setWinVal(x, winName=win)
}
#-----------------------------------.updateHistory


#.updateHistoryButtons------------------2012-12-19
.updateHistoryButtons <- function( hisname )
{
	tget(PBS.history)
	i <- PBS.history[[hisname]][[1]]$index
	n <- length(PBS.history[[hisname]])-1

	first_name <- PBS.history[[hisname]][[1]]$buttonnames[["first"]]
	back_name <- PBS.history[[hisname]][[1]]$buttonnames[["back"]]
	next_name <- PBS.history[[hisname]][[1]]$buttonnames[["next"]]
	last_name <- PBS.history[[hisname]][[1]]$buttonnames[["last"]]

	if( i == n ) {
		#last position
		setWidgetState( next_name, "disabled" )
		setWidgetState( last_name, "disabled" )
	}
	if( i <= 1 ) {
		#first position
		setWidgetState( back_name, "disabled" )
		setWidgetState( first_name, "disabled" )
	}
	#enabled buttons
	if( i < n ) {
		setWidgetState( next_name, "normal" )
		setWidgetState( last_name, "normal" )
	}
	if( i > 1 ) {
		setWidgetState( back_name, "normal" )
		setWidgetState( first_name, "normal" )
	}
	#if( i == ( length(PBS.history[[hisname]]) ) - 1 ) {}
}
#----------------------------.updateHistoryButtons


#=================================================
# Helper functions for sortHistory
#=================================================


#.sortHelper----------------------------2012-12-20
.sortHelper <- function()
{
	act <- getWinAct()[1]
	if (act=="active") {
		hisname <- getWinVal("hisname")$hisname
		.sortHelperActive(hisname)
	} else if (act=="file") {
		openfile <- getWinVal("openfile")$openfile
		savefile <- getWinVal("savefile")$savefile
		.sortHelperFile(openfile, savefile)
	}
	sortHistory()
}
#--------------------------------------.sortHelper


#.sortHelperActive----------------------2012-12-20
.sortHelperActive <- function(hisname)
{
	#convert history into a data.frame (which is used for sorting)
	tget(PBS.history)
	x <- PBS.history[[hisname]]
	if( length( x ) < 2 )
		stop( "History does not contain any items - unable to sort" )
	i <- 1
	hist <- NULL
	for( h in x[-1] ) {
		if( is.null(hist) ) {
			hist <- list() #matrix( nrow=length(x[-1]), ncol=length(h) )
			if( length( h ) )
				for( j in 1:length( h ) )
					hist[ j ] <- list(c())
		}
		if( length( h ) )
			for( j in 1:length( h ) )
				hist[[j]][i] <- as.vector( h[[j]] )[ 1 ]
		i <- i + 1
	}
	names( hist ) <- names( x[[2]] )
	#display sort widget, .done_sorting() takes care of saving the data
	hist <- as.data.frame( hist, stringsAsFactors = FALSE )

	#trim down long strings (if multiple lines take first non empty line)
	MAX_STRING_LEN <- 15
	.shortenStrings <- function(x)
	{
		x <- strsplit(x,"\n")[[1]]
		needs_dots <- FALSE
		#grab first non empty element
		x <- x[x!=""]
		if( length( x ) == 0 )
			return( "" )
		if( length( x ) > 1 ) {
			needs_dots <- TRUE
			x <- x[1]
		}
		if( nchar( x ) > MAX_STRING_LEN ) {
			needs_dots <- TRUE
			x <- strtrim( x, MAX_STRING_LEN - 3 )
		}
		if( needs_dots == TRUE )
			x <- paste( x, "...", sep="" )
		return( x )
	}
	#eval(parse(text="tmp.before <<- hist"))
	tmp.before <- hist; tput(tmp.before)
	if( ncol( hist ) > 0 ) {
		for( i in 1:ncol( hist ) ) {
			if( is.character( hist[,i] ) )
				hist[,i] <- unlist( lapply( hist[,i], .shortenStrings ) )
		}
	}
	.sortWidget( hist, hisname )
}
#--------------------------------.sortHelperActive


#.sortHelperFile------------------------2012-12-20
.sortHelperFile <- function(openfile, savefile)
{
	inHis <- readList(openfile)
	len <- length(inHis) - 1
	if (len < 1) stop("unable to sort empty history")
	x <- data.frame(new = 1:len)
	x <- fix(x); xnew <- order(x$new, na.last=NA);
	
	outHis <- list(inHis[[1]])
	j <- 2
	for (i in xnew) {
		if (!is.na(i)) {
			outHis[[j]] <- inHis[[i + 1]]
			j <- j + 1
		}
	}
	writeList(outHis, savefile)
}
#----------------------------------.sortHelperFile


#.updateFile----------------------------2012-12-20
.updateFile <- function()
{
	act <- getWinAct()[1]
	if (act=="open") {
		f <- selectFile( mode="open" )
		s <- getWinVal("savefile")$savefile
		if (s=="")
			setWinVal(list(savefile=f))
		setWinVal(list(openfile=f))
	} else if (act=="save") {
		f <- selectFile( mode="save" )
		setWinVal(list(savefile=f))
	}
}
#--------------------------------------.updateFile


#===== THE END ===================================

