################################
# GUI for dudi.nsc function
################################
"dialog.dudi.nsc" <- function(show, history)
{
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"dudi.nsc")	
#
# Variables for text fields
#
	dudivar <- tclVar()

	dfvar <- tclVar()
	nfvar <- tclVar()
#
# Checkboxes
#
	scannfvar <- tclVar(1)
#
# Title
#
	TFrame <- tkframe(tt, relief="groove")
	labh <- tklabel(TFrame, bitmap="questhead")
	tkgrid(tklabel(TFrame,text="Non symmetric correspondence analysis", font="Times 18", foreground="red"), labh)
	tkbind(labh, "<Button-1>", function() print(help("dudi.nsc")))
	tkgrid(TFrame)
#
# Dataframe and dudi
#
	IOFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(IOFrame,text="- Input & output -", foreground="blue"), columnspan=5)
	dudi.entry <- tkentry(IOFrame, textvariable=dudivar)
	df.entry <- tkentry(IOFrame, textvariable=dfvar)
	dfnr.label <- tklabel(IOFrame, width=4)
	dfnc.label <- tklabel(IOFrame, width=4)
	choosedf.but <- tkbutton(IOFrame, text="Set", command=function() choosedf(df.entry, dfnr.label, dfnc.label))
	tkgrid(tklabel(IOFrame,text="Input data frame : "), df.entry, choosedf.but, dfnr.label, dfnc.label, sticky="w")
	tkgrid(tklabel(IOFrame,text="Output dudi name : "), dudi.entry, sticky="w")
	tkgrid(IOFrame)
#
# Number of axes
#
	NAFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(NAFrame,text="- Number of principal axes -", foreground="blue"), columnspan=2)
	nf.entry <- tkentry(NAFrame, textvariable=nfvar, width=4, state="disabled")
	scannf.cbut <- tkcheckbutton(NAFrame,text="Ask number of axes interactively", variable=scannfvar,
		command=function() if (!as.logical(tclObj(scannfvar))) tkconfigure(nf.entry, state="normal") else tkconfigure(nf.entry, state="disabled") )
	tkgrid(scannf.cbut, sticky="w", columnspan=2)
	tkgrid(tklabel(NAFrame,text="Number of axes : "), nf.entry, sticky="w")
	tkgrid(NAFrame)
#
# Local variables
#
	vnr=NULL			# Vector of dataframes row numbers
	vnc=NULL			# Vector of dataframes column numbers
	numi=1				# Number of choosed element
	done <- tclVar(0)	# To terminate the dialog

################################
# Function to build the command line from dialog widgets
################################
	"build" <- function()
	{
	#
	# Check that the data frame is not empty and get its name
	#
		if (tclvalue(dfvar) != "") {
			df  <- parse(text=tclvalue(dfvar))[[1]]
		} else {
			return(0)
		}
	#
	# If scannf is false, check that nf is not empty and get it
	#
		scannf <- as.logical(tclObj(scannfvar))
		if (!scannf) {
			if (tclvalue(nfvar) != "")
				nf <- parse(text=tclvalue(nfvar))[[1]]
				else nf <- 2
		} else nf <- 2
	#
	# Make the command line
	#
		substitute(dudi.nsc(df = df, scannf = scannf, nf = nf))
	}
		
################################
# Function to reset all dialog elements to default values
################################
	"reset" <- function()
	{
		tclvalue(dfvar)<-""
		tclvalue(nfvar)<-""
		tclvalue(scannfvar)<-"1"
		tkconfigure(dfnr.label, text="")
		tkconfigure(dfnc.label, text="")
	
	}
	
################################
# Function to launch computations
################################
	"execcomp" <- function()
	{
	#
	# Check that the analysis name is not empty and get it
	#
		if (tclvalue(dudivar) == "") tkinsert(dudi.entry, "end", "untitled1")
		dudiname <- parse(text=paste("\"",tclvalue(dudivar)[[1]],"\"",sep=""))
	#
	# Build and display the command line so that the user can check it
	#
		cmd <- build()
		if (cmd == 0) return(0)
		if (show) {
			#
			# Echoe the command line to the console
			#
			pr1 <- substr(options("prompt")$prompt, 1,2)
			cat(eval(dudiname), " <- ", deparse(cmd, width.cutoff = 256), "\n", pr1, sep="")
		}
	#
	# Execute the command
	#
		assign("ade4TkGUIFlag", 1, envir=env_ade4tkgui)
		mydudi <- eval.parent(cmd)
		assign(eval(dudiname), mydudi, envir=env_ade4tkgui)
		dialog.dudi.display(show, history, eval(dudiname))
		rm("ade4TkGUIFlag", envir=env_ade4tkgui)
		if (history) {
			commande = paste(eval(dudiname), " <- ", deparse(cmd, width.cutoff = 500), sep = "")
			rewriteHistory(commande)
		}
	}

#
# Reset and Submit buttons
#
	RCSFrame <- tkframe(tt, relief="groove")
	reset.but <- tkbutton(RCSFrame, text="Reset", command=reset)
	cancel.but <- tkbutton(RCSFrame, text="Dismiss", command=function() tkdestroy(tt))
	submit.but <- tkbutton(RCSFrame, text="Submit", default="active", command=function() execcomp())
	tkgrid(cancel.but, submit.but, reset.but, ipadx=20)	
	tkgrid(RCSFrame)
#
# If window is closed by user, terminate the dialog
#
	tkbind(tt, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tt, "<KeyPress-Return>", function() execcomp())
	tkbind(tt, "<KeyPress-Escape>", function() tkdestroy(tt))
#
# User closed the window
#
	if(tclvalue(done)=="2") return()
}

