################################
# GUI for pcaiv function
################################
"dialog.pcaiv" <- function(show, history)
{
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"pcaiv")	
#
# Variables for text fields
#
	dudi1var <- tclVar()
	df2var <- tclVar()
	dudioutvar <- tclVar()
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
	tkgrid(tklabel(TFrame,text="Principal components analysis w/r to instrumental variables", font="Times 18", foreground="red"), labh)
	tkbind(labh, "<Button-1>", function() print(help("pcaiv")))
	tkgrid(TFrame)
#
# Dataframe and dudi
#	
	IOFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(IOFrame,text="- Input & output -", foreground="blue"), columnspan=5)
	dudi1.entry <- tkentry(IOFrame, textvariable=dudi1var)
	dfnr1.label <- tklabel(IOFrame, width=4)
	dfnc1.label <- tklabel(IOFrame, width=4)
	df2.entry <- tkentry(IOFrame, textvariable=df2var)
	dfnr2.label <- tklabel(IOFrame, width=4)
	dfnc2.label <- tklabel(IOFrame, width=4)
	dudiout.entry <- tkentry(IOFrame, textvariable=dudioutvar)
	choose1.but <- tkbutton(IOFrame, text="Set", command=function() choosedudirc(dudi1.entry, dfnr1.label, dfnc1.label))
	choose2.but <- tkbutton(IOFrame, text="Set", command=function() choosedf(df2.entry, dfnr2.label, dfnc2.label))
	tkgrid(tklabel(IOFrame,text="Input dudi name : "), dudi1.entry, choose1.but, dfnr1.label, dfnc1.label)
	tkgrid(tklabel(IOFrame,text="Input dataframe name : "), df2.entry, choose2.but, dfnr2.label, dfnc2.label)
	tkgrid(tklabel(IOFrame,text="Output dudi name : "), dudiout.entry)
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
	done <- tclVar(0)	# To terminate the dialog
	
################################
# Function to build the command line from dialog widgets
################################
	"build" <- function()
	{
	#
	# Check that the dudi1 is not empty and get its name
	#
		if (tclvalue(dudi1var) != "") {
			dudi1  <- parse(text=tclvalue(dudi1var))[[1]]
		} else {
			return(0)
		}
	#
	# Check that the df2 is not empty and get its name
	#
		if (tclvalue(df2var) != "") {
			df2  <- parse(text=tclvalue(df2var))[[1]]
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
		substitute(pcaiv(dudi = dudi1, df = df2, scannf = scannf, nf = nf))
	}
		
################################
# Function to reset all dialog elements to default values
################################
	"reset" <- function()
	{
		tclvalue(dudi1var) <- ""
		tclvalue(df2var) <- ""
		tclvalue(dudioutvar) <- ""
		tclvalue(nfvar) <- ""
		tkconfigure(dfnr1.label, text="")
		tkconfigure(dfnc1.label, text="")
		tkconfigure(dfnr2.label, text="")
		tkconfigure(dfnc2.label, text="")
	}
	
################################
# Function to launch computations
################################
	"execcomp" <- function()
	{
	#
	# Check that the output dudi name is not empty and set it
	#
		if (tclvalue(dudioutvar) == "") tkinsert(dudiout.entry, "end", "untitled1")
		dudiout <- parse(text=paste("\"",tclvalue(dudioutvar)[[1]],"\"",sep=""))
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
			cat(eval(dudiout), " <- ", deparse(cmd, width.cutoff = 256), "\n", pr1, sep="")
		}
	#
	# Execute the command
	#
		assign("ade4TkGUIFlag", 1, envir=env_ade4tkgui)
		myObject <- eval.parent(cmd)
		assign(eval(dudiout), myObject, envir=env_ade4tkgui)
		dialog.dudi.display(show, history, eval(dudiout))
		rm("ade4TkGUIFlag", envir=env_ade4tkgui)
		if (history) {
			commande <- paste(eval(dudiout), " <- ", deparse(cmd, width.cutoff = 500), sep = "")
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
