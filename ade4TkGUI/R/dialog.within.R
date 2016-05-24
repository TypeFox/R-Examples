################################
# GUI for within function
################################
"dialog.within" <- function(show, history)
{
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"within")	
#
# Variables for text fields
#
	dudioutvar <- tclVar()
	dudivar <- tclVar()
	facvar <- tclVar()
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
	tkgrid(tklabel(TFrame,text="Within analysis", font="Times 18", foreground="red"), labh)
	tkbind(labh, "<Button-1>", function() print(help("within")))
	tkgrid(TFrame)
#
# Input and output dudis
#
	IOFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(IOFrame,text="- Input & output -", foreground="blue"), columnspan=3)
	dudiout.entry <- tkentry(IOFrame, textvariable=dudioutvar)
	dudi.entry <- tkentry(IOFrame, textvariable=dudivar)
	choosedudi.but <- tkbutton(IOFrame, text="Set", command=function() choosedudi(dudi.entry))
	tkgrid(tklabel(IOFrame,text="Input dudi name : "), dudi.entry, choosedudi.but)
	tkgrid(tklabel(IOFrame,text="Output dudi name : "), dudiout.entry)
	tkgrid(IOFrame)
#
# Grouping factor
#
	GFFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(GFFrame,text="- Groups -", foreground="blue"), columnspan=3)
	fac.entry <- tkentry(GFFrame, textvariable=facvar)
	choosefac.but <- tkbutton(GFFrame, text="Set", command=function() choosefac2(fac.entry, dudi.entry))
	tkgrid(tklabel(GFFrame,text="Grouping factor : "), fac.entry, choosefac.but)
	tkgrid(GFFrame)
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
	# Check that the input dudi is not empty and get its name
	#
		if (tclvalue(dudivar) != "") {
			dudi  <- parse(text=tclvalue(dudivar))[[1]]
		} else {
			return(0)
		}
	#
	# Check that the factor is not empty and get its name
	#
		if (tclvalue(facvar) != "") {
			fac  <- parse(text=tclvalue(facvar))[[1]]
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
		substitute(wca(x = dudi, fac = fac, scannf = scannf, nf = nf))
	}
		
################################
# Function to reset all dialog elements to default values
################################
	"reset" <- function()
	{
		tclvalue(dudivar)<-""
		tclvalue(nfvar)<-""
		tclvalue(scannfvar)<-"1"
	
	}
	
################################
# Function to launch computations
################################
	"execcomp" <- function()
	{
	#
	# Check that the analysis name is not empty and get it
	#
		if (tclvalue(dudioutvar) == "") tkinsert(dudiout.entry, "end", "untitled1")
		dudiname <- parse(text=paste("\"",tclvalue(dudioutvar)[[1]],"\"",sep=""))
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
