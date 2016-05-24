################################
# GUI for randtest functions
################################
"dialog.MCTests" <- function(show, history)
{
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"randtest")	
#
# Variables for text fields
#
	randtestvar <- tclVar()
	objectvar <- tclVar()
	npermvar <- tclVar()
	fixedvar <- tclVar(0)
#
# Checkboxes
#
	scannfvar <- tclVar(1)
#
# Title
#
	TFrame <- tkframe(tt, relief="groove")
	labh <- tklabel(TFrame, bitmap="questhead")
	tkgrid(tklabel(TFrame,text="Monte-Carlo test", font="Times 18", foreground="red"), labh)
	tkbind(labh, "<Button-1>", function() print(help("randtest")))
	tkgrid(TFrame)
#
# Frame 1 : object
#
	IOFrame <- tkframe(tt, relief="groove", borderwidth=2)
	randtest.entry <- tkentry(IOFrame, textvariable=randtestvar)
	object.entry <- tkentry(IOFrame, textvariable=objectvar)
	chooseobject.but <- tkbutton(IOFrame, text="Set", command=function() chooseduditest(object.entry))
	tkgrid(tklabel(IOFrame,text="- Input & output -", foreground="blue"), columnspan=5)
	tkgrid(tklabel(IOFrame,text="Input object name : "), object.entry, chooseobject.but)
	tkgrid(tklabel(IOFrame,text="Output object name : "), randtest.entry)
	tkgrid(IOFrame)
#
# Number of permutations
#
	NPFrame <- tkframe(tt, relief="groove", borderwidth=2)
	nperm.entry <- tkentry(NPFrame, textvariable=npermvar)
	tkgrid(tklabel(NPFrame,text="- Permutations -", foreground="blue"), columnspan=2)
	tkgrid(tklabel(NPFrame,text="Number of permutations : "), nperm.entry)
	tkgrid(NPFrame)
#
# Fixed table
#
	FTFrame <- tkframe(tt, relief="groove", borderwidth=2)
	fixed.entry <- tkentry(FTFrame, textvariable=fixedvar)
	tkgrid(tklabel(FTFrame,text="- Fixed table -", foreground="blue"), columnspan=2)
	tkgrid(tklabel(FTFrame,text="Fixed table : "), fixed.entry)
	tkgrid(FTFrame)
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
	# Check that the object is not empty and get its name
	#
		if (tclvalue(objectvar) != "") {
			object  <- parse(text=tclvalue(objectvar))[[1]]
		} else object <- NULL
	#
	# Check that the number of permutations is not empty and get it
	#
		if (tclvalue(npermvar) != "") {
			nperm  <- parse(text=tclvalue(npermvar))[[1]]
		} else nperm <- 999
	#
	# fixed
	#
		if (tclvalue(fixedvar) != "") {
			fixed  <- parse(text=tclvalue(fixedvar))[[1]]
		} else fixed <- 0
	#
	# Make the command line
	#
		substitute(randtest(xtest = object, nperm, fixed))
	}
		
################################
# Function to reset all dialog elements to default values
################################
	"reset" <- function()
	{
		tclvalue(randtestvar) <- "untitled"
		tclvalue(objectvar) <- ""
		tclvalue(npermvar) <- "999"
		tclvalue(fixedvar) <- "0"
	}
	
################################
# Function to launch computations
################################
	"execcomp" <- function()
	{
	#
	# Check that the object name is not empty and set it
	#
		if (tclvalue(randtestvar) == "") tkinsert(randtest.entry, "end", "untitled1")
		randtestname <- parse(text=paste("\"",tclvalue(randtestvar)[[1]],"\"",sep=""))
		if (tclvalue(npermvar) == "") tkinsert(nperm.entry, "end", "999")
		nperm <- parse(text=paste("\"",tclvalue(npermvar)[[1]],"\"",sep=""))
	#
	# Build and display the command line so that the user can check it
	#
		cmd <- build()
		if (show) {
			#
			# Echoe the command line to the console
			#
			pr1 <- substr(options("prompt")$prompt, 1,2)
			cat("plot(", eval(randtestname), " <- ", deparse(cmd, width.cutoff = 256), ")\n", pr1, sep="")
		}
	#
	# Execute the command
	#
		assign("ade4TkGUIFlag", 1, envir=env_ade4tkgui)
		myObject <- eval.parent(cmd)
		assign(eval(randtestname), myObject, envir=env_ade4tkgui)
		plot(myObject)
		rm("ade4TkGUIFlag", envir=env_ade4tkgui)
		if (history) {
			commande = paste("plot(", eval(randtestname), " <- ", deparse(cmd, width.cutoff = 500), ")", sep = "")
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
