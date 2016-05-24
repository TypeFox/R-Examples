# MSeasy Tcl/Tk GUI for some basic functions in the MSeasy package.
# Part of the code is based on the ade4TkGUI package by Jean Thioulouse <jthioulouse@biomserv.univ-lyon1.fr>, Stephane
#       Dray <dray@biomserv.univ-lyon1.fr>
#

dialog.MS.nist <-
function()
{
	#done <- tclVar(0)
	typesvar <- tclVar()
	#tabnamevar <- tclVar("")
	#filenamevar <- tclVar("clipboard")
	fictrt <- tclVar()
	#varnames <- tclVar(1)
	#sepVar <- tclVar(1)
	#otherSepVar <- tclVar("")
	#decSepVar <- tclVar(1)
	#colClassVar <- tclVar(4)

	"choosefile" <- function()
	{
		tclvalue(typesvar)="{{NIST MSP files} {.msp}}"
		fictrt <- tkgetOpenFile(filetypes=tclvalue(typesvar))
		fpath <- tclvalue(fictrt)
		tkdelete(file.entry, 0, "end")
		tkinsert(file.entry, "end", fpath)
	}
	
if (exists("DemoFlag", envir=as.environment("e1"))) {
	#
	tkmessageBox(title="Launch Demonstration",message="To launch the demonstration as is, \n just click on the Submit button \n in the next window",icon="info",type="ok")

	#eval(data(Data_testclust), envir=as.environment("e1"))
	op=options()
	options(warn=-1)
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"SearchNIST ")	
#
# Variables for text fields
#

	
	msclustvar<-tclVar("Out_NIST-demo")
	mfile<-system.file("doc/Output_examples/ForNIST.msp",package="MSeasy")
	dfvar <- tclVar(mfile)
	clustvart <- tclVar()
	


	}
	else{
	op=options()
	options(warn=-1)
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"SearchNIST ")	
#
# Variables for text fields
#

	
	msclustvar<-tclVar()
	dfvar <- tclVar()
	clustvart <- tclVar()

	

	}
	
#
# Title 
#
	TFrame <- tkframe(tt, relief="groove")
	labh <- tklabel(TFrame, bitmap="questhead")
	tkgrid(tklabel(TFrame,text="Search mass spectra with NIST ms search", font="Times 14", foreground="red"), labh)
	tkbind(labh, "<Button-1>", function() print(help("SearchNIST")))
	tkgrid(TFrame, sticky="we")
#
# MSeasyToARISTO

	IOFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(IOFrame,text="- SearchNIST Input & Output files-", foreground="blue"), columnspan=5)
	file.entry <- tkentry(IOFrame, textvariable=dfvar)
	choosefile.but <- tkbutton(IOFrame, text="Choose file", command=function() choosefile())
	tkgrid(tklabel(IOFrame,text="Info: Only *.MSP files"),labh, columnspan=2)
	tkgrid(tklabel(IOFrame,text="File to read (*.msp): "), file.entry, choosefile.but)
	msclust.entry <- tkentry(IOFrame, textvariable=msclustvar)
	tkgrid(tklabel(IOFrame,text="Output name : "), msclust.entry, sticky="n")
	tkgrid(IOFrame, sticky="we")
	
#
# Local variables
#

	numi=1				# Number of choosed element
	done <- tclVar(0)	# To terminate the dialog
	

		
################################
# Function to reset all dialog elements to default values
################################
	"reset" <- function()
	{
		
		
		tclvalue(dfvar)<-""
		tclvalue(msclustvar)<-"FromNIST"
	
	
	}
	
		"build" <- function()
	{
		if (tclvalue(dfvar) != "") {
			#df  <- parse(text=tclvalue(dfvar))[[1]]
			df<-tclvalue(dfvar)
		} else {
			return(0)
		}
		#
	# Make the command line
	#
		
		
		substitute(SearchNIST(mspfile=paste(df)))
		
		
		
	}	
################################
# Function to launch computations
################################
	"execcomp" <- function()
	{
		
		
	# Check that the analysis name is not empty and get it
	#
		if (tclvalue(msclustvar) == "") tkinsert(msclust.entry, "end", "untitled1")
		dudiname <- parse(text=paste("\"",tclvalue(msclustvar)[[1]],"\"",sep=""))
	#
	# Build and display the command line so that the user can check it
	#
		cmd <- build()
		if (cmd == 0) return(0)
		
	#
	# Execute the command
	#
		mbox <- tkProgressBar(title = "R progress bar", label = "Processing ...", min = 0, max = 1, initial = 0, width = 300)
		Sys.sleep(0.5)
		tkdestroy(tt) 
		print("Processing...")
		setTkProgressBar(mbox, .05, title = "R progress bar", label = "Processing ...")
		Sys.sleep(0.5)
		setTkProgressBar(mbox, .15, title = "R progress bar", label = "Please wait ...")
		mydudi <- eval.parent(cmd)
		assign(eval(dudiname), mydudi, envir=as.environment("e1"))
		setTkProgressBar(mbox, 1, title = "100 % Done", label = "100% Done")
		Sys.sleep(0.5)
		close(mbox)
		tkmessageBox(title="check",message= "Done",icon="info",type="ok")
		
		
	}
#
# Reset Cancel and Submit buttons
#
	RCSFrame <- tkframe(tt, relief="groove")
	reset.but <- tkbutton(RCSFrame, text="Reset", command=reset)
	cancel.but <- tkbutton(RCSFrame, text="Cancel", command=function() tkdestroy(tt))
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

