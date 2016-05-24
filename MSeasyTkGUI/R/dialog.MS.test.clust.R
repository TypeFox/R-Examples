# MSeasy Tcl/Tk GUI for some basic functions in the MSeasy package.
# Part of the code is based on the ade4TkGUI package by Jean Thioulouse <jthioulouse@biomserv.univ-lyon1.fr>, Stephane
#       Dray <dray@biomserv.univ-lyon1.fr>
#

dialog.MS.test.clust <-
function()
{

if (exists("DemoFlag", envir=as.environment("e1"))) {
	#
	tkmessageBox(title="Launch Demonstration",message="To launch the demonstration as is, \n just click on the Submit button \n in the next window",icon="info",type="ok")

	eval(Data_testclust, envir=as.environment("e1"))
	op=options()
	options(warn=-1)
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"MS.test.clust ")	
#
# Variables for text fields
#

	dirvar <- tclVar()
	msclustvar<-tclVar("Out_test.clust-demo")
	dfvar <- tclVar("Data_testclust")
	clustvart <- tclVar()
#
# Variable for number fields
#

	nctotvar <-tclVar(10)
	

	}
	else{
	op=options()
	options(warn=-1)
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"MS.test.clust ")	
#
# Variables for text fields
#

	dirvar <- tclVar()
	msclustvar<-tclVar()
	dfvar <- tclVar()
	clustvart <- tclVar()
#
# Variable for number fields
#

	nctotvar <-tclVar()
	

	}
	
#
# Title 
#
	TFrame <- tkframe(tt, relief="groove")
	labh <- tklabel(TFrame, bitmap="questhead")
	tkgrid(tklabel(TFrame,text="Test for the best clustering method", font="Times 14", foreground="red"), labh)
	tkbind(labh, "<Button-1>", function() print(help("MS.test.clust ")))
	tkgrid(TFrame, sticky="we")
#
# MS.test.clust
#	

	IOFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(IOFrame,text="- MS.test.clust Input & Output files-", foreground="blue"), columnspan=5)
	msclust.entry <- tkentry(IOFrame, textvariable=msclustvar)
	df.entry <- tkentry(IOFrame, textvariable=dfvar)
	dfnr.label <- tklabel(IOFrame, width=4)
	dfnc.label <- tklabel(IOFrame, width=4)
	choosedf.but <- tkbutton(IOFrame, text="Set", command=function() choosedf(df.entry, dfnr.label, dfnc.label))
	tkgrid(tklabel(IOFrame,text="Input data frame : "), df.entry, choosedf.but, dfnr.label, dfnc.label, sticky="w")
	tkgrid(tklabel(IOFrame,text="Output name : "), msclust.entry, sticky="w")
	tkgrid(IOFrame, sticky="we")
	


	
#
# Number of clusters
#
	clVFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(clVFrame,text="                                    - Parameter -", foreground="blue"), sticky="n")
	nctot.entry <- tkentry(clVFrame, textvariable=nctotvar, width=4, state="normal")
	tkgrid(tklabel(clVFrame,text="Number of molecules in the dataset : "), nctot.entry)
	tkgrid(clVFrame, sticky="we")

	

#
# Local variables
#
    vnr=NULL			# Vector of dataframes row numbers
	vnc=NULL			# Vector of dataframes column numbers
	numi=1				# Number of choosed element
	done <- tclVar(0)	# To terminate the dialog
	

		
################################
# Function to reset all dialog elements to default values
################################
	"reset" <- function()
	{
		tclvalue(dirvar)<-""
		tclvalue(nctotvar)<-10
		tclvalue(dfvar)<-""
		tclvalue(msclustvar)<-""
		tkconfigure(dfnr.label, text="")
		tkconfigure(dfnc.label, text="")
	
	}
	
		"build" <- function()
	{
		if (tclvalue(dfvar) != "") {
			df  <- parse(text=tclvalue(dfvar))[[1]]
		} else {
			return(0)
		}
		#
	# Make the command line
	#
		substitute(MS.test.clust(df,nclust=as.integer(tclvalue(nctotvar))))
		
		
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

