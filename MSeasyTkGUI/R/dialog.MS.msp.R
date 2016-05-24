# MSeasy Tcl/Tk GUI for some basic functions in the MSeasy package.
# Part of the code is based on the ade4TkGUI package by Jean Thioulouse <jthioulouse@biomserv.univ-lyon1.fr>, Stephane
#       Dray <dray@biomserv.univ-lyon1.fr>
#

dialog.MS.msp <-
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
		tclvalue(typesvar)="{{Text files} {.txt}}"
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
	tkwm.title(tt,"MSeasyToMSP ")	
#
# Variables for text fields
#

	
	msclustvar<-tclVar("Out_msp-demo")
	mfile<-system.file("doc/Output_examples/output_cluster21.txt",package="MSeasy")
	dfvar <- tclVar(mfile)
	clustvart <- tclVar()
	

#
# Variable for number fields
#

	nctotvar <-tclVar("all")
#
# Checkboxes
#
	autonfvar <- tclVar(0)	

	}
	else{
	op=options()
	options(warn=-1)
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"MSeasyToMSP ")	
#
# Variables for text fields
#

	
	msclustvar<-tclVar()
	dfvar <- tclVar()
	clustvart <- tclVar()
#
# Variable for number fields
#

	nctotvar <-tclVar("all")
	
#
# Checkboxes
#
	autonfvar <- tclVar(0)
	
	}
	
#
# Title 
#
	TFrame <- tkframe(tt, relief="groove")
	labh <- tklabel(TFrame, bitmap="questhead")
	tkgrid(tklabel(TFrame,text="Export mass spectra to MSP format for NIST search", font="Times 14", foreground="red"), labh)
	tkbind(labh, "<Button-1>", function() print(help("MSeasyToMSP")))
	tkgrid(TFrame, sticky="we")
#
# MSeasyToARISTO

	IOFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(IOFrame,text="- MSeasyToMSP Input & Output files-", foreground="blue"), columnspan=5)
	file.entry <- tkentry(IOFrame, textvariable=dfvar)
	choosefile.but <- tkbutton(IOFrame, text="Choose file", command=function() choosefile())
	tkgrid(tklabel(IOFrame,text="Info: input files are output files from MS.clust (output_cluster or output_peak)"),labh, columnspan=2)
	tkgrid(tklabel(IOFrame,text="File to read (txt, csv..): "), file.entry, choosefile.but)
	msclust.entry <- tkentry(IOFrame, textvariable=msclustvar)
	tkgrid(tklabel(IOFrame,text="Output name : "), msclust.entry, sticky="n")
	tkgrid(IOFrame, sticky="we")
	


	
#
# Number of clusters
#
	clVFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(clVFrame,text="                                    - Option cluster=X -", foreground="blue"), sticky="n")
	nctot.entry <- tkentry(clVFrame, textvariable=nctotvar, width=8, state="normal")
	tkgrid(tklabel(clVFrame,text="Selected cluster(s) to export (e.g all or 1,5,60) : "), nctot.entry)
	tkgrid(clVFrame, sticky="we")

#
# autosearch
#	
	quantFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(quantFrame,text="                       - Automatic Search on NIST ms search -", foreground="blue"), columnspan=2)
	quantnf.cbut <- tkcheckbutton(quantFrame,text="Check for automatic launch of NIST ms search", variable=autonfvar  )
	tkgrid(quantnf.cbut, sticky="we", columnspan=2)
	tkgrid(quantFrame, sticky="we", columnspan=2)	

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
		
		tclvalue(nctotvar)<-"all"
		tclvalue(dfvar)<-""
		tclvalue(msclustvar)<-"ForMSP"
	
	
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
	if(tclvalue(autonfvar )==0){	
		if(tclvalue(nctotvar)=="all"){
		substitute(MSeasyToMSP(filename=paste(df)))
		}
		else
		{
				substitute(MSeasyToMSP(filename=paste(df),cluster=as.integer(as.vector(as.numeric(unlist(strsplit(tclvalue(nctotvar),",")))))))
		}
	}
	else
	{
		if(tclvalue(nctotvar)=="all"){
		substitute(MSeasyToMSP(filename=paste(df), autosearch=TRUE))
		}
		else
		{
		substitute(MSeasyToMSP(filename=paste(df),cluster=as.integer(as.vector(as.numeric(unlist(strsplit(tclvalue(nctotvar),","))))), autosearch=TRUE))
		}
	}	
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

