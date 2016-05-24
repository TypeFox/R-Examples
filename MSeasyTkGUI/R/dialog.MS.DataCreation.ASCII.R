# MSeasy Tcl/Tk GUI for some basic functions in the MSeasy package.
# Part of the code is based on the ade4TkGUI package by Jean Thioulouse <jthioulouse@biomserv.univ-lyon1.fr>, Stephane
#       Dray <dray@biomserv.univ-lyon1.fr>
#
dialog.MS.DataCreation.ASCII <-
function()
{

	if (exists("DemoFlag", envir=as.environment("e1"))) {
	#
		tkmessageBox(title="Launch Demonstration",message="To launch the demonstration as is, \n just click on the Submit button \n in the next window",icon="info",type="ok")

	op=options()
	options(warn=-1)
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"MS.DataCreation From ASCII ")	
#
# Variables for text fields
#
	pathdemo<-paste(system.file("doc/ASCII_TransASCII",package="MSeasy"),"/", sep="")
	rbValue <- tclVar()
	dirvar <- tclVar(pathdemo)
	outvar <- tclVar("Out_ASCIIdemo")
	
	
#
# Variable for number fields
#
	mzminvar <- tclVar(30)
	mzmaxvar <- tclVar(250)
	nfiltvar <- tclVar(3)

	}
	else{
	op=options()
	options(warn=-1)
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"MS.DataCreation From ASCII ")	
#
# Variables for text fields
#
	rbValue <- tclVar()
	dirvar <- tclVar()
	outvar <- tclVar()
	
#
# Variable for number fields
#
	mzminvar <- tclVar(30)
	mzmaxvar <- tclVar(250)
	nfiltvar <- tclVar(3)

	}
	
	
#
# Title
#
	TFrame <- tkframe(tt, relief="groove")
	labh <- tklabel(TFrame, bitmap="questhead")
	tkgrid(tklabel(TFrame,text="Data Creation from ASCII", font="Times 18", foreground="red"), labh)
	tkbind(labh, "<Button-1>", function() print(help("MS.DataCreation")))
	tkgrid(TFrame)
#
# DataCreation
#	
	IOFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(IOFrame,text="- MS.DataCreation Input & Output files -", foreground="blue"), columnspan=5, sticky="e")
	MS.DataCreation.entry <- tkentry(IOFrame, textvariable=outvar)
	dir.entry <- tkentry(IOFrame, textvariable=dirvar)
	choosedir.but <- tkbutton(IOFrame, text="Set path", command=function() tkinsert(dir.entry, "end", tclvalue(tkchooseDirectory())))
	tkgrid(tklabel(IOFrame,text="Data path : "), dir.entry, choosedir.but, sticky="w")
	tkgrid(tklabel(IOFrame,text="Output name : "), MS.DataCreation.entry, sticky="w")
	tkgrid(IOFrame, sticky="we")
#
# Number of mz
#
	MZFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(MZFrame,text="- Range of mass fragments -", foreground="blue"), columnspan=2)
	mzmin.entry <- tkentry(MZFrame, textvariable=mzminvar, width=4, state="normal")
	mzmax.entry <- tkentry(MZFrame, textvariable=mzmaxvar, width=4, state="normal")
	tkgrid(tklabel(MZFrame,text="mz (min) : "), mzmin.entry,tklabel(MZFrame,text="mz (max) : "),mzmax.entry)
	tkgrid(MZFrame, sticky="we")

#
# N_Filt value for smoothing
#
	NFFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(NFFrame,text="- Peak detection & smoothing parameter -", foreground="blue"), columnspan=2)
	nfilt.entry <- tkentry(NFFrame, textvariable=nfiltvar, width=4, state="normal")
	tkgrid(tklabel(NFFrame,text="N_filt : "), nfilt.entry, sticky="n")
	tkgrid(NFFrame, sticky="we")	
	
#
# apex type
#	
	ApexFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(ApexFrame,text="- Mass spectrum selection -", foreground="blue"), columnspan=3)
	rb1 <- tkradiobutton(ApexFrame)
	rb2 <- tkradiobutton(ApexFrame)
	rbValue <- tclVar(TRUE)
	tkconfigure(rb1,variable=rbValue,value=TRUE)
	tkconfigure(rb2,variable=rbValue,value=FALSE)
	tkgrid(tklabel(ApexFrame,text="Apex only (apex=T)"),rb1)
	tkgrid(tklabel(ApexFrame,text="Mean of 5% mass spectra around the apex (apex=F)"),rb2)
	tkgrid(ApexFrame, sticky="we")
#
# Local variables
#
	
	done <- tclVar(0)	# To terminate the dialog
	

		
################################
# Function to reset all dialog elements to default values
################################
	"reset" <- function()
	{
		tclvalue(dirvar)<-""
		tclvalue(mzminvar)<-30
		tclvalue(mzmaxvar)<-250
		tclvalue(nfiltvar)<-3
		tclvalue(rbValue)<-TRUE
	
	
	}
	
################################
# Function to launch computations
################################
	"execcomp" <- function()
	{
	path<-as.character(tclvalue(dirvar))
	mz<-tclvalue(mzminvar):tclvalue(mzmaxvar)
	apex<-tclvalue(rbValue)
	N_filt<-tclvalue(nfiltvar)
	
    mbox <- tkProgressBar(title = "R progress bar", label = "Processing ...", min = 0, max = 1, initial = 0, width = 300)
	Sys.sleep(0.5)
	tkdestroy(tt)
	
	print("Processing...")
    setTkProgressBar(mbox, 0.15, title = "Processing", label = "Please wait ...")	
	if (apex=="0"){
	      setTkProgressBar(mbox, 0.25, title = "R progress bar", label = "Processing ...")
		    Sys.sleep(0.5)
		outtrans<-trans.ASCII(path,mz)
		path2<-paste("output",lapply(strsplit(outtrans, "output"),"[",2), sep="")
		
		assign(tclvalue(outvar), MS.DataCreation(path=path2,mz=mz,DataType="ASCII", N_filt=N_filt, apex=TRUE), envir=as.environment("e1"))
#Progress bar
		    setTkProgressBar(mbox, 1, title = "100 % Done", label = "100 % Done")
			Sys.sleep(0.10)
			close(mbox)
			tkmessageBox(title="check",message= "Done",icon="info",type="ok")
		
	}else{
		setTkProgressBar(mbox, 0.25, title = "R progress bar", label = "Processing ...")
		Sys.sleep(0.5)
		outtrans<-trans.ASCII(path,mz)
		setTkProgressBar(mbox, 0.5, title = "R progress bar", label = "Processing ...")
		path2<-paste("output",lapply(strsplit(outtrans, "output"),"[",2), sep="")
		assign(tclvalue(outvar), MS.DataCreation(path=path2,mz=mz,DataType="ASCII", N_filt=N_filt, apex=FALSE), envir=as.environment("e1"))
#Progress bar
			setTkProgressBar(mbox, 1, title = "100 % Done", label = "100 % Done")
			Sys.sleep(0.10)
			close(mbox)
			tkmessageBox(title="Done",message= "Done",icon="info",type="ok")
	}
	
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

