# MSeasy Tcl/Tk GUI for some basic functions in the MSeasy package.
# Part of the code is based on the ade4TkGUI package by Jean Thioulouse <jthioulouse@biomserv.univ-lyon1.fr>, Stephane
#       Dray <dray@biomserv.univ-lyon1.fr>
#
dialog.MS.DataCreation.Agilent <-
function()
{
	if (exists("DemoFlag", envir=as.environment("e1"))) {
	#
	tkmessageBox(title="Launch Demonstration",message="To launch the demonstration as is, \n just click on the Submit button \n in the next window \n If needed examples files will be downloaded from internet \n this may take a minute",icon="info",type="ok")
	op=options()
	options(warn=-1)
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"MS.DataCreation with CDF + rteres.txt (Agilent Technologies) ")	
#
# Variables for text fields
#
if(file.exists(file.path(getwd(),"Agilent_rteres"))==FALSE){
url1<-"http://sites.google.com/site/rpackagemseasy/downloads/Agilent_example.zip"
download.file(url=url1, destfile="AgilentCDF.zip")
unzip(zipfile="AgilentCDF.zip", exdir=".") 
print(paste("Two example folders have been downloaded in ", getwd(), sep=" "))
unlink("AgilentCDF.zip") 
}
	pathdemo<-paste(file.path(getwd(),"Agilent_rteres"),sep="")
	pathcdfdemo<-paste(file.path(getwd(),"Agilent_CDF"),sep="")
	rbValue <- tclVar()
	dirvar <- tclVar(pathdemo)
	dirvar1 <- tclVar(pathcdfdemo)
	outvar <- tclVar("Out_Agilentdemo")
	outvar1 <- tclVar()
	
	
#
# Variable for number fields
#
	mzminvar <- tclVar(30)
	mzmaxvar <- tclVar(250)

#
# Checkboxes
#
	mznfvar <- tclVar(0)
	quantnfvar <- tclVar(1)  #for quant=TRUE/FALSE option in MS.Clust with Agilent files

	}
	else{
	#
	op=options()
	options(warn=-1)
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"MS.DataCreation with CDF + rteres.txt (Agilent Technologies) ")	
#
# Variables for text fields
#
	rbValue <- tclVar()
	dirvar <- tclVar()
	dirvar1 <- tclVar()
	outvar <- tclVar()
	outvar1 <- tclVar()
	
#
# Variable for number fields
#
	mzminvar <- tclVar(30)
	mzmaxvar <- tclVar(250)
	
#
# Checkboxes
#
	mznfvar <- tclVar(0)
	quantnfvar <- tclVar(1)  #for quant=TRUE/FALSE option in MS.Clust with Agilent files
}

#
# Title
#
	TFrame <- tkframe(tt, relief="groove")
	labh <- tklabel(TFrame, bitmap="questhead")
	tkgrid(tklabel(TFrame,text="Data Creation with CDF + rteres.txt (Agilent Technologies) ", font="Times 18", foreground="red"), labh)
	tkbind(labh, "<Button-1>", function() print(help("MS.DataCreation")))
	tkgrid(TFrame)
#
# DataCreation
#	
	IOFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(IOFrame,text="- MS.DataCreation Input & Output  files -", foreground="blue"), columnspan=5)
	MS.DataCreation.entry <- tkentry(IOFrame, textvariable=outvar)
	dir.entry <- tkentry(IOFrame, textvariable=dirvar)
	choosedir.but <- tkbutton(IOFrame, text="Set path to the directory containing sample folders", command=function() tkinsert(dir.entry,"end", tclvalue(tkchooseDirectory())))
	tkgrid(tklabel(IOFrame,text="Data path : "), dir.entry, choosedir.but, sticky="w")
	MS.DataCreation.entry <- tkentry(IOFrame, textvariable=outvar1)
	dir1.entry <- tkentry(IOFrame, textvariable=dirvar1)
	choosedir1.but <- tkbutton(IOFrame, text="Set path to the directory containing CDF files", command=function() tkinsert(dir1.entry, "end", tclvalue(tkchooseDirectory())))
	tkgrid(tklabel(IOFrame,text="CDF path : "), dir1.entry, choosedir1.but, sticky="w")
	tkgrid(tklabel(IOFrame,text="Output name : "), MS.DataCreation.entry, sticky="w")
	tkgrid(IOFrame, sticky="we")
#
# Number of mz
#
	MZFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(MZFrame,text="- Range of mass fragments -", foreground="blue"), columnspan=2)
	scannf.cbut <- tkcheckbutton(MZFrame,text="Check for all mz", variable=mznfvar,
		command=function() if (tclvalue(mznfvar)==1) {tkconfigure(mzmin.entry, state="disabled") ; tkconfigure(mzmax.entry, state="disabled")}else {tkconfigure(mzmin.entry, state="normal") ; tkconfigure(mzmax.entry, state="normal")})
	tkgrid(scannf.cbut, sticky="we", columnspan=2)
	mzmin.entry <- tkentry(MZFrame, textvariable=mzminvar, width=4, state="normal")
	mzmax.entry <- tkentry(MZFrame, textvariable=mzmaxvar, width=4, state="normal")
	tkgrid(tklabel(MZFrame,text="mz (min) : "), mzmin.entry,tklabel(MZFrame,text="mz (max) : "),mzmax.entry, sticky="w")
	tkgrid(MZFrame, sticky="we")

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
	tkgrid(tklabel(ApexFrame,text="Mean of 5% scan around the apex (apex=F)"),rb2)
	tkgrid(ApexFrame, sticky="we")

#
# quant variable
#	
	quantFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(quantFrame,text="                       - Additional output -", foreground="blue"), columnspan=2)
	quantnf.cbut <- tkcheckbutton(quantFrame,text="Check box to add quantification measures of peak size (-> profiling matrix)", variable=quantnfvar )
	tkgrid(quantnf.cbut, sticky="we", columnspan=2)
	tkgrid(quantFrame, sticky="we", columnspan=2)
	
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
		tclvalue(dirvar1)<-""
		tclvalue(mzminvar)<-30
		tclvalue(mzmaxvar)<-250
		tclvalue(rbValue)<-TRUE
		tclvalue(quantnfvar)<-"1"
	
	}
	
################################
# Function to launch computations
################################
	"execcomp" <- function()
	{
	require(xcms)
	path<-as.character(tclvalue(dirvar))
	pathCDF<-as.character(tclvalue(dirvar1))
	mz<-tclvalue(mzminvar):tclvalue(mzmaxvar)
	apex<-tclvalue(rbValue)
   # Check that the analysis name is not empty and get it
	#
		if (tclvalue(outvar1) == "") tkinsert(MS.DataCreation.entry, "end", "initial_DATA")
		dudiname <- parse(text=paste("\"",tclvalue(outvar)[[1]],"\"",sep=""))
		
	print("Processing...")
	Sys.sleep(0.5)
	tkdestroy(tt)
	if (tclvalue(quantnfvar)==1) {
		if (apex=="0"){
			
			assign(tclvalue(outvar1), MS.DataCreation(path=path,pathCDF=pathCDF,mz=mz,DataType="Agilent",apex=TRUE, quant=TRUE), envir=as.environment("e1"))
			tkmessageBox(title="check",message= "Done",icon="info",type="ok")
		}else{
			
			assign(tclvalue(outvar1), MS.DataCreation(path=path,pathCDF=pathCDF,mz=mz,DataType="Agilent",apex=FALSE, quant=TRUE), envir=as.environment("e1"))
			tkmessageBox(title="Done",message= "Done",icon="info",type="ok")
		}
		}
		else{
		if (apex=="0"){
			
			assign(tclvalue(outvar1), MS.DataCreation(path=path,pathCDF=pathCDF,mz=mz,DataType="Agilent",apex=TRUE, quant=FALSE), envir=as.environment("e1"))
			tkmessageBox(title="check",message= "Done",icon="info",type="ok")
		}else{
			
			assign(tclvalue(outvar1), MS.DataCreation(path=path,pathCDF=pathCDF,mz=mz,DataType="Agilent",apex=FALSE, quant=FALSE), envir=as.environment("e1"))
			tkmessageBox(title="Done",message= "Done",icon="info",type="ok")
		}
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

