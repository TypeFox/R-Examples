# MSeasy Tcl/Tk GUI for some basic functions in the MSeasy package.
# Part of the code is based on the ade4TkGUI package by Jean Thioulouse <jthioulouse@biomserv.univ-lyon1.fr>, Stephane
#       Dray <dray@biomserv.univ-lyon1.fr>
#
dialog.MS.clust <-
function()
{
if (exists("DemoFlag", envir=as.environment("e1"))) {
	#
	tkmessageBox(title="Launch Demonstration",message="To launch the demonstration, select: \n clusMeth=hierarchical \n disMeth=euclidean \n linkMeth=ward \n in the next window \n Then click on the Submit button",icon="info",type="ok")
	eval(Agilent_quantF_MSclust, envir=as.environment("e1"))
	op=options()
	options(warn=-1)
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"MS.clust ")	
#
# Variables for text fields
#

	dirvar <- tclVar()
	msclustvar<-tclVar("Out_MS.clustdemo")
	dfvar <- tclVar("Agilent_quantF_MSclust")
	clustvar1<-tclVar()
	
#
# Variable for number fields
#
	ncminvar <- tclVar(10)
	ncmaxvar <- tclVar(50)
	nctotvar <-tclVar(21)
	varRTvar <-tclVar(0.1)
	
#
# Checkboxes
#

	scannfvar <- tclVar(0)
	quantnfvar <- tclVar(0)  #for quant=TRUE/FALSE option in MS.Clust with Agilent files
#
# liste
#	
	
	
	distancevar <- tclVar()
	clustvar <- tclVar()
	agglovar<-tclVar()
	}
	else{
	op=options()
	options(warn=-1)
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"MS.clust ")	
#
# Variables for text fields
#

	dirvar <- tclVar()
	msclustvar<-tclVar()
	dfvar <- tclVar()
	clustvar1<-tclVar()
	
#
# Variable for number fields
#
	ncminvar <- tclVar(10)
	ncmaxvar <- tclVar(50)
	nctotvar <-tclVar(21)
	varRTvar <-tclVar(0.1)
	
#
# Checkboxes
#

	scannfvar <- tclVar(1)
	quantnfvar <- tclVar(1)  #for quant=TRUE/FALSE option in MS.Clust with Agilent files
#
# liste
#	
	distancevar <- tclVar()
	clustvar <- tclVar()
	agglovar<-tclVar()
	}
	
#
# Title 
#
	TFrame <- tkframe(tt, relief="groove")
	labh <- tklabel(TFrame, bitmap="questhead")
	tkgrid(tklabel(TFrame,text="Clustering of Mass Spectra", font="Times 18", foreground="red"), labh)
	tkbind(labh, "<Button-1>", function() print(help("MS.clust ")))
	tkgrid(TFrame)
#
# MS.clust
#	

	IOFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(IOFrame,text="- MS.clust Input & Output files-", foreground="blue"), columnspan=5)
	msclust.entry <- tkentry(IOFrame, textvariable=msclustvar)
	df.entry <- tkentry(IOFrame, textvariable=dfvar)
	dfnr.label <- tklabel(IOFrame, width=4)
	dfnc.label <- tklabel(IOFrame, width=4)
	choosedf.but <- tkbutton(IOFrame, text="Set", command=function() choosedf(df.entry, dfnr.label, dfnc.label))
	tkgrid(tklabel(IOFrame,text="Input data frame : "), df.entry, choosedf.but, dfnr.label, dfnc.label, sticky="w")
	tkgrid(tklabel(IOFrame,text="Output name : "), msclust.entry, sticky="w")
	tkgrid(IOFrame, sticky="we")
	
#
# Cluster validation 
#	
	clVFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(clVFrame,text="                      - Find optimal number(s) of clusters ?-", foreground="blue"), columnspan=2)
	scannf.cbut <- tkcheckbutton(clVFrame,text="Check for Yes or Uncheck for No", variable=scannfvar,
		command=function() if (tclvalue(scannfvar)==1) {tkconfigure(ncmin.entry, state="normal") ; tkconfigure(ncmax.entry, state="normal") ; tkconfigure(nctot.entry, state="disabled")}else {tkconfigure(ncmin.entry, state="disabled") ; tkconfigure(ncmax.entry, state="disabled") ; tkconfigure(nctot.entry, state="normal")})
	tkgrid(scannf.cbut, sticky="we", columnspan=2)
	
	
#
# Number of clusters
#
	listFrame <- tkframe(clVFrame, relief="groove", borderwidth=2)
	tkgrid(tklabel(clVFrame,text="IF Yes", foreground="red"), columnspan=2, sticky="w")
	ncmin.entry <- tkentry(clVFrame, textvariable=ncminvar, width=4, state="normal")
	ncmax.entry <- tkentry(clVFrame, textvariable=ncmaxvar, width=4, state="normal")
	nctot.entry <- tkentry(clVFrame, textvariable=nctotvar, width=4, state="disabled")
	tkgrid(tklabel(clVFrame,text="ncmin : "), ncmin.entry,tklabel(clVFrame,text="     ncmax : "),ncmax.entry,sticky="e")
	tkgrid(tklabel(clVFrame,text="IF No", foreground="red"), columnspan=2, sticky="w")
	
	tkgrid(tklabel(clVFrame,text="Number of clusters (Nbc) : "), nctot.entry, sticky="e")
	tkgrid(clVFrame, sticky="we")
	tkgrid(listFrame, sticky="we")

	
#
#list
#	

modify_command1 <- function() {
		clustvar <- clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1]
		if (clustvar=="hierarchical"){
				tkconfigure(comboBox2, state="normal") 
				tkconfigure(comboBox, state="normal")
			}
			else {
				tkconfigure(comboBox2, state="disabled")
				tkconfigure(comboBox, state="disabled")
			}
		
		
		if ((clustvar=="pam")||(clustvar=="diana")||(clustvar=="hierarchical")){
				tkconfigure(comboBox, state="normal") 
				if ((clustvar=="pam")||(clustvar=="diana")){
					tkmessageBox(title="Warning",message="Only Euclidean or Manhattan as disMeth",icon="info",type="ok")
				}
			} 
			else { 
				tkconfigure(comboBox, state="disabled")
				
				}
		if ((clustvar=="kmeans"))
			tkconfigure(comboBox4, state="normal") else tkconfigure(comboBox4, state="disabled")
	}

	
	
	tclRequire("BWidget")
	listFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(listFrame,text="                            - Clustering method parameters-", foreground="blue"), columnspan=2, sticky="n")

#clusMeth	
#clustvar<-"pam"
	clustls <- c("hierarchical","diana","kmeans","pam")
	comboBoxclust <- tkwidget(listFrame,text="Method (clustMeth)","ComboBox",editable=FALSE,"-modifycmd", 
                    modify_command1,values=clustls)
	clustvar <- clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1]
	tkgrid(comboBoxclust,sticky="e")
	
#disMeth	
	distancels <- c("euclidean","correlation","manhattan")
	comboBox <- tkwidget(listFrame,text="Distance (disMeth)","ComboBox",editable=FALSE,values=distancels, state="disabled")
	distancevar <- distancels[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]
	tkgrid(comboBox, sticky="e")
	
#linkMeth	
	agglols <- c("ward","single","complet","average", "centroid")
	comboBox2 <- tkwidget(listFrame,text="Linkage (linkMeth)","ComboBox",editable=FALSE,values=agglols,state="disabled")
	agglovar <- agglols[as.numeric(tclvalue(tcl(comboBox2,"getvalue")))+1]
	tkgrid(comboBox2, sticky="e")
	tkgrid(listFrame, sticky="we")
	
#
# varRT (accepted shift in RT)
#	
	varRTFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(varRTFrame,text="                       - Quality controls of putative molecules -", foreground="blue"), columnspan=2)
	varRT.entry <- tkentry(varRTFrame, textvariable=varRTvar, width=4, state="normal")
	tkgrid(tklabel(varRTFrame,text="Accepted shift in RT/RI (min) : "), varRT.entry, sticky="w")
	
#DisMeth if kmeans or pam	
	qualdiss <- c("euclidean","correlation","manhattan")
	comboBox4 <- tkwidget(varRTFrame,text="Distance for Silhouette","ComboBox",editable=FALSE,values=qualdiss,state="disabled")
	qualdissvar <- qualdiss[as.numeric(tclvalue(tcl(comboBox4,"getvalue")))+1]
	tkgrid(comboBox4, sticky="e")
	tkgrid(varRTFrame, sticky="we")	
#
# quant variable
#	
	quantFrame <- tkframe(tt, relief="groove", borderwidth=2)
	tkgrid(tklabel(quantFrame,text="                       - Additional output -", foreground="blue"), columnspan=2)
	quantnf.cbut <- tkcheckbutton(quantFrame,text="Check box to add quantification measures of peak size (->profiling matrix)", variable=quantnfvar )
	tkgrid(quantnf.cbut, sticky="we", columnspan=2)
	tkgrid(quantFrame, sticky="we", columnspan=2)

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
		tclvalue(ncminvar)<-10
		tclvalue(ncmaxvar)<-50
		tclvalue(nctotvar)<-21
		tclvalue(dfvar)<-""
		tclvalue(scannfvar)<-"1"
		tclvalue(quantnfvar)<-"1"
		tclvalue(varRTvar)<-0.1
		tclvalue(msclustvar)<-""
		tkconfigure(comboBox2, state="disabled")
		tkconfigure(comboBox4, state="disabled")
		tkconfigure(comboBox, state="disabled")
		tkconfigure(dfnr.label, text="")
		tkconfigure(dfnc.label, text="")
		tclvalue(clustvar1)<-""
	
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
	clustvar1<-clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1]
	if (tclvalue(quantnfvar)==1) {
	
		if (clustvar1=="hierarchical"){
		
			if (tclvalue(scannfvar)==1) substitute(MS.clust(df,clV=T,ncmin=tclvalue(ncminvar),ncmax=tclvalue(ncmaxvar), varRT=tclvalue(varRTvar), disMeth=distancels[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1],linkMeth=agglols[as.numeric(tclvalue(tcl(comboBox2,"getvalue")))+1],clustMeth=clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1], quant=T))
			else  substitute(MS.clust(df,clV=F,Nbc=tclvalue(nctotvar), varRT=tclvalue(varRTvar), disMeth=distancels[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1],linkMeth=agglols[as.numeric(tclvalue(tcl(comboBox2,"getvalue")))+1],clustMeth=clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1], quant=T))
		}
		else{
			if(clustvar1=="kmeans"){
				if (tclvalue(scannfvar)==1) substitute(MS.clust(df,clV=T,ncmin=tclvalue(ncminvar),ncmax=tclvalue(ncmaxvar), varRT=tclvalue(varRTvar), disMeth=qualdiss[as.numeric(tclvalue(tcl(comboBox4,"getvalue")))+1],linkMeth=NULL,clustMeth=clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1], quant=T))
				else  substitute(MS.clust(df,clV=F,Nbc=tclvalue(nctotvar), varRT=tclvalue(varRTvar), disMeth=qualdiss[as.numeric(tclvalue(tcl(comboBox4,"getvalue")))+1],linkMeth=NULL,clustMeth=clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1], quant=T))
			}
			else{
				if (tclvalue(scannfvar)==1) substitute(MS.clust(df,clV=T,ncmin=tclvalue(ncminvar),ncmax=tclvalue(ncmaxvar), varRT=tclvalue(varRTvar), disMeth=distancels[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1],linkMeth=NULL,clustMeth=clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1], quant=T))
				else  substitute(MS.clust(df,clV=F,Nbc=tclvalue(nctotvar), varRT=tclvalue(varRTvar), disMeth=distancels[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1],linkMeth=NULL,clustMeth=clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1], quant=T))
	
			}
		} 
		
	}	
	 else{
	 if (clustvar1=="hierarchical"){
		
			if (tclvalue(scannfvar)==1) substitute(MS.clust(df,clV=T,ncmin=tclvalue(ncminvar),ncmax=tclvalue(ncmaxvar), varRT=tclvalue(varRTvar), disMeth=distancels[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1],linkMeth=agglols[as.numeric(tclvalue(tcl(comboBox2,"getvalue")))+1],clustMeth=clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1], quant=F))
			else  substitute(MS.clust(df,clV=F,Nbc=tclvalue(nctotvar), varRT=tclvalue(varRTvar), disMeth=distancels[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1],linkMeth=agglols[as.numeric(tclvalue(tcl(comboBox2,"getvalue")))+1],clustMeth=clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1], quant=F))
		}
		else{
			if(clustvar1=="kmeans"){
				if (tclvalue(scannfvar)==1) substitute(MS.clust(df,clV=T,ncmin=tclvalue(ncminvar),ncmax=tclvalue(ncmaxvar), varRT=tclvalue(varRTvar), disMeth=qualdiss[as.numeric(tclvalue(tcl(comboBox4,"getvalue")))+1],linkMeth=NULL,clustMeth=clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1], quant=F))
				else  substitute(MS.clust(df,clV=F,Nbc=tclvalue(nctotvar), varRT=tclvalue(varRTvar), disMeth=qualdiss[as.numeric(tclvalue(tcl(comboBox4,"getvalue")))+1],linkMeth=NULL,clustMeth=clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1], quant=F))
			}
			else{
				if (tclvalue(scannfvar)==1) substitute(MS.clust(df,clV=T,ncmin=tclvalue(ncminvar),ncmax=tclvalue(ncmaxvar), varRT=tclvalue(varRTvar), disMeth=distancels[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1],linkMeth=NULL,clustMeth=clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1], quant=F))
				else  substitute(MS.clust(df,clV=F,Nbc=tclvalue(nctotvar), varRT=tclvalue(varRTvar), disMeth=distancels[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1],linkMeth=NULL,clustMeth=clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1], quant=F))
	
			}
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
			pr1 <- substr(options("prompt")$prompt, 1,2)
			#cat(eval(dudiname), " <- ", deparse(cmd, width = 256), cat("\n"), pr1, sep="")
	#
	# Execute the command
	#
		
		mydudi <- eval.parent(cmd)
        #bringToTop(which=-1)
		assign(eval(dudiname), mydudi, envir=as.environment("e1"))
	
	}

"execcomp1" <- function()
	{		
	print(paste(cat("\n")," Processing..."))
	execcomp()
	#tkmessageBox(title="check",message= paste("disMeth=",distancels[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]," \n linkMeth=" , agglols[as.numeric(tclvalue(tcl(comboBox2,"getvalue")))+1]," \n clustMeth=" ,clustls[as.numeric(tclvalue(tcl(comboBoxclust,"getvalue")))+1]),icon="info",type="ok")
		tkconfigure(comboBox2, state="disabled")
		tkconfigure(comboBox4, state="disabled")
		tkconfigure(comboBox, state="disabled")
		tkconfigure(comboBoxclust, state="disabled")
		tkconfigure(ncmin.entry, state="disabled")
		tkconfigure(ncmax.entry, state="disabled")
		tkconfigure(nctot.entry, state="disabled")
		tkconfigure(submit.but, state="disabled")
		tkconfigure(reset.but, state="disabled")
		
		mbox <- tkProgressBar(title = "R progress bar", label = "Processing ...", min = 0, max = 1, initial = 0, width = 300)
		
		Sys.sleep(0.5)
		setTkProgressBar(mbox, .5, title ="R progress bar", label = "Please wait ...")
		Sys.sleep(0.5)
		print("...")
		setTkProgressBar(mbox, .15, title ="R progress bar", label = "Processing ...")	
		tkfocus(tt)
		tkdestroy(tt)
		
		setTkProgressBar(mbox, 1, title = "100% Done", label = "100% Done")
		Sys.sleep(0.5)
		close(mbox)
		tkmessageBox(title="Done",message= "Done",icon="info",type="ok")
		tkdestroy(tt)
		#if (tclvalue(scannfvar)==1)return(0)
	}
#
# Reset Cancel and Submit buttons
#
	RCSFrame <- tkframe(tt, relief="groove")
	reset.but <- tkbutton(RCSFrame, text="Reset", command=reset)
	cancel.but <- tkbutton(RCSFrame, text="Cancel", command=function() tkdestroy(tt))
	submit.but <- tkbutton(RCSFrame, text="Submit", default="active", command=function() execcomp1())
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

