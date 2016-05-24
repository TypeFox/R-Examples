# MSeasy Tcl/Tk GUI for some basic functions in the MSeasy package.
# Part of the code is based on the ade4TkGUI package by Jean Thioulouse <jthioulouse@biomserv.univ-lyon1.fr>, Stephane
#       Dray <dray@biomserv.univ-lyon1.fr>
#
choosepackage <-
function()
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose MSeasy data set")
	
	done <- tclVar(0)
	
	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	

	cancel.but <- tkbutton(frame1, text="Cancel", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")
	
	lmseasy <- data(package="MSeasy")
	for (i in seq(1, length(lmseasy$result)/4)) {
		dsname <- lmseasy$result[i,3]
		tkinsert(tlb, "end", dsname)
	}

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tlb, "<Button-3>", function() if (tclvalue(tkcurselection(tlb)) != "") print(help(tclvalue(tkget(tlb, tclvalue(tkcurselection(tlb)))))))
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
	
	numc <- tclvalue(tkcurselection(tlb))
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
	#data(list=choix)
	#data(list=choix, package="MSeasy", envir=parent.frame())
	
	tkdestroy(tf)
	q <- tkmessageBox(icon="info", title="Data set loaded", type="ok", 
		message=paste("The \"",choix,"\" data set has been successfully loaded.", sep=""))
}

