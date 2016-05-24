# MSeasy Tcl/Tk GUI for some basic functions in the MSeasy package.
# Part of the code is based on the ade4TkGUI package by Jean Thioulouse <jthioulouse@biomserv.univ-lyon1.fr>, Stephane
#       Dray <dray@biomserv.univ-lyon1.fr>
#
dialog.MS.step3 <-
function()
{
	tf <- tktoplevel()
	tkwm.title(tf,"Step3- Export results to... ")
	
	done <- tclVar(0)
	sepVar <- tclVar(1)
	

	
	
	"choosefic" <- function()
	{	
		
		sep <- tclvalue(sepVar)
		if (sep == 1) sepch <- "aristo"
		if (sep == 2) sepch <- "msp"
		if (sep == 3) sepch <- "nist"
		
	
		

			rdcom <- paste("dialog.MS.", sepch, "()",sep="")
		
			eval(parse(text=rdcom), envir=as.environment("e1"))
		    tkdestroy(tf)
		
	}
	
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	labh <- tklabel(frame1, bitmap="questhead")
	tkbind(labh, "<Button-1>", function() print(help("MSeasy")))
	tkgrid(tklabel(frame1,text="Step3- Export results to...", font="Times 18", foreground="red"), labh, columnspan=2)
	tkpack(frame1, fill = "x")

		
	sepFrame <- tkframe(tf, relief="groove", borderwidth=2)
    tkgrid(tklabel(sepFrame, text="Options:", foreground="blue"))
    tkgrid(tkradiobutton(sepFrame, text="ARISTO format", value=1, variable=sepVar), sticky="w")
    tkgrid(tkradiobutton(sepFrame, text="MSP format (for NIST)", value=2, variable=sepVar), sticky="w")
    tkgrid(tkradiobutton(sepFrame, text="Search NIST (Windows only!)", value=3, variable=sepVar), sticky="w")
	tkpack(sepFrame, fill = "x")
	
    


	ok.but <- tkbutton(tf, text="Submit", command=function() choosefic())
	cancel.but <- tkbutton(tf, text="Cancel", command=function() tkdestroy(tf))
	tkpack(cancel.but, ok.but, side="left", fill="x", expand=1)

	tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
	tkbind(tf, "<KeyPress-Return>", function() choosefic())
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	tkwait.variable(done)
	if(tclvalue(done) == "2") return(0)
	tkdestroy(tf)
}

