# MSeasy Tcl/Tk GUI for some basic functions in the MSeasy package.
# Part of the code is based on the ade4TkGUI package by Jean Thioulouse <jthioulouse@biomserv.univ-lyon1.fr>, Stephane
#       Dray <dray@biomserv.univ-lyon1.fr>
#
dialog.MS.demon <-
function()
{
	tf <- tktoplevel()
	tkwm.title(tf,"Demonstration ")
	
	done <- tclVar(0)
	sepVar <- tclVar(1)
	

	
	
	"choosefic" <- function()
	{	
		
		sep <- tclvalue(sepVar)
		if (sep == 1) sepch <- "choose"
		if (sep == 2) sepch <- "MS.DataCreation.Agilent"
		if (sep == 3) sepch <- "MS.DataCreation.ASCII"
		if (sep == 4) sepch <- "MS.test.clust"
		if (sep == 5) sepch <- "MS.clust"
		if (sep == 6) sepch <- "MS.aristo"
		if (sep == 7) sepch <- "MS.msp"
		if (sep == 8) sepch <- "MS.nist"
		if (sep == 9) sepch <- "MS.DataCreation.CDF"
		
			if (sep==1){
				eval(choosepackage())
			}
			else{
			rdcom <- paste("dialog.", sepch, "()",sep="")
		 #data(list=c("Data_testclust","Agilent_quantF_MSclust","Agilent_quantT_MSclust","ASCII_MSclust"), package="MSeasy" ) #load demo data sets from MSeasy
			eval(parse(text=rdcom), envir=as.environment("e1"))
			}
		tkdestroy(tf)
		
	}
	
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	labh <- tklabel(frame1, bitmap="questhead")
	tkbind(labh, "<Button-1>", function() print(help("MSeasy")))
	tkgrid(tklabel(frame1,text="Demonstration", font="Times 18", foreground="red"), labh, columnspan=2)
	tkpack(frame1, fill = "x")

		
	sepFrame <- tkframe(tf, relief="groove", borderwidth=2)
    tkgrid(tklabel(sepFrame, text="Options:", foreground="blue"))
	tkgrid(tkradiobutton(sepFrame, text="Load an MSeasy demo Dataset                          ", value=1, variable=sepVar), sticky="w")
	tkgrid(tkradiobutton(sepFrame, text="MS.DataCreation (CDF+peaklist)                               ", value=9, variable=sepVar), sticky="w")
    tkgrid(tkradiobutton(sepFrame, text="MS.DataCreation (CDF+rteres.txt (Agilent))                   ", value=2, variable=sepVar), sticky="w")
    tkgrid(tkradiobutton(sepFrame, text="MS.DataCreation (ASCII)         ", value=3, variable=sepVar), sticky="w")
    tkgrid(tkradiobutton(sepFrame, text="MS.test.clust                   ", value=4, variable=sepVar), sticky="w")
	tkgrid(tkradiobutton(sepFrame, text="MS.clust                        ", value=5, variable=sepVar), sticky="w")
	tkgrid(tkradiobutton(sepFrame, text="MSeasyToARISTO                        ", value=6, variable=sepVar), sticky="w")
	tkgrid(tkradiobutton(sepFrame, text="MSeasyToMSP                        ", value=7, variable=sepVar), sticky="w")
	tkgrid(tkradiobutton(sepFrame, text="SearchNIST (Only for Windows)                       ", value=8, variable=sepVar), sticky="w")
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

