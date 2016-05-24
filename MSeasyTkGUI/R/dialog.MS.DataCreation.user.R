# MSeasy Tcl/Tk GUI for some basic functions in the MSeasy package.
# Part of the code is based on the ade4TkGUI package by Jean Thioulouse <jthioulouse@biomserv.univ-lyon1.fr>, Stephane
#       Dray <dray@biomserv.univ-lyon1.fr>
#
dialog.MS.DataCreation.user <-
function()
{
	op=options()
	options(warn=-1)
	
	tt <- tktoplevel()
	tkwm.title(tt,"Load a user made file")
	
	done <- tclVar(0)
	typesvar <- tclVar()
	tabnamevar <- tclVar("")
	filenamevar <- tclVar("clipboard")
	fictrt <- tclVar()
	varnames <- tclVar(1)
	sepVar <- tclVar(1)
	otherSepVar <- tclVar("")
	decSepVar <- tclVar(1)
	colClassVar <- tclVar(4)

	"choosefile" <- function()
	{
		tclvalue(typesvar)="{{Text files} {.txt}}"
		fictrt <- tkgetOpenFile(filetypes=tclvalue(typesvar))
		fpath <- tclvalue(fictrt)
		tkdelete(file.entry, 0, "end")
		tkinsert(file.entry, "end", fpath)
	}
	
	"choosefic" <- function()
	{	
		if (tclvalue(tabnamevar) != "") {
			tabname  <- parse(text=tclvalue(tabnamevar))[[1]]
		} else tabname <- "initial_DATA"

		if (tclvalue(filenamevar) != "") {
			filename  <- tclvalue(filenamevar)
		} else return()
		
		varn <- as.logical(tclObj(varnames))
		sep <- tclvalue(sepVar)
		if (sep == 1) sepch <- ""
		if (sep == 2) sepch <- ","
		if (sep == 3) sepch <- ";"
		if (sep == 4) {
			if (tclvalue(otherSepVar) != "") {
				otherSep <- tclvalue(otherSepVar)
			} else otherSep <- ""
			sepch <- otherSep
		}
		decSep <- tclvalue(decSepVar)
		if (decSep == 1) decsepch <- "."
		if (decSep == 2) decsepch <- ","

		colClass <- tclvalue(colClassVar)
		if (colClass == 1) colClass <- "numeric"
		if (colClass == 2) colClass <- "character"
		if (colClass == 3) colClass <- "factor"
		if (colClass == 4) colClass <- "auto"

		if (colClass == "auto") {
			rdcom <- paste(tabname," <<- read.table(file='", filename, "', header=",varn,
				", sep='",sepch,"', dec='",decsepch,"')", sep="", check.names=FALSE)
		} else {
			rdcom <- paste(tabname," <<- read.table(file='", filename, "', header=",varn,
				", sep='",sepch,"', dec='",decsepch,"', colClasses='",colClass,"')", sep="", check.names=FALSE)
		}
	
		

		eval(parse(text=rdcom))
		tkdestroy(tt)
		#displaytable(tabname)
		eval(parse(text=paste("edit(", tabname, ")", sep="")))
	}
	
	frame1 <- tkframe(tt, relief="groove", borderwidth=2)
	labh <- tklabel(frame1, bitmap="questhead")
	tkbind(labh, "<Button-1>", function() print(help("read.table")))
	tkgrid(tklabel(frame1,text="Load a user made file", font="Times 18", foreground="red"), labh, columnspan=2)
	tkpack(frame1, fill = "x")

	frame2 <- tkframe(tt, relief="groove", borderwidth=2)	
	tab.entry <- tkentry(frame2, textvariable=tabnamevar)
	file.entry <- tkentry(frame2, textvariable=filenamevar)
	choosefile.but <- tkbutton(frame2, text="Set", command=function() choosefile())
	tkgrid(tklabel(frame2,text="Info: files  should be formatted with: col1=sample name, col2=RT/RI, following columns=mz intensity"),labh, columnspan=2)
	tkgrid(tklabel(frame2,text="File to read (txt, csv..): "), file.entry, choosefile.but)
	tkgrid(tklabel(frame2,text="Name of the project : "), tab.entry)
	varnames.cbut <- tkcheckbutton(frame2,text="Columns headings in the first row of data file", variable=varnames)
	tkgrid(varnames.cbut, columnspan=2, sticky="w")
	
	sepFrame <- tkframe(frame2, relief="groove", borderwidth=2)
	sep.entry <- tkentry(sepFrame, textvariable=otherSepVar, width=10)
    tkgrid(tklabel(sepFrame, text="Field separator:", foreground="blue"))
    tkgrid(tkradiobutton(sepFrame, text="Default", value=1, variable=sepVar), sticky="w")
    tkgrid(tkradiobutton(sepFrame, text="Commas", value=2, variable=sepVar), sticky="w")
    tkgrid(tkradiobutton(sepFrame, text="Semicolon", value=3, variable=sepVar), sticky="w")
    tkgrid(tkradiobutton(sepFrame, text="Other", value=4, variable=sepVar), sep.entry, sticky="w")

    decSepFrame <- tkframe(frame2, relief="groove", borderwidth=2)
    tkgrid(tklabel(decSepFrame, text="Decimal separator:", foreground="blue"))
    tkgrid(tkradiobutton(decSepFrame, text="Period [.]", value=1, variable=decSepVar), sticky="w")
    tkgrid(tkradiobutton(decSepFrame, text="Comma [,]", value=2, variable=decSepVar), sticky="w")
 	#tkgrid(sepFrame, decSepFrame, sticky="n")

    colClassFrame <- tkframe(frame2, relief="groove", borderwidth=2)
    tkgrid(tklabel(colClassFrame, text="Column types:", foreground="blue"))
    tkgrid(tkradiobutton(colClassFrame, text="Default", value=4, variable=colClassVar), sticky="w")
    tkgrid(tkradiobutton(colClassFrame, text="Numeric", value=1, variable=colClassVar), sticky="w")
    tkgrid(tkradiobutton(colClassFrame, text="Character", value=2, variable=colClassVar), sticky="w")
    tkgrid(tkradiobutton(colClassFrame, text="Factor", value=3, variable=colClassVar), sticky="w")
 	tkgrid(sepFrame, decSepFrame, colClassFrame, sticky="n")

	tkpack(frame2, fill = "x")

	ok.but <- tkbutton(tt, text="Submit", command=function() choosefic())
	cancel.but <- tkbutton(tt, text="Cancel", command=function() tkdestroy(tt))
	tkpack(cancel.but, ok.but, side="left", fill="x", expand=1)

	tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
	tkbind(tt, "<KeyPress-Return>", function() choosefic())
	tkbind(tt, "<KeyPress-Escape>", function() tkdestroy(tt))
	tkwait.variable(done)
	
	#
# User closed the window
#
	if(tclvalue(done)=="2") return()
	
	if(tclvalue(done) == "2") return(0)
	tkdestroy(tt)
}

