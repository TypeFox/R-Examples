################################
# GUI for s.corcircle function
################################
"dialog.s.corcircle" <- function(show, history)
{
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"s.corcircle")
  
  frame1 <- tkframe(tt, relief="groove", borderwidth=2)
  frame2 <- tkframe(tt, relief="groove", borderwidth=2)
  frame3 <- tkframe(tt, relief="groove", borderwidth=2)
  xyframe <- tkframe(frame1, relief="groove", borderwidth=2)
  labframe <- tkframe(frame1, relief="groove", borderwidth=2)
  posframe <- tkframe(frame2, relief="groove", borderwidth=2)
  legframe <- tkframe(frame2, relief="groove", borderwidth=2)
  optframe <- tkframe(frame3, relief="groove", borderwidth=2)
  origframe <- tkframe(frame3, relief="groove", borderwidth=2)
  gridframe <- tkframe(frame3, relief="groove", borderwidth=2)
#
# Variables for text fields
#
	xyvar <- tclVar()
	nxvar <- tclVar(1)
	nyvar <- tclVar(2)
	labvar <- tclVar()
	clabvar <- tclVar(1)
	cgrvar <- tclVar(1)
	subvar <- tclVar()
	csubvar <- tclVar(1)
#
# Checkboxes variables
#
	gridvar <- tclVar(1)
	posvar <- tclVar(3)
	addvar <- tclVar(0)
	fullcircvar <- tclVar(1)
	boxvar <- tclVar(1)
  drawBoxvar <- tclVar(1)
#
# Title
#
	TFrame <- tkframe(tt, relief="groove")
	labh <- tklabel(TFrame, bitmap="questhead")
	tkgrid(tklabel(TFrame,text="Correlation circle", font="Times 18", foreground="red"), labh)
	tkbind(labh, "<Button-1>", function() print(help("s.corcircle", package = "adegraphics")))
	tkpack(TFrame)
#
# Coordinates frame
#
	xy.entry <- tkentry(xyframe, textvariable=xyvar, width=10)
	nx.entry <- tkentry(xyframe, textvariable=nxvar, width=3)
	ny.entry <- tkentry(xyframe, textvariable=nyvar, width=3)
	dfnr.label <- tklabel(xyframe, width=4)
	dfnc.label <- tklabel(xyframe, width=4)
	choosexy.but <- tkbutton(xyframe, text="Set", command=function() choosedf(xy.entry, dfnr.label, dfnc.label))
	tkgrid(tklabel(xyframe, text="- Coordinates -", foreground="blue"), columnspan=5)
	tkgrid(tklabel(xyframe,text="XY coordinates"), xy.entry, choosexy.but, dfnr.label, dfnc.label)
	tkgrid(tklabel(xyframe,text="X axis col. #"), nx.entry)
	tkgrid(tklabel(xyframe,text="Y axis col. #"), ny.entry)
#
# Labels frame
#
	lab.entry <- tkentry(labframe, textvariable=labvar, width=10)
	clab.entry <- tkentry(labframe, textvariable=clabvar, width=10)
	chooselab.but <- tkbutton(labframe, text="Set", command=function() chooselab(tt, dfnr.label, lab.entry))
	drawBox.cbut <- tkcheckbutton(labframe,text="Boxes", variable=drawBoxvar)
  
  tkgrid(tklabel(labframe, text="- Labels -", foreground="blue"), columnspan=3)
	tkgrid(tklabel(labframe,text="Labels"), lab.entry, chooselab.but)
	tkgrid(tklabel(labframe,text="Label size"), clab.entry, drawBox.cbut)

	tkpack(xyframe, labframe, side="left")
	tkpack(frame1)
#
# Legend frame
#
  tkpack(tklabel(posframe, text="- Sub-title position -", foreground="blue"), anchor="w")
  tkpack(tkradiobutton(posframe, text="Top left", value=1, variable=posvar), anchor="w")
  tkpack(tkradiobutton(posframe, text="Top right", value=2, variable=posvar), anchor="w")
  tkpack(tkradiobutton(posframe, text="Bottom left", value=3, variable=posvar), anchor="w")
  tkpack(tkradiobutton(posframe, text="Bottom right", value=4, variable=posvar), anchor="w")
  
	sub.entry <- tkentry(legframe, textvariable=subvar)
	csub.entry <- tkentry(legframe, textvariable=csubvar, width=10)
	tkgrid(tklabel(legframe, text="- Sub-title -", foreground="blue"), columnspan=2)
	tkgrid(tklabel(legframe,text="Sub-title string"), sub.entry)
	tkgrid(tklabel(legframe,text="Sub-title size"), csub.entry)
  
  tkpack(legframe, posframe, side="left", expand=1)
  tkpack(frame2, fill="x")
#
# Options frame
#
	tkgrid(tklabel(optframe, text="- Draw options -", foreground="blue"))
	add.cbut <- tkcheckbutton(optframe,text="Add to plot", variable=addvar)
	tkgrid(add.cbut)
	
	fullcirc.cbut <- tkcheckbutton(optframe,text="Full circle", variable=fullcircvar)
	tkgrid(fullcirc.cbut)
	
	box.cbut <- tkcheckbutton(optframe,text="Draw box", variable=boxvar)
	tkgrid(box.cbut)
	
	tkgrid(tklabel(gridframe, text="- Grid -", foreground="blue"), columnspan=2)
	grid.cbut <- tkcheckbutton(gridframe,text="Draw grid", variable=gridvar)
	cgr.entry <- tkentry(gridframe, textvariable=cgrvar, width=10)
	tkgrid(grid.cbut)
	tkgrid(tklabel(gridframe,text="Grid legend size"), cgr.entry)
	
	tkpack(optframe, gridframe, origframe, side="left", expand=1)

	tkpack(frame3, fill="x")
#
# Local variables
#
	vnr=NULL			# Vector of dataframes row numbers
	vnc=NULL			# Vector of dataframes column numbers
	numi=1				# Number of choosed element
	done <- tclVar(0)	# To terminate the dialog
	
################################
# Function to build the command line from dialog widgets
################################
"build" <- function()	{
  l <- list(dfxy = .test1value(tclvalue(xyvar), ""),
            xax = .test1value(tclvalue(nxvar), ""),
            yax = .test1value(tclvalue(nyvar), ""),
            labels = .test1value(tclvalue(labvar), ""),
            plabels.cex = .test1value(tclvalue(clabvar), ""),
            plabels.boxes.draw = as.logical(tclObj(drawBoxvar)),
            psub.text = tclvalue(subvar),
            psub.cex = .test1value(tclvalue(csubvar), ""),
            add = as.logical(tclObj(addvar)),
            pgrid.draw = as.logical(tclObj(gridvar)),
            pgrid.text.cex = .test1value(tclvalue(cgrvar), ""),
            fullcircle = as.logical(tclObj(fullcircvar)),
		        pbackground.box = as.logical(tclObj(boxvar)),
            plot = FALSE
          	)
  
	if (tclvalue(posvar) == 1) l$psub.position <- "topleft"
	if (tclvalue(posvar) == 2) l$psub.position <- "topright"
	if (tclvalue(posvar) == 3) l$psub.position <- "bottomleft"
	if (tclvalue(posvar) == 4) l$psub.position <- "bottomright"
  
  l <- l[which(l != "")]
  return(do.call("s.corcircle", l))
	}

################################
# Function to reset all dialog elements to default values
################################
	"reset" <- function()
	{
		tclvalue(xyvar) <- ""
		tclvalue(nxvar) <- "1"
		tclvalue(nyvar) <- "2"
		tclvalue(labvar) <- ""
		tclvalue(clabvar) <- "1"
		tclvalue(cgrvar) <- "1"
		tclvalue(subvar) <- ""
		tclvalue(csubvar) <- "1"
		tkconfigure(dfnr.label, text="")
		tkconfigure(dfnc.label, text="")
		tclvalue(gridvar) <- "1"
    tclvalue(posvar) <- "3"
		tclvalue(addvar) <- "0"
    tclvalue(fullcircvar) <- "1"
	  tclvalue(boxvar) <- "1"
    tclvalue(drawBoxvar) <- "1"
	}
	
################################
# Function to draw the graphic
################################
	"drawgraph" <- function()
	{
		#
		# Build and display the command line so that the user can check it
		#
		cmd <- print(build())
		if (show) {
			#
			# Echoe the command line to the console
			#
			pr1 <- substr(options("prompt")$prompt, 1, 2)
			cat(deparse(cmd@Call, width.cutoff = 500), "\n", pr1, sep="")
		}
		#
		# Execute the command
		#
		eval.parent(cmd)
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), cmd@Call), envir=env_ade4tkgui)
		if (history) rewriteHistory(deparse(cmd@Call, width.cutoff = 500))
	}
#
# Place the three buttons
#
	RCSFrame <- tkframe(tt, relief="groove")
	reset.but <- tkbutton(RCSFrame, text="Reset", command=reset)
	cancel.but <- tkbutton(RCSFrame, text="Dismiss", command=function() tkdestroy(tt))
	submit.but <- tkbutton(RCSFrame, text="Submit", default="active", command=function() drawgraph())
	tkgrid(cancel.but, submit.but, reset.but, ipadx=20)	
	tkpack(RCSFrame)
#
# If window is closed by user, terminate the dialog
#
	tkbind(tt, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tt, "<KeyPress-Return>", function() drawgraph())
	tkbind(tt, "<KeyPress-Escape>", function() tkdestroy(tt))
#
# User closed the window
#
	if(tclvalue(done)=="2") return()
}
