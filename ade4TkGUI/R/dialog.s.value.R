################################
# GUI for s.value function
################################
"dialog.s.value" <- function(show, history)
{
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"s.value")
  
	frame1 <- tkframe(tt, relief="groove", borderwidth=2)
	frame2 <- tkframe(tt, relief="groove", borderwidth=2)
	frame3 <- tkframe(tt, relief="groove", borderwidth=2)
	frame4 <- tkframe(tt, relief="groove", borderwidth=2)
	xyframe <- tkframe(frame1, relief="groove", borderwidth=2)
	symframe <- tkframe(frame1, relief="groove", borderwidth=2)
	methframe <- tkframe(frame1, relief="groove", borderwidth=2)
	limframe <- tkframe(frame2, relief="groove", borderwidth=2)
	posframe <- tkframe(frame2, relief="groove", borderwidth=2)
	legframe <- tkframe(frame2, relief="groove", borderwidth=2)
	optframe <- tkframe(frame3, relief="groove", borderwidth=2)
	origframe <- tkframe(frame3, relief="groove", borderwidth=2)
	gridframe <- tkframe(frame3, relief="groove", borderwidth=2)
	miscframe <- tkframe(frame4, relief="groove", borderwidth=2)
	valframe <- tkframe(frame4, relief="groove", borderwidth=2)
#
# Variables for text fields
#
	xyvar <- tclVar()
	nxvar <- tclVar(1)
	nyvar <- tclVar(2)
	valvar <- tclVar()
	pchvar <- tclVar(1)
	cpvar <- tclVar(1)
	clegendvar <- tclVar(1)
	xl1var <- tclVar()
	xl2var <- tclVar()
	yl1var <- tclVar()
	yl2var <- tclVar()
	cgrvar <- tclVar(1)
	orxvar <- tclVar(0)
	oryvar <- tclVar(0)
	subvar <- tclVar()
	csubvar <- tclVar(1)
	neigvar <- tclVar()
	cneigvar <- tclVar(1)
	#pmvar <- tclVar()
	spvar <- tclVar()
#
# Checkboxes variables
#
	methvar <- tclVar(1)
	gridvar <- tclVar(1)
	axesvar <- tclVar(0)
	origvar <- tclVar(1)
	posvar <- tclVar(1)
	addvar <- tclVar(0)
#
# Title
#
	TFrame <- tkframe(tt, relief="groove")
	labh <- tklabel(TFrame, bitmap="questhead")
	tkgrid(tklabel(TFrame,text="Values", font="Times 18", foreground="red"), labh)
	tkbind(labh, "<Button-1>", function() print(help("s.value", package = "adegraphics")))
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
# Symbols frame
#
  tkpack(tklabel(symframe, text="- Symbols -", foreground="blue"))
  tkpack(tkradiobutton(symframe, text="Square", value=1, variable=pchvar), anchor="w")
  tkpack(tkradiobutton(symframe, text="Circle", value=2, variable=pchvar), anchor="w")
  tkpack(tkradiobutton(symframe, text="Diamond", value=3, variable=pchvar), anchor="w")
  tkpack(tkradiobutton(symframe, text="Up triangle", value=4, variable=pchvar), anchor="w")
  tkpack(tkradiobutton(symframe, text="Down triangle", value=5, variable=pchvar), anchor="w")
#
# Method frame
#
  tkpack(tklabel(methframe, text="- Method -", foreground="blue"))
  tkpack(tkradiobutton(methframe, text="B&W squares", value=1, variable=methvar), anchor="w")
  tkpack(tkradiobutton(methframe, text="Grey levels", value=2, variable=methvar), anchor="w")
	tkpack(xyframe, symframe, methframe, side="left")
	tkpack(frame1, fill="x")
#
# Values frame
#
	tkgrid(tklabel(valframe, text="- Values -", foreground="blue"), columnspan=3)
	chooseval.but <- tkbutton(valframe, text="Set", command=function() chooseval(tt, dfnr.label, val.entry))
	val.entry <- tkentry(valframe, textvariable=valvar, width=10)
  cp.entry <- tkentry(valframe, textvariable=cpvar, width=10)
  clegend.entry <- tkentry(valframe, textvariable=clegendvar, width=10)
	tkgrid(tklabel(valframe,text="Values"), val.entry, chooseval.but)
  tkgrid(tklabel(valframe,text="Char. size"), cp.entry)
  tkgrid(tklabel(valframe,text="Legend size"), clegend.entry)

#
# Misc frame
#
	neig.entry <- tkentry(miscframe, textvariable=neigvar, width=10)
	cneig.entry <- tkentry(miscframe, textvariable=cneigvar, width=3)
	#pm.entry <- tkentry(miscframe, textvariable=pmvar, width=10, state="disabled")
	sp.entry <- tkentry(miscframe, textvariable=spvar, width=10)

	chooseneig.but <- tkbutton(miscframe, text="Set", command=function() chooseneig(neig.entry))
	#choosepm.but <- tkbutton(miscframe, text="Set", state="disabled", command=function() choosepm(pm.entry))
	choosesp.but <- tkbutton(miscframe, text="Set", command=function() choosesp(sp.entry))

  tkgrid(tklabel(miscframe, text="- Misc. options -", foreground="blue"), columnspan=3)
	tkgrid(tklabel(miscframe,text="Neighbouring relation"), neig.entry, chooseneig.but)
	tkgrid(tklabel(miscframe,text="Neighbouring size"), cneig.entry)
	#tkgrid(tklabel(miscframe,text="Pixmap"), pm.entry, choosepm.but)
	tkgrid(tklabel(miscframe,text="Spatial object"), sp.entry, choosesp.but)

	tkpack(valframe, miscframe, side="left", expand=1)
	tkpack(frame4, fill="x")
#
# Limits frame
#
	xl1.entry <- tkentry(limframe, textvariable=xl1var, width=10)
	xl2.entry <- tkentry(limframe, textvariable=xl2var, width=10)
	yl1.entry <- tkentry(limframe, textvariable=yl1var, width=10)
	yl2.entry <- tkentry(limframe, textvariable=yl2var, width=10)
	tkgrid(tklabel(limframe, text="- Limits -", foreground="blue"), columnspan=2)
	tkgrid(tklabel(limframe,text="X min"), xl1.entry)
	tkgrid(tklabel(limframe,text="X max"), xl2.entry)
	tkgrid(tklabel(limframe,text="Y min"), yl1.entry)
	tkgrid(tklabel(limframe,text="Y max"), yl2.entry)
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
  
  tkpack(limframe, legframe, posframe, side="left", expand=1)
  tkpack(frame2, fill="x")
#
# Options frame
#
	axes.cbut <- tkcheckbutton(optframe,text="Draw axes", variable=axesvar)
	add.cbut <- tkcheckbutton(optframe,text="Add to plot", variable=addvar)
	tkgrid(tklabel(optframe, text="- Draw options -", foreground="blue"))
	tkgrid(axes.cbut)
	tkgrid(add.cbut)
	
	orig.cbut <- tkcheckbutton(origframe,text="Include origin", variable=origvar)
	orx.entry <- tkentry(origframe, textvariable=orxvar, width=10)
	ory.entry <- tkentry(origframe, textvariable=oryvar, width=10)
	tkgrid(tklabel(origframe, text="- Origin -", foreground="blue"), columnspan=2)
	tkgrid(orig.cbut)
	tkgrid(tklabel(origframe,text="X Origin"), orx.entry)
	tkgrid(tklabel(origframe,text="Y Origin"), ory.entry)

	grid.cbut <- tkcheckbutton(gridframe,text="Draw grid", variable=gridvar)
	cgr.entry <- tkentry(gridframe, textvariable=cgrvar, width=10)
	tkgrid(tklabel(gridframe, text="- Grid -", foreground="blue"), columnspan=2)
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
"build" <- function() {
  l <- list(dfxy = .test1value(tclvalue(xyvar), ""),
    				xax = .test1value(tclvalue(nxvar), ""),
            yax = .test1value(tclvalue(nyvar), ""),
            ppoints.cex = .test1value(tclvalue(cpvar), ""),
            xlim = .test2values(tclvalue(xl1var), tclvalue(xl2var), ""),
            ylim = .test2values(tclvalue(yl1var), tclvalue(yl2var), ""),
            psub.text = tclvalue(subvar),
            psub.cex = .test1value(tclvalue(csubvar), ""),
            paxes.draw = as.logical(tclObj(axesvar)),
            add = as.logical(tclObj(addvar)),
            pgrid.draw = as.logical(tclObj(gridvar)),
            pgrid.text.cex = .test1value(tclvalue(cgrvar), ""),
            porigin.include = as.logical(tclObj(origvar)),
            porigin.origin = .test2values(tclvalue(orxvar), tclvalue(oryvar), ""),
            nbobject = .test1value(tclvalue(neigvar), ""),
            pnb.edge.lwd = .test1value(tclvalue(cneigvar), ""),
            Sp = .test1value(tclvalue(spvar), ""),
            plegend.size = .test1value(tclvalue(clegendvar), ""),
            z = .test1value(tclvalue(valvar), rep(1, as.numeric(tkcget(dfnr.label, "-text")))),
            plot = FALSE
          	)
  
	if (tclvalue(posvar) == 1) l$psub.position <- "topleft"
	if (tclvalue(posvar) == 2) l$psub.position <- "topright"
	if (tclvalue(posvar) == 3) l$psub.position <- "bottomleft"
	if (tclvalue(posvar) == 4) l$psub.position <- "bottomright"

	if (tclvalue(methvar) == 1) l$method <- "size"
	if (tclvalue(methvar) == 2) l$method <- "color"

  if (tclvalue(pchvar) == 1) l$symbol <- "square"
  if (tclvalue(pchvar) == 2) l$symbol <- "circle"
  if (tclvalue(pchvar) == 3) l$symbol <- "diamond"
  if (tclvalue(pchvar) == 4) l$symbol <- "uptriangle"
  if (tclvalue(pchvar) == 5) l$symbol <- "downtriangle"

  l <- l[which(l != "")]
  return(do.call("s.value", l))
}

################################
# Function to reset all dialog elements to default values
################################
	"reset" <- function()
	{
		tclvalue(xyvar) <- ""
		tclvalue(nxvar) <- "1"
		tclvalue(nyvar) <- "2"
		tclvalue(valvar) <- ""
		tclvalue(csizevar) <- "1"
    tclvalue(pchvar) <- "1"
		tclvalue(cpvar) <- "1"
    tclvalue(clegendvar) <- "1"
		tclvalue(xl1var) <- ""
		tclvalue(xl2var) <- ""
		tclvalue(yl1var) <- ""
		tclvalue(yl2var) <- ""
		tclvalue(cgrvar) <- "1"
		tclvalue(orxvar) <- "0"
		tclvalue(oryvar) <- "0"
		tclvalue(subvar) <- ""
		tclvalue(csubvar) <- "1"
		tclvalue(neigvar) <- ""
		tclvalue(cneigvar) <- "1"
		tclvalue(pmvar) <- ""
		tclvalue(spvar) <- ""
		tkconfigure(dfnr.label, text="")
		tkconfigure(dfnc.label, text="")
		tclvalue(methvar) <- "1"
		tclvalue(gridvar) <- "1"
		tclvalue(axesvar) <- "0"
		tclvalue(origvar) <- "1"
		tclvalue(posvar) <- "1"
		tclvalue(addvar) <- "0"
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
		dcmd <- deparse(cmd@Call, width.cutoff = 500)
		tcmd <- paste(cmd@Call, sep="", collapse="; ")
		
		if (show) {
			#
			# Echoe the command line to the console
			#
			pr1 <- substr(options("prompt")$prompt, 1, 2)
			if (length(grep("expression", dcmd, fixed=TRUE)) == 0)
				cat(dcmd, "\n", pr1, sep="")
			else
				cat(tcmd, "\n", pr1, sep="")
		}
		#
		# Execute the command
		#
		eval.parent(cmd)
		if (length(grep("expression", dcmd, fixed=TRUE)) == 0) {
		  assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), cmd@Call), envir=env_ade4tkgui)
			if (history) rewriteHistory(deparse(cmd@Call, width.cutoff = 500))
		} else {
			if (history) rewriteHistory(tcmd)
		}
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
