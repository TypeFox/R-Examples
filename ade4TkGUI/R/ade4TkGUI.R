################################
# GUI for ade4 functions 
################################
env_ade4tkgui <- .GlobalEnv

"ade4TkGUI" <- function(show = FALSE, history = FALSE)
{
  
	assign("cmdlist", "cmdlist", envir=env_ade4tkgui)
	assign("winlist", 1, envir=env_ade4tkgui)
	
	if (exists("ade4TkGUIFlag")) rm("ade4TkGUIFlag", envir=env_ade4tkgui)
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"ade4TkGUI")
#
# Checkboxes
#
#	if (show) showvar <- tclVar(1) else showvar <- tclVar(0)
#	if (history) histvar <- tclVar(1) else histvar <- tclVar(0)
#
# Menu setup
#
	frame0 <- tkframe(tt, relief="groove", borderwidth=2, background="white")

	topMenuFile <- tkmenubutton(frame0, text="File", background="white")
	fileMenu <- tkmenu(topMenuFile, tearoff=TRUE)
	tkconfigure(topMenuFile, menu=fileMenu)
	tkadd(fileMenu,"command",label="Read text file",command=function() readtable(show, history))
	tkadd(fileMenu,"command",label="Load data set",command=function() choosepackage(show, history))
	tkadd(fileMenu,"command",label="Edit dataframe",command=function() editdf())
	tkadd(fileMenu,"separator")
	tkadd(fileMenu,"command",label="Quit R",command=function() askQuit())

	topMenuWin <- tkmenubutton(frame0, text="Windows", background="white")
	WinMenu <- tkmenu(topMenuWin, tearoff=TRUE)
	tkconfigure(topMenuWin, menu=WinMenu)
	tkadd(WinMenu,"command",label="New graphic window",command=function() newGr())
	tkadd(WinMenu,"command",label="Change graphic window",command=function() selectGr())
	tkadd(WinMenu,"command",label="Save graphic window",command=function() outgraph())

	topMenu1Tab <- tkmenubutton(frame0, text="1table", background="white")
	oneTabMenu <- tkmenu(topMenu1Tab, tearoff=TRUE)
	tkconfigure(topMenu1Tab, menu=oneTabMenu)
	tkadd(oneTabMenu,"command",label="Principal components analysis",command=function() dialog.dudi.pca(show, history))
	tkadd(oneTabMenu,"command",label="Correspondence analysis",command=function() dialog.dudi.coa(show, history))
	tkadd(oneTabMenu,"command",label="Multiple correspondence analysis",command=function() dialog.dudi.acm(show, history))
	tkadd(oneTabMenu,"command",label="Principal coordinates analysis",command=function() dialog.dudi.pco(show, history))
	tkadd(oneTabMenu,"command",label="Fuzzy correspondence analysis",command=function() dialog.dudi.fca(show, history))
	tkadd(oneTabMenu,"command",label="Fuzzy principal components analysis",command=function() dialog.dudi.fpca(show, history))
	tkadd(oneTabMenu,"command",label="Mixed (quant + qual) variables analysis",command=function() dialog.dudi.mix(show, history))
	tkadd(oneTabMenu,"command",label="Non symmetric correspondence analysis",command=function() dialog.dudi.nsc(show, history))
	tkadd(oneTabMenu,"command",label="Decentered correspondence analysis",command=function() dialog.dudi.dec(show, history))

	topMenu1TabG <- tkmenubutton(frame0, text="1table+groups", background="white")
	oneTabMenuG <- tkmenu(topMenu1TabG, tearoff=TRUE)
	tkconfigure(topMenu1TabG, menu=oneTabMenuG)
	tkadd(oneTabMenuG,"command",label="Between groups analysis",command=function() dialog.between(show, history))
	tkadd(oneTabMenuG,"command",label="Within groups analysis",command=function() dialog.within(show, history))
	tkadd(oneTabMenuG,"command",label="Discriminant analysis",command=function() dialog.discrimin(show, history))

	topMenu2Tab <- tkmenubutton(frame0, text="2tables", background="white")
	twoTabMenu <- tkmenu(topMenu2Tab, tearoff=TRUE)
	tkconfigure(topMenu2Tab, menu=twoTabMenu)
	tkadd(twoTabMenu,"command",label="Coinertia analysis",command=function() dialog.coinertia(show, history))
	tkadd(twoTabMenu,"command",label="Canonical correspondence analysis",command=function() dialog.cca(show, history))
	tkadd(twoTabMenu,"command",label="PCA w/r to instrumental variables",command=function() dialog.pcaiv(show, history))
	tkadd(twoTabMenu,"command",label="PCA w/r to orthogonal I.V.",command=function() dialog.pcaivortho(show, history))
	tkadd(twoTabMenu,"command",label="Double principal coordinate analysis",command=function() dialog.dpcoa(show, history))

	topMenuGraph <- tkmenubutton(frame0, text="Graphics", background="white")
	graphMenu <- tkmenu(topMenuGraph, tearoff=TRUE)
	tkconfigure(topMenuGraph, menu=graphMenu)
	tkadd(graphMenu,"command",label="Labels",command=function() dialog.s.label(show, history))
	tkadd(graphMenu,"command",label="Arrows",command=function() dialog.s.arrow(show, history))
	tkadd(graphMenu,"command",label="Classes",command=function() dialog.s.class(show, history))
	tkadd(graphMenu,"command",label="Values",command=function() dialog.s.value(show, history))
	tkadd(graphMenu,"command",label="Correlation circle",command=function() dialog.s.corcircle(show, history))
	tkadd(graphMenu,"command",label="Convex hulls",command=function() dialog.s.chull(show, history))
	tkadd(graphMenu,"command",label="Match",command=function() dialog.s.match(show, history))
	tkadd(graphMenu,"separator")
	tkadd(graphMenu,"command",label="Display dudi",command=function() dudisp(show, history))
	tkadd(graphMenu,"command",label="Monte-Carlo test",command=function() dialog.MCTests(show, history))
	# tkadd(graphMenu,"command",label="Explore graph",command=function() exploregraph(show, history))
	tkadd(graphMenu,"command",label="ordiClust",command=function() ordiClust())
	tkadd(graphMenu,"command",label="Reset graph list",command=function() resetgraph())

	tkpack(topMenuFile, topMenuWin, topMenu1Tab, topMenu1TabG, topMenu2Tab, topMenuGraph, side="left")
	tkpack(frame0, expand="TRUE", fill="x")
#
# title and icons
#
	frame1 <- tkframe(tt, relief="groove", borderwidth=2, background="white")

	for (i in 1:length(.libPaths())) {
		icnfnameR <- file.path(paste(.libPaths()[i],"/ade4TkGUI/Rlogo.gif",sep=""))
		if (file.exists(icnfnameR)) {
			icn <- tkimage.create("photo", file=icnfnameR)
			Rlabel <- tklabel(frame1, image=icn, background="white")
			oldopt <- options("htmlhelp")$htmlhelp
			tkbind(Rlabel, "<Button-1>", function() help.start())
			options("htmlhelp" = oldopt)
		}
		icnfnameTk <- file.path(paste(.libPaths()[i],"/ade4TkGUI/tcltk.gif",sep=""))
		if (file.exists(icnfnameTk)) {
			icn <- tkimage.create("photo", file=icnfnameTk)
			TclTklabel <- tklabel(frame1, image=icn, background="white")
			tkbind(TclTklabel, "<Button-1>", function() print(help("tcltk")))
		}
	}

	labh <- tklabel(frame1, bitmap="questhead", background="white")
	tkbind(labh, "<Button-1>", function() print(help("ade4TkGUI")))
	if (show) if (history) titre <- tklabel(frame1,text="ade4TkGUI (T, T)", font="Times 20", foreground="red", background="white")
		else titre <- tklabel(frame1,text="ade4TkGUI (T, F)", font="Times 20", foreground="red", background="white")
		else if (history) titre <- tklabel(frame1,text="ade4TkGUI (F, T)", font="Times 20", foreground="red", background="white")
		else titre <- tklabel(frame1,text="ade4TkGUI (F, F)", font="Times 20", foreground="red", background="white")
	
	helplab <- tklabel(frame1,text="- Right click buttons for help - Double click in lists to select -", font="Times 12", foreground="dark green", background="white")

#	frameCheck <- tkframe(frame1, relief="flat", borderwidth=0, background="white")
#	if (show) show.cbut <- tklabel(frameCheck,text="T", background="white", font="system 10") else
#		show.cbut <- tklabel(frameCheck,text="F", background="white", font="system 10")
#	if (history) hist.cbut <- tklabel(frameCheck,text="T", background="white", font="system 10") else
#		hist.cbut <- tklabel(frameCheck,text="F", background="white", font="system 10")
#	tkgrid(show.cbut, sticky = "w")
#	tkgrid(hist.cbut, sticky = "w")
#
#	tkgrid(frameCheck, Rlabel, titre, labh, TclTklabel, padx=10, sticky = "w")
	tkgrid(Rlabel, titre, labh, TclTklabel, padx=10, sticky = "w")
	tkgrid(helplab, columnspan=4)
	tkpack(frame1, expand="TRUE", fill="x")
#
# read text files and load ade4 data set
#
	frame1b <- tkframe(tt, relief="groove", borderwidth=2, background="white")
	tkpack(tklabel(frame1b,text="- Data sets -", font="Times 14", foreground="blue", background="white"))	
	readtable.but <- tkbutton(frame1b, text="Read a data file", command=function() readtable(show, history))
	getdata.but <- tkbutton(frame1b, text="Load a data set", command=function() choosepackage(show, history))
	tkpack(readtable.but, getdata.but, ipadx=20, side="left", expand="TRUE", fill="x")
	tkpack(frame1b, expand="TRUE", fill="x")

	tkbind(readtable.but, "<Button-3>", function() print(help("read.table")))
	tkbind(getdata.but, "<Button-3>", function() print(help("data")))

#
# One table analyses
#
	frame2 <- tkframe(tt, relief="groove", borderwidth=2, background="white")
	tkpack(tklabel(frame2,text="- One table analyses -", font="Times 14", foreground="blue", background="white"))
	dudi.pca.but <- tkbutton(frame2, text="PCA", command=function() dialog.dudi.pca(show, history))
	dudi.coa.but <- tkbutton(frame2, text="COA", command=function() dialog.dudi.coa(show, history))
	dudi.acm.but <- tkbutton(frame2, text="MCA", command=function() dialog.dudi.acm(show, history))
	dudi.pco.but <- tkbutton(frame2, text="PCO", command=function() dialog.dudi.pco(show, history))
	tkpack(dudi.pca.but, dudi.coa.but, dudi.acm.but, dudi.pco.but, ipadx=20, side="left", expand="TRUE", fill="x")
	tkpack(frame2, expand="TRUE", fill="x")

	tkbind(dudi.pca.but, "<Button-3>", function() print(help("dudi.pca")))
	tkbind(dudi.coa.but, "<Button-3>", function() print(help("dudi.coa")))
	tkbind(dudi.acm.but, "<Button-3>", function() print(help("dudi.acm")))
	tkbind(dudi.pco.but, "<Button-3>", function() print(help("dudi.pco")))
#
# One table with groups
#
	frame2b <- tkframe(tt, relief="groove", borderwidth=2, background="white")
	tkpack(tklabel(frame2b,text="- One table with groups -", font="Times 14", foreground="blue", background="white"))
	between.but <- tkbutton(frame2b, text="BGA", command=function() dialog.between(show, history))
	within.but <- tkbutton(frame2b, text="WGA", command=function() dialog.within(show, history))
	discrimin.but <- tkbutton(frame2b, text="DA", command=function() dialog.discrimin(show, history))
	tkpack(between.but, within.but, discrimin.but, ipadx=20, side="left", expand="TRUE", fill="x")
	tkpack(frame2b, expand="TRUE", fill="x")

	tkbind(between.but, "<Button-3>", function() print(help("bca")))
	tkbind(within.but, "<Button-3>", function() print(help("within")))
	tkbind(discrimin.but, "<Button-3>", function() print(help("discrimin")))
#
# Two tables analyses
#
	frame2c <- tkframe(tt, relief="groove", borderwidth=2, background="white")
	tkpack(tklabel(frame2c,text="- Two tables analyses -", font="Times 14", foreground="blue", background="white"))
	coinertia.but <- tkbutton(frame2c, text="Coinertia", command=function() dialog.coinertia(show, history))
	cca.but <- tkbutton(frame2c, text="CCA", command=function() dialog.cca(show, history))
	pcaiv.but <- tkbutton(frame2c, text="PCAIV", command=function() dialog.pcaiv(show, history))
	tkpack(coinertia.but, cca.but, pcaiv.but, ipadx=20, side="left", expand="TRUE", fill="x")
	tkpack(frame2c, expand="TRUE", fill="x")

	tkbind(coinertia.but, "<Button-3>", function() print(help("coinertia")))
	tkbind(cca.but, "<Button-3>", function() print(help("cca")))
	tkbind(pcaiv.but, "<Button-3>", function() print(help("pcaiv")))
#
# Graphics
#
	frame3 <- tkframe(tt, relief="groove", borderwidth=2, background="white")
	tkpack(tklabel(frame3,text="- Graphic functions -", font="Times 14", foreground="blue", background="white"))
	s.label.but <- tkbutton(frame3, text="Labels", command=function() dialog.s.label(show, history))
	s.class.but <- tkbutton(frame3, text="Classes", command=function() dialog.s.class(show, history))
	s.value.but <- tkbutton(frame3, text="Values", command=function() dialog.s.value(show, history))
	tkpack(s.label.but, s.class.but, s.value.but, ipadx=20, side="left", expand="TRUE", fill="x")
	tkpack(frame3, expand="TRUE", fill="x")

	tkbind(s.label.but, "<Button-3>", function() print(help("s.label", package = "adegraphics")))
	tkbind(s.value.but, "<Button-3>", function() print(help("s.value", package = "adegraphics")))
	tkbind(s.class.but, "<Button-3>", function() print(help("s.class", package = "adegraphics")))
#
# Advanced graphics
#
	frame4 <- tkframe(tt, relief="groove", borderwidth=2, background="white")
	tkpack(tklabel(frame4,text="- Advanced graphics -", font="Times 14", foreground="blue", background="white"))
	dudisp.but <- tkbutton(frame4, text="dudi display", command=function() dudisp(show, history))
	#explore.but <- tkbutton(frame4, text="Explore", state="disabled", command=function() exploregraph(show, history)) ## bouton pour explore
  MC.but <- tkbutton(frame4, text="MCTests", command=function() dialog.MCTests(show, history))
	OC.but <- tkbutton(frame4, text="ordiClust", command=function() ordiClust())
	#tkpack(dudisp.but, MC.but, explore.but, OC.but, ipadx=5, side="left", expand="TRUE", fill="x")
  tkpack(dudisp.but, MC.but, OC.but, ipadx=5, side="left", expand="TRUE", fill="x")
	tkpack(frame4, expand="TRUE", fill="x")

	tkbind(dudisp.but, "<Button-3>", function() print(help("dudi")))
	tkbind(MC.but, "<Button-3>", function() print(help("randtest")))
	#tkbind(explore.but, "<Button-3>", function() print(help("explore")))
	tkbind(OC.but, "<Button-3>", function() print(help("ordiClust")))
#
# Quit
#
	frame5 <- tkframe(tt, relief="groove", borderwidth=2, background="white")
	cancel.but <- tkbutton(frame5, text="Dismiss", command=function() tkdestroy(tt), foreground="darkgreen", background="white", font="Times 14")
	quity.but <- tkbutton(frame5, text="Quit R (save)", command=function() q("yes"), foreground="darkgreen", background="white", font="Times 14")
	quitn.but <- tkbutton(frame5, text="Quit R (don't save)", command=function() q("no"), foreground="darkgreen", background="white", font="Times 14")
	tkpack(quity.but, cancel.but, quitn.but, side="left", expand="TRUE", fill="x")
	tkpack(frame5, expand="TRUE", fill="x")
	tkfocus(tt)
	return(invisible())
}
