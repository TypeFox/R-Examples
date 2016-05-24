# MSeasy Tcl/Tk GUI for some basic functions in the MSeasy package.
# Part of the code is based on the ade4TkGUI package by Jean Thioulouse <jthioulouse@biomserv.univ-lyon1.fr>, Stephane
#       Dray <dray@biomserv.univ-lyon1.fr>
#

MSeasyTkGUI <-
function()
{
	require(tcltk) || stop("tcltk support is absent")
	require(MSeasy) || stop("MSeasy support is absent")
	require(grDevices) || stop("grDevices support is absent")
	
	 attach(what=NULL, name="e1")
#	cmdlist <<- "cmdlist"
	assign("cmdlist", "cmdlist", envir=as.environment("e1"))
#	winlist <<- 1
	assign("winlist", 1, envir=as.environment("e1"))
	
	if (exists("MSeasyTkGUIFlag", envir=as.environment("e1"))) rm("MSeasyTkGUIFlag", envir=as.environment("e1"))
	if (exists("DemoFlag", envir=as.environment("e1"))) rm("DemoFlag", envir=as.environment("e1"))
#
# Main dialog window with title
#
	tt <- tktoplevel()
	tkwm.title(tt,"MSeasyTkGUI")
	
#
#Demo functions
#

demochoose<-function(){
assign("DemoFlag", 1, envir=as.environment("e1"))
choosepackage()
}

demodialog.MS.DataCreation.CDF<-function(){
assign("DemoFlag", 1, envir=as.environment("e1"))
dialog.MS.DataCreation.CDF()
}

demodialog.MS.DataCreation.Agilent<-function(){
assign("DemoFlag", 1, envir=as.environment("e1"))
dialog.MS.DataCreation.Agilent()
}

demodialog.MS.DataCreation.ASCII<-function(){
assign("DemoFlag", 1, envir=as.environment("e1"))
dialog.MS.DataCreation.ASCII()
}

demodialog.MS.test.clust<-function(){
assign("DemoFlag", 1, envir=as.environment("e1"))
dialog.MS.test.clust()
}

demodialog.MS.clust<-function(){
assign("DemoFlag", 1, envir=as.environment("e1"))
dialog.MS.clust()
}
demodialog.MS.aristo<-function(){
assign("DemoFlag", 1, envir=as.environment("e1"))
dialog.MS.aristo()
}
demodialog.MS.msp<-function(){
assign("DemoFlag", 1, envir=as.environment("e1"))
dialog.MS.msp()
}
demodialog.MS.nist<-function(){
assign("DemoFlag", 1, envir=as.environment("e1"))
dialog.MS.nist()
}
demobutton<-function(){
assign("DemoFlag", 1, envir=as.environment("e1"))
dialog.MS.demon()
}

#
#Usal functions
#

nodialog.MS.step1<-function(){
if (exists("DemoFlag", envir=as.environment("e1"))) rm("DemoFlag", envir=as.environment("e1"))
dialog.MS.step1()
}
nodialog.MS.step2<-function(){
if (exists("DemoFlag", envir=as.environment("e1"))) rm("DemoFlag", envir=as.environment("e1"))
dialog.MS.step2()
}
nodialog.MS.step3<-function(){
if (exists("DemoFlag", envir=as.environment("e1"))) rm("DemoFlag", envir=as.environment("e1"))
dialog.MS.step3()
}
nodialog.MS.DataCreation.user<-function(){
if (exists("DemoFlag", envir=as.environment("e1"))) rm("DemoFlag", envir=as.environment("e1"))
dialog.MS.DataCreation.user()
}

nodialog.MS.DataCreation.Agilent<-function(){
if (exists("DemoFlag", envir=as.environment("e1"))) rm("DemoFlag", envir=as.environment("e1"))
dialog.MS.DataCreation.Agilent()
}

nodialog.MS.DataCreation.CDF<-function(){
if (exists("DemoFlag", envir=as.environment("e1"))) rm("DemoFlag", envir=as.environment("e1"))
dialog.MS.DataCreation.CDF()
}

nodialog.MS.DataCreation.ASCII<-function(){
if (exists("DemoFlag", envir=as.environment("e1"))) rm("DemoFlag", envir=as.environment("e1"))
dialog.MS.DataCreation.ASCII()
}

nodialog.MS.test.clust<-function(){
if (exists("DemoFlag", envir=as.environment("e1"))) rm("DemoFlag", envir=as.environment("e1"))
dialog.MS.test.clust()
}

nodialog.MS.clust<-function(){
if (exists("DemoFlag", envir=as.environment("e1"))) rm("DemoFlag", envir=as.environment("e1"))
dialog.MS.clust()
}



#
# Menu setup
#
	frame0 <- tkframe(tt, relief="groove", borderwidth=2, background="grey")

	topMenuMSeasy <- tkmenubutton(frame0, text="Menu", background="grey")
	MSeasyMenu <- tkmenu(topMenuMSeasy, tearoff=TRUE)
	tkconfigure(topMenuMSeasy, menu=MSeasyMenu)
	
	openRecentMenu <- tkmenu(topMenuMSeasy,tearoff=FALSE)
		tkadd(openRecentMenu,"command",label="CDF or mzXML files and peaklist.txt",	command=function() nodialog.MS.DataCreation.CDF())
		tkadd(openRecentMenu,"command",label="CDF and rteres.txt (Agilent Technologies)",    command=function() nodialog.MS.DataCreation.Agilent())
		tkadd(openRecentMenu,"command",label="ASCII Files",	command=function() nodialog.MS.DataCreation.ASCII())
		tkadd(openRecentMenu,"command",label="User made Files",	command=function() nodialog.MS.DataCreation.user())
		
	
	openRecentMenu1 <- tkmenu(topMenuMSeasy,tearoff=FALSE)
		tkadd(openRecentMenu1,"command",label="MS.clust",    command=function() nodialog.MS.clust())
		tkadd(openRecentMenu1,"command",label="MS.test.clust",	command=function() nodialog.MS.test.clust())
	
	openRecentMenu2 <- tkmenu(topMenuMSeasy,tearoff=FALSE)
		tkadd(openRecentMenu2,"command",label="Load an MSeasy demo Dataset",    command=function() demochoose())
		tkadd(openRecentMenu2,"command",label="MS.DataCreation (CDF)",    command=function() demodialog.MS.DataCreation.CDF())
		tkadd(openRecentMenu2,"command",label="MS.DataCreation (Agilent)",    command=function() demodialog.MS.DataCreation.Agilent())
		tkadd(openRecentMenu2,"command",label="MS.DataCreation (ASCII)",	command=function() demodialog.MS.DataCreation.ASCII())
		tkadd(openRecentMenu2,"command",label="MS.test.clust",	command=function() demodialog.MS.test.clust())
		tkadd(openRecentMenu2,"command",label="MS.clust",	command=function() demodialog.MS.clust())
		tkadd(openRecentMenu2,"command",label="MSeasyToARISTO",command=function()demodialog.MS.aristo())
		tkadd(openRecentMenu2,"command",label="MSeasyToMSP",command=function()demodialog.MS.msp())
		tkadd(openRecentMenu2,"command",label="SearchNIST",command=function()demodialog.MS.nist())
	
	
	openRecentMenu3 <- tkmenu(topMenuMSeasy,tearoff=FALSE)
		tkadd(openRecentMenu3,"command",label="MSeasyToARISTO",    command=function() dialog.MS.aristo())
		tkadd(openRecentMenu3,"command",label="MSeasyToMSP",    command=function() dialog.MS.msp())
		tkadd(openRecentMenu3,"command",label="SearchNIST",   command=function() dialog.MS.nist())
		
	tkadd(MSeasyMenu,"cascade",label="Step1-GC-MS Data formatting...",menu=openRecentMenu)
	tkadd(MSeasyMenu,"separator")	
	tkadd(MSeasyMenu,"cascade",label="Step2-Mass spectra clustering...",menu=openRecentMenu1)
	tkadd(MSeasyMenu,"separator")	
	tkadd(MSeasyMenu,"cascade",label="Step3-Export to NIST or ARISTO...",menu=openRecentMenu3)
	tkadd(MSeasyMenu,"separator")	
	tkadd(MSeasyMenu,"cascade",label="Demonstration",menu=openRecentMenu2)
	tkadd(MSeasyMenu,"separator")
	tkadd(MSeasyMenu,"command",label="Quit R",command=function() q())
	
	topMenuHelp <- tkmenubutton(frame0, text="Help", background="grey")
	MSeasyMenu <- tkmenu(topMenuHelp, tearoff=TRUE)
	tkconfigure(topMenuHelp, menu=MSeasyMenu)
	tkadd(MSeasyMenu,"command",label="MSeasy(Help)",command=function()print(help("MSeasy")) )
	tkadd(MSeasyMenu,"command",label="trans.ASCII (Help)",command=function()print(help("trans.ASCII")) )
	tkadd(MSeasyMenu,"command",label="MS.DataCreation (Help)",command=function()print(help("MS.DataCreation")) )
	tkadd(MSeasyMenu,"command",label="MS.test.clust (Help)",command=function()print(help("MS.test.clust")) )
	tkadd(MSeasyMenu,"command",label="MS.clust (Help)",command=function()print(help("MS.clust")) )
	tkadd(MSeasyMenu,"command",label="MSeasyToARISTO (Help)",command=function()print(help("MSeasyToARISTO")) )
	tkadd(MSeasyMenu,"command",label="MSeasyToMSP (Help)",command=function()print(help("MSeasyToMSP")) )
	tkadd(MSeasyMenu,"command",label="SearchNIST (Help)",command=function()print(help("SearchNIST")) )
	
	topMenuURL <- tkmenubutton(frame0, text="Links", background="grey")
	MSeasyMenu <- tkmenu(topMenuURL, tearoff=TRUE)
	tkconfigure(topMenuURL, menu=MSeasyMenu)
	tkadd(MSeasyMenu,"command",label="MSeasy web site",command=function()shell.exec("http://sites.google.com/site/rpackagemseasy/") )
	tkadd(MSeasyMenu,"command",label="NIST ms search web site",command=function()shell.exec("http://chemdata.nist.gov/mass-spc/ms-search/") )
	tkadd(MSeasyMenu,"command",label="ARISTO web site",command=function()shell.exec("http://www.ionspectra.org/aristo/") )
	
	tkpack(topMenuMSeasy,topMenuHelp,topMenuURL,side="left")
	tkpack(frame0, expand="TRUE", fill="x")
#
# title and icons
#
	frame1 <- tkframe(tt, relief="groove", borderwidth=2, background="white")

		labh <- tklabel(frame1, bitmap="questhead", background="white")
	tkbind(labh, "<Button-1>", function() print(help("MSeasyTkGUI")))
	titre <- tklabel(frame1,text="MSeasy even more easy", font="Times 14", foreground="red", background="white")
	
	helplab <- tklabel(frame1,text="- Right click buttons for help - Double click in lists to select -", font="Times 11", foreground="dark green", background="white")


	tkgrid(titre, labh, padx=10, sticky = "w")
	tkgrid(helplab, columnspan=4)
	tkpack(frame1, expand="TRUE", fill="x")
#
# MSeasy step 1 button
#
	frame1b <- tkframe(tt,  borderwidth=2, background="white")
	tkpack(tklabel(frame1b,text="- MSeasy -", font="Times 12", foreground="blue", background="white"))	
	step1.but <- tkbutton(frame1b, text="Step1- GC-MS Data formatting",width=20,  command=function() nodialog.MS.step1())
	
	tkpack(step1.but,ipadx=25, side="bottom", expand="TRUE")
	tkpack(frame1b, expand="TRUE", fill="x")

	tkbind(step1.but, "<Button-3>", function() print(help("MS.DataCreation")))
	

#
# MSeasy Step 2 button
#
	frame2 <- tkframe(tt,  borderwidth=2, background="white")
	step2.but <- tkbutton(frame2, text="Step2- Mass spectra clustering", width=20, command=function() nodialog.MS.step2())
	
	tkpack(step2.but,ipadx=25, side="top", expand="TRUE")
	tkpack(frame2, expand="TRUE", fill="x")
	
	tkbind(step2.but, "<Button-3>", function() print(help("MS.clust")))
#
# MSeasy Step 3 button
#
	frame2 <- tkframe(tt,  borderwidth=2, background="white")
	step2.but <- tkbutton(frame2, text="Step3- Export to ARISTO or NIST", width=20,command=function() nodialog.MS.step3())
	
	tkpack(step2.but,ipadx=25, side="top", expand="TRUE")
	tkpack(frame2, expand="TRUE", fill="x")
	
	tkbind(step2.but, "<Button-3>", function() print(help("MS.clust")))	
#
# MSeasy Demo button
#
	frame3 <- tkframe(tt, relief="groove", borderwidth=2, background="white")
	tkpack(tklabel(frame3,text="- Demonstration -", font="Times 12", foreground="blue", background="white"))
	demon.but <- tkbutton(frame3, text="Demonstration", width=20, command=function() demobutton())
	
	tkpack(demon.but,ipadx=25, side="top", expand="TRUE")
	tkpack(frame3, expand="TRUE", fill="x")
	
	tkbind(demon.but, "<Button-3>", function() print(help("MSeasyTkGUI")))	
	
#
# Quit
#
	frame5 <- tkframe(tt, relief="groove", borderwidth=2, background="white")
	cancel.but <- tkbutton(frame5, text="Cancel", command=function() tkdestroy(tt), font="Times 14")
	quity.but <- tkbutton(frame5, text="Quit R (save)", command=function() q("yes"), font="Times 14")
	quitn.but <- tkbutton(frame5, text="Quit R (don't save)", command=function() q("no"),  font="Times 14")
	tkpack(quity.but, cancel.but, quitn.but, side="left", expand="TRUE", fill="x")
	tkpack(frame5, expand="TRUE", fill="x")
	tkfocus(tt)
	return(invisible())
##on.exit
on.exit(detach(e1))
}

