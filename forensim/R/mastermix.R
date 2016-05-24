#MasterMix for forensim
#Hinda Haned,
#July 2009
mastermix<-function()
{
	# if(!require(tcltk)) stop("package tcltk is required")
	#tclRequire("Tktable")
	#fonts definition
	font0 <- tkfont.create(family="times",size=35,weight="bold",slant="italic")
	font1<-tkfont.create(family="times",size=14,weight="bold")#,slant="italic")
	font2<-tkfont.create(family="times",size=16,weight="bold",slant="italic")
	font3<-tkfont.create(family="times",size=10,weight="bold")#,slant="italic")
	font4<-tkfont.create(family="times",size=10)#,slant="italic")
	#
	tt <- tktoplevel()
	tkwm.title(tt,"forensic DNA mixtures resolution")	
	dudioutvar <- tclVar()
	dudivar <- tclVar()
	facvar <- tclVar()
	nfvar <- tclVar()

	scannfvar <- tclVar(1)
#
# Title
#
	TFrame <- tkframe(tt, relief="groove")
	labh <- tklabel(TFrame)
	font0 <- tkfont.create(family="times",size=35,weight="bold",slant="italic")
	font1<-tkfont.create(family="times",size=14,weight="bold")#,slant="italic")
	font2<-tkfont.create(family="times",size=16,weight="bold",slant="italic")
	
	tkgrid(tklabel(TFrame,text="MasterMix", font=font0, foreground="red"), labh)
	tkgrid(tklabel(TFrame,text="Two-person DNA mixtures resolution using allele peak height/ or area information", font=font2, foreground="red"), labh)
	#tkbind(labh, "<Button-1>", function() print(help("between")))
	tkgrid(TFrame)

	



#____________________________________________
	
RCSFrame <- tkframe(tt, relief="groove")
A2.but <- tkbutton(RCSFrame, text="Two-allele model", font=font1, command=A2.simu)
tkbind(A2.but, "<Button-3>", function() print(help("A2.simu")))
A3.but <- tkbutton(RCSFrame, text="Three-allele model",font=font1, command=A3.simu)
tkbind(A3.but, "<Button-3>", function() print(help("A3.simu")))
A4.but <- tkbutton(RCSFrame, text="Four-allele model", font=font1,command=A4.simu)
tkbind(A4.but, "<Button-3>", function() print(help("A4.simu")))
tkgrid(A2.but, A3.but,A4.but,ipadx=20)	
tkgrid(RCSFrame)

}




	# tkbind(getdata.but, "<Button-3>", function() print(help("data")))








