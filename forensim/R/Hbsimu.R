# Hinda Haned, May 2010- Lyon
"Hbsimu" <- function()
{
	# if(!require(tcltk)) stop("package tcltk is required")
	# if(!require(tkrplot)) stop("package tkrplot is required")
	tclRequire("Tktable")
	# tclRequire("Tktable")
	font0 <- tkfont.create(family="courrier",size=35,weight="bold",slant="italic")
	font1<-tkfont.create(family="times",size=14,weight="bold")#,slant="italic")
	font2<-tkfont.create(family="times",size=16,weight="bold",slant="italic")
	font3<-tkfont.create(family="times",size=12)#,slant="italic")
	font4<-tkfont.create(family="courrier",size=14)#,slant="italic")
	font5<-tkfont.create(family="courrier",size=13,weight="bold")#,slant="italic")
	font6<-tkfont.create(family="times",size=13)#tkframe entries labels
	tf <- tktoplevel()
	tkwm.title(tf,"Hbdemo: heterozygote balance simulation")
	
	
	
	
	#function Hb

	done <- tclVar(0)
		
	frame1 <- tkframe(tf, relief="groove", borderwidth=4)
	labh <- tklabel(tf)
	#labh <- tklabel(frame1)
	#tkbind(labh, "<Button-1>", function() 'hh')
	tkgrid(tklabel(frame1,text="   Heterozygote balance simulator  ", font="courrier 22", 
	foreground="darkblue"),labh)
	#tkbind(frame1)
	tkpack(frame1)#,padx=10)#,pady=10)

	frame2 <- tkframe(tf, relief="flat", borderwidth=2)#,bg="white")	
	frame3 <- tkframe(tf, relief="groove", borderwidth=2)#,bg="white")	
	# tab.entry <- tkentry(frame2, textvariable=tabnamevar)
	# file.entry <- tkentry(frame2, textvariable=filenamevar)
	# choosefile.but <- tkbutton(frame2, text="Set", command=function() print('hh'))
	#tkgrid(tklabel(frame1, text="- Input parameters -",font=font3,foreground="blue"), columnspan=6)
	# tkgrid(tklabel(frame2,text=" "), file.entry)
	# tkgrid(tklabel(frame2,text="Dataframe to receive the data : "), tab.entry)
	# varnames.cbut <- tkcheckbutton(frame2,text="Variables names on the first row of data file", variable=varnames)
	# tkgrid(varnames.cbut, columnspan=2, sticky="w")
	
	inFrame <- tkframe(frame2, relief="groove", borderwidth=2)
	inFrame2 <- tkframe(frame2, relief="groove", borderwidth=2)
	#inFrame3 <- tkframe(frame2, relief="groove", borderwidth=2)

	#define new entries for input parameters
	#
	# Variables for text fields
	# Input parameters
	#PCR FRAME
	ncells<-tclVar(5)
	varncells<-tclVar(1)
	Pextrac <- tclVar(0.6)
	Paliquot<-tclVar(0.30)
	PPCR<-tclVar(0.80)
	Tdrop<-tclVar(2e+07)#threshold of detection for each allele signal
	Tcyc<-tclVar(28)#of PCR cycles
	#ncells.entry <- tkentry(inFrame, textvariable=ncells, width=4,relief="solid",justify="center")
	ncells.entry<-tkentry(inFrame,textvariabl=ncells,width=4,relief="solid",justify="center")
	varncells.entry<-tkentry(inFrame,textvariabl=varncells,width=4,relief="solid",justify="center")

	Pextrac.entry <-tkentry(inFrame, textvariable=Pextrac, width=4,relief="solid",justify="center")
	Paliquot.entry<-tkentry(inFrame, textvariable=Paliquot, width=4,relief="solid",justify="center")
	#repl.entry<-tkentry(frame2, textvariable=repl, width=5)
	PPCR.entry<-tkentry(inFrame2, textvariable=PPCR, width=4,relief="solid",justify="center")
	Tcyc.entry<-tkentry(inFrame2, textvariable=Tcyc, width=4,relief="solid",justify="center")
	Tdrop.entry<-tkentry(inFrame2, textvariable=Tdrop, width=8,relief="solid",justify="center")
	
	#tkgrid(tklabel(inFrame, text="________________",font=font3,foreground="blue"), columnspan=9)
	tkgrid(tklabel(inFrame, text="pre-PCR parameters",font=font4,foreground="blue"), columnspan=9)
	#tkgrid(tklabel(inFrame, text="Number of cells",font=font3,foreground="blue"), columnspan=9)
	#tkgrid(tklabel(inFrame,text="#Cells",font=font6), ncells.entry,sticky="we")#title varibale ncells.entry
	tkgrid(tklabel(inFrame,text="Mean #cells",font=font6), ncells.entry,sticky="we")#title varibale ncells.entry
	tkgrid(tklabel(inFrame,text="Variance #cells",font=font6), varncells.entry,sticky="we")#title varibale ncells.entry
	tkgrid(tklabel(inFrame,text="Prob. extraction",font=font6), Pextrac.entry,sticky="we")
	tkgrid(tklabel(inFrame,text="Prob. surviving aliquot",font=font6), Paliquot.entry,sticky="we")
	#tkgrid(tklabel(inFrame, text="________________",font=font3,foreground="blue"), columnspan=9)
	
	tkgrid(tklabel(inFrame2, text="PCR parameters",font=font4,foreground="blue"), columnspan=9)
	tkgrid(tklabel(inFrame2,text="Prob. PCR efficiency",font=font6), PPCR.entry,sticky="n")
	tkgrid(tklabel(inFrame2,text="#PCR cycles",font=font6), Tcyc.entry,sticky="n")#title varibale ncells.entry
	tkgrid(tklabel(inFrame2,text="Allele detection threshold",font=font6), Tdrop.entry,sticky="n")
	#tkgrid(tklabel(inFrame2, text="________________",font=font3,foreground="blue"), columnspan=9)
	
	# tkgrid(tklabel(inFrame3, text="Cells ploidy",font=font4,foreground="blue"), columnspan=9)
	# dip<-tclVar(TRUE)#ploidy variable, TRUE cells are diploid
	# hap<-tclVar(0.50)#probability of encountering allele of type A in haploid cells
	# dip.entry<-tkentry(inFrame3, textvariable=hap, width=4,state="disabled",relief="solid",justify="center")
	
	# tkgrid(tkcheckbutton(inFrame3, text="Diploid", variable=dip,font=font6, 
	# command=function() if (!as.logical(tclObj(dip))) tkconfigure(dip.entry, state="normal") else tkconfigure(dip.entry, state="disabled") ))
	# tkgrid(tklabel(inFrame3,text="prob. allele A in haploid cells",font=font6), dip.entry,sticky="n")

    repFrame<-tkframe(frame3, relief="flat", borderwidth=2)
    repl<-tclVar(100)
	repl.entry<-tkentry(repFrame,textvariable=repl,width=10,highlightthickness=2,relief="solid",justify="center")
	tkgrid(tklabel(repFrame,text="# Replicate simulations",font=font6),repl.entry,sticky="we")
    tkgrid(inFrame, inFrame2, padx=12,pady=8 )
	tkgrid(repFrame,padx=8)
	tkpack(frame2,frame3, pady=18, padx=28)#,side="top")
	
	#now, the following functions check that the user entered the correct values
	
	
	ncells.verif<-function()
	{
		
		q<-tclvalue(ncells)
		if(q=='NULL'){return(NULL)}
		else{
			if(!is.null(q) & q<0){tkmessageBox(message="Invalid value for the number of cells",icon="error",type="ok")}
			else{return(as.numeric(q))}
		}
	}
	
	varncells.verif<-function()
	{
		
		q<-tclvalue(varncells)
		if(q=='NULL'){return(NULL)}
		else{
			if(!is.null(q) & q<0){tkmessageBox(message="Invalid value for the variance of the number of cells",icon="error",type="ok")}
			else{return(as.numeric(q))}
		}
	}
	Pextrac.verif<-function()
	{
		p2<-as.numeric(tclvalue(Pextrac))
		if(p2<0 || p2>1){tkmessageBox(message="Invalid value for the probability of extraction",icon="error",type="ok")
		}
		else{return(p2)}
	}
	
	Paliquot.verif<-function()
	{
		p3<-as.numeric(tclvalue(Paliquot))
		if(p3<0 || p3>1){tkmessageBox(message="Invalid value, please enter a probability ([0,1])",icon="error",type="ok")
		}
		else{return(p3)}
	}
	
		
	PPCR.verif<-function()
	{
		p4<-as.numeric(tclvalue(PPCR))
		if(p4<0 || p4>1){tkmessageBox(message="Invalid value for the probability of PCR efficiency",icon="error",type="ok")
		}
		else{return(p4)}
	}
	Tcyc.verif<-function()
	{
		p5<-as.numeric(tclvalue(Tcyc))
		if(p5<=0) {tkmessageBox(message="Invalid value for the number of PCR cycles",icon="error",type="ok")
		}
		else{return(p5)}
	}
	repl.verif<-function()
	{
		p6<-as.numeric(tclvalue(repl))
		if(p6<1) {tkmessageBox(message="At least one replicate is needed",icon="error",type="ok")
		}
		else{return(p6)}
	}
	Tdrop.verif<-function()
	{
		p6<-as.numeric(tclvalue(Tdrop))
		if(p6<10^7) {tkmessageBox(message="Detection threshold, in number of molecules, must be at last 10^7",icon="error",type="ok")
		}
		else{return(p6)}
	}
	########### function Hb.local in tcltk 
	# main frame
	Hb.local<-function()
	{
		#first, get the parameters, check their validity
		
		ncells<-ncells.verif()
		varncells<-varncells.verif()
		#n<-ncells.verif()
		p1<-Pextrac.verif()
		p2<-Paliquot.verif()
		p3<-PPCR.verif()
		cyc<-Tcyc.verif()
		repl<-repl.verif()
		Tdrop<-Tdrop.verif()#allele detection threshold
		# an then get the results from the pcrsim function
		# all replicates are concatenated in a singla table
		
		
		#diploid case
		nsim<-as.integer(pmax(rnorm(repl,mean=ncells,sd=sqrt(varncells)),1))
		tab1<-t(sapply(nsim,function(i)simPCR2(i,probEx=p1,probAlq=p2,probPCR=p3,cyc=cyc,dip=TRUE) ))
		Hbsim<-na.omit(unlist(tab1[,5]))
		#hist(p1,col="gray",prob=T,xlab="Hb",ylab="dF(x)",cex.lab=1.3,main="Heterozygous balance, mu=1,sd=1\nProbability density function")
		
		dd <- tktoplevel()
		frameC <- tkframe(dd, relief="flat", borderwidth=4)
		tkwm.title(dd,"Hb distribution")
		Myhscale <- 1    # Horizontal scaling
		Myvscale <- 1

		# Vertical scaling
		plotHb<-function()
		{
			params <- par(bg="white")
			hist(Hbsim,col="gray",prob=TRUE,xlab="Hb",ylab="f(Hb)",cex.lab=1.3,main="Heterozygote balance\n Probability density function",xlim=c(0,1))
			par(params)
		}	

		img <- tkrplot(frameC,fun= plotHb,hscale=Myhscale,vscale=Myvscale)
		CopyToClip <- function()
		{
			tkrreplot(img)
		}
		copy.but <- tkbutton(frameC,text="Copy to Clipboard",font="courrier 10",fg="darkblue",command=CopyToClip)
		tkgrid(img)
		tkgrid(copy.but)
			
		tkpack(frameC)
	}
		
		
		
	

	#Bottom buttons #
	ok.but <- tkbutton(tf, text="Simulate!", font=font5,fg="darkblue",command=function() Hb.local())
	cancel.but <- tkbutton(tf, text="Dismiss", font=font5,fg="darkblue", command=function() tkdestroy(tf))
	tkpack(cancel.but, ok.but, side="left", fill="x", expand=1)

}
