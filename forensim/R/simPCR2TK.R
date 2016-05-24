"simPCR2TK" <- function()
{
	# if(!require(tcltk)) stop("package tcltk is required")
	# if(!require(tkrplot)) stop("package tkrplot is required")
	# tclRequire("Tktable")
	tclRequire("Tktable")
	font0 <- tkfont.create(family="courrier",size=35,weight="bold",slant="italic")
	font1<-tkfont.create(family="times",size=14,weight="bold")#,slant="italic")
	font2<-tkfont.create(family="times",size=16,weight="bold",slant="italic")
	font3<-tkfont.create(family="times",size=12)#,slant="italic")
	font4<-tkfont.create(family="courrier",size=14)#,slant="italic")
	font5<-tkfont.create(family="courrier",size=13,weight="bold")#,slant="italic")
	font6<-tkfont.create(family="times",size=13)#tkframe entries labels
	tf <- tktoplevel()
	tkwm.title(tf,"simPCR2TK: a graphical simulation interface the PCR")
	
	done <- tclVar(0)
		
	frame1 <- tkframe(tf, relief="groove", borderwidth=4)
	#icn <- tkimage.create("photo", file=system.file("files/test.GIF", package = "forensim"))#"test.GIF")
	#TclTklabel <- tklabel(frame1, image=icn, background="white")
	labh <- tklabel(tf)#, image=icn)
	#labh <- tklabel(frame1)
	#tkbind(labh, "<Button-1>", function() 'hh')
	tkgrid(tklabel(frame1,text="   simPCR2TK: a graphical interface for PCR simulation   ", font="courrier 22", 
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
	inFrame3 <- tkframe(frame2, relief="groove", borderwidth=2)

	#define new entries for input parameters
	#
	# Variables for text fields
	# Input parameters
	#PCR FRAME
	ncells<-tclVar(5)
	Pextrac <- tclVar(0.6)
	Paliquot<-tclVar(0.30)
	PPCR<-tclVar(0.80)
	Tdrop<-tclVar(2e+07)#threshold of detection for each allele signal
	Tcyc<-tclVar(28)#of PCR cycles
	#ncells.entry <- tkentry(inFrame, textvariable=ncells, width=4,relief="solid",justify="center")
	ncells.entry<-tkentry(inFrame,textvariabl=ncells,width=4,relief="solid",justify="center")
	Pextrac.entry <-tkentry(inFrame, textvariable=Pextrac, width=4,relief="solid",justify="center")
	Paliquot.entry<-tkentry(inFrame, textvariable=Paliquot, width=4,relief="solid",justify="center")
	#repl.entry<-tkentry(frame2, textvariable=repl, width=5)
	PPCR.entry<-tkentry(inFrame2, textvariable=PPCR, width=4,relief="solid",justify="center")
	Tcyc.entry<-tkentry(inFrame2, textvariable=Tcyc, width=4,relief="solid",justify="center")
	Tdrop.entry<-tkentry(inFrame2, textvariable=Tdrop, width=8,relief="solid",justify="center")
	
	#tkgrid(tklabel(inFrame, text="________________",font=font3,foreground="blue"), columnspan=9)
	tkgrid(tklabel(inFrame, text="pre-PCR parameters",font=font4,foreground="blue"), columnspan=9)
	#tkgrid(tklabel(inFrame,text="#Cells",font=font6), ncells.entry,sticky="we")#title varibale ncells.entry
	tkgrid(tklabel(inFrame,text="number of cells",font=font6), ncells.entry,sticky="we")#title varibale ncells.entry
	tkgrid(tklabel(inFrame,text="Prob. extraction",font=font6), Pextrac.entry,sticky="we")
	tkgrid(tklabel(inFrame,text="Prob. surviving aliquot",font=font6), Paliquot.entry,sticky="we")
	#tkgrid(tklabel(inFrame, text="________________",font=font3,foreground="blue"), columnspan=9)
	
	tkgrid(tklabel(inFrame2, text="PCR parameters",font=font4,foreground="blue"), columnspan=9)
	tkgrid(tklabel(inFrame2,text="Prob. PCR efficiency",font=font6), PPCR.entry,sticky="n")
	tkgrid(tklabel(inFrame2,text="#PCR cycles",font=font6), Tcyc.entry,sticky="n")#title varibale ncells.entry
	tkgrid(tklabel(inFrame2,text="Allele detection threshold",font=font6), Tdrop.entry,sticky="n")
	#tkgrid(tklabel(inFrame2, text="________________",font=font3,foreground="blue"), columnspan=9)
	
	tkgrid(tklabel(inFrame3, text="Cells ploidy",font=font4,foreground="blue"), columnspan=9)
	dip<-tclVar(TRUE)#ploidy variable, TRUE cells are diploid
	hap<-tclVar(0.50)#probability of encountering allele of type A in haploid cells
	dip.entry<-tkentry(inFrame3, textvariable=hap, width=4,state="disabled",relief="solid",justify="center")
	
	tkgrid(tkcheckbutton(inFrame3, text="Diploid", variable=dip,font=font6, 
	command=function() if (!as.logical(tclObj(dip))) tkconfigure(dip.entry, state="normal") else tkconfigure(dip.entry, state="disabled") ))
	tkgrid(tklabel(inFrame3,text="prob. allele A in haploid cells",font=font6), dip.entry,sticky="n")

    repFrame<-tkframe(frame3, relief="flat", borderwidth=2)
    repl<-tclVar(10)
	repl.entry<-tkentry(repFrame,textvariable=repl,width=3,highlightthickness=2,relief="solid",justify="center")
	tkgrid(tklabel(repFrame,text="# Replicate simulations",font=font6),repl.entry,sticky="we")
    tkgrid(inFrame, inFrame2, inFrame3, padx=12,pady=8 )
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
	########### function simPCR2TK.local in tcltk 
	# main frame
	simPCR2TK.local<-function()
	{
		#first, get the parameters, check their validity
		dip<-(tclvalue(dip))#ploidy: 1= dip, 0= hap
		hap<-as.numeric(tclvalue(hap))#probability of encountering allele of type A in haploid cells
		ncells<-ncells.verif()
		#n<-ncells.verif()
		p1<-Pextrac.verif()
		p2<-Paliquot.verif()
		p3<-PPCR.verif()
		cyc<-Tcyc.verif()
		repl<-repl.verif()
		Tdrop<-Tdrop.verif()#allele detection threshold
		# an then get the results from the pcrsim function
		# all replicates are concatenated in a singla table
		tmp<-NULL
		
		#diploid case
		if(dip==1)
		{
			for(i in 1:repl)
			{
				tmp<-rbind(tmp,simPCR2(ncells=ncells,p1,p2,p3,cyc=cyc,Tdrop=Tdrop,dip=TRUE))
			}
		}
		else{
			for(i in 1:repl)
			{
				tmp<-rbind(tmp,simPCR2(ncells=ncells,p1,p2,p3,cyc=cyc,Tdrop=Tdrop,probSperm=hap,dip=FALSE,KH=55))
			}
		}
		
		save<-function()
		{
			tclarray <- tclArray()
			myRarray<- c("HeightA",tmp[,1],"DropoutA",tmp[,2],"HeightB",tmp[,3],"DropoutB",tmp[,4])
			dim(myRarray)<-c(repl+1,4)
			for(i in 0:(repl)){
			for (j in (0:3)){
			tclarray[[i,j]] <- myRarray[i+1,j+1]}}
			tt <- tktoplevel()
			tkwm.title(tt,"simPCR2: Edit simulation results")
			frame0 <- tkframe(tt, relief="flat", borderwidth=4)
			#frameHb<-tkframe(tt,relief="flat",borderwidth=4)
			tkgrid(tklabel(frame0,text="Edit simulation results", font="courrier 18", foreground="darkblue"))

			framet<-tkframe(tt, relief="groove", borderwidth=2)

			table1 <- tkwidget(framet,"table",variable=tclarray,rows=(repl+1),cols=4,titlerows=1,titlecols=0,
			yscrollcommand=function(...) tkset(yscr,...))
			#xscr <-tkscrollbar(tt,orient="horizontal", command=function(...)tkxview(table1,...))
			yscr <- tkscrollbar(framet,repeatinterval=5, command=function(...)tkyview(table1,...))
			
			#tkgrid.configure(yscr,sticky="nsw")
			tkpack(frame0,framet)
			tkgrid(table1,yscr,sticky="new", padx=20,pady=18)
			tkgrid.configure(yscr,sticky="nsw")
			#tkconfigure(table1,variable=tclarray,background="white",selectmode="extended")
			tkconfigure(table1,variable=tclarray,background="white",selectmode="extended",rowseparator="\"\n\"",colseparator="\"\t\"")

			#tkpack(table1,frame0)
			# function to save the simulation reuslts into a text file
			saveFile<-function()
			{
				
				ff<-tktoplevel()
				Fframe<- tkframe(ff, relief="groove")
				tkwm.title(ff,"Filenames")
				tkgrid(tklabel(Fframe, text="- Enter filename -",font=font3,foreground="darkblue"), columnspan=9)
				filtervar<- tclVar(paste('simutable',repl,'.txt',sep=""))
				filtervar.entry <- tkentry(Fframe, textvariable=filtervar, width=12)
				saveF.butt<-tkbutton(Fframe, text="Enter",fg="darkblue")
				tkgrid(filtervar.entry, saveF.butt)		
				tkpack(Fframe,padx=12,pady=18,side="left")
				filen<-tclvalue(filtervar)
				write.table(tmp,file=filen,row.names=FALSE)
			}
		
			#filename variable and entry
			#save.entry<-tkentry(framet, textvariable=filen, width=18)
			save.button<-tkbutton(tt, text="Save table",font="courrier 10", fg="darkblue",command=function() saveFile())
			cancel.but1 <- tkbutton(tt, text="Dismiss",font="courrier 10", fg="darkblue", command=function() tkdestroy(tt))
			tkpack(cancel.but1,save.button,side="left", fill="x", expand=1)
		}
		#frame where the simulations table will be displayed
		logis.frames<-function(res)
		{
			ll<-tktoplevel()
			tkwm.title(ll,"simPCR2TK: Logistic regression")
			#icn.log<-tkimage.create("photo",file="test.GIF")
			
			logFramet<-tkframe(ll,relief="groove",borderwidth=2)
			frameunder<-tkframe(ll,relief="groove",borderwidth=2)#,bg="white")
			logFrame<-tkframe(frameunder,relief="groove",borderwidth=2)
			dataFrame<-tkframe(frameunder,relief="groove",borderwidth=2)
			outFrame<-tkframe(frameunder,relief="groove",borderwidth=2)
			goFrame<-tkframe(frameunder,relief="flat",borderwidth=2)
			
			lablog<-tklabel(logFramet)#,image=icn.log)
			tkgrid(tklabel(logFramet,text="Logistic regression on dropout data",font="courrier 18",foreground="darkblue"),lablog)
			tkgrid(tklabel(logFrame, text="Model",font=font4,foreground="blue"), columnspan=9)
			
			intercept<-tclVar(1)#initial value for intercept
			outplot<-tclVar(1)
			samp1<-tclVar(TRUE)
			dev<-tclVar(1)
			dataselec<-tclVar(repl)
			sampdata.entry<-tkentry(dataFrame, textvariable=dataselec, state="disabled",width=4,relief="solid",justify="center")
			#- intercept- #
			tkgrid(tkcheckbutton(logFrame, text="Include intercept?", variable=intercept,font=font6))
				
			tkgrid(tklabel(dataFrame, text="Data",font=font4,foreground="blue"), columnspan=9)
					
			tkgrid(tkcheckbutton(dataFrame, text="All data", variable=samp1,font=font6, 
			command=function() if (!as.logical(tclObj(samp1))) tkconfigure(sampdata.entry, state="normal") else tkconfigure(sampdata.entry, state="disabled") ))
			tkgrid(tklabel(dataFrame,text="# points",font="times 12"), sampdata.entry,sticky="n")

			#data sekection
						
			tkgrid(tklabel(outFrame, text="Output",font=font4,foreground="blue"), columnspan=9)
			tkgrid(tkcheckbutton(outFrame, text="Residuals diagnosis", variable=outplot,font=font6), sticky="we")
			#tkgrid(tkcheckbutton(outFrame, text="Deviance tests", variable=dev,font=1), sticky="we")
		
			tkpack(logFramet)
			tkpack(logFrame,dataFrame,outFrame,padx=18,pady=18)
			
			#function to verify the entry in the # of selected points
			dataselec.verif<-function()
			{
				p<-as.numeric(tclvalue(dataselec))
				if(p <0 || p > repl){tkmessageBox(message="Invalid value for the number of selected points",icon="error",type="ok")
				}
				else{return(p)}
			}
			#perfroms the logitic regression according to the chosen parameters
			logi.regression<-function()
			{
				
				#print(tclvalue(outplot))
				
				# logitic model #
				#intercept
				intercept<-as.numeric(tclvalue(intercept))
				samp1<-as.numeric(tclvalue(samp1))
				#print('samp:');print(samp1)
				dataselec<-as.numeric(tclvalue(dataselec))
				outplot<-as.numeric(tclvalue(outplot))
				dev<-as.numeric(tclvalue(dev))
				#subset of data
				if(samp1==1)
				{
					tab1<-tmp[,c(2,3)]
					colnames(tab1)<-c('Dropout','PHeights')
					tab2<-tmp[,c(4,1)]
					colnames(tab2)<-c('Dropout','PHeights')
					res<-rbind(tab1,tab2)
					if(intercept==1){ glm0<-glm(res[,1]~res[,2],binomial)}
					else{glm0<-glm(res[,1]~res[,2]-1,binomial)}
					#if(outplot)	plot(glm100,which=1)
				}
				else
				{	
					s<-1:nrow(tmp)#select the rows to include in the calculation
					tmp2<-tmp[sample(s,dataselec.verif(),replace=FALSE),]#sample data 
					tab1<-tmp2[,c(2,3)]
					colnames(tab1)<-c('Dropout','PHeights')
					tab2<-tmp2[,c(4,1)]
					colnames(tab2)<-c('Dropout','PHeights')
					res<-rbind(tab1,tab2)
					if(intercept){ glm0<-glm(res[,1]~res[,2],binomial)}
					else{glm0<-glm(res[,1]~res[,2]-1,binomial)}
					#if(outplot)	plot(glm100,which=1)
				}
				
				# regression details
				
				
				dd <- tktoplevel()
				tkwm.title(dd,"Logistic regession results")
				frameC <- tkframe(dd, relief="flat", borderwidth=4)
				tclcoeff <- tclArray()	#tcltk array that wil receive the results
				tclcoeff2<-tclArray() #tcltk array that wil receive the results
				resume<-summary(glm0)$coeff#summary function applied to a glm objcet: deviance tests and estimates
				
				#print(resume)
				#if intrecept is included in the calculations...
				if(intercept==1)
				{
					resume<-matrix(c("Coefficient","B0","B1","Estimate",signif(resume[,1],5),"Std.error",signif(resume[,2],5),"z.value", signif(resume[,3],5),
					"Pr(>|z|)", signif(resume[,4],5)),ncol=5)
					a0<-signif(glm0$coeff[1],5)
					b0<-signif(glm0$coeff[2],5)
					
					tab.coeff<-matrix(c("B0",a0,"B1",b0),ncol=2)
					#print(tab.coeff)
					#tab.coeff containing the data is copied to tclcoeff (the arry objcet)
					
					for(j in 0:2){
					for(i in 0:4){
					tclcoeff2[[j,i]] <- resume[j+1,i+1]}}
					#print('ok')
					
					for(j in 0:1){
					for(i in 0:1){
					tclcoeff[[j,i]] <- tab.coeff[j+1,i+1]}}
					#if the intercept is included
					tableC2<-tkwidget(frameC,"table",background="white",variable=tclcoeff2,colwidth=18,rows=3,cols=5,titlerows=1,titlecols=1)

				}
				#if no intercept
				else{
					resume<-matrix(c("Coefficient","B1","Estimate",signif(resume[,1],5),"Std.error",signif(resume[,2],5),"z.value", signif(resume[,3],5),
					"Pr(>|z|)", signif(resume[,4],5)),ncol=5)
					for(j in 0:1){
					for(i in 0:4){
					tclcoeff2[[j,i]] <- resume[j+1,i+1]}}
					#if the intercept is not included
					tableC2<-tkwidget(frameC,"table",background="white",variable=tclcoeff2,colwidth=18,rows=2,cols=5,titlerows=1,titlecols=1)

				}
					
					
				#tkgrid(tklabel(frameC,text="Estimates", font=font4, foreground="darkblue"))
				tableC<-tkwidget(frameC,"table",background="white",variable=tclcoeff,colwidth=18,rows=2,cols=(1+intercept),titlerows=1,titlecols=0)

				#tkgrid(tableC)
				tkgrid(tklabel(frameC,text="Estimates, deviance test", font="times 18", foreground="darkblue"))

				tkgrid(tableC2,pady=12)
				#function to coopy deviance test and estimates results to filetext
				saveLog<-function(filen)
				{
			
					tl<-tktoplevel()
					Lframe<- tkframe(tl, relief="groove")
					tkwm.title(tl,"Filenames")
					tkgrid(tklabel(Lframe, text="- Enter filename -",font=font3,foreground="darkblue"), columnspan=9)
					filtervar<- tclVar(paste('logres',repl,'.txt',sep=""))
					filtervar.entry <- tkentry(Lframe, textvariable=filtervar, width=12)
					saveF.butt<-tkbutton(Lframe, text="Enter",fg="darkblue")
					tkgrid(filtervar.entry, saveF.butt)		
					tkpack(Lframe,padx=12,pady=18,side="left")
					filen<-tclvalue(filtervar)
					write.table(resume,file=filen,col.names=FALSE,row.names=FALSE)
			
				}
				but<-tkbutton(frameC,text="Save",fg="darkblue",font="courrier 10",command=function() saveLog())
				tkgrid(but)
				#tkpack(frameC)
				#print(summary(glm0))
				
			
				
			if(outplot==1)
			{
				
				
				# dd <- tktoplevel()
				# tkwm.title(dd,"Residuals diagnosis")
				 Myhscale <- 1.55    # Horizontal scaling
				 Myvscale <- 1.55	 # Vertical scaling

				plotglm0<-function()
				{
					params <- par(bg="white")
					#plot(glm0,which=1)
					plot(predict(glm0,link="response"),residuals(glm0,type="pearson"),las=1,xlab="Predicted values",ylab="Residuals",cex.lab=1.3)
					abline(h=0,lty=3,col="red")
					title(main=list('Residuals diagnosis',cex=1.5))
					par(params)
				}
				img <- tkrplot(frameC,fun= plotglm0,hscale=Myhscale,vscale=Myvscale)
				CopyToClip <- function()
				{
				  tkrreplot(img)
				}
				copy.but <- tkbutton(frameC,text="Copy to Clipboard",font="courrier 10",fg="darkblue",command=CopyToClip)
				tkgrid(img)
				tkgrid(copy.but)
				#tkpack(frameC)
		
			}
			
			#if the results or the plots are to be displayed, the frameC is packed, otherwise, nothing happens
			tkpack(frameC)
			#else{tkdestroy(dd)}
			}#this button will dislay the results for the logitic regression with the set parameters
			go.but <- tkbutton(ll, text="Go!", font=font5,fg="darkblue",command=logi.regression)
			cancel.log <- tkbutton(ll, text="Dismiss", font=font5,fg="darkblue", command=function() tkdestroy(ll))
			#tkpack(goFrame)
			tkpack(frameunder, padx=12,pady=12)
			tkpack(cancel.log, go.but, side="left", fill="x", expand=1)
		
		}
		hh<-tktoplevel()
		tkwm.title(hh,"simPCR2TK.local: Simulation results")
		#icn2 <- tkimage.create("photo", file="test.GIF")
		sFrame<- tkframe(hh, relief="groove",borderwidth=4)
		labh <- tklabel(sFrame,bitmap="questhead")#, image=icn2)
		tkgrid(tklabel(sFrame,text="    Simulation results    ", font="courrier 18", foreground="darkblue"), labh)
		editFrame<-tkframe(hh,relief="sunken",borderwidth=2,bg="white")
		tabFrame<-tkframe(hh, relief="sunken",borderwidth=2,bg="white")
		butFrame<-tkframe(hh,relief="flat")
		
		#tkgrid(tklabel(editFrame,text="Simulations", font="courrier 14", foreground="blue"))
		tab.but <- tkbutton(editFrame, text="Edit simulations", highlightthickness=1,fg="blue",font=font5, command= save)
		
		
		#tkgrid(tklabel(tabFrame,text="Logistic regression", font="courrier 14", foreground="blue"))
		glm.but<-tkbutton(tabFrame,text="Logistic regression",font=font5,highlightthickness=1,highlightbackground="black",fg="blue",command=function() logis.frames(tmp))
		
		glm.cancel<- tkbutton(butFrame, text="Dismiss", fg="darkblue",font=font3,command=function() tkdestroy(hh))
		tkgrid(tab.but,padx=18,pady=9)
		tkgrid(glm.but, padx=8,pady=19)
		tkgrid(glm.cancel)
		tkpack(sFrame)
		tkpack(editFrame, tabFrame,side="top",padx=15,pady=19)
		tkpack(butFrame,side="bottom")
	}

	#Bottom buttons #
	ok.but <- tkbutton(tf, text="Simulate!", font=font5,fg="darkblue",command=simPCR2TK.local)
	cancel.but <- tkbutton(tf, text="Dismiss", font=font5,fg="darkblue", command=function() tkdestroy(tf))
	tkpack(cancel.but, ok.but, side="left", fill="x", expand=1)

}
#simPCR2TK()