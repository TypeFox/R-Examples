# Hinda, Hague, june 2011
# Hinda, Oslo Dec 2011
# Hinda, Hague May 2012
#GUI for LRmix
# with Help from Oyvind Bleka
#patch 3.2.1 (NFI) # 17/01/2013
# patch 4.1, Delft
LRmixTK <-function()
{
	# if(!require(tcltk)) stop("package tcltk is required")
	# if(!require(tcltk2)) stop("package tcltk2 is required")
	# if(!require(tkrplot)) stop("package tkrplot is required")
	tclRequire("Tktable")
	font0 <- tkfont.create(family="courrier",size=35,weight="bold",slant="italic")
	font1<-tkfont.create(family="times",size=14,weight="bold")#,slant="italic")
	font2<-tkfont.create(family="times",size=16,weight="bold",slant="italic")
	font3<-tkfont.create(family="times",size=12)#,slant="italic")
	font4<-tkfont.create(family="courrier",size=14)#,slant="italic")
	font5<-tkfont.create(family="courrier",size=13,weight="bold")#,slant="italic")
	font6<-tkfont.create(family="times",size=8,weight='bold')#tkframe entries labels
	font7<-tkfont.create(family="courrier",size=10,weight='bold')#tkframe entries labels
	verifAF<-function(a,b)
	{
		if(tclvalue(a)=='txt')
		freqFile<-read.table(tclvalue(b),header=TRUE,as.is=TRUE,sep='\t',na.strings='')#strings as strings, avoid converting to factors
		else{
		freqFile<-read.csv(tclvalue(b),header=TRUE,as.is=TRUE,na.strings='')#strings as strings, avoid converting to factors
		}
		return(freqFile)
	}

	#------------- Begin< formatting functions

	# function to convert tables of data, as exported from Genemapper to a list
		#------------- Begin> formatting functions
	# function to convert tables of data, as exported from Genemapper to a list
	ConvertSamp.old<-function(tab)
	{
			whichsamp<-unique(tab[,1])
			whichmark<-unique(tab$Marker)

			if(any('AMEL' %in% whichmark)) whichmark<-whichmark[-which(whichmark=='AMEL')]
			replist<-vector('list',length(whichmark))
			names(replist)<-whichmark
			for(i in 1:length(whichmark))
			{
				reptmp<-NULL
				for(k in whichsamp)
				{		
					allele<-tab[tab$Marker==whichmark[i] & tab$SampleName==k,-c(1,2)]
					tmpo<-allele[which(!is.na(allele))]
					if(length(tmpo)!=0) allele2<-as.numeric(tmpo)
					else{allele2<-NULL}
				
					list0<-list(allele2)
					names(list0)<-k
					reptmp<-c(reptmp, list0)
				}
				

				replist[[i]]<-reptmp
		}	
			
			replist
	}
	ConvertSamp<-function(tab)
	{
		# sample names (ex. suspect or replicate names)
		whichsamp<-unique(tab[,1])
		# get marker names
		whichmark<-unique(tab$Marker)
		#remove AMEL from the data
		if(any('AMEL' %in% whichmark)) whichmark<-whichmark[-which(whichmark=='AMEL')]
		# intialize 
		replist<-vector('list',length(whichmark))
		# initialize 
		repVec<-vector('list',length(whichmark))

		names(replist)<-whichmark
		for(i in 1:length(whichmark))
		{
			reptmp<-NULL
			for(k in whichsamp)
			{		
				# get the alleles for a given marker, for all samples
				allele<-tab[tab$Marker==whichmark[i] & tab$SampleName==k,-c(1,2)]# get only the alleles, so unsekect columns 1 & 2
				# only select alleles, no NASs (empty cells)
				tmpo<-allele[which(!is.na(allele))]
				if(length(tmpo)!=0) allele2<-as.numeric(tmpo)
				else{allele2<-0}
				# allele 2 is the vector of alleles at marker i and sample k
				# list0<-list(allele2)
				# names(list0)<-k
				reptmp<-c(reptmp, allele2,0)
			}
			

			replist[[i]]<-reptmp
		}	
			
			replist
	}
	
	# fromatting function for reference profiles
	ConvertSamp2<-function(tab)
	{	
		# sample names (ex. suspect or replicate names)
		whichsamp<-unique(tab[,1])
		# get marker names
		whichmark<-unique(tab$Marker)
		#remove AMEL from the data
		if(any('AMEL' %in% whichmark)) whichmark<-whichmark[-which(whichmark=='AMEL')]
		# intialize 
		replist<-vector('list',length(whichmark))
		# initialize 
		repVec<-vector('list',length(whichmark))

		names(replist)<-whichmark
		for(i in 1:length(whichmark))
		{
			reptmp<-NULL
			for(k in whichsamp)# exemple victm 1, victim 2, victim 3 etc
			{		
				# get the alleles for a given marker, for all samples
				allele<-tab[tab$Marker==whichmark[i] & tab$SampleName==k,-c(1,2)]# get only the alleles, so unsekect columns 1 & 2
				# only select alleles, no NASs (empty cells)
				tmpo<-allele[which(!is.na(allele))]
				if(length(tmpo)!=0) allele2<-as.numeric(tmpo)
				else{allele2<-0}
				# allele 2 is the vector of alleles at marker i and sample k
				# list0<-list(allele2)
				# names(list0)<-k
				reptmp<-c(reptmp, allele2)
			}
			

			replist[[i]]<-reptmp
		}	
			
			replist
	}
	#-------------End> formatting functions

	#Analyse function 
	infoStain<-NULL
	analyse <-function(loc,repl, file1,file2,file3,ext1,ext2,ext3,infoStain)
	{
		

		
		# import allele frequencies and tranform into tabfreq objcet

		
		# verify that the user explored the profiles and selected the loci 
		if(setequal(tclvalue(loc),''))
		{
			stop(tkmessageBox(message="First select the loci",icon="error",type="ok"))
			
		}
		else
		{
		loc0<- strsplit(tclvalue(loc),' ')[[1]]

		}
		# check replicates
		if(setequal(tclvalue(repl),''))
		{
			stop(tkmessageBox(message="First select the replicates",icon="error",type="ok"))
			
		}
		else
		{
			repl0<- strsplit(tclvalue(repl),' ')[[1]]

		}
		#check wether the user uploaded the files for the sample and the refrenecs profiles
		veriFile<-function(filename,ext,error=TRUE,txt='')
		{
			if(tclvalue(filename)=='')
			{	
				if(error){					
				stop(tkmessageBox(message=paste("First load the",txt,sep="","profile"),icon="error",type="ok"))
				}
				else {return(0)}
				# else
					# tkmessageBox(message="First load the reference profile",icon="info",type="ok")
			}
			else
			{
				if(tclvalue(ext)=='txt')
				{
					tab<-read.table(tclvalue(filename),header=TRUE,as.is=TRUE,sep='\t',na.strings='')
					if(any('AMEL' %in% tab$Marker)) tab<-tab[-which(tab$Marker=='AMEL'),]

					#strings as strings, avoid converting to factors
				}
				else
				{
					tab<-read.csv(tclvalue(filename),header=TRUE,as.is=TRUE,na.strings='')#strings as strings, avoid converting to factors
					if(any('AMEL' %in% tab$Marker)) tab<-tab[-which(tab$Marker=='AMEL'),]

				}# rm(file
				return(tab)
			}
				
		}

		#suspect samplename to be inserted into the listbox of the contributors under Hp, or eventually the non contributors under Hd
		
		verifFormat<-function(tab)
		{
			if(infoStain[1]!='SampleName')
			{
				tkmessageBox(message="Format error, pleae check your file",icon="error",type="ok")
			}
			if(infoStain[2]!='Marker')
			{
				tkmessageBox(message="Format error, please check your file",icon="error",type="ok")
			}
		}
			
		
		main <- tktoplevel()
		main2 <- tkframe(main)
		# tkgrid(main2)
		frame.tit<-tkframe(main2,relief="groove", borderwidth=4)
		tkgrid(frame.tit)
		tkgrid(tklabel(frame.tit,text="LR calculation", font="courrier 22", 	foreground="darkblue"),sticky='n')
		tkwm.title(main,paste('Analyse the profiles, forensim v.',versionNum,sep=''))

	#readFile button
		firstFrame<-tkframe(main, relief="groove", borderwidth=4)
		tkgrid(tklabel(firstFrame,text='Parameters',font='courrier 14', fg='blue'))
		
		ncFrame<-tkframe(firstFrame,relief="flat", borderwidth=4)
		secondFrame<-tkframe(main,relief="groove", borderwidth=4)

		#---first frame 
		probFrame <- tkframe(firstFrame, relief="flat", borderwidth=4)
		repFrame<- tkframe(firstFrame, relief="groove", borderwidth=4)
		#--- second frame
		# extraFrame<-tkframe(secondFrame,relief='groove', borderwidth=4)
		hypoFrame<-tkframe(secondFrame, relief="groove", borderwidth=1)
		
		#third frame for uploading allele frequencies
		thirdFrame<-tkframe(main,relief="groove", borderwidth=4)
		#bottom frame
		bottomFrame<-tkframe(main,relief='flat',borderwidth=4)
		#-------------------- number of contributors ----------------- #
		tkgrid(tklabel(ncFrame,text="   Unknown contributors     ",font='courrier 10 bold'),sticky="w")
		ncHd<-tclVar(1)
		ncHp<-tclVar(0)
		ncHd.entry<-tkentry(ncFrame,textvariable=ncHd,width=4,highlightthickness=1,relief="solid",justify="center")
		ncHp.entry<-tkentry(ncFrame,textvariable=ncHp,width=4,highlightthickness=1,relief="solid",justify="center")
		tkgrid(tklabel(ncFrame,text=" Under Hp   ",font='courrier 10'),ncHp.entry,sticky="w")
		tkgrid(tklabel(ncFrame,text=" Under Hd   ",font='courrier 10'),ncHd.entry,sticky="w")
		#-------------------Hypotheses
		
		hypoFrame1<-tkframe(secondFrame)
		titre<-tkframe(hypoFrame1)

		tkgrid(titre,pady=10)
		tkgrid(tklabel(titre,text="Hypotheses",font='courrier 14', fg='blue'))
		
		
		
		
		#--- first frame
		hypoFrame1.tit<-tklabel(hypoFrame1,text='Contributors under Hp',font='courrier 10 bold')
		hypoFrame11<-tkframe(secondFrame,relief='groove')
		# hypoFrame11.tit<-tklabel(secondFrame,text='Non-contributors under Hp')


		hypoFrame2<-tkframe(secondFrame)
		hypoFrame2.tit<-tklabel(hypoFrame2,text='Contributors under Hd',font='courrier 10 bold')
		# hypoFrame22.tit<-tklabel(secondFrame,text='Non-contributors under Hd')

		#-------------------------PROSECUTION-------------------------------------------------------------#
		#suspects variables
			
		csp<-veriFile(file1,ext1,error=TRUE,'crime scene profile')
		suspect<-veriFile(file2,ext2,error=TRUE,'suspect profile')
		suspectID<-unique(suspect$SampleName)
		# contributors under Hp

		#victim and other profiles not necesarily loaded, depends on the case
		victim<-veriFile(file3,ext3,error=FALSE,'')
		if(is.data.frame(victim))
		{
			vicID<-unique(victim$SampleName)
			for(v in 1:length(vicID))
			{
				#choice of contributors under Hp/Hd: victim profiles (if applicable)
				assign(paste('vicHp',v,sep=''), tclVar(1))
				assign(paste('vicHd',v,sep=''), tclVar(1))
				#dorpout for victims, if applicable

			}
		}
		tkgrid(hypoFrame1.tit)
		for(i in 1:length(suspectID))
		{
			# Tk entry for suspect(s) under Hp
			tkgrid(tklabel(hypoFrame1, text=suspectID[i],font='courrier 10'),sticky='w')
			# Tk entry for suspect(s) under Hd
		}
		
		
		tkgrid(hypoFrame1,padx=10)
		tkgrid(hypoFrame2.tit,pady=2)
		# for(i in 1:length(suspectID))
		# {
			# tkgrid(tkcheckbutton(hypoFrame2, text=suspectID[i], variable=get(paste('suspectHd',i,sep='')),font='courrier 8'),sticky='w')
		# }
		tkgrid(hypoFrame2)
		
		if(is.data.frame(victim)){
		for(i in 1:length(vicID))
		{
				tkgrid(tkcheckbutton(hypoFrame1, text=vicID[i], variable=get(paste('vicHp',i,sep='')),font='courrier 10'),sticky='w')
		
				tkgrid(tkcheckbutton(hypoFrame2, text=vicID[i], variable=get(paste('vicHd',i,sep='')),font='courrier 10'),sticky='w')
						
		}
		}	
		
	#--------- Performance plots ------------------#
      #Simulation subframe:
        thirdFrame <- tkframe(main, relief = "groove", borderwidth = 4)
        simFrame1 <- tkframe(thirdFrame)
        titre <- tkframe(simFrame1)
        tkgrid(titre, pady = 10)
        tkgrid(tklabel(titre, text = "Performance test", font = "courrier 14", fg = "blue"))
        simFrame11 <- tkframe(thirdFrame, relief = "groove")
        tkgrid(simFrame1, padx = 10)
        #init. checkbuttons: suspects need not to be checked, they have to be given
        for (v in 1:length(suspectID))
		{
			assign(paste("simSus", v, sep = ""), tclVar(0))
        }
        #List suspects (but user has to give at least one suspect) that may be subsitited in the Tippet calculations
		
		tkgrid(tklabel(simFrame1,text="Choose suspect",font='courrier 10 bold'),sticky="w")
		# tkgrid(tklabel(simFrame1,text="(suspect to replace by random man)",font='courrier 7 bold'),sticky="w")
      	for (i in 1:length(suspectID))
		{
			
			tkgrid(tklabel(simFrame1,text=suspectID[i],font = "courrier 10",fg='darkred'),     tkcheckbutton(simFrame1, text = "",variable = get(paste("simSus", i, sep = "")),font = "courrier 10"), sticky = "w")   
		}  
        # }
		
		
        simval <- tclVar(100)
        simval.entry <- tkentry(simFrame1, textvariable = simval, width = 4, highlightthickness = 1, relief = "solid", justify = "center")
		tkgrid(tklabel(simFrame1, text = "number of iterations", font = "courrier 10"), simval.entry, sticky = "w")
		# here put function simulation
		sim.loc<-function()
		{
			# check that only one suspect is substituted at the time
			# simSus were initilized to ero, thei values change if user check a box
			nBoxChecked <- 0
			for(v in 1:length(suspectID))
			{   
			  nBoxChecked <- nBoxChecked + as.numeric(tclvalue(get(paste("simSus",v,sep=""))))
			}
			if (nBoxChecked!=1)	# warning if more than one suspect is selected
			{
				stop(tkmessageBox(message="You must select exactly one suspect!",icon="warning"))
			}
			
			#-- get the ID of the suspect that we want to substitute
			id.tmp<-rep(0,length(suspectID))
			for(j in 1:length(suspectID))
			{
					id.tmp[j]<-as.numeric(tclvalue(get(paste('simSus',j,sep=''))))
			}
			# at this stage, suspect is a dataframe, with possibly several suspest, so nee to select the data for the relevant suspects, and remove the one that is going to be replaced by a random man in the Tippet analysis
			
			#--- get the number of simulations
			M <- as.numeric(tclvalue(simval)) 
			# create the table of random men
			data1<-tabfreq(verifAF(extens22,filePath22))
			data0<-data1$tab
			#extens22 and filePath22 are updated by importfreq, these are local variables to analyse() function
			
			#----------VICTIM PROFILES  get users choice: victim
			#default values, eventually modified by what follows
			Tp.vic<-Td.vic<-0
			if(is.data.frame(victim))#if a file is given for the victim
			{
				
				selecVicHp<-rep(0,length(vicID))#vic ID is of length the number of victims provides in the sample
				selecVicHd<-rep(0,length(vicID))

				for(j in 1:length(vicID))
				{
					selecVicHp[j]<-as.numeric(tclvalue(get(paste('vicHp',j,sep=''))))
					selecVicHd[j]<-as.numeric(tclvalue(get(paste('vicHd',j,sep=''))))
				}
				if(length(which(selecVicHp==1))!=0)#checks if victim ID is selected if so then add under Hp
				{
					Tp.vic<-victim[victim$SampleName %in% vicID[which(selecVicHp==1)],]
				}
				# else
				# {
					# Tp.vic<-0
				# }
				if(length(which(selecVicHd==1))!=0)# same treatement under Hd
				{
					Td.vic<-victim[victim$SampleName %in% vicID[which(selecVicHd==1)],]#if user selects victim, then add
					
				}
				# else#otherwise no contribution from a victim
				# {
					# Td.vic<-0
				# }
			}
			# else{victim<-0}
			#------------EXTRA PROFILES
			Vp<-0#default values, modified if users upload extra profiles
			
			#------------ SUSPECT PROFILES
			#the suspect become known non-contributor under Hd, under Hp, suspects that are still evaluated
			Tp.sus<-suspect[suspect$SampleName %in% suspectID[which(id.tmp==0)],]
			# other potential suspects who are contributors under Hp, if the matrix is empty (nrow=0), then the ConvertSamp will yield a NULL, there is always one contributor under Hp, if unique suspect + no victims, that will be the random man, so tp!=NULL
			# the replaced suspect becomes contributor under Hd
			Vd.sus<-suspect[suspect$SampleName %in% suspectID[which(id.tmp!=0)],]

			#-------------------------------------------
			# assign the contributors under Hp: Tp
			if(is.data.frame(Tp.vic)) #if victim is providedm add to suspect
			{
				Tp<-rbind.data.frame(Tp.sus, Tp.vic)
			}
			else{#not really needed
			Tp<-Tp.sus}#if not do nothing
			
			

			#-------------------------------------------
			cspFinal<-ConvertSamp(csp[csp$Marker %in% loc0 & csp$SampleName %in% repl0, ])
			TpFinal<-NULL#and changed if applicable
			if(is.data.frame(Tp)){
			if(nrow(Tp)!=0) {
			TpFinal<-ConvertSamp2(Tp[Tp$Marker %in% loc0, ])}
			}
			
			VdFinal<-ConvertSamp2(Vd.sus[Vd.sus$Marker %in% loc0, ])
		

			# ------ select subset with the markers in loc0, difficulty for the victim is that it is no necessarily --------- #
			if(is.data.frame(Td.vic) )#if a victim is given
			{
				if(nrow(Td.vic)!=0){
				TdFinal<-ConvertSamp2(Td.vic[Td.vic$Marker %in% loc0, ])}
				# if(is.character(TdFinal)) TdFinal <-matrix(TdFinal,ncol=length(loc0))
			}
			else{
			TdFinal<-Td.vic}#else its null, for code clarity only
			
			
			xp<-as.numeric(tclvalue(ncHp))
			xd<-as.numeric(tclvalue(ncHd))
			theta0<-as.numeric(tclvalue(theta))
			# init. vector for the storage of LRs for each simulated profile
			lr0<-rep(1,M)
			print('========= Performance plot ============')
      	   	drop0<-as.numeric(tclvalue(prD)) #get dropout-prob from GUI
			for(jj in loc0) {
                   popfreq = data0[[jj]] #get allele-frequencies
			 rep0<-cspFinal[[jj]] #get evidence

			 #denumerator Pr(E|Hd)
  	 		 if(is.list(TdFinal )){ tmpTd<-unlist(TdFinal[jj])}
			 else{ tmpTd<-0}
                   Vd<-unlist(VdFinal[[jj]])
        	       hd_val = likEvid(Repliste=(rep0),T=tmpTd,V=Vd,x=xd,theta=theta0,prDHet=rep(drop0,5),prDHom=rep(drop0^2,5),prC=as.numeric(tclvalue(prC)), freq=data0[[jj]]) # V does not contribute to replicate probability)
 
			 #numerator Pr(E|Hp). G is all possible genotypes:
                   hp_val <- rep(NA,M) #vector for hp_values of random man
                   G <- t(as.matrix(expand.grid(rep(list(as.numeric(names(popfreq)),as.numeric(names(popfreq)) )))))
                   keep <- G[2,]>=G[1,] #unique genotypes
                   G <- G[,keep]  #store genotypes
                   tmpP <- t(as.matrix(expand.grid(rep(list(as.numeric(popfreq),as.numeric(popfreq) )))))
                   Gprob <- exp(colSums(log(tmpP[,keep]))) #get genotype probs
                   ishet <- G[1,]!=G[2,]
  			 Gprob[ishet] <- 2*Gprob[ishet] #multiply with two to get heterozygote probs
                   Gsampled <- sample(1:length(Gprob),size=M,prob=Gprob,replace=TRUE)
                   unGsampled <- unique(Gsampled) #get unique sampled
                   for(uu in 1:length(unGsampled)) {
                    randoman <- as.numeric(G[,unGsampled[uu]]) #get allele-frequence of unique random man
        		  if(!is.null(TpFinal)){tp<-c(unlist(TpFinal[[jj]]), randoman)}
			  else{tp<-randoman}
                    hp_val[ Gsampled==unGsampled[uu] ] <-likEvid(Repliste=(rep0),T=tp,V=0,x=xp,theta=theta0,prDHet=rep(drop0,5),prDHom=rep(drop0^2,5),prC=as.numeric(tclvalue(prC)),freq=data0[[jj]]) #calculate for unique random man
                   } #end for each unique calculation
		 	 lr0 <- lr0*hp_val/hd_val #multiply with LR-value for current locus for all randoms
     			 cat(paste(jj, 'completed','\n'))
			} #end 'jj' for each locus
			
			print('===================================')
		#--- sensitivity analysis
			distriLR<-log(lr0, 10) 
			# # plot the empirical cumultaive distribution of the log10  LR, using function ecdf 
			# plot(ecdf(log(distriLR,10)),xlab='log10 LR')
			qvals <- c(0.01,0.05,0.5,0.95,0.99)
			quantiles<-quantile(distriLR,qvals)
			minmax <-range(distriLR)
			tab<-cbind(c("min",as.character(qvals),"max"),round(c(minmax[1],quantiles,minmax[2]),4))
			colnames(tab) <- c("quantile","value")
			write.table(tab,paste("LRdistrQuantiles",M,".txt",sep=""),row.names=FALSE)
			tkmessageBox(message=paste('Performance test percentiles saved to',
			paste("LRdistrQuantiles",M,".txt",sep="")), icon='info',type='ok')	
			
			if('Inf' %in% range(distriLR) | NaN %in% range(distriLR) )
			{
				stop(tkmessageBox(message="infinite LR values, please change the model parameters",icon="error",type="ok"))
			}
			Myhscale <- 1   # Horizontal scaling
			Myvscale <- 1
			dd <- tktoplevel()
			frameC<-tkframe(dd)
			tkwm.title(dd,paste('Performance plot, forensim v.',versionNum,sep=''))
			# tkwm.title(res,paste('LRmix: Results, forensim v.',versionNum,sep=''))

			Dplot.loc<-function()
			{
				# tkconfigure(dd,cursor="watch")
				params <- par(bg="white")
				plot(ecdf(distriLR),xlab='log10 LR',main='Empirical distribution function')
				grid()
				par(params)
			}
			img <- tkrplot(frameC,fun= Dplot.loc,hscale=Myhscale,vscale=Myvscale)
			CopyToClip <- function()
			{
				tkrreplot(img)
			}
			copy.but <- tkbutton(frameC,text="Copy to Clipboard",font="courrier 10",fg="darkblue",command=CopyToClip)
			# LRtab2<-LRres
			# excel.but2<-tkbutton(frameC, text="Export results",fg="blue", font="courrier 10",command=function() exportFile(LRtab2))#,command=function() openFile())
			tkgrid(img)
			tkgrid(copy.but)#,excel.but2,rowspan=10,sticky='ew')
			# tkgrid(plot.but, excel.but,rowspan=10,sticky='ew')
			tkpack(frameC)
		}
		
		sim.butt <- tkbutton(simFrame1, text = "OK", font = "courrier 10",fg = "black", command =sim.loc)
		tkgrid(sim.butt, columnspan = 20,pady=10)

		#---------prob of dropout and dropoin
		prD<-tclVar(0.10)
		prC<-tclVar(0.05)
		theta<-tclVar(0)
		#distribution
		prD.entry<-tkentry(probFrame,textvariable=prD,width=4,highlightthickness=1,relief="solid",justify="center")
		prC.entry<-tkentry(probFrame,textvariable=prC,width=4,highlightthickness=1,relief="solid",justify="center")
		theta.entry<-tkentry(probFrame,textvariable=theta,width=4,highlightthickness=1,relief="solid",justify="center")
		# merde=tkbutton(main)
		tkgrid(tklabel(probFrame,text="   Pr(D), Pr(C), theta     ",font='courrier 10 bold'))#,sticky="w")
		tkgrid(tklabel(probFrame,text="   Probability of Dropout   Pr(D)     ",font='courrier 10'),prD.entry,sticky="w")
		tkgrid(tklabel(probFrame,text="   Probability of Contamination   Pr(C)     ",font='courrier 10'),prC.entry,sticky="w")
		tkgrid(tklabel(probFrame,text="   Theta Correction (Fst)   ",font='courrier 10'),theta.entry,sticky="w")

			
			
		# function to import allele frequencies from JFS-like file format
		# tclvalue(ext) <- strsplit(foo3[length(foo3)], "\\.")[[1]][2] 
		
		
		#-- function to read the AF
		
		#---------------------------------------------------------------#
		
		#-------- ANALYSIS OF RESULTS ---------- #
		
		
			

		
		analyse.loc<-function()
		{
			data0<-tabfreq(verifAF(extens22,filePath22))$tab
			#extens22 and filePath22 are updated by importfreq, these are local variables to analyse() function
			
			#----------VICTIM PROFILES  get users choice: victim
			#default values, eventually modified by what follows
			Tp.vic<-Td.vic<-0
			if(is.data.frame(victim))#if a file is given for the victim
			{
				
				selecVicHp<-rep(0,length(vicID))#vic ID is of length the number of victims provides in the sample
				selecVicHd<-rep(0,length(vicID))

				for(j in 1:length(vicID))
				{
					selecVicHp[j]<-as.numeric(tclvalue(get(paste('vicHp',j,sep=''))))
					selecVicHd[j]<-as.numeric(tclvalue(get(paste('vicHd',j,sep=''))))
				}
				if(length(which(selecVicHp==1))!=0)#checks if victim ID is selected if so then add under Hp
				{
					Tp.vic<-victim[victim$SampleName %in% vicID[which(selecVicHp==1)],]
				}
				# else
				# {
					# Tp.vic<-0
				# }
						
				if(length(which(selecVicHd==1))!=0)# same treatement under Hd
				{
					Td.vic<-victim[victim$SampleName %in% vicID[which(selecVicHd==1)],]#if user selects victim, then add
					
				}
				# else#otherwise no contribution from a victim
				# {
					# Td.vic<-0
				# }
			}
			# else{victim<-0}
			#------------EXTRA PROFILES
			Vp<-0#default values, modified if users upload extra profiles
			
			#------------ SUSPECT PROFILES
		
			Tp.sus<-suspect#the suspect become known non-contributor under Hd
			Vd.sus<-suspect
			#-------------------------------------------
			# assign the contributors under Hp: Tp
			if(is.data.frame(Tp.vic)) #if victim is providedm add to suspect
			{
			
				indVicHp<-unique(Tp.vic$SampleName)

				Tp<-rbind.data.frame(Tp.sus, Tp.vic)
			}
			else{#not really needed
			Tp<-Tp.sus
			indVicHp<-NULL
			}#if not do nothing
		
			
			#-------------------------------------------
			
			
			
			cspFinal<-ConvertSamp(csp[csp$Marker %in% loc0 & csp$SampleName %in% repl0, ])
			cspFinal2<-ConvertSamp.old(csp[csp$Marker %in% loc0 & csp$SampleName %in% repl0, ])
			nbAll<-rep(0,length(repl0))
			for(mm in 1:length(nbAll))
			{
				repmm<-repl0[mm]
				nbAll[mm]<-sum(sapply(cspFinal2,function(k) length(unique(k[[repmm]]))))
			}
			# nbAll<-sum(sapply(cspFinal,function(k) length(unique(k))))
			nbAll.mean<-round(mean(nbAll))
			# locosimuHp<-function(d=vecD,contri=TpFinal,x=xp,freq=data0,nrep=100,nbAll=nbAll.mean)

			#contributors under Hp
			TpFinal<-ConvertSamp2(Tp[Tp$Marker %in% loc0, ])
			TpFinal2<-ConvertSamp.old(Tp[Tp$Marker %in% loc0, ])


			VdFinal<-ConvertSamp2(Vd.sus[Vd.sus$Marker %in% loc0, ])
			VdFinal2<-ConvertSamp.old(Vd.sus[Vd.sus$Marker %in% loc0, ])

			# ------ select subset with the markers in loc0, difficulty for the victim is that it is no necessarily --------- #
			if(is.data.frame(Td.vic)   )#if a victim is given
			{
				if(nrow(Td.vic)!=0){
					indVicHd<-unique(Td.vic$SampleName)
					TdFinal<-ConvertSamp2(Td.vic[Td.vic$Marker %in% loc0, ])
					TdFinal2<-ConvertSamp.old(Td.vic[Td.vic$Marker %in% loc0, ])
					# if(is.character(TdFinal)) TdFinal <-matrix(TdFinal,ncol=length(loc0))
				}
			}
			else{
			indVicHd<-NULL
			TdFinal<-Td.vic
			TdFinal2<-Td.vic
			}#else its null, for code clarity only
			tdfinal<-TdFinal
			locus.hp<-rep(0,length(loc0))
			locus.hd<-rep(0,length(loc0))

			names(locus.hd)<-loc0
			names(locus.hp)<-loc0
			xp<-as.numeric(tclvalue(ncHp))
			xd<-as.numeric(tclvalue(ncHd))
			theta0<-as.numeric(tclvalue(theta))

			for(jj in loc0)
			{
				print(jj)
				rep0<-cspFinal[[jj]]
				if(is.list(TdFinal )){ tmpTd<-unlist(TdFinal[jj])}
				else{ tmpTd<-0}
				#numerator Pr(E|Hp)
				drop0<-as.numeric(tclvalue(prD))
				tp<-unlist(TpFinal[[jj]])
				locus.hp[jj]<-likEvid(Repliste=(rep0),T=tp,V=0,x=xp,theta=theta0,prDHet=rep(drop0,length(tp)/2 + xp),prDHom=rep(drop0^2,length(tp)/2 + xp),prC=as.numeric(tclvalue(prC)),freq=data0[[jj]])
				locus.hd[jj]<-likEvid(Repliste=(rep0),T=tmpTd,V=unlist(VdFinal[[jj]]),x=xd,theta=theta0,prDHet=rep(drop0,length(tmpTd)/2 + xd),prDHom=rep(drop0^2,length(tmpTd)/2 + xd),prC=as.numeric(tclvalue(prC)), freq=data0[[jj]])# V does not contribute to replicate probability
			}
			# create a table of the results that will be exported in an Excel file
LRtab<-cbind.data.frame('Locus'=c(loc0,'product'),
'Pr(E|Hp)'=signif(c(locus.hp,prod(locus.hp)),4),
'Pr(E|Hd)'=signif(c(locus.hd,prod(locus.hd)),4), 
'LR'=signif(c(locus.hp/locus.hd,prod(locus.hp/locus.hd)),4),
'log(LR)'=signif(c(log10(locus.hp/locus.hd),log10(prod(locus.hp/locus.hd))),4))


	
			#display the results
			res<-tktoplevel()
			tkwm.title(res,paste('LRmix: Results, forensim v.',versionNum,sep=''))

			f1<-tkframe(res)
			tkgrid(tklabel(f1,text="   Results     ",font='courrier 14', foreground="blue"),sticky="w")
			
			#
			array1 <- tclArray()
			array2 <- tclArray()

			myRarray<-c("LR per Locus",loc0,"LR",signif(locus.hp/locus.hd,4))
			myRarray2<-c("Overall LR",signif(prod(locus.hp/locus.hd),4))

			dim(myRarray)<-c(length(loc0)+1,2)
			dim(myRarray2)<-c(2,1)

			for(i in 0:(length(loc0))){
			array1[[i,0]] <- myRarray[i+1,1]
			array1[[i,1]] <- myRarray[i+1,2]

			}
			#table2
			array2[[0,0]] <- myRarray2[1,1]
			array2[[1,0]] <- myRarray2[2,1]

			#if no known non-contributors under Hp     
			table1 <- tkwidget(f1,"table",variable=array1,rows=(length(loc0))+1 ,cols=2,titlerows=1,titlecols=0, colwidth=25)
			
			table2 <- tkwidget(f1,"table",variable=array2,rows=2,cols=1,titlerows=1,titlecols=0, colwidth=25)
			
			
			#xscr <-tkscrollbar(tt,orient="horizontal", command=function(...)tkxview(table1,...))
			# yscr <- tkscrollbar(f1,repeatinterval=5, command=function(...)tkyview(table1,...))
			# tkgrid(table1, table2,yscr,sticky="new", padx=20,pady=18)
			tkgrid(table1, table2,sticky="new", padx=20,pady=18)
			# tkgrid.configure(yscr,sticky="nsw")
			#tkconfigure(table1,variable=array1,background="white",selectmode="extended")
			tkconfigure(table1,variable=array1,background="white",selectmode="extended",rowseparator="\"\n\"",colseparator="\"\t\"")
			tkconfigure(table2,variable=array2,background="white",selectmode="extended",rowseparator="\"\n\"",colseparator="\"\t\"")
			
			# button to display the plot of the LR vs the PrD
			plot.but<-tkbutton(f1, text="Plot LR vs PrD",fg="blue", font="courrier 10",command=function() Dplot())#,command=function() openFile())
			#button for the fisrt phase of the analysis: point estimate
			excel.but<-tkbutton(f1, text="Export results",fg="blue", font="courrier 10",command=function() exportFile(LRtab))#,command=function() openFile())
			#---export LR
			# Export the results from the LR calculations, and let the user decide the name
			exportFile<-function(tmp)
			{
				ff<-tktoplevel()
				Fframe<- tkframe(ff, relief="groove")
				tkwm.title(ff,"Filenames")
				tkgrid(tklabel(Fframe, text="===== Enter filename ====",font='courrier 12',foreground="darkblue"), columnspan=9)
				filtervar<- tclVar('LRs.txt')
				filtervar.entry <- tkentry(Fframe, textvariable=filtervar, width=12)
				saveF.butt<-tkbutton(Fframe, text="Enter",fg="darkblue",font='courrier 8',command=				  function() functionMAJ() )
				
				
				
				functionMAJ<-function(){
				filen<-tclvalue(filtervar)
				# write.table('____________________________________________________',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('======= Log file: LRs per locus ==========',file=filen,append=FALSE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('============= Input files ================',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				aaaa1<-paste('SAMPLE FILE:',(tclvalue(filePath)),sep='\n')
				aaaa2<-paste('SUSPECT FILE:',(tclvalue(filePath2)),sep='\n')
				aaaa3<-paste('VICTIM FILE:',(tclvalue(filePath3)),sep='\n')
				aaaa4<-paste('SELECTED LOCI:',tclvalue(locus), sep='\n')
				aaaa5<-paste('SELECTED REPLICATES:',tclvalue(repl), sep='\n')
				# options(warn=-1)
				# write.table(x=header,sep=',',file=fileName)
				write.table(aaaa1,file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(aaaa2,file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(aaaa3,file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(aaaa4,file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(aaaa5,file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				# write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)

				write.table('============ User parameters =============',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(paste('Drop-out value:',tclvalue(prD),sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(paste('Drop-in value:',tclvalue(prC),sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				# theta
				write.table(paste('Theta value:',tclvalue(theta),sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				# xp under Hp
				write.table(paste('Unknowns under Hp:',xp,sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				# xd under Hd
				write.table(paste('Unknowns under Hd:',xd,sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				

				#-----------------------------------------------------#

				
write.table('=============== Hypotheses ===============',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
# print('ici');print(suspectID)
if(!is.null(indVicHp)){
write.table(paste('Under Hp:', paste(suspectID,collapse=' + '), '+',xp,'unknown(s)','+',paste(indVicHp,collapse=' + '),sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)}
else{write.table(paste('Under Hp:', paste(suspectID,collapse=' + '), '+',xp,'unknown(s)',sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)}

if(!is.null(indVicHd)){
write.table(paste('Under Hd:',xd,'unknown(s)','+',paste(indVicHd,collapse=' + '),sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)}
else{write.table(paste('Under Hd:',xd,'unknown(s)',sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)}



write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				
# write.table('-----------------------------------------',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE)
				write.table('=========== Log10(LR) vs. Pr(D) ===========',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(tmp,file=filen,row.names=FALSE,append=TRUE,quote=FALSE)
				# write.table('____________________________________________________',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)

				}
				
				tkgrid(filtervar.entry, saveF.butt)		
				tkpack(Fframe,padx=12,pady=18,side="left")
			}
			exportFile2<-function(tmp,r0,r1)
			{
				options(warn=-1)
				ff<-tktoplevel()
				Fframe<- tkframe(ff, relief="groove")
				tkwm.title(ff,"Filenames")
				tkgrid(tklabel(Fframe, text="===== Enter filename ====",font='courrier 12',foreground="darkblue"), columnspan=9)
				filtervar<- tclVar('sensitivity.txt')
				filtervar.entry <- tkentry(Fframe, textvariable=filtervar, width=12)
				saveF.butt<-tkbutton(Fframe, text="Enter",fg="darkblue",font='courrier 8',command=function() functionMAJ() )
				# function called when exporting, function of tclvalue (updated for each run)
				functionMAJ<-function(){
				# get the file chosen by the user
				filen<-tclvalue(filtervar)
			    # files
				# write.table('____________________________________________________',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('===== Log file: Sensitivity analysis =====',file=filen,append=FALSE,row.names=FALSE,col.names=FALSE,quote=FALSE)

				write.table('============= Input files ================',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				
				aaaa1<-paste('SAMPLE FILE:',(tclvalue(filePath)),sep='\n')
				aaaa2<-paste('SUSPECT FILE:',(tclvalue(filePath2)),sep='\n')
				aaaa3<-paste('VICTIM FILE:',(tclvalue(filePath3)),sep='\n')
				aaaa4<-paste('SELECTED LOCI:',tclvalue(locus), sep='\n')
				aaaa5<-paste('SELECTED REPLICATES:',tclvalue(repl), sep='\n')
				# options(warn=-1)
				# write.table(x=header,sep=',',file=fileName)
				write.table(aaaa1,file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(aaaa2,file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(aaaa3,file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(aaaa4,file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(aaaa5,file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				# write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)

				write.table('============ User parameters =============',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(paste('Drop-in value:',tclvalue(prC),sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				# theta
				write.table(paste('Theta value:',tclvalue(theta),sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				# xp under Hp
				write.table(paste('Unknowns under Hp:',xp,sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				# xd under Hd
				write.table(paste('Unknowns under Hd:',xd,sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				#---- Hypotheses
write.table('=============== Hypotheses ===============',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)


if(!is.null(indVicHp)){
write.table(paste('Under Hp:', paste(suspectID,collapse=' + '), '+',xp,'unknown(s)','+',paste(indVicHp,collapse=' + '),sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)}
else{write.table(paste('Under Hp:', paste(suspectID,collapse=' + '), '+',xp,'unknown(s)',sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)}

if(!is.null(indVicHd)){
write.table(paste('Under Hd:',xd,'unknown(s)','+',paste(indVicHd,collapse=' + '),sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)}
else{write.table(paste('Under Hd:',xd,'unknown(s)',sep=' '),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)}



write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
	
	
				#------------------- drop-out ranges
				#-----------------------------------------------------#
				write.table('======== Drop-out ranges: under Hp ========',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(paste('5% percentile',r0[1],sep=" "),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(paste('95% percentile',r0[2],sep=" "),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)

				write.table('======== Drop-out ranges: under Hd ========',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(paste('5% percentile',r1[1],sep=" "),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(paste('95% percentile',r1[2],sep=" "),file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table('\n',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				# write.table('____________________________________________________',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				#-----------------------------------------------------#
				# LRs	
				# write.table('-----------------------------------------',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE)
				write.table('==== Likelihoods & likelihood ratios =====',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
				write.table(tmp,file=filen,row.names=FALSE,append=TRUE,quote=FALSE)
				# write.table('-----------------------------------------',file=filen,append=TRUE,row.names=FALSE,col.names=FALSE)
	

				}
				
				tkgrid(filtervar.entry, saveF.butt)		
				tkpack(Fframe,padx=12,pady=18,side="left")
			}
			# function which displays hep for sensitivity analysis
			infoSP<-function()
			{
				tkmessageBox(message="The plot displays the sensitivity analysis of the LRs to variations of the drop-out probability.\nThe same drop-out probability is applied to all contributors under Hp and under Hd.\nThe red arrow report the ranges of the probabilities of drop-out, estimated via the\n Monte-Carlo simulation procedure described in Gill et al (Forensic Sci. Int. 2007). ",icon="info",type="ok")

			}
			#--- function which draws the LR vs PrD
			#--- new from 3.1: ranges of drop-put are also reported
			Dplot<-function()
			{
				
				vecD<-signif(seq(0.01,0.99,length=20),2)
				# vecD<-seq(0.01,1,by=0.02)
				# Hinda, Delft May 28th 2012
				# simulates the number of alleles x in the whole profile, for a range of PrD values
				#-- under Hd
				locosimu<-function(d=vecD,prC=cc,contri=TdFinal2,x=xd,freq=data0,nrep=100,nb=nbAll.mean,locnames=loc0)
				{

				# there may be no known contributors under Hd
				if(!setequal(contri,0))
				{
					locosimu.loc<-function(d,contri,x,freq)
					{
						#-- get the names of TpFinal to get the contributors
						index<-names(contri[[1]])
						for(i in 1:length(index)){assign(paste("contri",i,sep=''),
						sapply(contri,function(ee) ee[[index[i]]]))}
						#-------------- if x=0, no uknown contributors
						if(x==0)
						{
							locvec<-vector('list',length=length(locnames))
							names(locvec)<-locnames
							locvec2<-locvec
							for(l in (locnames))
							{
								z<-NULL
								for(rr in 1:length(d))# rth drop-out value
								{
									# sample from runif, in order to select the alleles that will drop-out
									# we don't know in advance the user choice and the number of contributors, so we need to make the code flexible
									#loop to go through all potential contributors under H
									z0<-NULL
									for(i in 1:length(index))#ith contributor
									{
										if(runif(1) <= d[rr]){
											a1<-NULL #remove allele from data
										}
										else{
											a1<-get(paste('contri',i,sep=''))[1,l]}
										
										if(runif(1) <= d[rr]){ a2<-NULL}#remove allele from data
										else{a2<-get(paste('contri',i,sep=''))[2,l]}
										z0<-c(z0,a1,a2)#list of alleles
										# z1<-c(z1,length(unique(c(a1,a2))))
									}
									#z1 contains the alleles of the all the contrubutors at locus l
									# add drop-in
									cc<-NULL
									if(prC!=0)
									{
										cc<-NULL
										if(runif(1)<=prC){ cc<-sample(names(freq[[l]]),1,prob=freq[[l]])}
									}
									z<-c(z,length(unique(c(z0,cc))))
									# now treat unknown contr
								}
								locvec[[l]]<-z#sum the #of alleles among the unknown and the 

							}
							tmp0<-apply(sapply(locvec,rbind),1,sum)
							return(tmp0)
						}
						
						# if there are unknown contributors
						#----------------------------------------------------#
						if(x!=0)
						{
							locvec<-vector('list',length=length(locnames))
							names(locvec)<-locnames
							locvec2=locvec
							for(l in (locnames))
							{
								zz<-NULL
								for(rr in 1:length(d))
								{
									# rth drop-out value
									# sample from runif, in order to select the alleles that will drop-out
									# we don't know in advance the user choice and the number of contributors, so we need to make the code flexible
									#loop to go trhourgh 
									zz1<-NULL
									for(i in 1:length(index))#ith contributor
									{
										#first allele of contri i
										# second allele of contri i
										A<-NULL
										B<-NULL
										if(runif(1) >= d[rr])
										{
											A<-get(paste('contri',i,sep=''))[1,l] #remove allele from data
										}
										if(runif(1) >= d[rr])
										{
											B<-get(paste('contri',i,sep=''))[2,l] #remove allele from data
										}
										
										zz1<-c(zz1,(c(A,B)))
									}
									# now treat unknown contr
									zz2<-NULL
									for(a in 1:x)#ith contributor
									{
										
										x1<-NULL
										x2<-NULL
										if(runif(1) >= d[rr]){ x1<-sample(names(freq[[l]]),1,prob=freq[[l]])} #remove allele from data
										if(runif(1) >= d[rr]) {x2<-sample(names(freq[[l]]),1,prob=freq[[l]])}#remove allele from data
										zz2<-c(zz2,c(x1,x2))
									}
									cc<-NULL
									if(prC!=0){
										din<-runif(1)
										if(din<=prC){cc<-sample(names(freq[[l]]),1,prob=freq[[l]])}
									}
									zz<-c(zz,length(unique(c(zz1,zz2,cc))))

								}
								locvec[[l]]<-zz#sum the #of alleles among the unknown and the 
								# knwon contributors
								# locvec2[[l]]<-hh
							}
							locus2<-locvec

							tmp0<-apply(sapply(locvec,rbind),1,sum)
							return(tmp0)
						}


				}}
				#--- now if no profiled individuals under Hd
				else
				{
					locosimu.loc<-function(d,contri,x,freq)
					{
						#-- get the names of TpFinal to get the contributors
						
						
						#-------------- if x=0, no uknown contributors
						if(x==0)
						{
							stop('at least one contributor must be stated under Hd')
						}
						
						# if there are unknown contributors
						
						if(x!=0)
						{
							locvec<-vector('list',length=length(locnames))
							names(locvec)<-locnames
							for(l in (locnames))
							{
								zz<-NULL
								for(rr in 1:length(d))
								{
									# now treat unknown contr
									yy1<-NULL
									for(a in 1:x)#ith contributor
									{
										
										x1<-NULL
										x2<-NULL
										if(runif(1) >= d[rr]){ x1<-sample(names(freq[[l]]),1,prob=freq[[l]])} #remove allele from data
										if(runif(1) >= d[rr]) {x2<-sample(names(freq[[l]]),1,prob=freq[[l]])}#remove allele from data
										yy1<-c(yy1,x1,x2)
									}
									if(prC==0){zz<-c(zz,(length(unique(yy1))))}
									else
									{
										cc<-NULL
										din<-runif(1)
										if(din<=prC){cc<-sample(names(freq[[l]]),1,prob=freq[[l]])}
										zz<-c(zz,(length(unique(c(yy1,cc)))))
									}

								}
								locvec[[l]]<-zz#sum the #of alleles among the unknown and the knwon contributors
							}
							locus2<-locvec
							tmp0<-apply(sapply(locvec,rbind),1,sum)
							return(tmp0)
						}


				}}
				tmp2<-replicate(nrep,locosimu.loc(d,contri,x,freq))
				tabHd<-d[which(tmp2==nb,arr.ind=TRUE)[,1]]
				signif(quantile(tabHd,c(5,95)/100,type=1,na.rm=TRUE),2)
				}
				
				cc<-as.numeric(tclvalue(prC))
				print('-- Determination of PrD ranges in progress --')
				r0<-locosimu(d=vecD,prC=cc,contri=TpFinal2,freq=data0,x=xp,nrep=100,nb=nbAll.mean,locnames=loc0)
				r1<-locosimu(d=vecD,prC=cc,contri=TdFinal2,freq=data0,x=xd,nrep=100,nb=nbAll.mean,locnames=loc0)

				print('-----------------  Done  --------------------')

				ranges0<-range(r0,r1)
				print('======== Sensitivity analysis ===========')
				xp<-as.numeric(tclvalue(ncHp))
				xd<-as.numeric(tclvalue(ncHd))
				theta0<-as.numeric(tclvalue(theta))
				
				
				LRres<-vector('list', length(loc0))
				Hdres<-vector('list',length(loc0))
				Hpres<-vector('list',length(loc0))
				names(Hpres)<-loc0
				names(Hdres)<-loc0
				names(LRres)<-loc0
				
				for(jjj in (loc0))
				{		
					# cat(paste(round(jjj*100/length(loc0)),'%',sep=''),'completed','\n')
					mark0<-jjj
					print(mark0)
					rep0<-cspFinal[[jjj]]
					if(is.list(TdFinal )){ tmpTd<-unlist(TdFinal[jjj])} else{ tmpTd<-0}
				
					tp<-unlist(TpFinal[jjj])
					tv<-unlist(VdFinal[[jjj]])
					# loop to evaluate different PrD probabilities in vecD
					tmp.hp<-NULL
					tmp.hd<-NULL
						 
						
						# print(data0[[mark0]])		
					for(kkk in 1:length(vecD))
					{
						
						#numerator Pr(E|Hp)
						d<-vecD[kkk]
						# print(rep(d,length(tp)/2 + xp))
						tmp.hp<-c(tmp.hp,likEvid(Repliste=(rep0),T=tp,V=0,x=xp,theta=theta0,prDHet=rep(d,length(tp)/2 + xp),prDHom=rep(d^2,length(tp)/2 + xp),prC=cc,freq=data0[[mark0]]))
						
									
						tmp.hd<-c(tmp.hd,likEvid(Repliste=(rep0),T=tmpTd,V=tv,x=xd,theta=theta0,prDHet=rep(d,length(tmpTd)/2 + xd),prDHom=rep(d^2,length(tmpTd)/2 + xd),prC=cc, freq=data0[[mark0]]))# V does not contribute to replicate probability
						# V does not contribute to replicate probability)
			
					}	
						
						LRres[[jjj]]<-tmp.hp/tmp.hd
						Hdres[[jjj]]<-tmp.hd
						Hpres[[jjj]]<-tmp.hp
				}
				
				print('================  Done   ================')
				tmp<-NULL
				for(i in LRres){
				tmp<-cbind(tmp,i)}

				tmp<-apply(tmp,1,prod)
				if('Inf' %in% range(tmp) | '-Inf' %in% range(tmp)| NaN %in% range(tmp) )
				{
					
					stop(tkmessageBox(message="infinite LR values, please change the model parameters",icon="error",type="ok"))
				}
				
				
				Myhscale <- 1   # Horizontal scaling
				Myvscale <- 1
				dd <- tktoplevel()
				# tkconfigure(dd,cursor="watch")

				tkwm.title(dd,paste('LR plot, forensim v.',versionNum,sep=''))
				frameC<-tkframe(dd)
				
				Dplot.loc<-function()
				{
					params <- par(bg="white")
					# plot(vecD,log(tmp,10),ylab='log10 LR',xlab='Probability of Dropout',cex.lab=1.3,xlim=c(0,1),pch=19,ylim=range(log(tmp,10),finite=TRUE))
					plot(vecD,log(tmp,10),ylab=expression(log[10](LR)),xlab=expression(Pr(D)),cex.lab=1.3,xlim=c(0.01,0.99),type='l',ylim=range(log(tmp,10),finite=TRUE),lty=1,xaxt='n')
					axis(1,at=c(0.01,0.1,0.2,0.4,0.6,0.8,0.99))
					# segments(ranges0[1],vecD)
					# quantiles function uses an algorithm that can yield intermediate values of d that are not in vecD, since limited numbers ar explored in [0.01,0.99], we need to make sure that the segments reported in the plot are accurate
					
					if(!any(is.na(ranges0)))#first chech that the estimation worked properly given the contributors chosen by the user
					{
						minLR<-min(log(tmp,10))
						if(minLR>0) minLR<-0 else minLR<- minLR-15
						x0<-log(tmp[which(vecD==ranges0[1])],10)
						y0<-log(tmp[which(vecD==ranges0[2])],10)
							
						arrows(ranges0[1],minLR,ranges0[1],x0,col='red',lwd=2)
						arrows(ranges0[2],minLR,ranges0[2],y0,col='red',lwd=2)
						
					}
					else
					{
						tkmessageBox(message=paste('Ranges of drop-out could not be determined with the chosen LR parameters.\n Please check that the hypothesised contributors are sufficient to explain the', nbAll.mean, 'observed alleles',sep=' '),icon='info')
					}
					
					# lines(vecD, log(tmp,10),lty=3,col='gray')
					grid()
					title('LR vs. probability of dropout', cex=1.3)
					par(params)
				}
				
				img <- tkrplot(frameC,fun= Dplot.loc,hscale=Myhscale,vscale=Myvscale)
				CopyToClip <- function()
				{
					tkrreplot(img)
				}
				
				# CopyToLog<-function()
				# {
				
				# }
				# log.but <- tkbutton(frameC,text="Generate log file",font="courrier 10",fg="darkblue",command=CopyToLog)
				
				copy.but <- tkbutton(frameC,text="Copy to Clipboard",font="courrier 10",fg="darkblue",command=CopyToClip)
				# export.but <- tkbutton(frameC,text="Export results",font="courrier 10",fg="darkblue",command=CopyToClip)
#				
			likHp<-apply(sapply(Hpres,rbind),1,prod)#get the prodcut of the likelihoods among loci for all Drop values
			likHd<-apply(sapply(Hdres,rbind),1,prod)

			LRtab2<-signif(cbind.data.frame(vecD,likHp,likHd,likHp/likHd,log10(likHp/likHd)),4)
			colnames(LRtab2)<-c('Pr(D)','Pr(E|Hp)','Pr(E|Hd)','LR','log10(LR)')

				excel.but2<-tkbutton(frameC, text="Export Log File",fg="darkblue", font="courrier 10",command=function() exportFile2(LRtab2,r0,r1))#,command=function() openFile())
				info.but<-tkbutton(frameC, text="Info?",fg="darkblue", font="courrier 10",command=function() infoSP())
				tkgrid(img)
				tkgrid(copy.but)
				# tkgrid(log.but)
				tkgrid(excel.but2,columnspan=13)
				tkgrid(info.but,columnspan=13)
				# tkgrid(plot.but, excel.but,rowspan=10,sticky='ew')
				tkpack(frameC)
			}
			
			
			
		# grid frame f1 contains the tables
		tkgrid(plot.but, excel.but,rowspan=10,sticky='ew')#,columnspan=10,pady=25,columnspan=10)
		# tkpack(plot.but, excel.but, side="left", fill="x", expand=1)

		# tkgrid(excel.but,pady=25,columnspan=10)

		tkgrid(f1,columnspan=10)
		}

		
		
			
		go.butt<-tkbutton(bottomFrame,text='OK!', font='courrier 14 bold',fg='blue', command= analyse.loc)#
		
		# 

		#----- grid and pack, analysis frame
			tkgrid(hypoFrame1, pady=12,padx=10)

		
		tkgrid(ncFrame,pady=12)#,pady=10, padx=10)
		tkgrid(probFrame,pady=12)
		tkgrid(secondFrame,firstFrame, thirdFrame,pady=10,padx=10)#, pady=10, padx=10)
		
		
		
		
		tkgrid(bottomFrame,pady=10,columnspan=45)
		tkgrid(go.butt,columnspan=45)

		

	}
	
	tf <- tktoplevel()
	versionNum<-'4.0'#sessionInfo()$otherPkgs$forensim$Version
	tkwm.title(tf,paste('LRmix: Likelihood Ratio Calculator',', forensim v.', versionNum,sep=''))
	# icn <- tkimage.create("photo", file=system.file("files/test.GIF", package = "forensim"))#"test.GIF")
	#TclTklabel <- tklabel(frame1, image=icn, background="white")
	done <- tclVar(0)
	filePath<-tclVar('')
	# filePath4<-tclVar('')
	extens1<-tclVar('')
	locus<-tclVar('')
	#replicates
	repl<-tclVar('')

	openFile<-function(file0,caselist,top,ext)
	{		
		fileName<-tclvalue(tkgetOpenFile(parent=top,initialdir=tclvalue(file0),multiple="true",
		filetypes="{{CSV Files} {.csv .txt}} ")) #tclvalue(tkgetOpenFile())
		if (!nchar(fileName))
		{
			tkmessageBox(message="No file was selected!")
		}
		else
		{
			tmp<-sub('\\}',fileName,replacement='')
			tmp2<-sub('\\{',tmp,replacement='')
			tclvalue(file0)<-tmp2
			foo3<-strsplit(tmp2,'/')[[1]]
			tclvalue(ext)<-strsplit(foo3[length(foo3)],'\\.')[[1]][2]
			tkinsert(caselist,0,paste(foo3[length(foo3)],sep=":"))
		}
	}
	
#---------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
'sampleProf'<-function()
{
	tprof <- tktoplevel()
	tkwm.title(tprof,"LRmix: import DNA sample profiles")
	filesFrame<-tkframe(tprof)
	tkgrid(tklabel(filesFrame, text="   DNA samples   ",font=font4,foreground="blue"), columnspan=9)
	bottomFrame<-tkframe(tprof)
	#white box will contain the lsist of the available files
	caselist <- tklistbox(filesFrame,height=10,selectmode="extended",background="white",width=25)
   #------------- Display profiles from crime stain, so that the user can choose the loci ---------#
   displayprofile<-function()
   {	
		prof<- tktoplevel()
		tkwm.title(prof,"Sample profile")
		done <- tclVar(0)
		f1 <- tkframe(prof, relief="groove", borderwidth=4)#,bg='white')
		repFrame1 <-tkframe(f1, relief="groove", borderwidth=4)#,bg='white')
		repFrame2 <-tkframe(f1, relief="groove", borderwidth=4)#,bg='white')
		repFrame3 <-tkframe(f1, relief="groove", borderwidth=4)#,bg='white')
		repFrame4 <-tkframe(f1, relief="groove", borderwidth=4)#,bg='white')

		f3 <-tkframe(f1, relief="groove", borderwidth=4)#,bg='white')
		f4 <-tkframe(prof, relief="flat", borderwidth=4)#,bg='white')
		tkgrid(tklabel(f1, text="      Replicate/Loci selection        ",font=font4,foreground="blue"), columnspan=9)
		#a max of four replicates
		

		if(tclvalue(extens1)=='txt')
		stainFile<-read.table(tclvalue(filePath),header=TRUE,as.is=TRUE,sep='\t',na.strings='')#strings as strings, avoid converting to factors

		else{
		stainFile<-read.csv(tclvalue(filePath),header=TRUE,as.is=TRUE,na.strings='')#strings as strings, avoid converting to factors
		}
		
		#remove the Amel Marker
		if("AMEL" %in% stainFile$Marker)
		{
			stainFile<-stainFile[-which(stainFile$Marker=='AMEL'),]
		}

		#handling replicates
		infoStain<-names(stainFile)
		# sample Infor;ation: replicates
		##### General verifications
		if(infoStain[1]!='SampleName')
		{
			tkmessageBox(message="Format error, pleae refer to the reference file",icon="error",type="ok")
		}
		if(infoStain[2]!='Marker')
		{
			tkmessageBox(message="Format error, pleae refer to the reference file",icon="error",type="ok")
		}
		
		
		# number of replicates
		sampname<-(unique(stainFile[,1]))
		
		if(length(sampname)>4)
		{
			tkmessageBox(message=paste("There are",length(sampname),"replicates, only the four first will be used"),icon="warning",	type="ok")
		}
		locusStain<-as.character(unique(stainFile$Marker))
		#okprof function selects the loci: these sould be selected only for the samples
		#and not the reference profiles
		
		okprof.but<-tkbutton(f1, text="OK!",fg="darkblue", font="courrier 12 bold",command=function() okprof())#,command=function() openFile())
		
		tkgrid(okprof.but,pady=2,columnspan=12)#sticky='s')

		okprof<-function()
		{
			selecLoci<-rep(0,length(locusStain))
			selecRep<-rep(0,length(sampname))
			for(k in 1:length(locusStain))
			{
				selecLoci[k]<-as.numeric(tclvalue(get(paste('loc',k,sep=''))))
			}
			
			for(k in 1:length(sampname))
			{
				selecRep[k]<-as.numeric(tclvalue(get(paste('tclRep',k,sep=''))))
				
			}
			#to get user's choice of replicates

			if(length(which(selecRep==1))==0)
			{
				(tkmessageBox(message="At least one replicate must be selected",icon="error",type="ok"))
			}
			else {
			
				tclvalue(repl)<-sampname[which(selecRep==1)]
			}
			
			# get user's choice of the loci
			if(length(which(selecLoci==1))==0)
			{
				(tkmessageBox(message="At least one locus must be selected",icon="error",type="ok"))
			
			}
			
			else
			{
				tclvalue(locus)<-locusStain[which(selecLoci==1)]
			}


			tkdestroy(prof)
			tkdestroy(tprof)
			# return(locusStain2)
		}
		
		
		
		#--------------Add replicates checkbuttons
		# first create the tcl variables
		for(m in 1:length(sampname)){
		assign(paste('tclRep',m,sep=''), tclVar(1))}
		for(i in 1:length(sampname))
		{
			tkgrid(tkcheckbutton(get(paste('repFrame',i,sep='')), text=sampname[i],fg='blue', variable=get(paste('tclRep',i,sep='')),font='courrier 14'),sticky='w',columnspan=9)
			
		}
		
		
		#create tclVar locus
		for(m in 1:length(locusStain)){
		assign(paste('loc',m,sep=''), tclVar(1))}
		j=1
		while(j<(length(sampname)+1))#replicates
		{
			for(i in 1:length(locusStain))
			{
				tmp0<-stainFile[stainFile$SampleName==sampname[j],]
				tmp1<-tmp0[tmp0$Marker==locusStain[i],-c(1,2)]
				# -- tmpProf: profile of the ith locus
				tmpProf<-unlist(tmp1[which(!is.na(tmp1))])
				names(tmpProf)<-NULL
				tkgrid(tkcheckbutton(get(paste('repFrame',j,sep='')), text=c(locusStain[i],"(",tmpProf,")"), variable=get(paste('loc',i,sep='')),font='courrier 8',padx=0),sticky='w')
				
			}
			j=j+1
		}
		
		#second create the checkbuttons
		
		if(length(sampname)>=4)
		{
			tkgrid(repFrame1,repFrame2,repFrame3,repFrame4,padx=12,pady=8 )
			
		}
		else
		{
		
			if(length(sampname)==1)
			{
				# tkgrid(tklabel(rep1, text=paste("Replicate 1",sampname,sep=':'),font=font4,foreground="blue"), columnspan=9)
				tkgrid(repFrame1,padx=12,pady=8 )
				
			}
			
			if(length(sampname)==2)
			{
				tkgrid(repFrame1,repFrame2,padx=12,pady=8 )
				# tkgrid(repFrame2,padx=12,pady=8 )

			
			}
			
			if(length(sampname)==3)
			{
			
				tkgrid(repFrame1,repFrame2,repFrame3,repFrame4,padx=12,pady=8 )
				# tkgrid(tklabel(rep1, text=paste("Replicate 1",sampname[1],sep=':'),font=font4,foreground="blue"), columnspan=9)
				# tkgrid(tklabel(rep2, text=paste("Replicate 2",sampname[2],sep=':'),font=font4,foreground="blue"), columnspan=9)
				# tkgrid(tklabel(rep3, text=paste("Replicate 1",sampname[3],sep=':'),font=font4,foreground="blue"), columnspan=9)
			
			}



		}
		
		tkgrid(f1,sticky='we')
   }

   
readFile.but<-tkbutton(bottomFrame, text="Import datafile",fg="darkblue", font="courrier 12",command=function() openFile(filePath,caselist,tf,extens1))
profile.but<-tkbutton(bottomFrame, text="Display profile",fg="darkblue", font="courrier 12",command=function() displayprofile())
tkgrid(filesFrame) 
tkgrid(caselist,padx=10,pady=10)
tkgrid(bottomFrame)
tkpack(readFile.but,  profile.but,side='left')
}		



#--------------------------------------------------------------------------------#
#--------- Reference profiles
#----------------------------------------------------------------------------------
filePath3<-tclVar(''); 	extens3<-tclVar('')
filePath2<-tclVar(''); 	extens2<-tclVar('')
'refProf'<-function()
{
	tref <- tktoplevel()
	tkwm.title(tref,"Reference profiles")
	#--------------------- refrence profiles---------------------#		
	refFrame<-tkframe(tref)
	buffer1<-tkframe(refFrame,relief='groove')
	buffer1.tit<-tkframe(refFrame,relief='groove')
	buffer2.tit<-tkframe(refFrame)
	buffer3.tit<-tkframe(refFrame)

	buffer2<-tkframe(refFrame)
	buffer3<-tkframe(refFrame,relief='groove')

	# tkgrid(tklabel(buffer1, text="   Suspect(s)   ",font=font4,foreground="blue"), columnspan=9)

	#------ suspect 
	suspectFrame<-tkframe(buffer1, relief='groove')

	#------- victim
	victimFrame<-tkframe(buffer2,relief='groove')
	elimFrame<-tkframe(refFrame)
	bottomFrame<-tkframe(refFrame)
	#-------- suspects frame scrollabr 
	scr2 <- tkscrollbar(buffer1, repeatinterval=10)#, command=function(...)tkyview(caselist,...))
	#white box will contain the lsist of the available files
	caselist2 <- tklistbox(buffer1,height=5,selectmode="extended",background="white",width=20)

	# tkgrid(refFrame)
	
	#-------victim(s) frame scrollabr 
	scr3 <- tkscrollbar(buffer2, repeatinterval=10)#, command=function(...)tkyview(caselist,...))
	#white box will contain the lsist of the available files
	caselist3 <- tklistbox(buffer2,height=5,selectmode="extended",background="white",width=20)
	
	
	#-------extra profiles frame scrollabr 
	# scr4 <- tkscrollbar(buffer3, repeatinterval=10)#, command=function(...)tkyview(caselist,...))
	#white box will contain the lsist of the available files
	# caselist4 <- tklistbox(buffer3,height=5,selectmode="extended",background="white",width=20)
	###---------Buttons: Open File and check?
	#openfile is called when button Openfile is pressed: values of filepath2,caselist2, etc are updated 
	ref.but<-tkbutton(suspectFrame, text="Open File",fg="darkblue", font="courrier 10", command= function() 
	openFile(filePath2,caselist2,tref,extens2))
	# check.but<-tkbutton(suspectFrame, text="Check?",fg="darkblue", font="courrier 10")#,command=function() openFile())
	# a lot of TCLvariables in order to 
	
	ref.but2<-tkbutton(victimFrame, text="Open File",fg="darkblue", font="courrier 10", command= function() 
	openFile(filePath3,caselist3,tref,extens3))
	# check.but2<-tkbutton(victimFrame, text="Check?",fg="darkblue", font="courrier 10")#,command=function() openFile())
	# a lot of TCLvariables in order to 
	
	# ref.but3<-tkbutton(elimFrame, text="Open File",fg="darkblue", font="courrier 10", command= function() 
	# openFile(filePath4,caselist4,tref,extens4))
	# check.but2<-tkbutton(victimFrame, text="Check?",fg="darkblue", font="courrier 10")#,command=function(
	
	
	ok.but<-tkbutton(bottomFrame,text=' OK ', font=font7,fg='blue',command=function()tkdestroy(tref))
	
	# checkInput<-function(){}
	#suspects
	tkgrid(caselist2)
	tkgrid(suspectFrame,padx=10,pady=10)
	tkgrid(buffer1.tit)
	tkgrid(tklabel(buffer1.tit, text="   Suspect(s)   ",font=font4,foreground="blue"), columnspan=9)
	tkgrid(buffer1,pady=10,padx=10)
	tkgrid(ref.but)#, check.but)
	# victims
	tkgrid(tklabel(buffer2.tit, text="   Victim(s)   ",font=font4,foreground="blue"),pady=3)
	tkgrid(caselist3)
	tkgrid(buffer2.tit)
	tkgrid(buffer2,pady=10,padx=10)
	tkgrid(victimFrame,pady=10,padx=10)
	# tkgrid(tklabel(buffer2, text="   Victim(s): known contributors   ",font=font4,foreground="blue"), columnspan=9)
	tkgrid(ref.but2)#, check.but2)
	tkgrid(ok.but,pady=2)
	tkgrid(bottomFrame)
	tkgrid(refFrame)	

}
	
	

	
	##############--------Main frame-----------------------------------------------
	
			
	main <- tkframe(tf, relief="groove", borderwidth=2)
	frame0<-tkframe(main, relief="groove", borderwidth=2)
	frame1<-tkframe(main)
	frame2<-tkframe(main)
	# frameButt <- tkframe(tf, relief="groove", borderwidth=4)
	#icn <- tkimage.create("photo", file=system.file("files/test.GIF", package = 		"forensim"))#"test.GIF")
	#TclTklabel <- tklabel(frame1, image=icn, background="white")
	labh <- tklabel(frame0,bitmap='questhead')#, image=icn)
	#labh <- tklabel(frame1)
	tkbind(labh, "<Button-1>", function() 'hh')
	# tkgrid(labh)
	
	
	#--------------
	# labh <- tklabel(tf)
	#labh <- tklabel(frame1)
	#tkbind(labh, "<Button-1>", function() 'hh')
	tkgrid(tklabel(frame0,text="   Evaluation of Likelihood Ratios  ", font="times 20", foreground="darkblue"),labh)
	
	#---------
	
	samp.butt<-	tkbutton(frame1, text="    Load Sample Profiles   ",fg="darkblue", font="courrier 12", command=sampleProf)
	ref.butt<-tkbutton(frame1, text="  Load Reference Profiles  ",fg="darkblue", font="courrier 12",command=refProf)
	
	
#-- define allele frequencies frame
	extens22<-tclVar('')
	filePath22<-tclVar('')
	
	freq.butt<-tkbutton(frame1,text=' Import allele frequencies', font='courrier 13',fg='darkblue', command= function() importAF(frame1,filePath22,extens22))
	

	#function to import the file that contains the AF
	# it will update the valueof filepath22 to its current value
	importAF<-function(fa,pa,ex,d0)#frame, pathfile, and extension var
	{
		file0<-tclvalue(tkgetOpenFile(parent=fa,initialdir=tclvalue(pa),multiple="true",filetypes="{{CSV Files} {.csv .txt}}"))
	
		if (!nchar(file0))
		{
			tkmessageBox(message="No file was selected!")
		}
		else
		{
			tmp<-sub('\\}',file0,replacement='')
			tmp2<-sub('\\{',tmp,replacement='')
			tclvalue(file0)<-tmp2
			foo3<-strsplit(tmp2,'/')[[1]]
			tclvalue(ex)<-strsplit(foo3[length(foo3)],'\\.')[[1]][2]
			tclvalue(pa)<-tmp2
		}
		#----------
		
		# return(dataFreq)
	}		
	
	#----- frame for allele frequencies
		
	
	# call analyse() with arguments locus and filepath, that are modified by 
	analyse.butt<-tkbutton(frame2, text=" Done! ",fg="blue", font="courrier 12 bold",command=function() analyse(locus,repl,filePath,filePath2,filePath3, extens1, extens2, extens3,infoStain),relief='raised')
	tkgrid(samp.butt, padx=10,pady=10)
	tkgrid(ref.butt,padx=10,pady=10)
	tkgrid(freq.butt, padx=10,pady=10)

	tkgrid(analyse.butt,padx=10,pady=3)
	tkgrid(frame0)
	tkgrid(frame1)
	tkgrid(frame2)
	tkgrid(main)

}
