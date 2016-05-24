#___________________________________________
#Four-allele model
#____________________________________________

"A4.simu" <- function()
{
	#four-allele model: numerical results
	# if(!require(tcltk)) stop("package tcltk is required")
	tclRequire("Tktable")
	#fonts definition
	font0 <- tkfont.create(family="times",size=35,weight="bold",slant="italic")
	font1<-tkfont.create(family="times",size=14,weight="bold")#,slant="italic")
	font2<-tkfont.create(family="times",size=16,weight="bold",slant="italic")
	font3<-tkfont.create(family="times",size=10,weight="bold")#,slant="italic")
	font4<-tkfont.create(family="times",size=12)#,slant="italic")
	fourAmod<- function(Mx,peak1,peak2,peak3,peak4)
	{
		expA <- c(Mx/2, Mx/2,Mx/2,(1-Mx)/2,(1-Mx)/2,(1-Mx)/2)
		expB <- c(Mx/2, (1-Mx)/2,(1-Mx)/2,Mx/2,Mx/2,(1-Mx)/2)
		expC<-c((1-Mx)/2,Mx/2,(1-Mx)/2,Mx/2,(1-Mx)/2,Mx/2)
		expD<-c((1-Mx)/2,(1-Mx)/2,Mx/2,(1-Mx)/2,Mx/2,Mx/2)
		
		obsA <- peak1/(peak1+peak2+peak3+peak4)
		obsB <- peak2/(peak1+peak2+peak3+peak4)
		obsC<-peak3/(peak1+peak2+peak3+peak4)
		obsD<-peak4/(peak1+peak2+peak3+peak4)
		resid <- (expA-obsA)^2+(expB-obsB)^2+(expC-obsC)^2+(expD-obsD)^2
		genotypes<-c("AB,CD","AC,BD","AD,BC","BC,AD","BD,AC","CD,AB")

		
		Mx.Conditioned<-abs(c((peak1+peak2)/(peak1+peak2+peak3+peak4),
		(peak1+peak3)/(peak1+peak2+peak3+peak4),
		(peak1+peak4)/(peak1+peak2+peak3+peak4),
		(peak2+peak3)/(peak1+peak2+peak3+peak4),
		(peak2+peak4)/(peak1+peak2+peak3+peak4),
		(peak3+peak4)/(peak1+peak2+peak3+peak4)))
		#to be completed
		#Hb1<-
		#Hb2<-
		result<-data.frame(genotypes,expA,expB,expC,expD,obsA,obsB,obsC,obsD,resid,Mx.Conditioned)
		return(result)
	}


	
	#four-allele model: graphical simulation
	fourAmodG<-function(obsMx,peak1,peak2,peak3,peak4)
	{
		Mx<-seq(0.1,0.9,by=0.05)
		res<-sapply(Mx, function(i) fourAmod(i,peak1,peak2,peak3,peak4)$resid)
		par(mar = c(4,4,2,3)+0.1,oma = c(1,2,0.3,5))
		#log="y"
		plot(Mx,(res[1,]),type='n',lab=c(6,5,1),xlab="Mx (mixture proportion)",ylab="Residuals",cex.lab=1.2,las=1,ylim=c(min(res),max(res)))
		title("Four-allele model simulations")
		pchref<-c(18,15,17,4,8,16)
		col<-c("lightgreen","magenta","purple","cyan","red","darkblue")
		for(i in 1:6)#7 genotype combinations
		{
			points(Mx,(res[i,]),pch=pchref[i],col=col[i])
			lines(Mx,(res[i,]),col=col[i])
		}
		abline(v=obsMx,col="gray",lty=2)
		#text(obsMx,median(res),"Observed Mx",   adj = c(0, -.1))
		legend(0.95,max(res),c("AB,CD","AC,BD","AD,BC","BC,AD","BD,AC","CD,AB",'Obs. Mx'),
		,pch =c(pchref,32),lty=2,col=c(col,'gray'),bg='white',xpd=NA,cex=1.2,
		 box.lwd=2)
	}

	#four-allele model:  likelihood estimation
	fourAmodT<-function(peak1,peak2,peak3,peak4)
	{
		Mx<-seq(0.1,0.9,by=0.05)
		geno<-c("AB,CD","AC,BD","AD,BC","BC,AD","BD,AC","CD,AB")
		res<-sapply(Mx, function(i) signif(fourAmod(i,peak1,peak2,peak3,peak4)$resid,digits=3))
		row.names(res)<-c(geno)
		colnames(res)<-paste('Mx=',Mx,sep='')
		minX<-which(res==min(res),arr.ind=TRUE)
		minC<-signif(Mx[minX[,2]],digits=2)
		minY<-row.names(minX)	
		return(list(minY,minC,res))
	}

	#main frame
	tt <- tktoplevel()
	tkwm.title(tt,"Four-allele model simulations")

    frame1 <- tkframe(tt, relief="groove", borderwidth=2)
    frame2 <- tkframe(tt, relief="groove", borderwidth=2)
   
    xyframe <- tkframe(frame1, relief="groove", borderwidth=2)
    labframe <- tkframe(frame1, relief="groove", borderwidth=2)
    limframe <- tkframe(frame2, relief="groove", borderwidth=2)
    posframe <- tkframe(frame2, relief="groove", borderwidth=2)
    legframe <- tkframe(frame2, relief="groove", borderwidth=2)
    

	xy1var <- tclVar(553)
	xy2var <- tclVar(472)
	xy3var<-  tclVar(623)
	xy4var<-  tclVar(738)
	mixvar <- tclVar(0.70)
	


	TFrame <- tkframe(tt, relief="groove")
	labh <- tklabel(TFrame)
	tkgrid(tklabel(TFrame,text="Four-allele model", font=font2, foreground="red"), labh)
	tkpack(TFrame)

	xy1.entry <- tkentry(xyframe, textvariable=xy1var, width=8)
	xy2.entry <- tkentry(xyframe, textvariable=xy2var, width=8)
	xy3.entry <- tkentry(xyframe, textvariable=xy3var, width=8)
	xy4.entry <- tkentry(xyframe, textvariable=xy4var, width=8)
	
	peak1<-function()
	{
		p1<-as.numeric(tclvalue(xy1var))
		if(p1<0){tkmessageBox(message="Invalid value for the peak height of allele #1",icon="error",type="ok")}
		else{return(p1)}
	}
	
	
	
	peak2<-function()
	{
		p2<-as.numeric(tclvalue(xy2var))
		if(p2<0){tkmessageBox(message="Invalid value for the peak height of allele #2",icon="error",type="ok")
		}
		else{return(p2)}
	}
	
	peak3<-function()
	{
		p3<-as.numeric(tclvalue(xy3var))
		if(p3<0){tkmessageBox(message="Invalid value for the peak height of allele #3",icon="error",type="ok")
		}
		else{return(p3)}
	}
	
	peak4<-function()
	{
		p4<-as.numeric(tclvalue(xy4var))
		if(p4<0){tkmessageBox(message="Invalid value for the peak height of allele #4",icon="error",type="ok")
		}
		else{return(p4)}
	}
	
	# choosexy1.but <- tkbutton(xyframe, text="Enter", command=function()peak1())
	# choosexy2.but <- tkbutton(xyframe, text="Enter", command=function() peak2())
	# choosexy3.but <- tkbutton(xyframe, text="Enter", command=function() peak3())
	# choosexy4.but <- tkbutton(xyframe, text="Enter", command=function() peak4())

	tkgrid(tklabel(xyframe, text="- Peak heights (rfu) -", font=font3,foreground="blue"), columnspan=5)
	tkgrid(tklabel(xyframe,text="Allele #1: "), xy1.entry)#, choosexy1.but)
	tkgrid(tklabel(xyframe,text="Allele #2: "), xy2.entry)#, choosexy2.but)
	tkgrid(tklabel(xyframe,text="Allele #3: "), xy3.entry)#, choosexy3.but)
	tkgrid(tklabel(xyframe,text="Allele #4: "), xy4.entry)#, choosexy4.but)
# Labels frame
#
	lab.entry <- tkentry(labframe, textvariable=mixvar, width=8)
	prop<-function() 
	{
		Mx<-as.numeric(tclvalue(mixvar))
		if(Mx>1 || Mx<0){
		tkmessageBox(message="Mx is the mixture proportion, it must be comprised in the interval [0,1]",icon="error",type="ok")

		}
		else{return(Mx)}
		#print(Mx)
		
	}
	#chooselab.but <- tkbutton(labframe, text="Set", command= function() prop())#function() print(tclvalue(mixvar)))
	tkgrid(tklabel(labframe, text="- Mixture proportion  -",font=font3, foreground="blue"), columnspan=3)
	#lab.entry pour dire que c'est l'entree lab.entry, et le label est celui de "set"
	tkgrid(tklabel(labframe,text="Mx : "), lab.entry)
	#tkgrid(tklabel(labframe,text="1-Mx : "), lab2.entry)
	tkpack(xyframe, labframe, side="left")
	tkpack(frame1)

# 
RCSFrame <- tkframe(tt, relief="groove")

plotFunction<-function()
{
	p1<-peak1()
	p2<-peak2()
	p3<-peak3()
	p4<-peak4()
	obsMx<-prop()
	fourAmodG(obsMx,p1,p2,p3,p4)
}

propFunction<-function()
{
	p1<-peak1()
	p2<-peak2()
	p3<-peak3()
	p4<-peak4()
	res<-fourAmodT(p1,p2,p3,p4)
	resdata<-res[[3]]
	resgeno<-res[[1]]
	resmix<-res[[2]]
	myRarray<- c("Genotype","Mx=0.1","Mx=0.15","Mx=0.20","Mx=0.25","Mx=0.30",
	"Mx=0.35","Mx=0.40","Mx=0.45","Mx=0.50","Mx=0.55","Mx=0.60","Mx=0.65",
	"Mx=0.70","Mx=0.75","Mx=0.80","Mx=0.85","Mx=0.90",
	"AB,CD",resdata[1,],
	"AC,BD",resdata[2,],
	"AD,BC",resdata[3,],
	"BC,AD",resdata[4,],
	"BD,AC",resdata[5,],
	"CD,AB",resdata[6,])

	dim(myRarray) <- c(18,7)
	tclarray <- tclArray()
	for(i in 0:17)
	{
		for (j in (0:6))
		{
		  tclarray[[i,j]] <- myRarray[i+1,j+1]
		}
	}
	 
	
	myRarray2<-resgeno
	dim(myRarray2)<-c(1,length(resgeno))
	myRarray3<-resmix
	dim(myRarray3)<-c(1,length(resmix))
	#myRarray3

	tclarray2 <- tclArray()
	for(k in 0:(length(resgeno)-1))
	{
		tclarray2[[0,k]] <- myRarray2[1,k+1]
	}		
	tclarray3 <- tclArray()
	for(h in 0:(length(resmix)-1))
	{
		tclarray3[[0,h]] <- myRarray3[1,h+1]
	}		
	
		
	save1<-function(filename1="simulation4.txt",filename2="likelihood4.txt")
		{
			#myRarray
			write.table(myRarray,file=filename1,row.names=FALSE,col.names=FALSE)
			loctab<-data.frame(resmix,resgeno)
			colnames(loctab)<-c('Mixture proportion','Genotype combination')
			#writeLines('\n',filename)
			write.table(loctab,file=filename2,row.names=FALSE)
			#write.table(myRarray3,"filter.txt",row.names=FALSE,col.names=FALSE)
		}
		saveFunction<-function()
		{
			
			ss<-tktoplevel()
			SSframe <- tkframe(ss, relief="groove",width=35)
			tkwm.title(ss,"")
			filevar1 <- tclVar("simulation4.txt")
			filevar1.entry <- tkentry(SSframe, textvariable=filevar1, width=12)
			filevar2 <- tclVar("likelihood4.txt")
			filevar2.entry <- tkentry(SSframe, textvariable=filevar2, width=12)
			#filevar2.entry <- tkentry(SSframe, textvariable=filevar, width=12)
			#tkgrid(tklabel(SSframe, text="- Enter filenames -",font=font1, foreground="blue"), columnspan=15)
			save1.butt<-tkbutton(ss, text="Enter", font=font3,command=function() save1(tclvalue(filevar1),tclvalue(filevar2)))
			tkgrid(tklabel(SSframe,text="Simulations results",font=font4), filevar1.entry)
			tkgrid(tklabel(SSframe,text="Maximum likelihood",font=font4), filevar2.entry)
			tkgrid(filevar2.entry, save1.butt)		
			tkpack(SSframe)
		}
		
		tclRequire("Tktable")
		tt<-tktoplevel()
		tkwm.title(tt,"Most likely genotypes combination")
		table1 <- tkwidget(tt,"table",variable=tclarray,rows=18,colwidth=8,cols=7,titlerows=1,background="white")
		table2 <- tkwidget(tt,"table",variable=tclarray2,
		cols=length(resgeno),selectmode="extended",colwidth=10,rows=1,background="lightblue")
		table3 <- tkwidget(tt,"table",variable=tclarray3,
		cols=length(resmix),selectmode="extended",colwidth=10,rows=1,background="lightblue")
			
		tit1<-tkwidget(tt,"label",text="Matrix of the residuals",font=font1,foreground="blue")	
		tit2<-tkwidget(tt,"label",text="Maximum likelihood estimation results",font=font1,foreground="blue")
		tit3<-tkwidget(tt,"label",text="Most likely genotype combinations",font=font1,foreground="blue")
		tit4<-tkwidget(tt,"label",text="Corresponding mixture proportions",font=font1,foreground="blue")
		
			#tkgrid(tklabel(tt,text="File name"), filevar.entry, save1.butt)
		filelab<-tkwidget(tt,"label",text="-Save the results-",font=font3,foreground="blue")
		save.butt<-tkbutton(tt, text="Save", font=font3,command=saveFunction)
		tkpack(tit1,table1,tit3,table2,tit4,table3,save.butt)
		
		
		
}
	

filterFunction<-function()
{
	p1<-peak1()
	p2<-peak2()
	p3<-peak3()
	p4<-peak4()
	obsMx<-prop()
	tmp1<-fourAmod(obsMx,p1,p2,p3,p4)
	Mx.C<-signif(tmp1$Mx.Conditioned,digits=2)
	genotypes<-as.character(tmp1$genotypes)
	tab1<-c("Genotype",genotypes,"Mx conditioned",Mx.C)
	dim(tab1)<-c(7,2)
	tab1array <- tclArray()

	for(m in 0:6){
	for (n in (0:1)){
	  tab1array[[m,n]] <- tab1[m+1,n+1]}}
	
	saveHH<-function(name1="filter4.txt")
		{
			#myRarray
			write.table(tab1,name1,row.names=FALSE,col.names=FALSE)
				
		}		
			
		saveFunction2<-function()
		{
			
			hh<-tktoplevel()
			HHframe<- tkframe(hh, relief="groove")
			tkwm.title(hh,"Filenames")
			filtervar<- tclVar("filter4.txt")
			filtervar.entry <- tkentry(HHframe, textvariable=filtervar, width=12)
			#filevar2.entry <- tkentry(SSframe, textvariable=filevar, width=12)
			#tkgrid(tklabel(HHframe, text="- Enter filename -",font=font1, foreground="blue"), columnspan=15)
			saveHH.butt<-tkbutton(hh, text="Enter", font=font3,command=function() saveHH(tclvalue(filtervar)))
			tkgrid(tklabel(HHframe,text="Genotypes filter",font=font4), filtervar.entry)
			#tkgrid(tklabel(SSframe,text="Maximum likelihood",font=font4), filevar2.entry)
			tkgrid(filtervar.entry, saveHH.butt)		
			tkpack(HHframe)
		}
		
		
		
		
			tclRequire("Tktable")
			tt2<-tktoplevel()
			tkwm.title(tt2,"Genotypes filter")
			#tab1.tit<-tkwidget(tt2,"label",text="Genotypes filter",font=font1,foreground="blue")
			save2.butt<-tkbutton(tt2, text="Save", font=font3,command=saveFunction2)
			tab1.tcl<-tkwidget(tt2,"table",variable=tab1array,rows=8,colwidth=18,cols=2,titlerows=1,background="white")
			tkpack(tab1.tcl,save2.butt)
			#tkpack(tab1.tit,tab1.tcl,save2.butt)
		
		
		
			
	}

	A1.but <- tkbutton(RCSFrame, text="Plot simulations",font=font3, command=plotFunction)#twoAmod())
	A2.but <- tkbutton(RCSFrame, text="Simulation details", font=font3,command=propFunction)
	A3.but <- tkbutton(RCSFrame, text="Genotype filter", font=font3,command=filterFunction)
	tkgrid(A1.but,A2.but,A3.but,ipadx=20)	
	tkpack(RCSFrame)
}


