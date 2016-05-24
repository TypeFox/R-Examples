sse<-function(h,...){
		try(({dat2Affy.m<-dat2Affy.m;datAgOne2.m<-datAgOne2.m;datAgTwo2.m<-datAgTwo2.m;
		datIllBA2.m2<-datIllBA2.m2;lumi_NQ.m<-lumi_NQ.m;data.matrix_Nimblegen2.m<-data.matrix_Nimblegen2.m;
		data.matrixNorm.m<-data.matrixNorm.m;data.matrix_onlineNorm.m<-data.matrix_onlineNorm.m;l<-l;tree<-tree;
		}),silent=TRUE)
	platforms=NULL
	aa=0;bb=0;cc=0;dd=0;ee=0;ff=0;gg=0;hh=0;
	try(({
		if(exists("dat2Affy.m"))aa=length(dat2Affy.m)
		if(exists("datAgOne2.m"))bb=length(datAgOne2.m)
		if(exists("datAgTwo2.m"))cc=length(datAgTwo2.m)
		if(exists("datIllBA2.m2"))dd=length(datIllBA2.m2)
		if(exists("lumi_NQ.m"))ee=length(lumi_NQ.m)
		if(exists("data.matrix_Nimblegen2.m"))ff=length(data.matrix_Nimblegen2.m)
		if(exists("data.matrixNorm.m"))gg=length(data.matrixNorm.m)
		if(exists("data.matrix_onlineNorm.m"))hh=length(data.matrix_onlineNorm.m)
		}),silent=TRUE)
	if(aa!=0)platforms=c(platforms,"Affymetrix")
	if(bb!=0)platforms=c(platforms,"Agilent_OneColor")
	if(cc!=0)platforms=c(platforms,"Agilent_TwoColor")
	if(dd!=0)platforms=c(platforms,"Illumina_Beadarray")
	if(ee!=0)platforms=c(platforms,"Illumina_Lumi")
	if(ff!=0)platforms=c(platforms,"Nimblegen")
	if(gg!=0)platforms=c(platforms,"Series_Matrix")
	if(hh!=0)platforms=c(platforms,"Online_Data")

	size_Affy=NULL;size_Ag1=NULL;size_Ag2=NULL;size_Il_B=NULL;size_Il_L=NULL;
	size_N=NULL;size_S=NULL;size_O=NULL;

	rm(size_Affy,size_Ag1,size_Ag2,size_Il_B,size_Il_L,size_N,size_S,size_O)
	
	x=NULL
	f<-function(h,...){
		x<<-svalue(h$obj)
		}
	z=NULL
	w_dge<-gwindow("Select your data",width=260,height=280,visible=FALSE,horizontal=FALSE)
	gp_dge<-ggroup(container=w_dge,horizontal=FALSE)
	cbg_dge<-gcheckboxgroup(platforms,container=gp_dge,handler=f)
	svalue(cbg_dge,index=FALSE)<-1:8
	
	gp2_dge<-ggroup(container=gp_dge,width=30,height=15,horizontal=TRUE)
	addSpring(gp2_dge)
	y<-gbutton("CANCEL",border=TRUE,handler=function(h,...){
		dispose(w_dge)
	},container=gp2_dge,anchor=c(1,-1))
	y2<-gbutton("OK",border=TRUE,handler=function(h,...){
		dispose(w_dge)
		if(length(x)!=0){
		if(length(which(x=="Affymetrix"))!=0){
			galert("Please wait while Sample size estimation",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			sds<-rowSds(dat2Affy.m)
			size_Affy<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
			ssize.plot(size_Affy,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
			if(length(size_Affy)!=0){
				visible(g1_1)<-FALSE
				l$Affymetrix$SSE_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			display()	
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
			}
		if(length(which(x=="Agilent_OneColor"))!=0){
			galert("Please wait while Sample size estimation",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			sds<-rowSds(datAgOne2.m)
			size_Ag1<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
			ssize.plot(size_Ag1,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
			if(length(size_Ag1)!=0){
				visible(g1_1)<-FALSE
				l$Agilent_OneColor$SSE_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
			display()	
			}
		if(length(which(x=="Agilent_TwoColor"))!=0){
			galert("Please wait while Sample size estimation",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			sds<-rowSds(datAgTwo2.m)
			size_Ag2<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
			ssize.plot(size_Ag2,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
			if(length(size_Ag2)!=0){
				visible(g1_1)<-FALSE
				l$Agilent_TwoColor$SSE_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
			display()	
			}
		if(length(which(x=="Illumina_Beadarray"))!=0){
			galert("Please wait while Sample size estimation",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			sds<-rowSds(datIllBA2.m2)
			size_Il_B<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
			ssize.plot(size_Il_B,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
			if(length(size_Il_B)!=0){
				visible(g1_1)<-FALSE
				l$Illumina_Beadarray$SSE_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
			display()	
			}
		if(length(which(x=="Illumina_Lumi"))!=0){
			galert("Please wait while Sample size estimation",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			sds<-rowSds(lumi_NQ.m)
			size_Il_L<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
			ssize.plot(size_Il_L,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
			if(length(size_Il_L)!=0){
				visible(g1_1)<-FALSE
				l$Illumina_Lumi$SSE_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			display()	
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
			}
		if(length(which(x=="Nimblegen"))!=0){
			galert("Please wait while Sample size estimation",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			sds<-rowSds(data.matrix_Nimblegen2.m)
			size_N<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
			ssize.plot(size_N,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
			if(length(size_N)!=0){
				visible(g1_1)<-FALSE
				l$Nimblegen$SSE_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			display()	
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
  			}
		if(length(which(x=="Series_Matrix"))!=0){
			galert("Please wait while Sample size estimation",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			sds<-rowSds(data.matrixNorm.m)
			size_S<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
			ssize.plot(size_S,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
			if(length(size_S)!=0){
				visible(g1_1)<-FALSE
				l$Series_Matrix$SSE_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			display()	
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
			}
		if(length(which(x=="Online_Data"))!=0){
			galert("Please wait while Sample size estimation",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			sds<-rowSds(data.matrix_onlineNorm.m)
			size_O<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
			ssize.plot(size_O,xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
			if(length(size_O)!=0){
				visible(g1_1)<-FALSE
				l$Online_Data$SSE_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			display()	
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
			}
		dispose(w_dge)
		}else{
			gmessage("Plz select the data for sample size estimation","Select Data")
			}
		},container=gp2_dge,anchor=c(1,-1))
	visible(w_dge)<-TRUE
	}

