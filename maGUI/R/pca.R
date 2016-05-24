pca<-function(h,...){
	f<-function(h,...){
		x<<-svalue(h$obj)
		}

	try(({
	dat2Affy.m<-dat2Affy.m;datAgOne2.m<-datAgOne2.m;datAgTwo2.m<-datAgTwo2.m;datIllBA2.m2<-datIllBA2.m2;
	lumi_NQ.m<-lumi_NQ.m;data.matrix_Nimblegen2.m<-data.matrix_Nimblegen2.m;
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

	pca_Affy=NULL;pca_Ag1=NULL;pca_Ag2=NULL;pca_Il_B=NULL;pca_Il_L=NULL;pca_N=NULL;pca_S=NULL;pca_O=NULL;
	rm(pca_Affy,pca_Ag1,pca_Ag2,pca_Il_B,pca_Il_L,pca_N,pca_S,pca_O)
	
	x=NULL
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
		if(length(x)!=0){
		if(length(which(x=="Affymetrix"))!=0){
			pca_Affy<<-prcomp(t(dat2Affy.m))
			}
		if(length(which(x=="Agilent_OneColor"))!=0){
			pca_Ag1<<-prcomp(t(datAgOne2.m))
			}
		if(length(which(x=="Agilent_TwoColor"))!=0){
			pca_Ag2<<-prcomp(t(datAgTwo2.m))
			}
		if(length(which(x=="Illumina_Beadarray"))!=0){
			pca_Il_B<<-prcomp(t(datIllBA2.m2))
			}
		if(length(which(x=="Illumina_Lumi"))!=0){
			pca_Il_L<<-prcomp(t(lumi_NQ.m))
			}
		if(length(which(x=="Nimblegen"))!=0){
			pca_N<<-prcomp(t(data.matrix_Nimblegen2.m))
  			}
		if(length(which(x=="Series_Matrix"))!=0){
			pca_S<<-prcomp(t(data.matrixNorm.m))
			}
		if(length(which(x=="Online_Data"))!=0){
			pca_O<<-prcomp(t(data.matrix_onlineNorm.m))
			}
		dispose(w_dge)
		}else{
			gmessage("Plz select the data for PCA","Select Data")
			}
		if(length(x)!=0){
		if(length(which(x=="Affymetrix"))!=0){
			plot(pca_Affy,main="Affymetrix PCA")
			if(length(pca_Affy)!=0){
				visible(g1_1)<-FALSE
				l$Affymetrix$PCA_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				}
				display()
			} else
		if(length(which(x=="Agilent_OneColor"))!=0){
			plot(pca_Ag1,main="Agilent_OneColor PCA")
			if(length(pca_Ag1)!=0){
				visible(g1_1)<-FALSE
				l$Agilent_OneColor$PCA_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				}
				display()
			} else
		if(length(which(x=="Agilent_TwoColor"))!=0){
			plot(pca_Ag2,main="Agilent_TwoColor PCA")
			if(length(pca_Ag2)!=0){
				visible(g1_1)<-FALSE
				l$Agilent_TwoColor$PCA_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				}
				display()
			} else
		if(length(which(x=="Illumina_Beadarray"))!=0){
			plot(pca_Il_B,main="Illumina_Beadarray PCA")
			if(length(pca_Il_B)!=0){
				visible(g1_1)<-FALSE
				l$Illumina_Beadarray$PCA_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				}
				display()
			} else
		if(length(which(x=="Illumina_Lumi"))!=0){
			plot(pca_Il_L,main="Illumina_Lumi PCA")
			if(length(pca_Il_L)!=0){
				visible(g1_1)<-FALSE
				l$Illumina_Lumi$PCA_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				}
				display()
			} else
		if(length(which(x=="Nimblegen"))!=0){
			plot(pca_N,main="Nimblegen PCA")
			if(length(pca_N)!=0){
				visible(g1_1)<-FALSE
				l$Nimblegen$PCA_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				}
				display()
			} else
		if(length(which(x=="Series_Matrix"))!=0){
			plot(pca_S,main="Series_Matrix PCA")
			if(length(pca_S)!=0){
				visible(g1_1)<-FALSE
				l$Series_Matrix$PCA_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				}
				display()
			} else
		if(length(which(x=="Online_Data"))!=0){
			plot(pca_O,main="Online_Data PCA")
			if(length(pca_O)!=0){
				visible(g1_1)<-FALSE
				l$Online_Data$PCA_Plot<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				}
				display()
			} else
		dispose(w_dge)
		}	
		},container=gp2_dge,anchor=c(1,-1))
	visible(w_dge)<-TRUE
	}
