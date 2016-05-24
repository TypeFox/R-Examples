ntwk<-function(h,...){
		try(({
		use.dat2Affy.m<-use.dat2Affy.m;use.datAgOne2.m<-use.datAgOne2.m;use.datAgTwo2.m<-use.datAgTwo2.m;
		use.datIllBA2.m2<-use.datIllBA2.m2;use.lumi_NQ.m<-use.lumi_NQ.m;
		use.data.matrix_Nimblegen2.m<-use.data.matrix_Nimblegen2.m;
		use.data.matrixNorm.m<-use.data.matrixNorm.m;use.data.matrix_onlineNorm.m<-use.data.matrix_onlineNorm.m;
		DE_Affy<-DE_Affy;DE_Ag1<-DE_Ag1;DE_Ag2<-DE_Ag2;DE_Il_B<-DE_Il_B;DE_Il_L<-DE_Il_L;
		DE_N<-DE_N;DE_S<-DE_S;DE_O<-DE_O;l<-l;tree<-tree;
			}),silent=TRUE)
	platforms=NULL;
	aa=0;bb=0;cc=0;dd=0;ee=0;ff=0;gg=0;hh=0;
	try(({
		if(exists("DE_Affy"))aa=length(DE_Affy)
		if(exists("DE_Ag1"))bb=length(DE_Ag1)
		if(exists("DE_Ag2"))cc=length(DE_Ag2)
		if(exists("DE_Il_B"))dd=length(DE_Il_B)
		if(exists("DE_Il_L"))ee=length(DE_Il_L)
		if(exists("DE_N"))ff=length(DE_N)
		if(exists("DE_S"))gg=length(DE_S)
		if(exists("DE_O"))hh=length(DE_O)
		}),silent=TRUE)
	if(aa!=0)platforms=c(platforms,"Affymetrix")
	if(bb!=0)platforms=c(platforms,"Agilent_OneColor")
	if(cc!=0)platforms=c(platforms,"Agilent_TwoColor")
	if(dd!=0)platforms=c(platforms,"Illumina_Beadarray")
	if(ee!=0)platforms=c(platforms,"Illumina_Lumi")
	if(ff!=0)platforms=c(platforms,"Nimblegen")
	if(gg!=0)platforms=c(platforms,"Series_Matrix")
	if(hh!=0)platforms=c(platforms,"Online_Data")

	adjMat_Affy=NULL;adjMat_Ag1=NULL;adjMat_Ag2=NULL;adjMat_Il_B=NULL;adjMat_Il_L=NULL;
	adjMat_N=NULL;adjMat_S=NULL;adjMat_O=NULL;
	myGraph_Affy=NULL;myGraph_Ag1=NULL;myGraph_Ag2=NULL;myGraph_Il_B=NULL;myGraph_Il_L=NULL;
	myGraph_N=NULL;myGraph_S=NULL;myGraph_O=NULL;

	rm(adjMat_Affy,adjMat_Ag1,adjMat_Ag2,adjMat_Il_B,adjMat_Il_L,
	adjMat_N,adjMat_S,adjMat_O,
	myGraph_Affy,myGraph_Ag1,myGraph_Ag2,myGraph_Il_B,myGraph_Il_L,
	myGraph_N,myGraph_S,myGraph_O)

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
		if(length(x)!=0){
			dispose(w_dge)
		if(length(which(x=="Affymetrix"))!=0){
			galert("Please wait while generating Co-expression network",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			myData_Sel<-use.dat2Affy.m[rownames(DE_Affy),]
			myData_Sel<-t(myData_Sel)
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
			adjMat_Affy<<-adjMat
			myGraph_Affy<<-as(adjMat_Affy,"graphNEL")
			plot(myGraph_Affy,nodeAttrs=makeNodeAttrs(myGraph_Affy,fontsize=18,fillcolor="grey"))
			if(length(myGraph_Affy)!=0){
				visible(g1_1)<-FALSE
				l$Affymetrix$Coexpression_Network<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			display()	
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
			}
		if(length(which(x=="Agilent_OneColor"))!=0){
			galert("Please wait while generating Co-expression network",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			myData_Sel<-use.datAgOne2.m[rownames(DE_Ag1),]
			myData_Sel<-t(myData_Sel)
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
			adjMat_Ag1<<-adjMat
			myGraph_Ag1<<-as(adjMat_Ag1,"graphNEL")
			plot(myGraph_Ag1,nodeAttrs=makeNodeAttrs(myGraph_Ag1,fontsize=18,fillcolor="grey"))
			if(length(myGraph_Ag1)!=0){
				visible(g1_1)<-FALSE
				l$Agilent_OneColor$Coexpression_Network<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			display()	
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
			}
		if(length(which(x=="Agilent_TwoColor"))!=0){
			galert("Please wait while generating Co-expression network",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			myData_Sel<-use.datAgTwo2.m[rownames(DE_Ag2),]
			myData_Sel<-t(myData_Sel)
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
			adjMat_Ag2<<-adjMat
			myGraph_Ag2<<-as(adjMat_Ag2,"graphNEL")
			plot(myGraph_Ag2,nodeAttrs=makeNodeAttrs(myGraph_Ag2,fontsize=18,fillcolor="grey"))
			if(length(myGraph_Ag2)!=0){
				visible(g1_1)<-FALSE
				l$Agilent_TwoColor$Coexpression_Network<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			display()	
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
			}
		if(length(which(x=="Illumina_Beadarray"))!=0){
			galert("Please wait while generating Co-expression network",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			myData_Sel<-use.datIllBA2.m2[rownames(DE_Il_B),]
			myData_Sel<-t(myData_Sel)
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
			adjMat_Il_B<<-adjMat
			myGraph_Il_B<<-as(adjMat_Il_B,"graphNEL")
			plot(myGraph_Il_B,nodeAttrs=makeNodeAttrs(myGraph_Il_B,fontsize=18,fillcolor="grey"))
			if(length(myGraph_Il_B)!=0){
				visible(g1_1)<-FALSE
				l$Illumina_Beadarray$Coexpression_Network<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			display()	
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
			}
		if(length(which(x=="Illumina_Lumi"))!=0){
			galert("Please wait while generating Co-expression network",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			myData_Sel<-use.lumi_NQ.m[rownames(DE_Il_L),]
			myData_Sel<-t(myData_Sel)
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
			adjMat_Il_L<<-adjMat
			myGraph_Il_L<<-as(adjMat_Il_L,"graphNEL")
			plot(myGraph_Il_L,nodeAttrs=makeNodeAttrs(myGraph_Il_L,fontsize=18,fillcolor="grey"))
			if(length(myGraph_Il_L)!=0){
				visible(g1_1)<-FALSE
				l$Illumina_Lumi$Coexpression_Network<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			display()	
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
			}
		if(length(which(x=="Nimblegen"))!=0){
			galert("Please wait while generating Co-expression network",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			myData_Sel<-use.data.matrix_Nimblegen2.m[rownames(DE_N),]
			myData_Sel<-t(myData_Sel)
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
			adjMat_N<<-adjMat
			myGraph_N<<-as(adjMat_N,"graphNEL")
			plot(myGraph_N,nodeAttrs=makeNodeAttrs(myGraph_N,fontsize=18,fillcolor="grey"))
			if(length(myGraph_N)!=0){
				visible(g1_1)<-FALSE
				l$Nimblegen$Coexpression_Network<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			display()	
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
  			}
		if(length(which(x=="Series_Matrix"))!=0){
			galert("Please wait while generating Co-expression network",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			myData_Sel<-use.data.matrixNorm.m[rownames(DE_S),]
			myData_Sel<-t(myData_Sel)
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
			adjMat_S<<-adjMat
			myGraph_S<<-as(adjMat_S,"graphNEL")
			plot(myGraph_S,nodeAttrs=makeNodeAttrs(myGraph_S,fontsize=18,fillcolor="grey"))
			if(length(myGraph_S)!=0){
				visible(g1_1)<-FALSE
				l$Series_Matrix$Coexpression_Network<<-list()
				tr<<-gtree(offspring=tree,container=g1_1)
				size(tr)<-c(300,400)
				visible(g1_1)<-TRUE
				
				}
			display()	
			galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Done"
			}
		if(length(which(x=="Online_Data"))!=0){
			galert("Please wait while generating Co-expression network",title="Miscellaneous",delay=5,parent=c(600,400))
			svalue(sb)<-"Working... Plz wait..."
			Sys.sleep(1)

			myData_Sel<-use.data.matrix_onlineNorm.m[rownames(DE_O),]
			myData_Sel<-t(myData_Sel)
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
			adjMat_O<<-adjMat
			myGraph_O<<-as(adjMat_O,"graphNEL")
			plot(myGraph_O,nodeAttrs=makeNodeAttrs(myGraph_O,fontsize=18,fillcolor="grey"))
			if(length(myGraph_O)!=0){
				visible(g1_1)<-FALSE
				l$Online_Data$Coexpression_Network<<-list()
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
			gmessage("Plz select the data for generating Co-expression network","Select Data")
			}
			dispose(w_dge)
		},container=gp2_dge,anchor=c(1,-1))
	visible(w_dge)<-TRUE
	}

