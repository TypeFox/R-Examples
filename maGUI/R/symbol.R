symbol<-function(h,...){
pre_rs<-c("GPL32","GPL33","GPL34","GPL71","GPL72","GPL74","GPL75","GPL76","GPL77","GPL78","GPL79","GPL80","GPL81","GPL82","GPL83","GPL85","GPL86","GPL87","GPL88","GPL89","GPL90","GPL91","GPL92","GPL93","GPL94","GPL95","GPL96","GPL97","GPL98","GPL99","GPL100","GPL101","GPL198","GPL199","GPL200","GPL201","GPL339","GPL340","GPL341","GPL342","GPL570","GPL571","GPL886","GPL887","GPL1261","GPL1318","GPL1319","GPL1322","GPL1352","GPL1355","GPL1708","GPL2112","GPL2529","GPL2891","GPL2898","GPL3154","GPL3213","GPL3533","GPL3738","GPL3921","GPL3979","GPL4032","GPL4191","GPL5689","GPL6097","GPL6102","GPL6244","GPL6947","GPL8300","GPL8490","GPL10558","GPL11532","GPL13497","GPL13534","GPL13667","GPL15380","GPL15396","GPL17897","mgu74a","mgu74b","mgu74c","ag","drosgenome1","hcg110","mu11ksuba","mu11ksubb","mu19ksuba","mu19ksubb","mu19ksubc","hu6800","mgu74av2","mgu74bv2","mgu74cv2","rgu34a","rgu34b","rgu34c","rnu34","rtu34","ygs98","hgu95av2","hgu95b","hgu95c","hgu95d","hgu95e","hgu133a","hgu133b","hu35ksuba","hu35ksubb","hu35ksubc","hu35ksubd","ath1121501","ecoli2","celegans","hgfocus","moe430a","mouse4302","rae230a","rae230b","hgu133plus2","hgu133a2","hgug4111a","hgug4110b","mouse430a2","xenopuslaevis","zebrafish","drosophila2","u133x3p","rat2302","hgug4112a","bovine","yeast2","h20kcod","adme16cod","ecoli2","chicken","porcine","canine2","hthgu133a","canine","","h10kcod","hgug4100a","illuminaHumanv1","illuminaHumanv2","hugene10sttranscriptcluster","illuminaHumanv3","hgu95av2","IlluminaHumanMethylation27k","illuminaHumanv4","hugene11sttranscriptcluster","HsAgilentDesign026652","IlluminaHumanMethylation450k","hgu219","GGHumanMethCancerPanelv1","hthgu133b","hthgu133a")

	rs=NULL;
	con=NULL;
	folder_Ann=NULL;
	rm(rs,con,folder_Ann)
	loc=gconfirm("Gene symbol needs annotation from GEOmetadb.sqlite database",icon="info")
	if(loc==TRUE)
	{
		choose_folder()
		folder_Ann<<-folderchoose
		folderchoose=NULL
		if(length(folder_Ann)!=0)
		{
			setwd(folder_Ann)
			err=NULL;
			try(con<-dbConnect(SQLite(),'GEOmetadb.sqlite'),silent=TRUE)
			try(geo_tables<-dbListTables(con),silent=TRUE)
			try(rs<<-dbGetQuery(con,'select gpl,bioc_package from gpl'),silent=TRUE)
			try(colnames(rs)<-c("gpl","bioc_package"),silent=TRUE)
		}
	}	
	else
	{
		loc2<-gconfirm("Annotation from inbuilt database",icon="question")
		if(loc2==TRUE)
		{
			rs<<-matrix(pre_rs,ncol=2)
			try(colnames(rs)<<-c("gpl","bioc_package"),silent=TRUE)
		}
	}
	if(length(rs)!=0)
	{
		try(({dat2Affy.m<-dat2Affy.m;datAgOne2.m<-datAgOne2.m;datAgTwo2.m<-datAgTwo2.m;
		datIllBA2.m2<-datIllBA2.m2;lumi_NQ.m<-lumi_NQ.m;data.matrix_Nimblegen2.m<-data.matrix_Nimblegen2.m;
		data.matrixNorm.m<-data.matrixNorm.m;data.matrix_onlineNorm.m<-data.matrix_onlineNorm.m;
		use.dat2Affy.m<-use.dat2Affy.m;use.datAgOne2.m<-use.datAgOne2.m;use.datAgTwo2.m<-use.datAgTwo2.m;
		use.datIllBA2.m2<-use.datIllBA2.m2;use.lumi_NQ.m<-use.lumi_NQ.m;
		use.data.matrix_Nimblegen2.m<-use.data.matrix_Nimblegen2.m;
		use.data.matrixNorm.m<-use.data.matrixNorm.m;use.data.matrix_onlineNorm.m<-use.data.matrix_onlineNorm.m;
		l<-l;tree<-tree;
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

	genes_Affy=NULL;genes_Ag1=NULL;genes_Ag2=NULL;genes_Il_B=NULL;genes_Il_L=NULL;
	genes_N=NULL;genes_S=NULL;genes_O=NULL;	x=NULL

	rm(genes_Affy,genes_Ag1,genes_Ag2,genes_Il_B,genes_Il_L,genes_N,genes_S,genes_O)
	
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
			if(length(ann_Affy)!=0 && is.na(ann_Affy)==FALSE)
			{
				
				try(({
				if(grep("GPL[0-9]",ann_Affy)==TRUE)
				{
					ann_Affy<<-rs[which(rs[,1]==ann_Affy),2]
					
					if(length(ann_Affy)==0)
					{
						gmessage("Annotation package not available for this platform",icon="warning")
						ann_Affy=NULL
					}
				}
				else
				{
					if(grep("GPL[0-9]",ann_Affy)==FALSE)
					{
						ann_Affy<<-ann_Affy
					}
				}
				}),silent=TRUE)
			}
			else
			{
				ginput("Please provide GPL name",icon="question",handler=function(h,...)
				{
					ann_Affy<<-h$input
					
					if(grep("GPL[0-9]",ann_Affy)==TRUE)
					{
						ann_Affy<<-rs[which(rs[,1]==ann_Affy),2]
						
						if(length(ann_Affy)==0)
						{
							gmessage("Annotation package not available for this platform",icon="warning")
							ann_Affy=NULL
							}
						}
					}
					)
				}	
			if(length(ann_Affy)!=0)
			{
				galert("Please wait while Identifing Symbols",title="Miscellaneous",delay=5,parent=c(600,400))
				svalue(sb)<-"Working... Plz wait..."
				Sys.sleep(1)
				Identifier<-rownames(dat2Affy.m)
				db<-annPkgName(ann_Affy)
				err<-try(library(db,character.only=TRUE),silent=TRUE)
				if(length(grep("Error",err))!=0)
				{
					source("http://bioconductor.org/biocLite.R")
					biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
					library(db,character.only=TRUE)
				}

				Symbol<-getSYMBOL(rownames(dat2Affy.m),db)
				genes_Affy<<-cbind(Identifier,Symbol)
				genes_Affy<<-as.data.frame(genes_Affy)
				if(length(genes_Affy)!=0){
					visible(g1_1)<-FALSE
					l$Affymetrix$Identifier_Symbol<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					
					}
				display()	
				galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
				svalue(sb)<-"Done"
				}
			}
		if(length(which(x=="Agilent_OneColor"))!=0){
			if(length(ann_Ag1)!=0 && is.na(ann_Ag1)==FALSE)
			{
				try(({
				if(grep("GPL[0-9]",ann_Ag1)==TRUE)
				{
					ann_Ag1<<-rs[which(rs[,1]==ann_Ag1),2]
					if(length(ann_Ag1)==0)
					{
						gmessage("Annotation package not available for this platform",icon="warning")
						ann_Ag1=NULL
					}
				}
				else
				{
					if(grep("GPL[0-9]",ann_Ag1)==FALSE)
					{
						ann_Ag1<<-ann_Ag1
					}
				}
				}),silent=TRUE)
			}
			else
			{
				ginput("Please provide GPL name",icon="question",handler=function(h,...)
				{
					x_a<-h$input
					ann_Ag1<<-x_a
					
					if(grep("GPL[0-9]",ann_Ag1)==TRUE)
					{
						ann_Ag1<<-rs[which(rs[,1]==ann_Ag1),2]
						
						if(length(ann_Ag1)==0)
						{
							gmessage("Annotation package not available for this platform",icon="warning")
							ann_Ag1<<-NULL
							}
						}
					}
					)
				}
			if(length(ann_Ag1)!=0)
			{
				galert("Please wait while Identifing Symbols",title="Miscellaneous",delay=5,parent=c(600,400))
				svalue(sb)<-"Working... Plz wait..."
				Sys.sleep(1)
				Identifier<-rownames(datAgOne2.m)
				db<-annPkgName(ann_Ag1)
				err<-try(library(db,character.only=TRUE),silent=TRUE)
				if(length(grep("Error",err))!=0)
				{
					source("http://bioconductor.org/biocLite.R")
					biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
					library(db,character.only=TRUE)
				}
								
				Symbol<-getSYMBOL(rownames(datAgOne2.m),db)
				genes_Ag1<<-cbind(Identifier,Symbol)
				genes_Ag1<<-as.data.frame(genes_Ag1)
				if(length(genes_Ag1)!=0){
					visible(g1_1)<-FALSE
					l$Agilent_OneColor$Identifier_Symbol<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					
					}
				display()	
				galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
				svalue(sb)<-"Done"
				}
			}
		if(length(which(x=="Agilent_TwoColor"))!=0){
			if(length(ann_Ag2)!=0 && is.na(ann_Ag2)==FALSE)
			{
				try(({
				if(grep("GPL[0-9]",ann_Ag2)==TRUE)
				{
					ann_Ag2<<-rs[which(rs[,1]==ann_Ag2),2]
					if(length(ann_Ag2)==0)
					{
						gmessage("Annotation package not available for this platform",icon="warning")
						ann_Ag2=NULL
					}
				}
				else
				{
					if(grep("GPL[0-9]",ann_Ag2)==FALSE)
					{
						ann_Ag2<<-ann_Ag2
					}
				}
				}),silent=TRUE)
			}
			else
			{
				ginput("Please provide GPL name",icon="question",handler=function(h,...)
				{
					ann_Ag2<<-h$input
					
					if(grep("GPL[0-9]",ann_Ag2)==TRUE)
					{
						ann_Ag2<<-rs[which(rs[,1]==ann_Ag2),2]
						
						if(length(ann_Ag2)==0)
						{
							gmessage("Annotation package not available for this platform",icon="warning")
							ann_Ag2=NULL
							}
						}
					}
					)
				}
			if(length(ann_Ag2)!=0)
			{
				Identifier<-rownames(datAgTwo2.m)
				db<-annPkgName(ann_Ag2)
				err<-try(library(db,character.only=TRUE),silent=TRUE)
				if(length(grep("Error",err))!=0)
				{
					source("http://bioconductor.org/biocLite.R")
					biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
					library(db,character.only=TRUE)
				}
								
				Symbol<-getSYMBOL(rownames(datAgTwo2.m),db)
				genes_Ag2<<-cbind(Identifier,Symbol)
				genes_Ag2<<-as.data.frame(genes_Ag2)
				if(length(genes_Ag2)!=0){
					visible(g1_1)<-FALSE
					l$Agilent_TwoColor$Identifier_Symbol<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					
					}
				display()	
				galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
				svalue(sb)<-"Done"
				}
			}
		if(length(which(x=="Illumina_Beadarray"))!=0){
			if(length(ann_Il_B)!=0 && is.na(ann_Il_B)==FALSE)
			{
				try(({
				if(grep("GPL[0-9]",ann_Il_B)==TRUE)
				{
					ann_Il_B<<-rs[which(rs[,1]==ann_Il_B),2]
					if(length(ann_Il_B)==0)
					{
						gmessage("Annotation package not available for this platform",icon="warning")
						ann_Il_B=NULL
					}
				}
				else
				{
					if(grep("GPL[0-9]",ann_Il_B)==FALSE)
					{
						ann_Il_B<<-ann_Il_B
					}
				}
				}),silent=TRUE)
			}
			else
			{
				ginput("Please provide GPL name",icon="question",handler=function(h,...)
				{
					ann_Il_B<<-h$input
					
					if(grep("GPL[0-9]",ann_Il_B)==TRUE)
					{
						ann_Il_B<<-rs[which(rs[,1]==ann_Il_B),2]
						
						if(length(ann_Il_B)==0)
						{
							gmessage("Annotation package not available for this platform",icon="warning")
							ann_Il_B=NULL
							}
						}
					}
					)
				}
			if(length(ann_Il_B)!=0)
			{
				Identifier<-rownames(datIllBA2.m2)
				db<-annPkgName(ann_Il_B)
				err<-try(library(db,character.only=TRUE),silent=TRUE)
				if(length(grep("Error",err))!=0)
				{
					source("http://bioconductor.org/biocLite.R")
					biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
					library(db,character.only=TRUE)
				}
								
				Symbol<-getSYMBOL(rownames(datIllBA2.m2),db)
				genes_Il_B<<-cbind(Identifier,Symbol)
				genes_Il_B<<-as.data.frame(genes_Il_B)
				if(length(genes_Il_B)!=0){
					visible(g1_1)<-FALSE
					l$Illumina_Beadarray$Identifier_Symbol<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					
					}
				display()	
				galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
				svalue(sb)<-"Done"
				}
			}
		if(length(which(x=="Illumina_Lumi"))!=0){
			if(length(ann_Il_L)!=0 && is.na(ann_Il_L)==FALSE)
			{
				try(({
				if(grep("GPL[0-9]",ann_Il_L)==TRUE)
				{
					ann_Il_L<<-rs[which(rs[,1]==ann_Il_L),2]
					if(length(ann_Il_L)==0)
					{
						gmessage("Annotation package not available for this platform",icon="warning")
						ann_Il_L=NULL
					}
				}
				else
				{
					if(grep("GPL[0-9]",ann_Il_L)==FALSE)
					{
						ann_Il_L<<-ann_Il_L
					}
				}
				}),silent=TRUE)
			}
			else
			{
				ginput("Please provide GPL name",icon="question",handler=function(h,...)
				{
					ann_Il_L<<-h$input
					
					if(grep("GPL[0-9]",ann_Il_L)==TRUE)
					{
						ann_Il_L<<-rs[which(rs[,1]==ann_Il_L),2]
						
						if(length(ann_Il_L)==0)
						{
							gmessage("Annotation package not available for this platform",icon="warning")
							ann_Il_L=NULL
							}
						}
					}
					)
				}
			if(length(ann_Il_L))
			{
				Identifier<-rownames(lumi_NQ.m)
				db<-annPkgName(ann_Il_L)
				err<-try(library(db,character.only=TRUE),silent=TRUE)
				if(length(grep("Error",err))!=0)
				{
					source("http://bioconductor.org/biocLite.R")
					biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
					library(db,character.only=TRUE)
				}
								
				Symbol<-getSYMBOL(rownames(lumi_NQ.m),db)
				genes_Il_L<<-cbind(Identifier,Symbol)
				genes_Il_L<<-as.data.frame(genes_Il_L)
				if(length(genes_Il_L)!=0){
					visible(g1_1)<-FALSE
					l$Illumina_Lumi$Identifier_Symbol<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					
					}
				display()	
				galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
				svalue(sb)<-"Done"
				}
			}
		if(length(which(x=="Nimblegen"))!=0){
			if(length(ann_N)!=0 && is.na(ann_N)==FALSE)
			{
				try(({
				if(grep("GPL[0-9]",ann_N)==TRUE)
				{
					ann_N<<-rs[which(rs[,1]==ann_N),2]
					if(length(ann_N)==0)
					{
						gmessage("Annotation package not available for this platform",icon="warning")
						ann_N=NULL
					}
				}
				else
				{
					if(grep("GPL[0-9]",ann_N)==FALSE)
					{
						ann_N<<-ann_N
					}
				}
				}),silent=TRUE)
			}
			else
			{
				ginput("Please provide GPL name",icon="question",handler=function(h,...)
				{
					ann_N<<-h$input
					
					if(grep("GPL[0-9]",ann_N)==TRUE)
					{
						ann_N<<-rs[which(rs[,1]==ann_N),2]
						
						if(length(ann_N)==0)
						{
							gmessage("Annotation package not available for this platform",icon="warning")
							ann_N=NULL
							}
						}
					}
					)
				}
			if(length(ann_N))
			{
				Identifier<-rownames(data.matrix_Nimblegen2.m)
				db<-annPkgName(ann_N)
				err<-try(library(db,character.only=TRUE),silent=TRUE)
				if(length(grep("Error",err))!=0)
				{
					source("http://bioconductor.org/biocLite.R")
					biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
					library(db,character.only=TRUE)
				}
								
				Symbol<-getSYMBOL(rownames(data.matrix_Nimblegen2.m),db)
				genes_N<<-cbind(Identifier,Symbol)
				genes_N<<-as.data.frame(genes_N)
				if(length(genes_N)!=0){
					visible(g1_1)<-FALSE
					l$Nimblegen$Identifier_Symbol<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					
					}
				display()	
				galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
				svalue(sb)<-"Done"
  				}
  			}
		if(length(which(x=="Series_Matrix"))!=0){
			if(length(ann_S)!=0 && is.na(ann_S)==FALSE)
			{
				try(({
				if(grep("GPL[0-9]",ann_S)==TRUE)
				{
					ann_S<<-rs[which(rs[,1]==ann_S),2]
					if(length(ann_S)==0)
					{
						gmessage("Annotation package not available for this platform",icon="warning")
						ann_S=NULL
					}
				}
				else
				{
					if(grep("GPL[0-9]",ann_S)==FALSE)
					{
						ann_S<<-ann_S
					}
				}
				}),silent=TRUE)
			}
			else
			{
				ginput("Please provide GPL name",icon="question",handler=function(h,...)
				{
					ann_S<<-h$input
					
					if(grep("GPL[0-9]",ann_S)==TRUE)
					{
						ann_S<<-rs[which(rs[,1]==ann_S),2]
						
						if(length(ann_S)==0)
						{
							gmessage("Annotation package not available for this platform",icon="warning")
							ann_S=NULL
							}
						}
					}
					)
				}
			if(length(ann_S))
			{
				Identifier<-rownames(data.matrixNorm.m)
				db<-annPkgName(ann_S)
				err<-try(library(db,character.only=TRUE),silent=TRUE)
				if(length(grep("Error",err))!=0)
				{
					source("http://bioconductor.org/biocLite.R")
					biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
					library(db,character.only=TRUE)
				}
								
				Symbol<-getSYMBOL(rownames(data.matrixNorm.m),db)
				genes_S<<-cbind(Identifier,Symbol)
				genes_S<<-as.data.frame(genes_S)
				if(length(genes_S)!=0){
					visible(g1_1)<-FALSE
					l$Series_Matrix$Identifier_Symbol<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					
					}
				display()	
				galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
				svalue(sb)<-"Done"
				}
			}
		if(length(which(x=="Online_Data"))!=0){
			if(length(ann_O)!=0 && is.na(ann_O)==FALSE)
			{
				try(({
				if(grep("GPL[0-9]",ann_O)==TRUE)
				{
					ann_O<<-rs[which(rs[,1]==ann_O),2]
					if(length(ann_O)==0)
					{
						gmessage("Annotation package not available for this platform",icon="warning")
						ann_O=NULL
					}
				}
				else
				{
					if(grep("GPL[0-9]",ann_O)==FALSE)
					{
						ann_O<<-ann_O
					}
				}
				}),silent=TRUE)
			}
			else
			{
				ginput("Please provide GPL name",icon="question",handler=function(h,...)
				{
					ann_O<<-h$input
					
					if(grep("GPL[0-9]",ann_O)==TRUE)
					{
						ann_O<<-rs[which(rs[,1]==ann_O),2]
						
						if(length(ann_O)==0)
						{
							gmessage("Annotation package not available for this platform",icon="warning")
							ann_O=NULL
							}
						}
					}
					)
				}
			if(length(ann_O))
			{
				Identifier<-rownames(data.matrix_onlineNorm.m)
				db<-annPkgName(ann_O)
				err<-try(library(db,character.only=TRUE),silent=TRUE)
				if(length(grep("Error",err))!=0)
				{
					source("http://bioconductor.org/biocLite.R")
					biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
					library(db,character.only=TRUE)
				}
								
				Symbol<-getSYMBOL(rownames(data.matrix_onlineNorm.m),db)
				genes_O<<-cbind(Identifier,Symbol)
				genes_O<<-as.data.frame(genes_O)
				if(length(genes_O)!=0){
					visible(g1_1)<-FALSE
					l$Online_Data$Identifier_Symbol<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,400)
					visible(g1_1)<-TRUE
					
					}
				display()	
				galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
				svalue(sb)<-"Done"
				}
			}
		dispose(w_dge)
		}else{
			gmessage("Plz select the data for annotation","Select Data")
			}
		},container=gp2_dge,anchor=c(1,-1))
	visible(w_dge)<-TRUE
	}
}
