gsta<-function(h,...){
pre_rs<-c("GPL32","GPL33","GPL34","GPL71","GPL72","GPL74","GPL75","GPL76","GPL77","GPL78","GPL79","GPL80","GPL81","GPL82","GPL83","GPL85","GPL86","GPL87","GPL88","GPL89","GPL90","GPL91","GPL92","GPL93","GPL94","GPL95","GPL96","GPL97","GPL98","GPL99","GPL100","GPL101","GPL198","GPL199","GPL200","GPL201","GPL339","GPL340","GPL341","GPL342","GPL570","GPL571","GPL886","GPL887","GPL1261","GPL1318","GPL1319","GPL1322","GPL1352","GPL1355","GPL1708","GPL2112","GPL2529","GPL2891","GPL2898","GPL3154","GPL3213","GPL3533","GPL3738","GPL3921","GPL3979","GPL4032","GPL4191","GPL5689","GPL6097","GPL6102","GPL6244","GPL6947","GPL8300","GPL8490","GPL10558","GPL11532","GPL13497","GPL13534","GPL13667","GPL15380","GPL15396","GPL17897","mgu74a","mgu74b","mgu74c","ag","drosgenome1","hcg110","mu11ksuba","mu11ksubb","mu19ksuba","mu19ksubb","mu19ksubc","hu6800","mgu74av2","mgu74bv2","mgu74cv2","rgu34a","rgu34b","rgu34c","rnu34","rtu34","ygs98","hgu95av2","hgu95b","hgu95c","hgu95d","hgu95e","hgu133a","hgu133b","hu35ksuba","hu35ksubb","hu35ksubc","hu35ksubd","ath1121501","ecoli2","celegans","hgfocus","moe430a","mouse4302","rae230a","rae230b","hgu133plus2","hgu133a2","hgug4111a","hgug4110b","mouse430a2","xenopuslaevis","zebrafish","drosophila2","u133x3p","rat2302","hgug4112a","bovine","yeast2","h20kcod","adme16cod","ecoli2","chicken","porcine","canine2","hthgu133a","canine","","h10kcod","hgug4100a","illuminaHumanv1","illuminaHumanv2","hugene10sttranscriptcluster","illuminaHumanv3","hgu95av2","IlluminaHumanMethylation27k","illuminaHumanv4","hugene11sttranscriptcluster","HsAgilentDesign026652","IlluminaHumanMethylation450k","hgu219","GGHumanMethCancerPanelv1","hthgu133b","hthgu133a")

	rs<-NULL;
	con<-NULL;
	folder_Ann<-NULL;
	rm(rs,con,folder_Ann)
	loc=gconfirm("Annotation needs GEOmetadb.sqlite database",icon="info")
	if(loc==TRUE)
	{
		choose_folder()
		folder_Ann<<-folderchoose
		folderchoose=NULL
		if(length(folder_Ann)!=0)
		{
			setwd(folder_Ann)
			err=NULL;
			try(con<<-dbConnect(SQLite(),'GEOmetadb.sqlite'),silent=TRUE)
			try(geo_tables<-dbListTables(con),silent=TRUE)
			try(rs<<-dbGetQuery(con,'select gpl,bioc_package from gpl'),silent=TRUE)
			try(colnames(rs)<<-c("gpl","bioc_package"),silent=TRUE)
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

	use.dat2Affy.m=NULL;use.datAgOne2.m=NULL;use.datAgTwo2.m=NULL;use.datIllBA2.m2=NULL;use.lumi_NQ.m=NULL;
	use.data.matrix_Nimblegen2.m=NULL;use.data.matrixNorm.m=NULL;use.data.matrix_onlineNorm.m=NULL;
	groups_go=NULL;groups_kegg=NULL;
	use.use.dat2Affy.m=NULL;use.use.datAgOne2.m=NULL;use.use.datAgTwo2.m=NULL;use.use.datIllBA2.m2=NULL;use.use.lumi_NQ.m=NULL;
	use.use.data.matrix_Nimblegen2.m=NULL;use.use.data.matrixNorm.m=NULL;use.use.data.matrix_onlineNorm.m=NULL;
	
	GOtable.outBP_Affy=NULL;GOtable.outBP_Ag1=NULL;GOtable.outBP_Ag2=NULL;GOtable.outBP_Il_B=NULL;GOtable.outBP_Il_L=NULL;
	GOtable.outBP_N=NULL;GOtable.outBP_S=NULL;GOtable.outBP_O=NULL;
	GOtable.outMF_Affy=NULL;GOtable.outMF_Ag1=NULL;GOtable.outMF_Ag2=NULL;GOtable.outMF_Il_B=NULL;GOtable.outMF_Il_L=NULL;
	GOtable.outMF_N=NULL;GOtable.outMF_S=NULL;GOtable.outMF_O=NULL;
	GOtable.outCC_Affy=NULL;GOtable.outCC_Ag1=NULL;GOtable.outCC_Ag2=NULL;GOtable.outCC_Il_B=NULL;GOtable.outCC_Il_L=NULL;
	GOtable.outCC_N=NULL;GOtable.outCC_S=NULL;GOtable.outCC_O=NULL;
	KEGGtable.out_Affy=NULL;KEGGtable.out_Ag1=NULL;KEGGtable.out_Ag2=NULL;KEGGtable.out_Il_B=NULL;KEGGtable.out_Il_L=NULL;
	KEGGtable.out_N=NULL;KEGGtable.out_S=NULL;KEGGtable.out_O=NULL;

	rm(use.dat2Affy.m,use.datAgOne2.m,use.datAgTwo2.m,use.datIllBA2.m2,use.lumi_NQ.m,
	use.data.matrix_Nimblegen2.m,use.data.matrixNorm.m,use.data.matrix_onlineNorm.m,
	groups_go,groups_kegg,
	use.use.dat2Affy.m,use.use.datAgOne2.m,use.use.datAgTwo2.m,use.use.datIllBA2.m2,use.use.lumi_NQ.m,
	use.use.data.matrix_Nimblegen2.m,use.use.data.matrixNorm.m,use.use.data.matrix_onlineNorm.m,
	GOtable.outBP_Affy,GOtable.outBP_Ag1,GOtable.outBP_Ag2,GOtable.outBP_Il_B,GOtable.outBP_Il_L,
	GOtable.outBP_N,GOtable.outBP_S,GOtable.outBP_O,
	GOtable.outMF_Affy,GOtable.outMF_Ag1,GOtable.outMF_Ag2,GOtable.outMF_Il_B,GOtable.outMF_Il_L,
	GOtable.outMF_N,GOtable.outMF_S,GOtable.outMF_O,
	GOtable.outCC_Affy,GOtable.outCC_Ag1,GOtable.outCC_Ag2,GOtable.outCC_Il_B,GOtable.outCC_Il_L,
	GOtable.outCC_N,GOtable.outCC_S,GOtable.outCC_O,
	KEGGtable.out_Affy,KEGGtable.out_Ag1,KEGGtable.out_Ag2,KEGGtable.out_Il_B,KEGGtable.out_Il_L,
	KEGGtable.out_N,KEGGtable.out_S,KEGGtable.out_O)

		genes=NULL;
		x=NULL
		f<-function(h,...){
			x<<-svalue(h$obj)
			}
	gsea_xx=NULL
	gsea_methods=c("GO categories","KEGG pathways")
	gsea_w<-gwindow("Select an analysis method",width=300,height=150)
	gsea_gp<-ggroup(container=gsea_w,horizontal=FALSE)
	gsea_x<-gtable(gsea_methods,chosencol=1,container=gsea_gp)
	size(gsea_x)=c(200,100)
	gsea_gp2<-ggroup(container=gsea_gp,width=30,height=15,horizontal=TRUE)
	addHandlerClicked(gsea_x,handler=function(h,...){
		gsea_x2<-svalue(h$obj)
		gsea_xx<<-gsea_x2
		}
	)
	addSpring(gsea_gp2)
		
	gsea_y<-gbutton("CANCEL",border=TRUE,handler=function(h,...){
		dispose(gsea_w)
		},container=gsea_gp2,anchor=c(1,-1)
	)
	gsea_y2<-gbutton("OK",border=TRUE,handler=function(h,...){
		dispose(gsea_w)
		if(gsea_xx=="GO categories")
		{
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
					dispose(w_dge);
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
						if(length(ann_Affy)!=0 && is.na(ann_Affy)==FALSE)
						{	
							
							galert("Please wait while GSTA GO analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_Affy)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}
							
							allprobes_id<-paste(ann_Affy,"GO2ALLPROBES",sep="")
							go2allprobes<-get(allprobes_id)
							go<-as.list(go2allprobes)
							design=NULL
							use.dat2Affy.m<-dat2Affy.m;
							xf=dim(use.dat2Affy.m)
							if(xf[2]%%2==0)
							{
								use.dat2Affy.m<<-dat2Affy.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.dat2Affy.m<<-dat2Affy.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_go<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.dat2Affy.m<<-dat2Affy.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							}

							galert("Performing GSTA GO Biological Process",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goBP<-gtGO(groups_go,t(use.dat2Affy.m),annotation=db,ontology="BP",sort=TRUE);
							BP_length<-length(test.goBP);
							test.goBP2<-test.goBP[1:BP_length,];
							table.outBP<-data.frame(pathway=names(test.goBP2),pvalue=p.value(test.goBP2));
							table.outBP2<-data.frame(table.outBP,Description=alias(test.goBP2))
							genes=NULL;
							for(i in 1:length(table.outBP2[,1])){genes[i]<-list(go[[rownames(table.outBP2)[i]]])}
							table.outBP2$genes<-genes
							table.outBP2$genes<-as.character(gsub("\n","",genes))
							GOtable.outBP_Affy<<-table.outBP2[order(table.outBP2$pvalue),]
							
							galert("Performing GSTA GO Molecular Function",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goMF<-gtGO(groups_go,t(use.dat2Affy.m),annotation=db,ontology="MF",sort=TRUE)
							MF_length<-length(test.goMF)
							test.goMF2<-test.goMF[1:MF_length,]
							table.outMF<-data.frame(pathway=names(test.goMF2),pvalue=p.value(test.goMF2))
							table.outMF2<-data.frame(table.outMF,Description=alias(test.goMF2))
							genes=NULL;
							for(i in 1:length(table.outMF2[,1])){genes[i]<-list(go[[rownames(table.outMF2)[i]]])}	
							table.outMF2$genes<-genes
							table.outMF2$genes<-as.character(gsub("\n","",genes))
							GOtable.outMF_Affy<<-table.outMF2[order(table.outMF2$pvalue),]
							
							galert("Performing GSTA GO Cellular Component",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goCC<-gtGO(groups_go,t(use.dat2Affy.m),annotation=db,ontology="CC",sort=TRUE)
							CC_length<-length(test.goCC)
							test.goCC2<-test.goCC[1:CC_length,]
							table.outCC<-data.frame(pathway=names(test.goCC2),pvalue=p.value(test.goCC2))
							table.outCC2<-data.frame(table.outCC,Description=alias(test.goCC2))
							genes=NULL;
							for(i in 1:length(table.outCC2[,1])){genes[i]<-list(go[[rownames(table.outCC2)[i]]])}
							table.outCC2$genes<-genes
							table.outCC2$genes<-as.character(gsub("\n","",genes))
							GOtable.outCC_Affy<<-table.outCC2[order(table.outCC2$pvalue),]

							if(length(GOtable.outBP_Affy)!=0){
								visible(g1_1)<-FALSE
								l$Affymetrix$GSTA_GO$BP<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outMF_Affy)!=0){
								visible(g1_1)<-FALSE
								l$Affymetrix$GSTA_GO$MF<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outCC_Affy)!=0){
								visible(g1_1)<-FALSE
								l$Affymetrix$GSTA_GO$CC<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
						}
					}
					if(length(which(x=="Agilent_OneColor"))!=0)
					{
						
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
								ann_Ag1<<-h$input
								
								if(grep("GPL[0-9]",ann_Ag1)==TRUE)
								{
									ann_Ag1<<-rs[which(rs[,1]==ann_Ag1),2]
									
									if(length(ann_Ag1)==0)
									{
										gmessage("Annotation package not available for this platform",icon="warning")
										ann_Ag1=NULL
									}
								}
							}
							)	
						}	
						if(length(ann_Ag1)!=0 && is.na(ann_Ag1)==FALSE)
						{	
							
							galert("Please wait while GSTA GO analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}
							
							allprobes_id<-paste(ann_Ag1,"GO2ALLPROBES",sep="")
							go2allprobes<-get(allprobes_id)
							go<-as.list(go2allprobes)
							design=NULL
							use.datAgOne2.m<-datAgOne2.m
							xf=dim(use.datAgOne2.m)
							if(xf[2]%%2==0)
							{
								use.datAgOne2.m<<-datAgOne2.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.datAgOne2.m<<-datAgOne2.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_go<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.datAgOne2.m<<-datAgOne2.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							}
	
							galert("Performing GSTA GO Biological Process",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goBP<-gtGO(groups_go,t(use.datAgOne2.m),annotation=db,ontology="BP",sort=TRUE)
							BP_length<-length(test.goBP)
							test.goBP2<-test.goBP[1:BP_length,]
							table.outBP<-data.frame(pathway=names(test.goBP2),pvalue=p.value(test.goBP2))
							table.outBP2<-data.frame(table.outBP,Description=alias(test.goBP2))
							genes=NULL;
							for(i in 1:length(table.outBP2[,1])){genes[i]<-list(go[[rownames(table.outBP2)[i]]])}
							table.outBP2$genes<-genes
							table.outBP2$genes<-as.character(gsub("\n","",genes))
							GOtable.outBP_Ag1<<-table.outBP2[order(table.outBP2$pvalue),]
							
							galert("Performing GSTA GO Molecular Function",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goMF<-gtGO(groups_go,t(use.datAgOne2.m),annotation=db,ontology="MF",sort=TRUE)
							MF_length<-length(test.goMF)
							test.goMF2<-test.goMF[1:MF_length,]
							table.outMF<-data.frame(pathway=names(test.goMF2),pvalue=p.value(test.goMF2))
							table.outMF2<-data.frame(table.outMF,Description=alias(test.goMF2))
							genes=NULL;
							for(i in 1:length(table.outMF2[,1])){genes[i]<-list(go[[rownames(table.outMF2)[i]]])}
							table.outMF2$genes<-genes
							table.outMF2$genes<-as.character(gsub("\n","",genes))
							GOtable.outMF_Ag1<<-table.outMF2[order(table.outMF2$pvalue),]
							
							galert("Performing GSTA GO Cellular Component",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goCC<-gtGO(groups_go,t(use.datAgOne2.m),annotation=db,ontology="CC",sort=TRUE)
							CC_length<-length(test.goCC)
							test.goCC2<-test.goCC[1:CC_length,]
							table.outCC<-data.frame(pathway=names(test.goCC2),pvalue=p.value(test.goCC2))
							table.outCC2<-data.frame(table.outCC,Description=alias(test.goCC2))
							genes=NULL;
							for(i in 1:length(table.outCC2[,1])){genes[i]<-list(go[[rownames(table.outCC2)[i]]])}
							table.outCC2$genes<-genes
							table.outCC2$genes<-as.character(gsub("\n","",genes))
							GOtable.outCC_Ag2<<-table.outCC2[order(table.outCC2$pvalue),]
	
							if(length(GOtable.outBP_Ag1)!=0){
								visible(g1_1)<-FALSE
								l$Agilent_OneColor$GSTA_GO$BP<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outMF_Ag1)!=0){
								visible(g1_1)<-FALSE
								l$Agilent_OneColor$GSTA_GO$MF<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outCC_Ag1)!=0){
								visible(g1_1)<-FALSE
								l$Agilent_OneColor$GSTA_GO$CC<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
						}
					}
					if(length(which(x=="Agilent_TwoColor"))!=0)
					{
						
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
						if(length(ann_Ag2)!=0 && is.na(ann_Ag2)==FALSE)
						{	
							
							galert("Please wait while GSTA GO analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)

							db<-annPkgName(ann_Ag2)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}

							allprobes_id<-paste(ann_Ag2,"GO2ALLPROBES",sep="")
							go2allprobes<-get(allprobes_id)
							go<-as.list(go2allprobes)
							design=NULL
							use.datAgTwo2.m<-datAgTwo2.m
							xf=dim(use.datAgTwo2.m)
							if(xf[2]%%2==0)
							{
								use.datAgTwo2.m<<-datAgTwo2.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.datAgTwo2.m<<-datAgTwo2.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_go<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.datAgTwo2.m<<-datAgTwo2.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							}
							
							galert("Performing GSTA GO Biological Process",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goBP<-gtGO(groups_go,t(use.datAgTwo2.m),annotation=db,ontology="BP",sort=TRUE)
							BP_length<-length(test.goBP)
							test.goBP2<-test.goBP[1:BP_length,]
							table.outBP<-data.frame(pathway=names(test.goBP2),pvalue=p.value(test.goBP2))
							table.outBP2<-data.frame(table.outBP,Description=alias(test.goBP2))
							genes=NULL;
							for(i in 1:length(table.outBP2[,1])){genes[i]<-list(go[[rownames(table.outBP2)[i]]])}
							table.outBP2$genes<-genes
							table.outBP2$genes<-as.character(gsub("\n","",genes))
							GOtable.outBP_Ag2<<-table.outBP2[order(table.outBP2$pvalue),]
							
							galert("Performing GSTA GO Molecular Function",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goMF<-gtGO(groups_go,t(use.datAgTwo2.m),annotation=db,ontology="MF",sort=TRUE)
							MF_length<-length(test.goMF)
							test.goMF2<-test.goMF[1:MF_length,]
							table.outMF<-data.frame(pathway=names(test.goMF2),pvalue=p.value(test.goMF2))
							table.outMF2<-data.frame(table.outMF,Description=alias(test.goMF2))
							genes=NULL;
							for(i in 1:length(table.outMF2[,1])){genes[i]<-list(go[[rownames(table.outMF2)[i]]])}
							table.outMF2$genes<-genes
							table.outMF2$genes<-as.character(gsub("\n","",genes))
							GOtable.outMF_Ag2<<-table.outMF2[order(table.outMF2$pvalue),]
							
							galert("Performing GSTA GO Cellular Component",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goCC<-gtGO(groups_go,t(use.datAgTwo2.m),annotation=db,ontology="CC",sort=TRUE)
							CC_length<-length(test.goCC)
							test.goCC2<-test.goCC[1:CC_length,]
							table.outCC<-data.frame(pathway=names(test.goCC2),pvalue=p.value(test.goCC2))
							table.outCC2<-data.frame(table.outCC,Description=alias(test.goCC2))
							genes=NULL;
							for(i in 1:length(table.outCC2[,1])){genes[i]<-list(go[[rownames(table.outCC2)[i]]])}
							table.outCC2$genes<-genes
							table.outCC2$genes<-as.character(gsub("\n","",genes))
							GOtable.outCC_Ag2<<-table.outCC2[order(table.outCC2$pvalue),]
	
							if(length(GOtable.outBP_Ag2)!=0){
								visible(g1_1)<-FALSE
								l$Agilent_TwoColor$GSTA_GO$BP<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outMF_Ag2)!=0){
								visible(g1_1)<-FALSE
								l$Agilent_TwoColor$GSTA_GO$MF<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outCC_Ag2)!=0){
								visible(g1_1)<-FALSE
								l$Agilent_TwoColor$GSTA_GO$CC<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
						}
					}
					if(length(which(x=="Illumina_Beadarray"))!=0)
					{
						
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
						if(length(ann_Il_B)!=0 && is.na(ann_Il_B)==FALSE)
						{	
							
							galert("Please wait while GSTA GO analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_Il_B)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}

								
							allprobes_id<-paste(ann_Il_B,"GO2ALLPROBES",sep="")
							go2allprobes<-get(allprobes_id)
							go<-as.list(go2allprobes)
							design=NULL
							use.datIllBA2.m2<-datIllBA2.m2
							xf=dim(use.datIllBA2.m2)
							if(xf[2]%%2==0)
							{
								use.datIllBA2.m2<<-datIllBA2.m2
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.datIllBA2.m2<<-datIllBA2.m2
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_go<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.datIllBA2.m2<<-datIllBA2.m2[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							}
							
							galert("Performing GSTA GO Biological Process",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goBP<-gtGO(groups_go,t(use.datIllBA2.m2),annotation=db,ontology="BP",sort=TRUE)
							BP_length<-length(test.goBP)
							test.goBP2<-test.goBP[1:BP_length,]
							table.outBP<-data.frame(pathway=names(test.goBP2),pvalue=p.value(test.goBP2))
							table.outBP2<-data.frame(table.outBP,Description=alias(test.goBP2))
							genes=NULL;
							for(i in 1:length(table.outBP2[,1])){genes[i]<-list(go[[rownames(table.outBP2)[i]]])}
							table.outBP2$genes<-genes
							table.outBP2$genes<-as.character(gsub("\n","",genes))
							GOtable.outBP_Il_B<<-table.outBP2[order(table.outBP2$pvalue),]
							
							galert("Performing GSTA GO Molecular Function",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goMF<-gtGO(groups_go,t(use.datIllBA2.m2),annotation=db,ontology="MF",sort=TRUE)
							MF_length<-length(test.goMF)
							test.goMF2<-test.goMF[1:MF_length,]
							table.outMF<-data.frame(pathway=names(test.goMF2),pvalue=p.value(test.goMF2))
							table.outMF2<-data.frame(table.outMF,Description=alias(test.goMF2))
							genes=NULL;
							for(i in 1:length(table.outMF2[,1])){genes[i]<-list(go[[rownames(table.outMF2)[i]]])}
							table.outMF2$genes<-genes
							table.outMF2$genes<-as.character(gsub("\n","",genes))
							GOtable.outMF_Il_B<<-table.outMF2[order(table.outMF2$pvalue),]
							
							galert("Performing GSTA GO Cellular Component",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goCC<-gtGO(groups_go,t(use.datIllBA2.m2),annotation=db,ontology="CC",sort=TRUE)
							CC_length<-length(test.goCC)
							test.goCC2<-test.goCC[1:CC_length,]
							table.outCC<-data.frame(pathway=names(test.goCC2),pvalue=p.value(test.goCC2))
							table.outCC2<-data.frame(table.outCC,Description=alias(test.goCC2))
							genes=NULL;
							for(i in 1:length(table.outCC2[,1])){genes[i]<-list(go[[rownames(table.outCC2)[i]]])}
							table.outCC2$genes<-genes
							table.outCC2$genes<-as.character(gsub("\n","",genes))
							GOtable.outCC_Il_B<<-table.outCC2[order(table.outCC2$pvalue),]

							if(length(GOtable.outBP_Il_B)!=0){
								visible(g1_1)<-FALSE
								l$Illumina_Beadarray$GSTA_GO$BP<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outMF_Il_B)!=0){
								visible(g1_1)<-FALSE
								l$Illumina_Beadarray$GSTA_GO$MF<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outCC_Il_B)!=0){
								visible(g1_1)<-FALSE
								l$Illumina_Beadarray$GSTA_GO$CC<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
						}
					}
					if(length(which(x=="Illumina_Lumi"))!=0)
					{
						
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
						if(length(ann_Il_L)!=0 && is.na(ann_Il_L)==FALSE)
						{	
							
							galert("Please wait while GSTA GO analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_Il_L)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}

							allprobes_id<-paste(ann_Il_L,"GO2ALLPROBES",sep="")
							go2allprobes<-get(allprobes_id)
							go<-as.list(go2allprobes)
							design=NULL
							use.lumi_NQ.m<-lumi_NQ.m
							xf=dim(use.lumi_NQ.m)
							if(xf[2]%%2==0)
							{
								use.lumi_NQ.m<<-lumi_NQ.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.lumi_NQ.m<<-lumi_NQ.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_go<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.lumi_NQ.m<<-lumi_NQ.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							}
							
							galert("Performing GSTA GO Biological Process",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goBP<-gtGO(groups_go,t(use.lumi_NQ.m),annotation=db,ontology="BP",sort=TRUE)
							BP_length<-length(test.goBP)
							test.goBP2<-test.goBP[1:BP_length,]
							table.outBP<-data.frame(pathway=names(test.goBP2),pvalue=p.value(test.goBP2))
							table.outBP2<-data.frame(table.outBP,Description=alias(test.goBP2))
							genes=NULL;
							for(i in 1:length(table.outBP2[,1])){genes[i]<-list(go[[rownames(table.outBP2)[i]]])}
							table.outBP2$genes<-genes
							table.outBP2$genes<-as.character(gsub("\n","",genes))
							GOtable.outBP_Il_L<<-table.outBP2[order(table.outBP2$pvalue),]
							
							galert("Performing GSTA GO Molecular Function",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goMF<-gtGO(groups_go,t(use.lumi_NQ.m),annotation=db,ontology="MF",sort=TRUE)
							MF_length<-length(test.goMF)
							test.goMF2<-test.goMF[1:MF_length,]
							table.outMF<-data.frame(pathway=names(test.goMF2),pvalue=p.value(test.goMF2))
							table.outMF2<-data.frame(table.outMF,Description=alias(test.goMF2))
							genes=NULL;
							for(i in 1:length(table.outMF2[,1])){genes[i]<-list(go[[rownames(table.outMF2)[i]]])}	
							table.outMF2$genes<-genes
							table.outMF2$genes<-as.character(gsub("\n","",genes))
							GOtable.outMF_Il_L<<-table.outMF2[order(table.outMF2$pvalue),]
							
							galert("Performing GSTA GO Cellular Component",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goCC<-gtGO(groups_go,t(use.lumi_NQ.m),annotation=db,ontology="CC",sort=TRUE)
							CC_length<-length(test.goCC)
							test.goCC2<-test.goCC[1:CC_length,]
							table.outCC<-data.frame(pathway=names(test.goCC2),pvalue=p.value(test.goCC2))
							table.outCC2<-data.frame(table.outCC,Description=alias(test.goCC2))
							genes=NULL;
							for(i in 1:length(table.outCC2[,1])){genes[i]<-list(go[[rownames(table.outCC2)[i]]])}
							table.outCC2$genes<-genes
							table.outCC2$genes<-as.character(gsub("\n","",genes))
							GOtable.outCC_Il_L<<-table.outCC2[order(table.outCC2$pvalue),]

							if(length(GOtable.outBP_Il_L)!=0){
								visible(g1_1)<-FALSE
								l$Illumina_Lumi$GSTA_GO$BP<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outMF_Il_L)!=0){
								visible(g1_1)<-FALSE
								l$Illumina_Lumi$GSTA_GO$MF<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outCC_Il_L)!=0){
								visible(g1_1)<-FALSE
								l$Illumina_Lumi$GSTA_GO$CC<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
						}
					}
					if(length(which(x=="Nimblegen"))!=0)
					{
						
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
						if(length(ann_N)!=0 && is.na(ann_N)==FALSE)
						{	
							
							galert("Please wait while GSTA GO analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_N)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}
								
							allprobes_id<-paste(ann_N,"GO2ALLPROBES",sep="")
							go2allprobes<-get(allprobes_id)
							go<-as.list(go2allprobes)
							design=NULL
							use.data.matrix_Nimblegen2.m<-data.matrix_Nimblegen2.m
							xf=dim(use.data.matrix_Nimblegen2.m)
							if(xf[2]%%2==0)
							{
								use.data.matrix_Nimblegen2.m<<-data.matrix_Nimblegen2.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.data.matrix_Nimblegen2.m<<-data.matrix_Nimblegen2.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_go<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.data.matrix_Nimblegen2.m<<-data.matrix_Nimblegen2.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							}
	
							galert("Performing GSTA GO Biological Process",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goBP<-gtGO(groups_go,t(use.data.matrix_Nimblegen2.m),annotation=db,ontology="BP",sort=TRUE)
							BP_length<-length(test.goBP)
							test.goBP2<-test.goBP[1:BP_length,]
							table.outBP<-data.frame(pathway=names(test.goBP2),pvalue=p.value(test.goBP2))
							table.outBP2<-data.frame(table.outBP,Description=alias(test.goBP2))
							genes=NULL;
							for(i in 1:length(table.outBP2[,1])){genes[i]<-list(go[[rownames(table.outBP2)[i]]])}
							table.outBP2$genes<-genes
							table.outBP2$genes<-as.character(gsub("\n","",genes))
							GOtable.outBP_N<<-table.outBP2[order(table.outBP2$pvalue),]
							
							galert("Performing GSTA GO Molecular Function",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goMF<-gtGO(groups_go,t(use.data.matrix_Nimblegen2.m),annotation=db,ontology="MF",sort=TRUE)
							MF_length<-length(test.goMF)
							test.goMF2<-test.goMF[1:MF_length,]
							table.outMF<-data.frame(pathway=names(test.goMF2),pvalue=p.value(test.goMF2))
							table.outMF2<-data.frame(table.outMF,Description=alias(test.goMF2))
							genes=NULL;
							for(i in 1:length(table.outMF2[,1])){genes[i]<-list(go[[rownames(table.outMF2)[i]]])}	
							table.outMF2$genes<-genes
							table.outMF2$genes<-as.character(gsub("\n","",genes))
							GOtable.outMF_N<<-table.outMF2[order(table.outMF2$pvalue),]
							
							galert("Performing GSTA GO Cellular Component",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goCC<-gtGO(groups_go,t(use.data.matrix_Nimblegen2.m),annotation=db,ontology="CC",sort=TRUE)
							CC_length<-length(test.goCC)
							test.goCC2<-test.goCC[1:CC_length,]
							table.outCC<-data.frame(pathway=names(test.goCC2),pvalue=p.value(test.goCC2))
							table.outCC2<-data.frame(table.outCC,Description=alias(test.goCC2))
							genes=NULL;
							for(i in 1:length(table.outCC2[,1])){genes[i]<-list(go[[rownames(table.outCC2)[i]]])}
							table.outCC2$genes<-genes
							table.outCC2$genes<-as.character(gsub("\n","",genes))
							GOtable.outCC_N<<-table.outCC2[order(table.outCC2$pvalue),]
	
							if(length(GOtable.outBP_N)!=0){
								visible(g1_1)<-FALSE
								l$Nimblegen$GSTA_GO$BP<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outMF_N)!=0){
								visible(g1_1)<-FALSE
								l$Nimblegen$GSTA_GO$MF<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outCC_N)!=0){
								visible(g1_1)<-FALSE
								l$Nimblegen$GSTA_GO$CC<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
  						}
  					}
					if(length(which(x=="Series_Matrix"))!=0)
					{
						
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
						if(length(ann_S)!=0 && is.na(ann_S)==FALSE)
						{	
							
							galert("Please wait while GSTA GO analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_S)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}
								
							allprobes_id<-paste(ann_S,"GO2ALLPROBES",sep="")
							go2allprobes<-get(allprobes_id)
							go<-as.list(go2allprobes)
							design=NULL
							use.data.matrixNorm.m<-data.matrixNorm.m
							xf=dim(use.data.matrixNorm.m)
							if(xf[2]%%2==0)
							{
								use.data.matrixNorm.m<<-data.matrixNorm.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.data.matrixNorm.m<<-data.matrixNorm.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_go<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.data.matrixNorm.m<<-data.matrixNorm.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							}
							
							galert("Performing GSTA GO Biological Process",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goBP<-gtGO(groups_go,t(use.data.matrixNorm.m),annotation=db,ontology="BP",sort=TRUE)
							BP_length<-length(test.goBP)
							test.goBP2<-test.goBP[1:BP_length,]
							table.outBP<-data.frame(pathway=names(test.goBP2),pvalue=p.value(test.goBP2))
							table.outBP2<-data.frame(table.outBP,Description=alias(test.goBP2))
							genes=NULL;
							for(i in 1:length(table.outBP2[,1])){genes[i]<-list(go[[rownames(table.outBP2)[i]]])}
							table.outBP2$genes<-genes
							table.outBP2$genes<-as.character(gsub("\n","",genes))
							GOtable.outBP_S<<-table.outBP2[order(table.outBP2$pvalue),]
							
							galert("Performing GSTA GO Molecular Function",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goMF<-gtGO(groups_go,t(use.data.matrixNorm.m),annotation=db,ontology="MF",sort=TRUE)
							MF_length<-length(test.goMF)
							test.goMF2<-test.goMF[1:MF_length,]
							table.outMF<-data.frame(pathway=names(test.goMF2),pvalue=p.value(test.goMF2))
							table.outMF2<-data.frame(table.outMF,Description=alias(test.goMF2))
							genes=NULL;
							for(i in 1:length(table.outMF2[,1])){genes[i]<-list(go[[rownames(table.outMF2)[i]]])}	
							table.outMF2$genes<-genes
							table.outMF2$genes<-as.character(gsub("\n","",genes))
							GOtable.outMF_S<<-table.outMF2[order(table.outMF2$pvalue),]
							
							galert("Performing GSTA GO Cellular Component",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goCC<-gtGO(groups_go,t(use.data.matrixNorm.m),annotation=db,ontology="CC",sort=TRUE)
							CC_length<-length(test.goCC)
							test.goCC2<-test.goCC[1:CC_length,]
							table.outCC<-data.frame(pathway=names(test.goCC2),pvalue=p.value(test.goCC2))
							table.outCC2<-data.frame(table.outCC,Description=alias(test.goCC2))
							genes=NULL;
							for(i in 1:length(table.outCC2[,1])){genes[i]<-list(go[[rownames(table.outCC2)[i]]])}
							table.outCC2$genes<-genes
							table.outCC2$genes<-as.character(gsub("\n","",genes))
							GOtable.outCC_S<<-table.outCC2[order(table.outCC2$pvalue),]
	
							if(length(GOtable.outBP_S)!=0){
								visible(g1_1)<-FALSE
								l$Series_Matrix$GSTA_GO$BP<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outMF_S)!=0){
								visible(g1_1)<-FALSE
								l$Series_Matrix$GSTA_GO$MF<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outCC_S)!=0){
								visible(g1_1)<-FALSE
								l$Series_Matrix$GSTA_GO$CC<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
						}
					}
					if(length(which(x=="Online_Data"))!=0)
					{
						
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
						if(length(ann_O)!=0 && is.na(ann_O)==FALSE)
						{	
							
							galert("Please wait while GSTA GO analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_O)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}

							allprobes_id<-paste(ann_O,"GO2ALLPROBES",sep="")
							go2allprobes<-get(allprobes_id)
							go<-as.list(go2allprobes)
							design=NULL
							use.data.matrix_onlineNorm.m<-data.matrix_onlineNorm.m
							xf=dim(use.data.matrix_onlineNorm.m)
							if(xf[2]%%2==0)
							{
								use.use.data.matrix_onlineNorm.m<<-use.data.matrix_onlineNorm.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.use.data.matrix_onlineNorm.m<<-use.data.matrix_onlineNorm.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_go<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.use.data.matrix_onlineNorm.m<<-use.data.matrix_onlineNorm.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_go<<-c(rep(group,yf))
								}
							}
							
							galert("Performing GSTA GO Biological Process",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goBP<-gtGO(groups_go,t(use.data.matrix_onlineNorm.m),annotation=db,ontology="BP",sort=TRUE)
							BP_length<-length(test.goBP)
							test.goBP2<-test.goBP[1:BP_length,]
							table.outBP<-data.frame(pathway=names(test.goBP2),pvalue=p.value(test.goBP2))
							table.outBP2<-data.frame(table.outBP,Description=alias(test.goBP2))
							genes=NULL;
							for(i in 1:length(table.outBP2[,1])){genes[i]<-list(go[[rownames(table.outBP2)[i]]])}
							table.outBP2$genes<-genes
							table.outBP2$genes<-as.character(gsub("\n","",genes))
							GOtable.outBP_O<<-table.outBP2[order(table.outBP2$pvalue),]
							
							galert("Performing GSTA GO Molecular Function",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goMF<-gtGO(groups_go,t(use.data.matrix_onlineNorm.m),annotation=db,ontology="MF",sort=TRUE)
							MF_length<-length(test.goMF)
							test.goMF2<-test.goMF[1:MF_length,]
							table.outMF<-data.frame(pathway=names(test.goMF2),pvalue=p.value(test.goMF2))
							table.outMF2<-data.frame(table.outMF,Description=alias(test.goMF2))
							genes=NULL;
							for(i in 1:length(table.outMF2[,1])){genes[i]<-list(go[[rownames(table.outMF2)[i]]])}	
							table.outMF2$genes<-genes
							table.outMF2$genes<-as.character(gsub("\n","",genes))
							GOtable.outMF_O<<-table.outMF2[order(table.outMF2$pvalue),]
							
							galert("Performing GSTA GO Cellular Component",title="GSTA Gene Ontology",delay=5,parent=c(600,400))
							test.goCC<-gtGO(groups_go,t(use.data.matrix_onlineNorm.m),annotation=db,ontology="CC",sort=TRUE)
							CC_length<-length(test.goCC)
							test.goCC2<-test.goCC[1:CC_length,]
							table.outCC<-data.frame(pathway=names(test.goCC2),pvalue=p.value(test.goCC2))
							table.outCC2<-data.frame(table.outCC,Description=alias(test.goCC2))
							genes=NULL;
							for(i in 1:length(table.outCC2[,1])){genes[i]<-list(go[[rownames(table.outCC2)[i]]])}
							table.outCC2$genes<-genes
							table.outCC2$genes<-as.character(gsub("\n","",genes))
							GOtable.outCC_O<<-table.outCC2[order(table.outCC2$pvalue),]

							if(length(GOtable.outBP_O)!=0){
								visible(g1_1)<-FALSE
								l$Online_Data$GSTA_GO$BP<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outMF_O)!=0){
								visible(g1_1)<-FALSE
								l$Online_Data$GSTA_GO$MF<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							if(length(GOtable.outCC_O)!=0){
								visible(g1_1)<-FALSE
								l$Online_Data$GSTA_GO$CC<<-list()
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
				}
				else{
					gmessage("Plz select the data for GSEA","Select Data")
					}
				
				},container=gp2_dge,anchor=c(1,-1)
				)
				visible(w_dge)<-TRUE
			}
		else if(gsea_xx=="KEGG pathways")
		{
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
					dispose(w_dge);
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
						if(length(ann_Affy)!=0 && is.na(ann_Affy)==FALSE)
						{	
							
							galert("Please wait while GSTA KEGG analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_Affy)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}
							path2probe_id<-paste(ann_Affy,"PATH2PROBE",sep="")
							pathway2probe<-get(path2probe_id)
							kegg<-as.list(pathway2probe)
							design=NULL
							use.dat2Affy.m<-dat2Affy.m
							xf=dim(use.dat2Affy.m)
							if(xf[2]%%2==0)
							{
								use.use.dat2Affy.m<<-use.dat2Affy.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.use.dat2Affy.m<<-use.dat2Affy.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_kegg<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.use.dat2Affy.m<<-use.dat2Affy.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							}
							test.kegg<-gtKEGG(groups_kegg,t(use.dat2Affy.m),annotation=db,sort=TRUE);
							table.out<-data.frame(pathway=names(test.kegg),pvalue=p.value(test.kegg));
							table.out<-data.frame(table.out,Description=alias(test.kegg));
							genes=NULL;
							for(i in 1:length(table.out[,1])){genes[i]<-list(kegg[[rownames(table.out)[i]]])}
							table.out$genes<-genes
							table.out$genes<-as.character(gsub("\n","",genes))	
							KEGGtable.out_Affy<<-table.out[order(table.out$pvalue),];
							
							if(length(KEGGtable.out_Affy)!=0){
								visible(g1_1)<-FALSE
								l$Affymetrix$GSTA_KEGG<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
						}
					}
					if(length(which(x=="Agilent_OneColor"))!=0)
					{
						
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
								ann_Ag1<<-h$input
								
								if(grep("GPL[0-9]",ann_Ag1)==TRUE)
								{
									ann_Ag1<<-rs[which(rs[,1]==ann_Ag1),2]
									
									if(length(ann_Ag1)==0)
									{
										gmessage("Annotation package not available for this platform",icon="warning")
										ann_Ag1=NULL
									}
								}
							}
							)	
						}	
						if(length(ann_Ag1)!=0 && is.na(ann_Ag1)==FALSE)
						{	
							
							galert("Please wait while GSTA KEGG analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_Ag1)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}

							path2probe_id<-paste(ann_Ag1,"PATH2PROBE",sep="")
							pathway2probe<-get(path2probe_id)
							kegg<-as.list(pathway2probe)
							design=NULL
							use.datAgOne2.m<-datAgOne2.m
							xf=dim(use.datAgOne2.m)
							if(xf[2]%%2==0)
							{
								use.use.datAgOne2.m<<-use.datAgOne2.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.use.datAgOne2.m<<-use.datAgOne2.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_kegg<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.use.datAgOne2.m<<-use.datAgOne2.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							}
							test.kegg<-gtKEGG(groups_kegg,t(use.datAgOne2.m),annotation=db,sort=TRUE)
							table.out<-data.frame(pathway=names(test.kegg),pvalue=p.value(test.kegg))
							table.out<-data.frame(table.out,Description=alias(test.kegg));
							genes=NULL;
							for(i in 1:length(table.out[,1])){genes[i]<-list(kegg[[rownames(table.out)[i]]])}
							table.out$genes<-genes
							table.out$genes<-as.character(gsub("\n","",genes))
							KEGGtable.out_Ag1<<-table.out[order(table.out$pvalue),]
							
							if(length(KEGGtable.out_Ag1)!=0){
								visible(g1_1)<-FALSE
								l$Agilent_OneColor$GSTA_KEGG<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
						}
					}	
					if(length(which(x=="Agilent_TwoColor"))!=0)
					{
						
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
						if(length(ann_Ag2)!=0 && is.na(ann_Ag2)==FALSE)
						{	
							
							galert("Please wait while GSTA KEGG analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_Ag2)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}
					
							path2probe_id<-paste(ann_Ag2,"PATH2PROBE",sep="")
							pathway2probe<-get(path2probe_id)
							kegg<-as.list(pathway2probe)
							design=NULL
							use.datAgTwo2.m<-datAgTwo2.m
							xf=dim(use.datAgTwo2.m)
							if(xf[2]%%2==0)
							{
								use.use.datAgTwo2.m<<-use.datAgTwo2.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.use.datAgTwo2.m<<-use.datAgTwo2.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_kegg<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.use.datAgTwo2.m<<-use.datAgTwo2.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							}
							test.kegg<-gtKEGG(groups_kegg,t(use.datAgTwo2.m),annotation=db,sort=TRUE);
							table.out<-data.frame(pathway=names(test.kegg),pvalue=p.value(test.kegg));
							table.out<-data.frame(table.out,Description=alias(test.kegg));
							genes=NULL;
							for(i in 1:length(table.out[,1])){genes[i]<-list(kegg[[rownames(table.out)[i]]])}
							table.out$genes<-genes
							table.out$genes<-as.character(gsub("\n","",genes))	
							KEGGtable.out_Ag2<<-table.out[order(table.out$pvalue),];
							
							if(length(KEGGtable.out_Ag2)!=0){
								visible(g1_1)<-FALSE
								l$Agilent_TwoColor$GSTA_KEGG<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
							}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
						}
					}
					if(length(which(x=="Illumina_Beadarray"))!=0)
					{
						
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
						if(length(ann_Il_B)!=0 && is.na(ann_Il_B)==FALSE)
						{	
							
							galert("Please wait while GSTA KEGG analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_Il_B)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}

							path2probe_id<-paste(ann_Il_B,"PATH2PROBE",sep="")
							pathway2probe<-get(path2probe_id)
							kegg<-as.list(pathway2probe)
							design=NULL
							use.datIllBA2.m2<-datIllBA2.m2
							xf=dim(use.datIllBA2.m2)
							if(xf[2]%%2==0)
							{
								use.use.datIllBA2.m2<<-use.datIllBA2.m2
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.use.datIllBA2.m2<<-use.datIllBA2.m2
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_kegg<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.use.datIllBA2.m2<<-use.datIllBA2.m2[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							}
							test.kegg<-gtKEGG(groups_kegg,t(use.datIllBA2.m2),annotation=db,sort=TRUE);
							table.out<-data.frame(pathway=names(test.kegg),pvalue=p.value(test.kegg));
							table.out<-data.frame(table.out,Description=alias(test.kegg));
							genes=NULL;
							for(i in 1:length(table.out[,1])){genes[i]<-list(kegg[[rownames(table.out)[i]]])}
							table.out$genes<-genes
							table.out$genes<-as.character(gsub("\n","",genes))	
							KEGGtable.out_Il_B<<-table.out[order(table.out$pvalue),];
							
							if(length(KEGGtable.out_Il_B)!=0){
								visible(g1_1)<-FALSE
								l$Illumina_Beadarray$GSTA_KEGG<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
						}
					}
					if(length(which(x=="Illumina_Lumi"))!=0)
					{
						
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
						if(length(ann_Il_L)!=0 && is.na(ann_Il_L)==FALSE)
						{	
							
							galert("Please wait while GSTA KEGG analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_Il_L)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}

							path2probe_id<-paste(ann_Il_L,"PATH2PROBE",sep="")
							pathway2probe<-get(path2probe_id)
							kegg<-as.list(pathway2probe)
							design=NULL
							use.lumi_NQ.m<-lumi_NQ.m
							xf=dim(use.lumi_NQ.m)
							if(xf[2]%%2==0)
							{
								use.use.lumi_NQ.m<<-use.lumi_NQ.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.use.lumi_NQ.m<<-use.lumi_NQ.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_kegg<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.use.lumi_NQ.m<<-use.lumi_NQ.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							}
							test.kegg<-gtKEGG(groups_kegg,t(use.lumi_NQ.m),annotation=db,sort=TRUE);
							table.out<-data.frame(pathway=names(test.kegg),pvalue=p.value(test.kegg));
							table.out<-data.frame(table.out,Description=alias(test.kegg));
							genes=NULL;
							for(i in 1:length(table.out[,1])){genes[i]<-list(kegg[[rownames(table.out)[i]]])}
							table.out$genes<-genes
							table.out$genes<-as.character(gsub("\n","",genes))	
							KEGGtable.out_Il_L<<-table.out[order(table.out$pvalue),];
						
							if(length(KEGGtable.out_Il_L)!=0){
								visible(g1_1)<-FALSE
								l$Illumina_Lumi$GSTA_KEGG<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
						}
					}
					if(length(which(x=="Nimblegen"))!=0)
					{
						
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
						if(length(ann_N)!=0 && is.na(ann_N)==FALSE)
						{	
							
							galert("Please wait while GSTA KEGG analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_N)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}

							path2probe_id<-paste(ann_N,"PATH2PROBE",sep="")
							pathway2probe<-get(path2probe_id)
							kegg<-as.list(pathway2probe)
							design=NULL
							use.data.matrix_Nimblegen2.m<-data.matrix_Nimblegen2.m
							xf=dim(use.data.matrix_Nimblegen2.m)
							if(xf[2]%%2==0)
							{
								use.use.data.matrix_Nimblegen2.m<<-use.data.matrix_Nimblegen2.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.use.data.matrix_Nimblegen2.m<<-use.data.matrix_Nimblegen2.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_kegg<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.use.data.matrix_Nimblegen2.m<<-use.data.matrix_Nimblegen2.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							}
							test.kegg<-gtKEGG(groups_kegg,t(use.data.matrix_Nimblegen2.m),annotation=db,sort=TRUE);
							table.out<-data.frame(pathway=names(test.kegg),pvalue=p.value(test.kegg));
							table.out<-data.frame(table.out,Description=alias(test.kegg));
							genes=NULL;
							for(i in 1:length(table.out[,1])){genes[i]<-list(kegg[[rownames(table.out)[i]]])}
							table.out$genes<-genes
							table.out$genes<-as.character(gsub("\n","",genes))	
							KEGGtable.out_N<<-table.out[order(table.out$pvalue),];
							
							if(length(KEGGtable.out_N)!=0){
								visible(g1_1)<-FALSE
								l$Nimblegen$GSTA_KEGG<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
  						}
  					}
					if(length(which(x=="Series_Matrix"))!=0)
					{
						
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
						if(length(ann_S)!=0 && is.na(ann_S)==FALSE)
						{	
							
							galert("Please wait while GSTA KEGG analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_S)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}

							path2probe_id<-paste(ann_S,"PATH2PROBE",sep="")
							pathway2probe<-get(path2probe_id)
							kegg<-as.list(pathway2probe)
							design=NULL
							use.data.matrixNorm.m<-data.matrixNorm.m
							xf=dim(use.data.matrixNorm.m)
							if(xf[2]%%2==0)
							{
								use.use.data.matrixNorm.m<<-use.data.matrixNorm.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.use.data.matrixNorm.m<<-use.data.matrixNorm.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_kegg<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.use.data.matrixNorm.m<<-use.data.matrixNorm.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							}
							test.kegg<-gtKEGG(groups_kegg,t(use.data.matrixNorm.m),annotation=db,sort=TRUE);
							table.out<-data.frame(pathway=names(test.kegg),pvalue=p.value(test.kegg));
							table.out<-data.frame(table.out,Description=alias(test.kegg));
							genes=NULL;
							for(i in 1:length(table.out[,1])){genes[i]<-list(kegg[[rownames(table.out)[i]]])}
							table.out$genes<-genes
							table.out$genes<-as.character(gsub("\n","",genes))	
							KEGGtable.out_S<<-table.out[order(table.out$pvalue),];
							
							if(length(KEGGtable.out_S)!=0){
								visible(g1_1)<-FALSE
								l$Series_Matrix$GSTA_KEGG<<-list()
								tr<<-gtree(offspring=tree,container=g1_1)
								size(tr)<-c(300,400)
								visible(g1_1)<-TRUE
								
								}
							display()	
							galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Done"
						}
					}
					if(length(which(x=="Online_Data"))!=0)
					{
						
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
						if(length(ann_O)!=0 && is.na(ann_O)==FALSE)
						{	
							
							galert("Please wait while GSTA KEGG analysis",title="Miscellaneous",delay=5,parent=c(600,400))
							svalue(sb)<-"Working... Plz wait..."
							Sys.sleep(1)
	
							db<-annPkgName(ann_O)
							err<-try(library(db,character.only=TRUE),silent=TRUE)
							if(length(grep("Error",err))!=0)
							{
								source("http://bioconductor.org/biocLite.R")
								biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
								library(db,character.only=TRUE)
							}

							path2probe_id<-paste(ann_O,"PATH2PROBE",sep="")
							pathway2probe<-get(path2probe_id)
							kegg<-as.list(pathway2probe)
							design=NULL
							use.data.matrix_onlineNorm.m<-data.matrix_onlineNorm.m
							xf=dim(use.data.matrix_onlineNorm.m)
							if(xf[2]%%2==0)
							{
								use.use.data.matrix_onlineNorm.m<<-use.data.matrix_onlineNorm.m
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							} else if(xf[2]%%3==0){
								use.use.data.matrix_onlineNorm.m<<-use.data.matrix_onlineNorm.m
								yf=xf[2]/3
								for(i in 1:yf)
								{
									group<-c(0,1,1)
									groups_kegg<<-c(rep(group,yf))
								}
								} else if(xf[2]%%3!=0 && xf[2]%%2!=0){
								use.use.data.matrix_onlineNorm.m<<-use.data.matrix_onlineNorm.m[,-xf[2]]
								yf=xf[2]/2
								for(i in 1:yf)
								{
									group<-c(0,1)
									groups_kegg<<-c(rep(group,yf))
								}
							}
							test.kegg<-gtKEGG(groups_kegg,t(use.data.matrix_onlineNorm.m),annotation=db,sort=TRUE);
							table.out<-data.frame(pathway=names(test.kegg),pvalue=p.value(test.kegg));
							table.out<-data.frame(table.out,Description=alias(test.kegg));
							genes=NULL;
							for(i in 1:length(table.out[,1])){genes[i]<-list(kegg[[rownames(table.out)[i]]])}
							table.out$genes<-genes
							table.out$genes<-as.character(gsub("\n","",genes))	
							KEGGtable.out_O<<-table.out[order(table.out$pvalue),];
							
							if(length(KEGGtable.out_O)!=0){
								visible(g1_1)<-FALSE
								l$Online_Data$GSTA_KEGG<<-list()
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
					galert("Done",title="Miscellaneous",delay=5,parent=c(600,400))
					svalue(sb)<-"Done"
				} else{
					gmessage("Plz select the data for GSEA","Select Data")
				}
					
			},container=gp2_dge,anchor=c(1,-1)
			)
			visible(w_dge)<-TRUE
			}
		},container=gsea_gp2,anchor=c(1,-1)
	)	
	}	
}		
