illlumi<-function(h,...){
	lumi_NQ=NULL;
	lumi_NQ.m=NULL;
	use.lumi_NQ.m=NULL;
	xf=NULL;
	lumi_NQ.f=NULL;
	design_Il_L=NULL;
	ttx=NULL;
	lumi_NQ.s=NULL;
	DE_Il_L=NULL;
	pca_Il_L=NULL;
	sample.dist_Il_L=NULL;
	sample.clust_Il_L=NULL;
	use_DE_Il_L=NULL;
	Clas_Il_L=NULL;
	design=NULL;

	try(({folder_Il_L<-folder_Il_L;lumi_data<-lumi_data;l<-l;tree<-tree;
	}),silent=TRUE)

	setwd(folder_Il_L)	
	galert("Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
#	rm(lumi_NQ)
#	lumi_NQ<<-lumiExpresso(lumi_data,QC.evaluation=TRUE)
#	summary(lumi_NQ,'QC')
#	rm(lumi_NQ.m)
#	lumi_NQ.m<<-lumi_NQ
#	lumi_NQ<<-lumiExpresso(lumi_data,QC.evaluation=TRUE)
	lumi_NQ<<-lumiExpresso(lumi_data,QC.evaluation=TRUE)
	summary(lumi_NQ,'QC')
	rm(lumi_NQ.m)
	lumi_NQ.m<<-lumiExpresso(lumi_data,QC.evaluation=TRUE)
	lumi_NQ.m<<-as.matrix(lumi_NQ.m)
	rownames(lumi_NQ.m)<<-as.character(rownames(lumi_NQ.m))
	try(if(length(rownames(lumi_NQ.m))==0)rownames(lumi_NQ.m)<<-1:length(lumi_NQ.m[,1]),silent=TRUE)
	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	plot(lumi_NQ.m,what="boxplot")
	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
	galert("Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(use.lumi_NQ.m)
	use.lumi_NQ.m<<-lumi_NQ.m
	try(({
	rm(xf)
	xf<<-dim(use.lumi_NQ.m)
	if(xf[2]%%3!=0 && xf[2]%%2!=0){
		use.lumi_NQ.m<<-lumi_NQ.m[,-xf[2]]
		xf<<-dim(use.lumi_NQ.m)
	}
	if(xf[2]%%2==0)
	{
		yf=xf[2]/2
		groups<-c(rep("C",yf),rep("T",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.lumi_NQ.m,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(lumi_NQ.f)
		lumi_NQ.f<<-fit2
		rm(design_Il_L)
		design_Il_L<<-design
#		print("Two-group Specific Filtering")
	}else
	{
		yf=xf[2]/3
		groups<-c(rep("C",yf),rep("T1",yf),rep("T2",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.lumi_NQ.m,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(lumi_NQ.f)
		lumi_NQ.f<<-fit2
		rm(design_Il_L)
		design_Il_L<<-design
#		print("Three-group Specific Filtering")
	}
	}),silent=TRUE)
	
	if(length(lumi_NQ.f)==0)
	{
		try(({
			rsd<-rowSds(use.lumi_NQ.m)
			i<-rsd>=2
			dat.f<-use.lumi_NQ.m[i,]
			fit<-lmFit(dat.f)
			yy<-try(toptable(fit,coef=2),silent=TRUE)
			if(length(grep("Error in",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
			rm(lumi_NQ.f)
			lumi_NQ.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
			
		if(length(lumi_NQ.f)==0)
		{
			try(({
				ff<-pOverA(A=1,p=0.5)
				i<-genefilter(use.lumi_NQ.m,ff)
				dat.fo<-use.lumi_NQ.m[i,]
				i<-genefilter(-use.lumi_NQ.m,ff)
				dat.fu<-use.lumi_NQ.m[i,]
				dat.f<-rbind(dat.fo,dat.fu)
				fit<-lmFit(dat.f)
				yy<-try(toptable(fit,coef=2),silent=TRUE)
				if(length(grep("Error in",yy))!=0){
					fit2<-eBayes(fit)
					}else{
						fit2<-eBayes(fit)
					}
				rm(lumi_NQ.f)
				lumi_NQ.f<<-fit2
#				print("Expression Filtering")
				}),silent=TRUE)
			}
		}

	if(length(lumi_NQ.f)!=0){
		rm(ttx)
		err<-try(ttx<<-toptable(lumi_NQ.f,coef=2,number=nrow(use.lumi_NQ.m)),silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(lumi_NQ.f,coef=3,number=nrow(use.lumi_NQ.m)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(lumi_NQ.f,coef=4,number=nrow(use.lumi_NQ.m)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)ttx<<-toptable(lumi_NQ.f,number=nrow(use.lumi_NQ.m))
		rn<-rownames(ttx)[ttx$P.Value<=0.01]
		rm(lumi_NQ.s)
		err=NULL;
		err<-try(lumi_NQ.s<<-use.lumi_NQ.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			lumi_NQ.s<<-use.lumi_NQ.m[rn,]
			}
		}
		rm(err)
		
	galert("Done",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Done"
	
	galert("Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(DE_Il_L)
	err<-try(DE_Il_L<<-toptable(lumi_NQ.f,coef=2,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p"),silent=TRUE)
	if(dim(DE_Il_L)[1]==0 || length(DE_Il_L)==0)
	{
		err<-try(DE_Il_L<<-toptable(lumi_NQ.f,coef=3,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p"),silent=TRUE)
	}
	if(dim(DE_Il_L)[1]==0 || length(DE_Il_L)==0)
	{
		err<-try(DE_Il_L<<-toptable(lumi_NQ.f,coef=4,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p"),silent=TRUE)
	}
	if(dim(DE_Il_L)[1]==0 || length(DE_Il_L)==0)
	{
		DE_Il_L<<-toptable(lumi_NQ.f,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p")
	}
	rm(err)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(pca_Il_L)
	pca_Il_L<<-prcomp(t(lumi_NQ.m))
	plot(pca_Il_L,main="Illumina_Lumi PCA")
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)

	sample.dist_c1<-dist(t(lumi_NQ.m))
	rm(sample.dist_Il_L) 
	sample.dist_Il_L<<-as.dist(1-cor(lumi_NQ.m,method="pearson"))
	rm(sample.clust_Il_L)
	sample.clust_Il_L<<-hclust(sample.dist_Il_L,method="complete")
	plot(sample.clust_Il_L)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)

	rm(use_DE_Il_L)
	if(dim(as.matrix(DE_Il_L))[1]>=20){
		er_x<-try(use_DE_Il_L<<-as.matrix(DE_Il_L)[1:20,],silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_Il_L)<<-as.character(DE_Il_L)
			use_DE_Il_L<<-as.matrix(DE_Il_L)[1:20,]
		}
	} else {
		er_x<-try(use_DE_Il_L<<-as.matrix(DE_Il_L),silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_Il_L)<<-as.character(DE_Il_L)
			use_DE_Il_L<<-as.matrix(DE_Il_L)
		}
	}
	
	rm(Clas_Il_L)
	err2<-try((Clas_Il_L<<-use.lumi_NQ.m[rownames(use_DE_Il_L),which(design_Il_L[,1]==1 | design_Il_L[,2]==1)]),silent=TRUE)
	if(length(grep("Error",err2))!=0)
	{
		err3<-try((Clas_Il_L<<-use.lumi_NQ.m[rownames(use_DE_Il_L),which(design_Il_L[,1]==1 | design_Il_L[,2]==1)] | design_Il_L[,3]==1),silent=TRUE)
		if(length(grep("Error",err3))!=0)
		{
			err4<-try(Clas_Il_L<<-use.lumi_NQ.m[rownames(use_DE_Il_L),],silent=TRUE)
			if(length(grep("Error",err4))!=0)
			{
				try(Clas_Il_L<<-use.lumi_NQ.m[as.numeric(rownames(use_DE_Il_L)),],silent=TRUE)
			}
		}
	}
	heatmap(Clas_Il_L,Rowv=NA)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
#		delete(g1,g1_1)
		visible(g1_1)<-FALSE
		l$Illumina_Lumi$Normalization<<-list()
		l$Illumina_Lumi$QC_Plot<<-list()
		l$Illumina_Lumi$Filtered<<-list()
		l$Illumina_Lumi$Stat_Significant<<-list()
		l$Illumina_Lumi$DGE<<-list()
		l$Illumina_Lumi$PCA_Plot<<-list()
		l$Illumina_Lumi$Cluster_Plot<<-list()
		l$Illumina_Lumi$Classification<<-list()
#		g1_1<<-ggroup(container=g1)
		tr<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,400)
#		tr<<-tr
		visible(g1_1)<-TRUE
		display()
		}),silent=TRUE)
#	display()	
}
