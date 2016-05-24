illbead<-function(h,...){
	datIllBA2=NULL;
	datIllBA2.m2=NULL;
	use.datIllBA2.m2=NULL;
	xf=NULL;
	datIllBA2.f=NULL;
	design_Il_B=NULL;
	ttx=NULL;
	datIllBA2.s=NULL;
	DE_Il_B=NULL;
	pca_Il_B=NULL;
	sample.dist_Il_B=NULL;
	sample.clust_Il_B=NULL;
	use_DE_Il_B=NULL;
	Clas_Il_B=NULL;
	design=NULL;

	try(({folder_Il_B<-folder_Il_B;datIllBA<-datIllBA;l<-l;tree<-tree;
	}),silent=TRUE)

	setwd(folder_Il_B)	
	galert("Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
#	rm(datIllBA2)
	datIllBA2<-normaliseIllumina(datIllBA,method="quantile")
	datIllBA2.m<-exprs(datIllBA2)
	rm(datIllBA2.m2)
	datIllBA2.m2<<-log2(datIllBA2.m)
	datIllBA2.m2<<-as.matrix(datIllBA2.m2)
	rownames(datIllBA2.m2)<<-as.character(rownames(datIllBA2.m2))
	try(if(length(rownames(datIllBA2.m2))==0)rownames(datIllBA2.m2)<<-1:length(datIllBA2.m2[,1]),silent=TRUE)
	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	boxplot(datIllBA2.m2)
	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
	galert("Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(use.datIllBA2.m2)
	use.datIllBA2.m2<<-datIllBA2.m2
	try(({
	rm(xf)
	xf<<-dim(use.datIllBA2.m2)
	if(xf[2]%%3!=0 && xf[2]%%2!=0){
		use.datIllBA2.m2<<-datIllBA2.m2[,-xf[2]]
		xf<<-dim(use.datIllBA2.m2)
	}
	if(xf[2]%%2==0)
	{
		yf=xf[2]/2
		groups<-c(rep("C",yf),rep("T",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.datIllBA2.m2,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(datIllBA2.f)
		datIllBA2.f<<-fit2
		rm(design_Il_B)
		design_Il_B<<-design
#		print("Two-group Specific Filtering")
	}else
	{
		yf=xf[2]/3
		groups<-c(rep("C",yf),rep("T1",yf),rep("T2",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.datIllBA2.m2,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(datIllBA2.f)
		datIllBA2.f<<-fit2
		rm(design_Il_B)
		design_Il_B<<-design
#		print("Three-group Specific Filtering")
	}
	}),silent=TRUE)
	
	if(length(datIllBA2.f)==0)
	{
		try(({
			rsd<-rowSds(use.datIllBA2.m2)
			i<-rsd>=2
			dat.f<-use.datIllBA2.m2[i,]
			fit<-lmFit(dat.f)
			yy<-try(toptable(fit,coef=2),silent=TRUE)
			if(length(grep("Error in",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
			rm(datIllBA2.f)
			datIllBA2.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
			
		if(length(datIllBA2.f)==0)
		{
			try(({
				ff<-pOverA(A=1,p=0.5)
				i<-genefilter(use.datIllBA2.m2,ff)
				dat.fo<-use.datIllBA2.m2[i,]
				i<-genefilter(-use.datIllBA2.m2,ff)
				dat.fu<-use.datIllBA2.m2[i,]
				dat.f<-rbind(dat.fo,dat.fu)
				fit<-lmFit(dat.f)
				yy<-try(toptable(fit,coef=2),silent=TRUE)
				if(length(grep("Error in",yy))!=0){
					fit2<-eBayes(fit)
					}else{
						fit2<-eBayes(fit)
					}
				rm(datIllBA2.f)
				datIllBA2.f<<-fit2
#				print("Expression Filtering")
				}),silent=TRUE)
			}
		}

	if(length(datIllBA2.f)!=0){
		rm(ttx)
		err<-try(ttx<<-toptable(datIllBA2.f,coef=2,number=nrow(use.datIllBA2.m2)),silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(datIllBA2.f,coef=3,number=nrow(use.datIllBA2.m2)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(datIllBA2.f,coef=4,number=nrow(use.datIllBA2.m2)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)toptable(datIllBA2.f,number=nrow(use.datIllBA2.m2))
		rn<-rownames(ttx)[ttx$P.Value<=0.01]
		rm(datIllBA2.s)
		err=NULL;
		err<-try(datIllBA2.s<<-use.datIllBA2.m2[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			datIllBA2.s<<-use.datIllBA2.m2[rn,]
			}
		}
		rm(err)
		
	galert("Done",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Done"
	
	galert("Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(DE_Il_B)
	err<-try(DE_Il_B<<-toptable(datIllBA2.f,coef=2,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p"),silent=TRUE)
	if(dim(DE_Il_B)[1]==0 || length(DE_Il_B)==0)
	{
		err<-try(DE_Il_B<<-toptable(datIllBA2.f,coef=3,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p"),silent=TRUE)
	}
	if(dim(DE_Il_B)[1]==0 || length(DE_Il_B)==0)
	{
		err<-try(DE_Il_B<<-toptable(datIllBA2.f,coef=4,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p"),silent=TRUE)
	}
	if(dim(DE_Il_B)[1]==0 || length(DE_Il_B)==0)
	{
		DE_Il_B<<-toptable(datIllBA2.f,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p")
	}
	rm(err)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(pca_Il_B)
	pca_Il_B<<-prcomp(t(datIllBA2.m2))
	plot(pca_Il_B,main="Illumina_Beadarray PCA")
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)

	sample.dist_c1<-dist(t(datIllBA2.m2))
	rm(sample.dist_Il_B) 
	sample.dist_Il_B<<-as.dist(1-cor(datIllBA2.m2,method="pearson"))
	rm(sample.clust_Il_B)
	sample.clust_Il_B<<-hclust(sample.dist_Il_B,method="complete")
	plot(sample.clust_Il_B)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)

	rm(use_DE_Il_B)
	if(dim(as.matrix(DE_Il_B))[1]>=20){
		er_x<-try(use_DE_Il_B<<-as.matrix(DE_Il_B)[1:20,],silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_Il_B)<<-as.character(DE_Il_B)
			use_DE_Il_B<<-as.matrix(DE_Il_B)[1:20,]
		}
	} else {
		er_x<-try(use_DE_Il_B<<-as.matrix(DE_Il_B),silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_Il_B)<<-as.character(DE_Il_B)
			use_DE_Il_B<<-as.matrix(DE_Il_B)
		}
	}
	
	rm(Clas_Il_B)
	err2<-try((Clas_Il_B<<-use.datIllBA2.m2[rownames(use_DE_Il_B),which(design_Il_B[,1]==1 | design_Il_B[,2]==1)]),silent=TRUE)
	if(length(grep("Error",err2))!=0)
	{
		err3<-try((Clas_Il_B<<-use.datIllBA2.m2[rownames(use_DE_Il_B),which(design_Il_B[,1]==1 | design_Il_B[,2]==1)] | design_Il_B[,3]==1),silent=TRUE)
		if(length(grep("Error",err3))!=0)
		{
			err4<-try(Clas_Il_B<<-use.datIllBA2.m2[rownames(use_DE_Il_B),],silent=TRUE)
			if(length(grep("Error",err4))!=0)
			{
				try(Clas_Il_B<<-use.datIllBA2.m2[as.numeric(rownames(use_DE_Il_B)),],silent=TRUE)
			}
		}
	}
	heatmap(Clas_Il_B,Rowv=NA)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
#		delete(g1,g1_1)
		visible(g1_1)<-FALSE
		l$Illumina_Beadarray$Normalization<<-list()
		l$Illumina_Beadarray$QC_Plot<<-list()
		l$Illumina_Beadarray$Filtered<<-list()
		l$Illumina_Beadarray$Stat_Significant<<-list()
		l$Illumina_Beadarray$DGE<<-list()
		l$Illumina_Beadarray$PCA_Plot<<-list()
		l$Illumina_Beadarray$Cluster_Plot<<-list()
		l$Illumina_Beadarray$Classification<<-list()
#		g1_1<<-ggroup(container=g1)
		tr<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,400)
#		tr<<-tr
		visible(g1_1)<-TRUE
		display()
		}),silent=TRUE)
#	display()	
}
