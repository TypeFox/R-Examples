agone<-function(h,...){
	datAgOne2=NULL;
	datAgOne2.m=NULL;
	use.datAgOne2.m=NULL;
	xf=NULL;
	datAgOne2.f=NULL;
	design_Ag1=NULL;
	ttx=NULL;
	datAgOne2.s=NULL;
	DE_Ag1=NULL;
	pca_Ag1=NULL;
	sample.dist_Ag1=NULL;
	sample.clust_Ag1=NULL;
	use_DE_Ag1=NULL;
	Clas_Ag1=NULL;
	design=NULL;

	try(({folder_Ag1<-folder_Ag1;datAgOne<-datAgOne;l<-l;tree<-tree;
	}),silent=TRUE)

	setwd(folder_Ag1)	
	galert("Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	try(datAgOne2<<-backgroundCorrect(datAgOne,"normexp"),silent=TRUE)
	if(length(datAgOne2)==0)datAgOne2<-datAgOne
	datAgOne2<-normalizeBetweenArrays(datAgOne2$R,method="quantile")
	datAgOne2<-log(datAgOne2)
	rm(datAgOne2.m)
	datAgOne2.m<<-datAgOne2
	datAgOne2.m<<-as.matrix(datAgOne2.m)
	rownames(datAgOne2.m)<<-as.character(rownames(datAgOne2.m))
	try(if(length(rownames(datAgOne2.m))==0)rownames(datAgOne2.m)<<-1:length(datAgOne2.m[,1]),silent=TRUE)
	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	boxplot(datAgOne2.m)
	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
	galert("Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(use.datAgOne2.m)
	use.datAgOne2.m<<-datAgOne2.m
	try(({
	rm(xf)
	xf<<-dim(use.datAgOne2.m)
	if(xf[2]%%3!=0 && xf[2]%%2!=0){
		use.datAgOne2.m<<-datAgOne2.m[,-xf[2]]
		xf<<-dim(use.datAgOne2.m)
	}
	if(xf[2]%%2==0)
	{
		yf=xf[2]/2
		groups<-c(rep("C",yf),rep("T",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.datAgOne2.m,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(datAgOne2.f)
		datAgOne2.f<<-fit2
		rm(design_Ag1)
		design_Ag1<<-design
#		print("Two-group Specific Filtering")
	}else
	{
		yf=xf[2]/3
		groups<-c(rep("C",yf),rep("T1",yf),rep("T2",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.datAgOne2.m,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(datAgOne2.f)
		datAgOne2.f<<-fit2
		rm(design_Ag1)
		design_Ag1<<-design
#		print("Three-group Specific Filtering")
	}
	}),silent=TRUE)
	
	if(length(datAgOne2.f)==0)
	{
		try(({
			rsd<-rowSds(use.datAgOne2.m)
			i<-rsd>=2
			dat.f<-use.datAgOne2.m[i,]
			fit<-lmFit(dat.f)
			yy<-try(toptable(fit,coef=2),silent=TRUE)
			if(length(grep("Error in",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
			rm(datAgOne2.f)
			datAgOne2.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
			
		if(length(datAgOne2.f)==0)
		{
			try(({
				ff<-pOverA(A=1,p=0.5)
				i<-genefilter(use.datAgOne2.m,ff)
				dat.fo<-use.datAgOne2.m[i,]
				i<-genefilter(-use.datAgOne2.m,ff)
				dat.fu<-use.datAgOne2.m[i,]
				dat.f<-rbind(dat.fo,dat.fu)
				fit<-lmFit(dat.f)
				yy<-try(toptable(fit,coef=2),silent=TRUE)
				if(length(grep("Error in",yy))!=0){
					fit2<-eBayes(fit)
					}else{
						fit2<-eBayes(fit)
					}
				rm(datAgOne2.f)
				datAgOne2.f<<-fit2
#				print("Expression Filtering")
				}),silent=TRUE)
			}
		}

	if(length(datAgOne2.f)!=0){
		rm(ttx)
		err<-try(ttx<<-toptable(datAgOne2.f,coef=2,number=nrow(use.datAgOne2.m)),silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(datAgOne2.f,coef=3,number=nrow(use.datAgOne2.m)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(datAgOne2.f,coef=4,number=nrow(use.datAgOne2.m)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)ttx<<-toptable(datAgOne2.f,number=nrow(use.datAgOne2.m))
		rn<-rownames(ttx)[ttx$P.Value<=0.01]
		rm(datAgOne2.s)
		err=NULL;
		err<-try(datAgOne2.s<<-use.datAgOne2.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			datAgOne2.s<<-use.datAgOne2.m[rn,]
			}
		}
		rm(err)
		
	galert("Done",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Done"
	
	galert("Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(DE_Ag1)
	err<-try(DE_Ag1<<-toptable(datAgOne2.f,coef=2,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p"),silent=TRUE)
	if(dim(DE_Ag1)[1]==0 || length(DE_Ag1)==0)
	{
		err<-try(DE_Ag1<<-toptable(datAgOne2.f,coef=3,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p"),silent=TRUE)
	}
	if(dim(DE_Ag1)[1]==0 || length(DE_Ag1)==0)
	{
		err<-try(DE_Ag1<<-toptable(datAgOne2.f,coef=4,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p"),silent=TRUE)
	}
	if(dim(DE_Ag1)[1]==0 || length(DE_Ag1)==0)
	{
		DE_Ag1<<-toptable(datAgOne2.f,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p")
	}
	rm(err)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(pca_Ag1)
	pca_Ag1<<-prcomp(t(datAgOne2.m))
	plot(pca_Ag1,main="Agilent_OneColor PCA")
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)

	sample.dist_c1<-dist(t(datAgOne2.m)) 
	rm(sample.dist_Ag1)
	sample.dist_Ag1<<-as.dist(1-cor(datAgOne2.m,method="pearson"))
	rm(sample.clust_Ag1)
	sample.clust_Ag1<<-hclust(sample.dist_Ag1,method="complete")
	plot(sample.clust_Ag1)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait...";
	Sys.sleep(1)
	
	rm(use_DE_Ag1)
	if(dim(as.matrix(DE_Ag1))[1]>=20){
		er_x<-try(use_DE_Ag1<<-as.matrix(DE_Ag1)[1:20,],silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_Ag1)<<-as.character(DE_Ag1)
			use_DE_Ag1<<-as.matrix(DE_Ag1)[1:20,]
		}
	} else {
		er_x<-try(use_DE_Ag1<<-as.matrix(DE_Ag1),silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_Ag1)<<-as.character(DE_Ag1)
			use_DE_Ag1<<-as.matrix(DE_Ag1)
		}
	}

	rm(Clas_Ag1)
	err2<-try((Clas_Ag1<<-use.datAgOne2.m[rownames(use_DE_Ag1),which(design_Ag1[,1]==1 | design_Ag1[,2]==1)]),silent=TRUE)
	if(length(grep("Error",err2))!=0)
	{
		err3<-try((Clas_Ag1<<-use.datAgOne2.m[rownames(use_DE_Ag1),which(design_Ag1[,1]==1 | design_Ag1[,2]==1)] | design_Ag1[,3]==1),silent=TRUE)
		if(length(grep("Error",err3))!=0)
		{
			err4<-try(Clas_Ag1<<-use.datAgOne2.m[rownames(use_DE_Ag1),],silent=TRUE)
			if(length(grep("Error",err4))!=0)
			{
				try(Clas_Ag1<<-use.datAgOne2.m[as.numeric(rownames(use_DE_Ag1)),],silent=TRUE)
			}
		}
	}
	heatmap(Clas_Ag1,Rowv=NA)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
		visible(g1_1)<-FALSE
		l$Agilent_OneColor$Normalization<<-list()
		l$Agilent_OneColor$QC_Plot<<-list()
		l$Agilent_OneColor$Filtered<<-list()
		l$Agilent_OneColor$Stat_Significant<<-list()
		l$Agilent_OneColor$DGE<<-list()
		l$Agilent_OneColor$PCA_Plot<<-list()
		l$Agilent_OneColor$Cluster_Plot<<-list()
		l$Agilent_OneColor$Classification<<-list()
		tr<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,400)
		visible(g1_1)<-TRUE
		display()
		}),silent=TRUE)
}
