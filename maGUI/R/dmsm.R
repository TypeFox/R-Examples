dmsm<-function(h,...){
	data.matrixNorm=NULL;
	data.matrixNorm.m=NULL;
	use.data.matrixNorm.m=NULL;
	xf=NULL;
	data.matrixNorm.f=NULL;
	design_S=NULL;
	ttx=NULL;
	data.matrixNorm.s=NULL;
	DE_S=NULL;
	pca_S=NULL;
	sample.dist_S=NULL;
	sample.clust_S=NULL;
	use_DE_S=NULL;
	Clas_S=NULL;
	design=NULL;

	try(({folder_S<-folder_S;data.matrixImp<-data.matrixImp;l<-l;tree<-tree;
	}),silent=TRUE)

	setwd(folder_S)	
	galert("Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	
	if(length(data.matrixImp)!=0){
	data.matrixNorm<-normalizeBetweenArrays(data.matrixImp,method="quantile")
	tmp<-aggregate(data.matrixNorm,list(rownames(data.matrixNorm)),median)
	rm(data.matrixNorm.m)
	data.matrixNorm.m<<-as.matrix(tmp[,-1])
	rownames(data.matrixNorm.m)<<-tmp[,1]
	try(if(length(rownames(data.matrixNorm.m))==0)rownames(data.matrixNorm.m)<<-1:length(data.matrixNorm.m[,1]),silent=TRUE)
	}
	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	boxplot(data.matrixNorm.m)
	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
	galert("Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(use.data.matrixNorm.m)
	use.data.matrixNorm.m<<-data.matrixNorm.m
	try(({
	rm(xf)
	xf<<-dim(use.data.matrixNorm.m)
	if(xf[2]%%3!=0 && xf[2]%%2!=0){
		use.data.matrixNorm.m<<-data.matrixNorm.m[,-xf[2]]
		xf<<-dim(use.data.matrixNorm.m)
	}
	if(xf[2]%%2==0)
	{
		yf=xf[2]/2
		groups<-c(rep("C",yf),rep("T",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.data.matrixNorm.m,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(data.matrixNorm.f)
		data.matrixNorm.f<<-fit2
		rm(design_S)
		design_S<<-design
#		print("Two-group Specific Filtering")
	}else
	{
		yf=xf[2]/3
		groups<-c(rep("C",yf),rep("T1",yf),rep("T2",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.data.matrixNorm.m,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(data.matrixNorm.f)
		data.matrixNorm.f<<-fit2
		rm(design_S)
		design_S<<-design
#		print("Three-group Specific Filtering")
	}
	}),silent=TRUE)
	
	if(length(data.matrixNorm.f)==0)
	{
		try(({
			rsd<-rowSds(use.data.matrixNorm.m)
			i<-rsd>=2
			dat.f<-use.data.matrixNorm.m[i,]
			fit<-lmFit(dat.f)
			yy<-try(toptable(fit,coef=2),silent=TRUE)
			if(length(grep("Error in",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
			rm(data.matrixNorm.f)
			data.matrixNorm.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
			
		if(length(data.matrixNorm.f)==0)
		{
			try(({
				ff<-pOverA(A=1,p=0.5)
				i<-genefilter(use.data.matrixNorm.m,ff)
				dat.fo<-use.data.matrixNorm.m[i,]
				i<-genefilter(-use.data.matrixNorm.m,ff)
				dat.fu<-use.data.matrixNorm.m[i,]
				dat.f<-rbind(dat.fo,dat.fu)
				fit<-lmFit(dat.f)
				yy<-try(toptable(fit,coef=2),silent=TRUE)
				if(length(grep("Error in",yy))!=0){
					fit2<-eBayes(fit)
					}else{
						fit2<-eBayes(fit)
					}
				rm(data.matrixNorm.f)
				data.matrixNorm.f<<-fit2
#				print("Expression Filtering")
				}),silent=TRUE)
			}
		}

	if(length(data.matrixNorm.f)!=0){
		rm(ttx)
		err<-try(ttx<<-toptable(data.matrixNorm.f,coef=2,number=nrow(use.data.matrixNorm.m)),silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(data.matrixNorm.f,coef=3,number=nrow(use.data.matrixNorm.m)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(data.matrixNorm.f,coef=4,number=nrow(use.data.matrixNorm.m)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)ttx<<-toptable(data.matrixNorm.f,number=nrow(use.data.matrixNorm.m))
		rn<-rownames(ttx)[ttx$P.Value<=0.01]
		rm(data.matrixNorm.s)
		err=NULL;
		err<-try(data.matrixNorm.s<<-use.data.matrixNorm.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			try(data.matrixNorm.s<<-use.data.matrixNorm.m[rn,],silent=TRUE)
			}
		}
		rm(err)
		
	galert("Done",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Done"
	
	galert("Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(DE_S)
	err<-try(DE_S<<-toptable(data.matrixNorm.f,coef=2,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p"),silent=TRUE)
	if(dim(DE_S)[1]==0 || length(DE_S)==0)
	{
		err<-try(DE_S<<-toptable(data.matrixNorm.f,coef=3,number=20,lfc=2,adjust.method="BH",p.value=0.01,
		sort.by="p"),silent=TRUE)
	}
	if(dim(DE_S)[1]==0 || length(DE_S)==0)
	{
		err<-try(DE_S<<-toptable(data.matrixNorm.f,coef=4,number=20,lfc=2,adjust.method="BH",p.value=0.01,
		sort.by="p"),silent=TRUE)
	}
	if(dim(DE_S)[1]==0 || length(DE_S)==0)
	{
		DE_S<<-toptable(data.matrixNorm.f,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p")
	}
	rm(err)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(pca_S)
	pca_S<<-prcomp(t(data.matrixNorm.m))
	plot(pca_S,main="Series_Matrix PCA")
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)

	sample.dist_c1<-dist(t(data.matrixNorm.m))
	rm(sample.dist_S)
	sample.dist_S<<-as.dist(1-cor(data.matrixNorm.m,method="pearson"))
	rm(sample.clust_S)
	sample.clust_S<<-hclust(sample.dist_S,method="complete")
	plot(sample.clust_S)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)

	rm(use_DE_S)
	if(dim(as.matrix(DE_S))[1]>=20){
		er_x<-try(use_DE_S<<-as.matrix(DE_S)[1:20,],silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_S)<<-as.character(DE_S)
			use_DE_S<<-as.matrix(DE_S)[1:20,]
		}
	} else {
		er_x<-try(use_DE_S<<-as.matrix(DE_S),silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_S)<<-as.character(DE_S)
			use_DE_S<<-as.matrix(DE_S)
		}
	}
	
	rm(Clas_S)
	err2<-try((Clas_S<<-use.data.matrixNorm.m[rownames(use_DE_S),which(design_S[,1]==1 | design_S[,2]==1)]),silent=TRUE)
	if(length(grep("Error",err2))!=0)
	{
		err3<-try((Clas_S<<-use.data.matrixNorm.m[rownames(use_DE_S),which(design_S[,1]==1 | design_S[,2]==1)] | design_S[,3]==1),silent=TRUE)
		if(length(grep("Error",err3))!=0)
		{
			err4<-try(Clas_S<<-use.data.matrixNorm.m[rownames(use_DE_S),],silent=TRUE)
			if(length(grep("Error",err4))!=0)
			{
				try(Clas_S<<-use.data.matrixNorm.m[as.numeric(rownames(use_DE_S)),],silent=TRUE)
			}
		}
	}
	heatmap(Clas_S,Rowv=NA)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
		visible(g1_1)<-FALSE
		l$Series_Matrix$Normalization<<-list()
		l$Series_Matrix$QC_Plot<<-list()
		l$Series_Matrix$Filtered<<-list()
		l$Series_Matrix$Stat_Significant<<-list()
		l$Series_Matrix$DGE<<-list()
		l$Series_Matrix$PCA_Plot<<-list()
		l$Series_Matrix$Cluster_Plot<<-list()
		l$Series_Matrix	$Classification<<-list()
		tr<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,400)
		visible(g1_1)<-TRUE
		display()
		}),silent=TRUE)
}
