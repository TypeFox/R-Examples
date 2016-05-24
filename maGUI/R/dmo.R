dmo<-function(h,...){
	data.matrix_onlineNorm=NULL;
	data.matrix_onlineNorm.m=NULL;
	use.data.matrix_onlineNorm.m=NULL;
	xf=NULL;
	data.matrix_onlineNorm.f=NULL;
	design_O=NULL;
	ttx=NULL;
	data.matrix_onlineNorm.s=NULL;
	DE_O=NULL;
	pca_O=NULL;
	sample.dist_O=NULL;
	sample.clust_O=NULL;
	use_DE_O=NULL;
	Clas_O=NULL;
	design=NULL;

	try(({data.matrix_onlineImp<-data.matrix_onlineImp;l<-l;tree<-tree;
	}),silent=TRUE)

	galert("Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	data.matrix_onlineNorm<-normalizeBetweenArrays(data.matrix_onlineImp,method="quantile")
	tmp<-aggregate(data.matrix_onlineNorm,list(rownames(data.matrix_onlineNorm)),median)
	rm(data.matrix_onlineNorm.m)
	data.matrix_onlineNorm.m<<-as.matrix(tmp[,-1])
	rownames(data.matrix_onlineNorm.m)<<-tmp[,1]
	rownames(data.matrix_onlineNorm.m)<<-as.character(rownames(data.matrix_onlineNorm.m))
	try(if(length(rownames(data.matrix_onlineNorm.m))==0)rownames(data.matrix_onlineNorm.m)<<-1:length(data.matrix_onlineNorm.m[,1]),silent=TRUE)
	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	boxplot(data.matrix_onlineNorm.m)
	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
	galert("Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(use.data.matrix_onlineNorm.m)
	use.data.matrix_onlineNorm.m<<-data.matrix_onlineNorm.m
	try(({
	rm(xf)
	xf<<-dim(use.data.matrix_onlineNorm.m)
	if(xf[2]%%3!=0 && xf[2]%%2!=0){
		use.data.matrix_onlineNorm.m<<-data.matrix_onlineNorm.m[,-xf[2]]
		xf<<-dim(use.data.matrix_onlineNorm.m)
	}
	if(xf[2]%%2==0)
	{
		yf=xf[2]/2
		groups<-c(rep("C",yf),rep("T",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.data.matrix_onlineNorm.m,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(data.matrix_onlineNorm.f)
		data.matrix_onlineNorm.f<<-fit2
		rm(design_O)
		design_O<<-design
#		print("Two-group Specific Filtering")
	}else
	{
		yf=xf[2]/3
		groups<-c(rep("C",yf),rep("T1",yf),rep("T2",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.data.matrix_onlineNorm.m,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(data.matrix_onlineNorm.f)
		data.matrix_onlineNorm.f<<-fit2
		rm(design_O)
		design_O<<-design
#		print("Three-group Specific Filtering")
	}
	}),silent=TRUE)
	
	if(length(data.matrix_onlineNorm.f)==0)
	{
		try(({
			rsd<-rowSds(use.data.matrix_onlineNorm.m)
			i<-rsd>=2
			dat.f<-use.data.matrix_onlineNorm.m[i,]
			fit<-lmFit(dat.f)
			yy<-try(toptable(fit,coef=2),silent=TRUE)
			if(length(grep("Error in",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
			rm(data.matrix_onlineNorm.f)
			data.matrix_onlineNorm.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
			
		if(length(data.matrix_onlineNorm.f)==0)
		{
			try(({
				ff<-pOverA(A=1,p=0.5)
				i<-genefilter(use.data.matrix_onlineNorm.m,ff)
				dat.fo<-use.data.matrix_onlineNorm.m[i,]
				i<-genefilter(-use.data.matrix_onlineNorm.m,ff)
				dat.fu<-use.data.matrix_onlineNorm.m[i,]
				dat.f<-rbind(dat.fo,dat.fu)
				fit<-lmFit(dat.f)
				yy<-try(toptable(fit,coef=2),silent=TRUE)
				if(length(grep("Error in",yy))!=0){
					fit2<-eBayes(fit)
					}else{
						fit2<-eBayes(fit)
					}
				rm(data.matrix_onlineNorm.f)
				data.matrix_onlineNorm.f<<-fit2
#				print("Expression Filtering")
				}),silent=TRUE)
			}
		}

	if(length(data.matrix_onlineNorm.f)!=0){
		rm(ttx)
		err<-try(ttx<<-toptable(data.matrix_onlineNorm.f,coef=2,number=nrow(use.data.matrix_onlineNorm.m)),silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(data.matrix_onlineNorm.f,coef=3,number=nrow(use.data.matrix_onlineNorm.m)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(data.matrix_onlineNorm.f,coef=4,number=nrow(use.data.matrix_onlineNorm.m)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)ttx<<-toptable(data.matrix_onlineNorm.f,number=nrow(use.data.matrix_onlineNorm.m))
		rn<-rownames(ttx)[ttx$P.Value<=0.01]
		rm(data.matrix_onlineNorm.s)
		err=NULL;
		err<-try(data.matrix_onlineNorm.s<<-use.data.matrix_onlineNorm.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			try(data.matrix_onlineNorm.s<<-use.data.matrix_onlineNorm.m[rn,],silent=TRUE)
			}
		}
		rm(err)
		
	galert("Done",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Done"
	
	galert("Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(DE_O)
	err<-try(DE_O<<-toptable(data.matrix_onlineNorm.f,coef=2,number=20,lfc=2,adjust.method="BH",p.value=0.01,
	sort.by="p"),silent=TRUE)
	if(dim(DE_O)[1]==0 || length(DE_O)==0)
	{
		err<-try(DE_O<<-toptable(data.matrix_onlineNorm.f,coef=3,number=20,lfc=2,adjust.method="BH",p.value=0.01,
		sort.by="p"),silent=TRUE)
	}
	if(dim(DE_O)[1]==0 || length(DE_O)==0)
	{
		err<-try(DE_O<<-toptable(data.matrix_onlineNorm.f,coef=4,number=20,lfc=2,adjust.method="BH",p.value=0.01,
		sort.by="p"),silent=TRUE)
	}
	if(dim(DE_O)[1]==0 || length(DE_O)==0)
	{
		DE_O<<-toptable(data.matrix_onlineNorm.f,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p")
	}
	rm(err)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(pca_O)
	pca_O<<-prcomp(t(data.matrix_onlineNorm.m))
	plot(pca_O,main="Online_Data PCA")
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)

	sample.dist_c1<-dist(t(data.matrix_onlineNorm.m)) 
	rm(sample.dist_O)
	sample.dist_O<<-as.dist(1-cor(use.data.matrix_onlineNorm.m,method="pearson"))
	rm(sample.clust_O)
	sample.clust_O<<-hclust(sample.dist_O,method="complete")
	plot(sample.clust_O)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)

	rm(use_DE_O)
	if(dim(as.matrix(DE_O))[1]>=20){
		er_x<-try(use_DE_O<<-as.matrix(DE_O)[1:20,],silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_O)<<-as.character(DE_O)
			use_DE_O<<-as.matrix(DE_O)[1:20,]
		}
	} else {
		er_x<-try(use_DE_O<<-as.matrix(DE_O),silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_O)<<-as.character(DE_O)
			use_DE_O<<-as.matrix(DE_O)
		}
	}
	
	rm(Clas_O)
	err2<-try((Clas_O<<-use.data.matrix_onlineNorm.m[rownames(use_DE_O),which(design_O[,1]==1 | design_O[,2]==1)]),silent=TRUE)
	if(length(grep("Error",err2))!=0)
	{
		err3<-try((Clas_O<<-use.data.matrix_onlineNorm.m[rownames(use_DE_O),which(design_O[,1]==1 | design_O[,2]==1)] | design_O[,3]==1),silent=TRUE)
		if(length(grep("Error",err3))!=0)
		{
			err4<-try(Clas_O<<-use.data.matrix_onlineNorm.m[rownames(use_DE_O),],silent=TRUE)
			if(length(grep("Error",err4))!=0)
			{
				try(Clas_O<<-use.data.matrix_onlineNorm.m[as.numeric(rownames(use_DE_O)),],silent=TRUE)
			}
		}
	}
	heatmap(Clas_O,Rowv=NA)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
		visible(g1_1)<-FALSE
		l$Online_Data$Normalization<<-list()
		l$Online_Data$QC_Plot<<-list()
		l$Online_Data$Filtered<<-list()
		l$Online_Data$Stat_Significant<<-list()
		l$Online_Data$DGE<<-list()
		l$Online_Data$PCA_Plot<<-list()
		l$Online_Data$Cluster_Plot<<-list()
		l$Online_Data	$Classification<<-list()
		tr<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,400)
		visible(g1_1)<-TRUE
		display()
		}),silent=TRUE)
}
