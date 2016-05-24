nimblg<-function(h,...){
	data.matrix_Nimblegen2=NULL;
	data.matrix_Nimblegen2.m=NULL;
	use.data.matrix_Nimblegen2.m=NULL;
	xf=NULL;
	data.matrix_Nimblegen2.f=NULL;
	design_N=NULL;
	ttx=NULL;
	data.matrix_Nimblegen2.s=NULL;
	DE_N=NULL;
	pca_N=NULL;
	sample.dist_N=NULL;
	sample.clust_N=NULL;
	use_DE_N=NULL;
	Clas_N=NULL;
	design=NULL;

	try(({folder_N<-folder_N;data.matrix_Nimblegen<-data.matrix_Nimblegen;l<-l;tree<-tree;
	}),silent=TRUE)

	setwd(folder_N)	
	galert("Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	data.matrix_Nimblegen2<<-normalizeBetweenArrays(data.matrix_Nimblegen,method="quantile")
	rm(data.matrix_Nimblegen2.m)
	data.matrix_Nimblegen2.m<<-normalizeBetweenArrays(data.matrix_Nimblegen,method="quantile")
	data.matrix_Nimblegen2.m<<-as.matrix(data.matrix_Nimblegen2.m)
	rownames(data.matrix_Nimblegen2.m)<<-as.character(rownames(data.matrix_Nimblegen2.m))
	try(if(length(rownames(data.matrix_Nimblegen2.m))==0)rownames(data.matrix_Nimblegen2.m)<<-1:length(data.matrix_Nimblegen2.m[,1]),silent=TRUE)
	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	boxplot(data.matrix_Nimblegen2.m)
	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
	galert("Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(use.data.matrix_Nimblegen2.m)
	use.data.matrix_Nimblegen2.m<<-data.matrix_Nimblegen2.m
	try(({
	rm(xf)
	xf<<-dim(use.data.matrix_Nimblegen2.m)
	if(xf[2]%%3!=0 && xf[2]%%2!=0){
		use.data.matrix_Nimblegen2.m<<-data.matrix_Nimblegen2.m[,-xf[2]]
		xf<<-dim(use.data.matrix_Nimblegen2.m)
	}
	if(xf[2]%%2==0)
	{
		yf=xf[2]/2
		groups<-c(rep("C",yf),rep("T",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.data.matrix_Nimblegen2.m,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(data.matrix_Nimblegen2.f)
		data.matrix_Nimblegen2.f<<-fit2
		rm(design_N)
		design_N<<-design
#		print("Two-group Specific Filtering")
	}else
	{
		yf=xf[2]/3
		groups<-c(rep("C",yf),rep("T1",yf),rep("T2",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.data.matrix_Nimblegen2.m,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(data.matrix_Nimblegen2.f)
		data.matrix_Nimblegen2.f<<-fit2
		rm(design_N)
		design_N<<-design
#		print("Three-group Specific Filtering")
	}
	}),silent=TRUE)
	
	if(length(data.matrix_Nimblegen2.f)==0)
	{
		try(({
			rsd<-rowSds(use.data.matrix_Nimblegen2.m)
			i<-rsd>=2
			dat.f<-use.data.matrix_Nimblegen2.m[i,]
			fit<-lmFit(dat.f)
			yy<-try(toptable(fit,coef=2),silent=TRUE)
			if(length(grep("Error in",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
			rm(data.matrix_Nimblegen2.f)
			data.matrix_Nimblegen2.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
			
		if(length(data.matrix_Nimblegen2.f)==0)
		{
			try(({
				ff<-pOverA(A=1,p=0.5)
				i<-genefilter(use.data.matrix_Nimblegen2.m,ff)
				dat.fo<-use.data.matrix_Nimblegen2.m[i,]
				i<-genefilter(-use.data.matrix_Nimblegen2.m,ff)
				dat.fu<-use.data.matrix_Nimblegen2.m[i,]
				dat.f<-rbind(dat.fo,dat.fu)
				fit<-lmFit(dat.f)
				yy<-try(toptable(fit,coef=2),silent=TRUE)
				if(length(grep("Error in",yy))!=0){
					fit2<-eBayes(fit)
					}else{
						fit2<-eBayes(fit)
					}
				rm(data.matrix_Nimblegen2.f)
				data.matrix_Nimblegen2.f<<-fit2
#				print("Expression Filtering")
				}),silent=TRUE)
			}
		}

	if(length(data.matrix_Nimblegen2.f)!=0){
		rm(ttx)
		err<-try(ttx<<-toptable(data.matrix_Nimblegen2.f,coef=2,number=nrow(use.data.matrix_Nimblegen2.m)),silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(data.matrix_Nimblegen2.f,coef=3,number=nrow(use.data.matrix_Nimblegen2.m)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(data.matrix_Nimblegen2.f,coef=4,number=nrow(use.data.matrix_Nimblegen2.m)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)ttx<<-toptable(data.matrix_Nimblegen2.f,number=nrow(use.data.matrix_Nimblegen2.m))
		rn<-rownames(ttx)[ttx$P.Value<=0.01]
		rm(data.matrix_Nimblegen2.s)
		err=NULL;
		err<-try(data.matrix_Nimblegen2.s<<-use.data.matrix_Nimblegen2.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			data.matrix_Nimblegen2.s<<-use.data.matrix_Nimblegen2.m[rn,]
			}
		}
		rm(err)
		
	galert("Done",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Done"
	
	galert("Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(DE_N)
	err<-try(DE_N<<-toptable(data.matrix_Nimblegen2.f,coef=2,number=20,lfc=2,adjust.method="BH",p.value=0.01,
	sort.by="p"),silent=TRUE)
	if(dim(DE_N)[1]==0 || length(DE_N)==0)
	{
		err<-try(DE_N<<-toptable(data.matrix_Nimblegen2.f,coef=3,number=20,lfc=2,adjust.method="BH",p.value=0.01,
		sort.by="p"),silent=TRUE)
	}
	if(dim(DE_N)[1]==0 || length(DE_N)==0)
	{
		err<-try(DE_N<<-toptable(data.matrix_Nimblegen2.f,coef=4,number=20,lfc=2,adjust.method="BH",p.value=0.01,
		sort.by="p"),silent=TRUE)
	}
	if(dim(DE_N)[1]==0 || length(DE_N)==0)
	{
		DE_N<<-toptable(data.matrix_Nimblegen2.f,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p")
	}
	rm(err)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(pca_N)
	pca_N<<-prcomp(t(data.matrix_Nimblegen2.m))
	plot(pca_N,main="Nimblegen PCA")
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)

	sample.dist_c1<-dist(t(data.matrix_Nimblegen2.m))
	rm(sample.dist_N)
	sample.dist_N<<-as.dist(1-cor(data.matrix_Nimblegen2.m,method="pearson"))
	rm(sample.clust_N)
	sample.clust_N<<-hclust(sample.dist_N,method="complete")
	plot(sample.clust_N)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)

	rm(use_DE_N)
	if(dim(as.matrix(DE_N))[1]>19){
		er_x<-try(use_DE_N<<-as.matrix(DE_N)[1:20,],silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_N)<<-as.character(DE_N)
			use_DE_N<<-as.matrix(DE_N)[1:20,]
		}
	} else {
		er_x<-try(use_DE_N<<-as.matrix(DE_N),silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_N)<<-as.character(DE_N)
			use_DE_N<<-as.matrix(DE_N)
		}
	}
	
	rm(Clas_N)
	err2<-try((Clas_N<<-use.data.matrix_Nimblegen2.m[rownames(use_DE_N),which(design_N[,1]==1 | design_N[,2]==1)]),silent=TRUE)
	if(length(grep("Error",err2))!=0)
	{
		err3<-try((Clas_N<<-use.data.matrix_Nimblegen2.m[rownames(use_DE_N),which(design_N[,1]==1 | design_N[,2]==1)] | design_N[,3]==1),silent=TRUE)
		if(length(grep("Error",err3))!=0)
		{
			err4<-try(Clas_N<<-use.data.matrix_Nimblegen2.m[rownames(use_DE_N),],silent=TRUE)
			if(length(grep("Error",err4))!=0)
			{
				try(Clas_N<<-use.data.matrix_Nimblegen2.m[as.numeric(rownames(use_DE_N)),],silent=TRUE)
			}
		}
	}
	heatmap(Clas_N,Rowv=NA)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
		visible(g1_1)<-FALSE
		l$Nimblegen$Normalization<<-list()
		l$Nimblegen$QC_Plot<<-list()
		l$Nimblegen$Filtered<<-list()
		l$Nimblegen$Stat_Significant<<-list()
		l$Nimblegen$DGE<<-list()
		l$Nimblegen$PCA_Plot<<-list()
		l$Nimblegen$Cluster_Plot<<-list()
		l$Nimblegen$Classification<<-list()
		tr<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,400)
		visible(g1_1)<-TRUE
		display()
		}),silent=TRUE)
}
