affym<-function(h,...){
	dat2Affy.m=NULL;
	aqc=NULL;
	use.dat2Affy.m=NULL;
	xf=NULL;
	dat2Affy.f=NULL;
	design_Affy=NULL;
	ttx=NULL;
	dat2Affy.s=NULL;
	DE_Affy=NULL;
	pca_Affy=NULL;
	sample.dist_Affy=NULL;
	sample.clust_Affy=NULL;
	use_DE_Affy=NULL;
	Clas_Affy=NULL;
	design=NULL;

	try(({folder_Affy<-folder_Affy;datAffy<-datAffy;l<-l;tree<-tree;
	}),silent=TRUE)
	setwd(folder_Affy)
	galert("Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	dat2Affy<-justRMA()
	dat2Affy.m<-exprs(dat2Affy)
	dat2Affy.m<<-as.matrix(dat2Affy.m)
	rownames(dat2Affy.m)<<-as.character(rownames(dat2Affy.m))
	try(if(length(rownames(dat2Affy.m))==0)rownames(dat2Affy.m)<<-1:length(dat2Affy.m[,1]),silent=TRUE)

	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(aqc)
	try(aqc<<-qc(datAffy),silent=TRUE)
	err<-try(plot(aqc),silent=TRUE)
	if(length(grep("Error in",err))!=0){
	plot(dat2Affy.m)
	}
	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(use.dat2Affy.m)
	use.dat2Affy.m<<-dat2Affy.m
	try(({
	rm(xf)
	xf<<-dim(use.dat2Affy.m)
	if(xf[2]%%3!=0 && xf[2]%%2!=0){
		use.dat2Affy.m<<-dat2Affy.m[,-xf[2]]
		xf<<-dim(use.dat2Affy.m)
	}
	if(xf[2]%%2==0)
	{
		yf=xf[2]/2
		groups<-c(rep("C",yf),rep("T",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.dat2Affy.m,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(dat2Affy.f)
		dat2Affy.f<<-fit2
		rm(design_Affy)
		design_Affy<<-design
#		print("Two-group Specific Filtering")
	}else
	{
		yf=xf[2]/3
		groups<-c(rep("C",yf),rep("T1",yf),rep("T2",yf))
		groups<-as.factor(groups)
		design<-model.matrix(~groups)
		fit<-lmFit(use.dat2Affy.m,design)
		yy<-try(toptable(fit,coef=2),silent=TRUE)
		if(length(grep("Error",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
		rm(dat2Affy.f)
		dat2Affy.f<<-fit2
		rm(design_Affy)
		design_Affy<<-design
#		print("Three-group Specific Filtering")
	}
	}),silent=TRUE)
	
	if(length(dat2Affy.f)==0)
	{
		try(({
			rsd<-rowSds(use.dat2Affy.m)
			i<-rsd>=2
			dat.f<-use.dat2Affy.m[i,]
			fit<-lmFit(dat.f)
			yy<-try(toptable(fit,coef=2),silent=TRUE)
			if(length(grep("Error in",yy))!=0){
			fit2<-eBayes(fit)
			}else{
				fit2<-eBayes(fit)
				}
			rm(dat2Affy.f)
			dat2Affy.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
			
		if(length(dat2Affy.f)==0)
		{
			try(({
				ff<-pOverA(A=1,p=0.5)
				i<-genefilter(use.dat2Affy.m,ff)
				dat.fo<-use.dat2Affy.m[i,]
				i<-genefilter(-use.dat2Affy.m,ff)
				dat.fu<-use.dat2Affy.m[i,]
				dat.f<-rbind(dat.fo,dat.fu)
				fit<-lmFit(dat.f)
				yy<-try(toptable(fit,coef=2),silent=TRUE)
				if(length(grep("Error in",yy))!=0){
					fit2<-eBayes(fit)
					}else{
						fit2<-eBayes(fit)
					}
				rm(dat2Affy.f)
				dat2Affy.f<<-fit2
#				print("Expression Filtering")
				}),silent=TRUE)
			}
		}

	if(length(dat2Affy.f)!=0){
		rm(ttx)
		err<-try(ttx<<-toptable(dat2Affy.f,coef=2,number=nrow(use.dat2Affy.m)),silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(dat2Affy.f,coef=3,number=nrow(use.dat2Affy.m)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)
		{
			err<-try(ttx<<-toptable(dat2Affy.f,coef=4,number=nrow(use.dat2Affy.m)),silent=TRUE)
		}
		if(length(grep("Error",err))!=0)ttx<<-toptable(dat2Affy.f,number=nrow(use.dat2Affy.m))
		rn<-row.names(ttx)[ttx$P.Value<=0.01]
		rm(dat2Affy.s)
		err=NULL;
		err<-try(dat2Affy.s<<-use.dat2Affy.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			dat2Affy.s<<-use.dat2Affy.m[rn,]
			}
		}
		rm(err)
		
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
	galert("Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(DE_Affy)
	err<-try(DE_Affy<<-toptable(dat2Affy.f,coef=2,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p"),silent=TRUE)
	if(dim(DE_Affy)[1]==0 || length(DE_Affy)==0)
	{
		err<-try(DE_Affy<<-toptable(dat2Affy.f,coef=3,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p"),silent=TRUE)
	}
	if(dim(DE_Affy)[1]==0 || length(DE_Affy)==0)
	{
		err<-try(DE_Affy<<-toptable(dat2Affy.f,coef=4,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p"),silent=TRUE)
	}
	if(dim(DE_Affy)[1]==0 || length(DE_Affy)==0)
	{
		DE_Affy<<-toptable(dat2Affy.f,number=20,lfc=2,adjust.method="BH",p.value=0.01,sort.by="p")
	}
	rm(err)
	
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)
	rm(pca_Affy)
	pca_Affy<<-prcomp(t(dat2Affy.m))
	plot(pca_Affy,main="Affymetrix PCA")
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)

	sample.dist_c1<-dist(t(dat2Affy.m)) 
	rm(sample.dist_Affy)
	sample.dist_Affy<<-as.dist(1-cor(dat2Affy.m,method="pearson"))
	rm(sample.clust_Affy)
	sample.clust_Affy<<-hclust(sample.dist_Affy,method="complete")
	plot(sample.clust_Affy)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	galert("Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"Working... Plz wait..."
	Sys.sleep(1)

	rm(use_DE_Affy)
	if(dim(as.matrix(DE_Affy))[1]>=20){
		er_x<-try(use_DE_Affy<<-as.matrix(DE_Affy)[1:20,],silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_Affy)<<-as.character(DE_Affy)
			use_DE_Affy<<-as.matrix(DE_Affy)[1:20,]
		}
	} else {
		er_x<-try(use_DE_Affy<<-as.matrix(DE_Affy),silent=TRUE)
		if(length(grep("Error",er_x))!=0)
		{
			rownames(DE_Affy)<<-as.character(DE_Affy)
			use_DE_Affy<<-as.matrix(DE_Affy)
		}
	}
	
	rm(Clas_Affy)		
	err2<-try((Clas_Affy<<-use.dat2Affy.m[rownames(use_DE_Affy),which(design_Affy[,1]==1 | design_Affy[,2]==1)]),silent=TRUE)
	if(length(grep("Error",err2))!=0)
	{
		err3<-try((Clas_Affy<<-use.dat2Affy.m[rownames(use_DE_Affy),which(design_Affy[,1]==1 | design_Affy[,2]==1)] | design_Affy[,3]==1),silent=TRUE)
		if(length(grep("Error",err3))!=0)
		{
			err4<-try(Clas_Affy<<-use.dat2Affy.m[rownames(use_DE_Affy),],silent=TRUE)
			if(length(grep("Error",err4))!=0)
			{
				try(Clas_Affy<<-use.dat2Affy.m[as.numeric(rownames(use_DE_Affy)),],silent=TRUE)
			}
		}
	}
	heatmap(Clas_Affy,Rowv=NA)
	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
		visible(g1_1)<-FALSE
		l$Affymetrix$Normalization<<-list()
		l$Affymetrix$QC_Plot<<-list()
		l$Affymetrix$Filtered<<-list()
		l$Affymetrix$Stat_Significant<<-list()
		l$Affymetrix$DGE<<-list()
		l$Affymetrix$PCA_Plot<<-list()
		l$Affymetrix$Cluster_Plot<<-list()
		l$Affymetrix$Classification<<-list()
		tr<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,400)
		visible(g1_1)<-TRUE
		display()
		}),silent=TRUE)
}
