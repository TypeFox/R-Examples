#----------------------
#' exploreBic function: calculates the mean, median, variance and the three quantiles
#' plots and prints either for a single or for all summary biclust results vs outside the bic.
#' plots and prints according to the user selection. either by gene or by condtion.
#' @ dset,bres, pfor and mname are parameters that must be fill by the user need.
#' @ gby, optional parameter, default is genes
#' @ dset and bres parameteres stands for dataset and biclust object respectively.
#' @ pfor is a parameter for summary selection(i.e. 'mean','median'.. etc)
#' @ mname, name of the biclust method
#' @ bnum, which biclust to be plotted? the first, second?... etc 
#' outcome is a user selection plot and 
#' a printed summary of the two groups(biclust and outbic)

#----------------------
exploreBic<-function(dset,bres,gby="genes",pfor="mean",mname="biclust",bnum=1,fabia.thresZ=0.5,fabia.thresL=NULL){
		
	if(any(!pfor %in% c("all","mean","variance","median","quant","mad"))) {
		stop("`pfor' must be one of `all', `mean', `variance','median', mad, or `quant'")
	}
	if(any(!mname %in% c("fabia","isa2","biclust","bicare"))){
		stop("`mname' must be one of `fabia',`isa2', 'biclust' or 'bicare'")
	} 
	if(any(!gby %in% c("genes","conditions"))){
		stop("`gby' must be one of `genes', or `conditions'")
	}
	ind.gc<-indexedBic(dset,bres,mname,bnum,fabia.thresZ=fabia.thresZ,fabia.thresL=fabia.thresL)
	indg<-ind.gc[[1]]
	indc<-ind.gc[[2]]
	#check the group
	if(gby=="conditions"){
		#group the genes in to two.
		rnams <- colnames(dset)
		grp <- rep(1, length(rnams))
		grp[indc] <- 2
		#calculate  for bic genes and conditions
		bic.mat<-dset[indg,indc]
		bic.mat<-t(bic.mat)
		sbic<-exploreCalc(bic.mat)
		sbic<-t(sbic)
		#run
		obic.mat<-dset[indg,-indc]
		obic.mat<-t(obic.mat)
		obic<-exploreCalc(obic.mat)
		obic<-t(obic)
		#call for plot
		explorePlot(sbic,obic,pfor)
	}
	if(gby=="genes"){
		#group the genes in to two.
		rnams <- rownames(dset)
		grp <- rep(1, length(rnams))
		grp[indg] <- 2
		#calculate  for bic genes and conditions
		bic.mat<-dset[indg,indc]
		#bic.mat<-t(bic.mat)
		sbic<-exploreCalc(bic.mat)
		#sbic<-t(sbic)
		#run
		obic.mat<-dset[-indg,indc]
		#obic.mat<-t(obic.mat)
		obic<-exploreCalc(obic.mat)
		#obic<-t(obic)
		#call for plot
		explorePlot(sbic,obic,pfor)
	}
	
}
