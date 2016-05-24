#-----------------

#' two way anova only for the biclusterd matrix
#' diagnostics of the biclustered matrix
#' @dset stands for the dataset,
#' @bres stands for biclust result, should be one of fabia, biclust or isa2
#' @fit a parameter takes string; 'all'(all plots in one frame), 'mean', 'median', 'variance' or 'mad'.
#' @gby, group by 'conditions' or 'genes'

#----------------
exploreOnlybic<-function(dset,bres,pfor="all",gby="genes",mname="biclust",bnum=1,fabia.thresZ=0.5,fabia.thresL=NULL){
	
	# Change of variable name (to have same name convention as 'explorebic')
	fit <- pfor

	if(any(!mname %in% c("fabia","isa2","biclust","bicare"))){
		stop("`mname' must be one of `fabia',`isa2', 'biclust' or 'bicare'")
	}
	if(any(!gby %in% c("genes","conditions"))){
		stop("`gby' must be one of Genes' or 'Condtions'")
	} 
	if(any(!fit %in% c("all","mean","median","variance","mad"))){
		stop("`fit' must be one of `all','mean','median','variance','mad'")
	} 
	
	indgc<-indexedBic(dset,bres,mname,bnum,fabia.thresZ=fabia.thresZ,fabia.thresL=fabia.thresL)# returns the required indecies based on thier method names
	indg<-indgc[[1]]
	indc<-indgc[[2]]
	bic<-dset[indg,indc]
	#two way anova
	#for all biclust genes
	if(gby=="genes"){
		bgmed<-bgmen<-bgvar<-bgmad<-vector()
		for(i in 1:nrow(bic)){
			bgmed<-c(bgmed,median(bic[i,]))
			bgmen<-c(bgmen,mean(bic[i,]))
			bgvar<-c(bgvar,var(bic[i,]))
			bgmad<-c(bgmad,mad(bic[i,]))
		}
		bgall<-list(bgmed,bgmen,bgvar,bgmad)
		plotOnlybic(bgall,fit,gby)
	}
	
	#for all conditions
	if(gby=="conditions"){
		bcmed<-bcmen<-bcvar<-bcmad<-vector()
		for(i in 1:ncol(bic)){
			bcmed<-c(bcmed,median(bic[,i]))
			bcmen<-c(bcmen,mean(bic[,i]))
			bcvar<-c(bcvar,var(bic[,i]))
			bcmad<-c(bcmad,mad(bic[,i]))
		}
		bcall<-list(bcmed,bcmen,bcvar,bcmad)
		plotOnlybic(bcall,fit,gby)
	}
}