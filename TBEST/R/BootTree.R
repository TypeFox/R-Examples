BootTree<-function(mydata, mymethod, mymetric, ntest, node)
{
    myinput <- mydata
    if(mymetric!="pearson"&mymetric!="kendall"&mymetric!="spearman"){
    	hc<-hclust(dist(myinput,method=mymetric),method=mymethod)
    }
    if(mymetric=="pearson"|mymetric=="kendall"|mymetric=="spearman"){
        hc<-hclust(as.dist(1-cor(t(myinput),method=mymetric,
        	use="pairwise.complete.obs")),method=mymethod)
    }
    mem<-LeafContent(hc,mynode=node)
    allcounts <- matrix(ncol = length(seq(0.5,1,by=0.5)), nrow = 1, data = 0)
    colnames(allcounts) <- paste("r",seq(0.5,1,by=0.5))
    for(j in 1:ncol(allcounts)){
    for (i in 1:ntest) {
        #myrdata<-myinput[sample(1:nrow(myinput),replace=T),]
	#myrdata<-myinput[unique(sample(1:nrow(myinput),replace=T)),]
        myrdata<-myinput[,sample(1:ncol(myinput),size=floor(seq(0.5,1,by=0.5)[j]*ncol(myinput)),replace=T)]
	#rindextable <- TreeStat(myrdata, mystat = mystat, method = mymethod, 
        #    metric = mymetric, metric.args = metric.args)
	if(mymetric!="pearson"&mymetric!="kendall"&mymetric!="spearman"){
        	hc<-hclust(dist(myrdata,method=mymetric),method=mymethod)
    	}
    	if(mymetric=="pearson"|mymetric=="kendall"|mymetric=="spearman"){
		if(sum(apply(myrdata,1,var)==0)>0){
			myrdata[which(apply(myrdata,1,var)==0),]<-
			jitter(as.numeric(myrdata[which(apply(myrdata,1,var)==0),]))
		}
        	hc<-hclust(as.dist(1-cor(t(myrdata),method=mymetric,
                	use="pairwise.complete.obs")),method=mymethod)
    	}
	pp<-0
        for(k in 1:nrow(hc$merge)){
		rmem<-LeafContent(hc,mynode=k)
		if(length(unlist(rmem))<=length(unlist(mem)))pp<-
			max(pp,sum(!is.na(match(unlist(mem),unlist(rmem))))/length(unlist(mem)))
	}
	allcounts[,j]<-allcounts[,j]+pp
	}
	}
	allcounts<-allcounts/1000
	return(allcounts)
}
