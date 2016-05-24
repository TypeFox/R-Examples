LeafContent<-function(myinput,mynode=NA){
	if(class(myinput)=="best"){
		indextable=myinput$indextable[,1:4]
		names<-dimnames(myinput$data)[[1]]
		if(any(is.na(mynode)))stop("Inappropriate node number")
	}
	if(class(myinput)=="partition"){
                indextable=myinput$best$indextable[,1:4];names<-dimnames(myinput$best$data)[[1]]
		mynode=as.numeric(unique(myinput$partition[,2]))
		if(any(mynode==nrow(indextable))){
			siglevel<-myinput$Call$siglevel
			if(is.null(siglevel))siglevel<-0.05
			if(myinput$sigvalue[mynode,2]>siglevel){
				stop("Most detailed partition does not exist! ")
			}
		}
        }
	if(class(myinput)=="hclust"){
		hc<-myinput
		if(any(is.na(mynode)))stop("Inappropriate node number")
		indextable<-cbind(hc$merge,hc$height)
        	dimnames(indextable)[[2]]<-c("index1","index2","height")
        	#cluster size
        	clustersize<-rep(NA,nrow(indextable))
        	csleft<-rep(NA,nrow(indextable))
        	csleft[indextable[,"index1"]<0]<-1
        	csright<-rep(NA,nrow(indextable))
        	csright[indextable[,"index2"]<0]<-1
        	while(is.na(sum(clustersize))){
                	clustersize<-csleft+csright
                	csleft[indextable[,"index1"]>0]<-
                        	clustersize[indextable[indextable[,"index1"]>0,"index1"]]
                	csright[indextable[,"index2"]>0]<-
                        	clustersize[indextable[indextable[,"index2"]>0,"index2"]]
        	}
		indextable<-cbind(indextable,clustersize)
	}
	if(class(mynode)!="numeric"&class(mynode)!="integer")stop("Inappropriate node number")
        if(max(mynode)>nrow(indextable))stop("Node number should be <= Sample Size - 1")
        if(min(mynode)< -(nrow(indextable)+1))stop("Node number should not be < -Sample Size")
	if(length(which(mynode==0))!=0)stop("Node number can not be 0")
	sigpartition<-mynode
        singleton<- -sigpartition[sigpartition<0]
        nodes<-sigpartition[sigpartition>0]
        membership<-vector("list",length(nodes))
	clusters<-vector("list",length(mynode))
	names(clusters)<-paste("branch",mynode)
	pos<-which(sigpartition>0)
	npos<-which(sigpartition<0)
	if(length(nodes)>=1){
        for(i in 1:length(nodes)){
        	myfamily<-nodes[i]
        	nmem<-0
                while(nmem<indextable[nodes[i],"clustersize"]){
                	membership[[i]]<-unique(c(myfamily,indextable[myfamily,"index1"],
                	indextable[myfamily,"index2"]))
                	myfamily<-membership[[i]][membership[[i]]>0]
                	nmem<-sum(membership[[i]]<0)
                }
                membership[[i]]<-(-membership[[i]][membership[[i]]<0])
        }
	for(i in 1:length(membership)){
		if(class(myinput)=="hclust")clusters[[pos[[i]]]]<-hc$labels[membership[[i]]]
                if(class(myinput)!="hclust")clusters[[pos[[i]]]]<-names[membership[[i]]]
        	}
	}
	if(length(npos)!=0){
	for(i in 1:length(npos)){
		if(class(myinput)=="hclust")clusters[[npos[[i]]]]<-hc$labels[singleton[i]]
		if(class(myinput)!="hclust")clusters[[npos[[i]]]]<-names[singleton[i]]
	}
	}
	return(clusters)
}
