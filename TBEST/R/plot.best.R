plot.best<-function(x,mystat="fldc",siglevel=0.05,sigtype=c("raw","corrected","fdr"),
	partition=NA,print.num=TRUE,print.lab=TRUE,float=0.01,
	col.best=c(2,3),cex.best=0.8,cex.leaf=0.8,font.best=NULL,main=NULL,sub=NULL,xlab=NULL,
	metric.args=list(),...){
    	clusterobj<-x
	method<-clusterobj$Call$mymethod
	metric<-clusterobj$Call$mymetric
    	if (is.null(sub)) 
        	sub = paste("Cluster method: ",method)
    	if (is.null(xlab)) 
        	xlab = paste("Distance: ", metric)
	myinput<-clusterobj$data
	if(data.class(myinput)=="dist"){
		hc<-hclust(myinput,method=method)
		myinput<-data.matrix(myinput)
	}
        if(data.class(myinput)=="matrix"){
		if(!(metric %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski",
                "pearson","kendall","spearman"))){
                        define.metric<-match.fun(metric)
			mymetric.args<-vector("list",length(metric.args)+1)
                        mymetric.args[[1]]<-myinput
                        if(length(mymetric.args)>1){mymetric.args[2:length(mymetric.args)]<-
                                metric.args}
                        mydis<-do.call(define.metric,mymetric.args)
                        #mydis<-define.metric(myinput,...)
                        mydis<-data.matrix(mydis)
			if(nrow(mydis)!=nrow(myinput)|ncol(mydis)!=nrow(myinput)){
                                stop("User defined metric returns wrong dimension")
                        }
                        hc<-hclust(as.dist(mydis),method=method)
                }
		else{
                	if(metric!="pearson"&metric!="kendall"&metric!="spearman"){
                        	hc<-hclust(dist(myinput,method=metric),method=method)
                	}
                	if(metric=="pearson"|metric=="kendall"|metric=="spearman"){
                        	hc<-hclust(as.dist(1-cor(t(myinput),method=metric,
                                	use="pairwise.complete.obs")),method=method)
                	}
		}
        }
        if(data.class(myinput)=="hclust")hc<-myinput
	if(length(partition)==1){
		sigtype<-match.arg(sigtype)
		if (is.null(main))main = switch(sigtype,
			raw="Dendrogram with P-values",
			corrected="Dendrogram with corrected P-values",
			fdr="Dendrogram with false discovery rate")
		if(print.lab){plot(hc, labels=hc$labels, main = main, sub = sub, xlab = xlab,cex=cex.leaf,...)}
        	if(!print.lab){plot(hc, labels=print.lab, main = main, sub = sub, xlab = xlab,cex=cex.leaf,...)}
		pval<-clusterobj$indextable[,paste("p",mystat,sep="")] 
		if(sigtype=="raw"){sigp<-siglevel}
        	if(sigtype=="fdr"){
                	qobj<-fdrtool::fdrtool(pval[!is.na(pval)], statistic="pvalue",plot=FALSE,verbose=FALSE)
                	sigp<-siglevel
                	pval<-qobj$qval
        	}
        	if(sigtype=="corrected"){
                	sigp<- siglevel
                	pval<- 1-(1-pval)^(nrow(myinput)-2)
        	}
		pval[pval>sigp]<-NA
        	treetext(hc,pval, col = col.best, cex = cex.best, float = float, 
			font=font.best,print.num = print.num)
	}
	if(length(partition)>1){
		sigtype<-partition$Call$sigtype
		if (is.null(main))main = switch(sigtype,
                        raw="Dendrogram with P-values",
                        corrected="Dendrogram with corrected P-values",
                        fdr="Dendrogram with false discovery rate")
		if(print.lab){plot(hc, labels=hc$labels, main = main, sub = sub, xlab = xlab,cex=cex.leaf,...)}
        	if(!print.lab){plot(hc, labels=print.lab, main = main, sub = sub, xlab = xlab,cex=cex.leaf,...)}
		pval<-partition$sigvalue[,2]
		sigp<-partition$Call$siglevel
		if(is.null(partition$Call$siglevel))sigp<-0.05
		if(length(unique(partition$partition[,2]))==1){
			pval[pval>sigp]<-NA
                	treetext(hc,pval, col = col.best, cex = cex.best, float = float, 
				font=font.best,print.num = print.num)
		}else{
		partitionp<-rep(NA,length(pval))
		partitionp[unique(partition$partition[,2])]<-pval[unique(partition$partition[,2])]
		treetext(hc,partitionp, col = col.best, cex = cex.best, float = float, 
			font=font.best,print.num = print.num)
		}
        }	
}
