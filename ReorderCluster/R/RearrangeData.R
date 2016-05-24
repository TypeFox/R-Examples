RearrangeData <- function(matr,class,clustering.method="complete",cDend=FALSE,dirpath="",cpp=TRUE)
{

  # to load the dataset the following command must be performed
  #sim <- read.delim("datafile.txt",header=T,sep='\t',row.names=1,as.is=T)
	matr=as.matrix(matr)
	dd=dim(matr)

	label=unique(class)
	weight=rep(1,length(label))
	names(weight)=as.character(label)

	Rowcolor=rainbow(length(label))
	rc=matrix(0,length(class),1)
	ds=matrix(0,length(class),1)

	for (j in 1:length(label))
  	{
		index=which(class==label[j])
		rc[index]=Rowcolor[j]
	}
	##algorithm insertion
	##label=rownames(sim)
##------------standardize
	matrS=scale(matr)
	matrS=matr
##---------------
	dist=dist(matrS)
	hc <- hclust(dist, method = clustering.method)
	##plot(hc,label=label)
	##barplot(rep(1,dd[1]),col=rc[hc$order],space=0,xpd=FALSE)

	dend=as.dendrogram(hc)
  
  if(cDend)
  {
	colorDendClass(dend,rc[hc$order])
	browser()
  }
  
	hv <-  heatmap(matrS,scale = "none",Rowv=dend, Colv=NA,
	               RowSideColors = rc, keep.dendro = TRUE)
  browser()

  if (dirpath!=""){
    pdf(file=paste(dirpath,"heatmap_raw.pdf",sep=""))
  }

	my_palette <- colorRampPalette(c("green", "black", "red"))(n = 399)
	hv <- heatmap.2(matrS,Rowv=dend,scale = "none",Colv=NA,
	                col=my_palette, RowSideColors = rc,trace="none",dendrogram="row",
	                symm=FALSE,symkey=FALSE,symbreaks=TRUE,labRow = rep("",dd[1]), labCol = rep("",dd[2]))#,key=FALSE)
  legend("topright",legend=label,col=Rowcolor,pch=15,cex=0.8)

  if (dirpath!=""){
    dev.off()
  }
	
	browser()

	
  res=RearrangeJoseph(hc,as.matrix(dist),class,cpp) 
  
	hcl=res$hcl
	
	#plot(hcl,label=label)

	#barplot(rep(1,dd[1]),col=rc[hcl$order],space=0,xpd=FALSE)
	
	dend=as.dendrogram(hcl)
  browser()
  
  if(cDend)
  {
	colorDendClass(dend,rc[hcl$order])
	browser()
  }
  
  hv <-  heatmap(matrS,scale = "none",Rowv=dend, Colv=NA,
               RowSideColors = rc, keep.dendro = TRUE)
  browser()

  if (dirpath!=""){
    pdf(file=paste(dirpath,"heatmap_raw.pdf",sep=""))
  }

	hv <- heatmap.2(matrS,Rowv=dend,scale = "none",Colv=NA,
	                col=my_palette, RowSideColors = rc,trace="none",dendrogram="row",
	                symm=FALSE,symkey=FALSE,symbreaks=TRUE,labRow = rep("",dd[1]), labCol = rep("",dd[2]))
  legend("topright",legend=label,col=Rowcolor,pch=15,cex=0.8)
  if (dirpath!=""){
    dev.off()
  }
  browser()
	
	return(list(order=hcl$order,A=res$A,class=class[hcl$order]))
	
}
