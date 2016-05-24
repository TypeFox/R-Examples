RearrangeJoseph <- function(hc, dis,class,cpp)
{

  coef=1.5
	num=dim(hc$merge)[1]
	A=array(1,c(num+1,num+1))
	r=array(1,c(num+1,num+1))

	## num - number of internal nodes

	maxI=array(0,c(num+1,num+1))
	maxJ=array(0,c(num+1,num+1))
	ind=num

	ptm <- Sys.time()
  ## function
	node=testBar(hc)
	## function for node selection
  flag=CalcMerge(hc,node,class)

## change matrix hc$merge
	hO<-hc
	nodeO<-node
	out=SubTree(ind,flag,node,hc,A,r,coef)
	hc=out$hc
	node=out$node
	A=out$A
	r=out$r

	vrem=hO$merge==hc$merge
	indexS=NULL
	for(i in 1:num)
	{
		if(all(vrem[i,])) indexS=c(indexS,i)
	}
  #browser()
  if(cpp)
  {
	res=OrderingJosephC(ind-1, hc$merge, node, A, r, maxI, maxJ, class, coef)
  }
  else
  {
    res=OrderingJoseph(ind, hc, node, A, r, maxI, maxJ, class, coef)
  }
	cat("Elapsed time =",format(Sys.time()-ptm), "\n")

	A=res$A
	maxI=res$maxI
	maxJ=res$maxJ
	r=res$r
  write.table(A,file="A.txt")
	write.table(maxI,file="minI.txt")
	write.table(maxJ,file="minJ.txt")
	
	col=which.max(apply(A[node[[ind]]$left,node[[ind]]$right,drop=FALSE],2,max))
	row=which.max(apply(A[node[[ind]]$left,node[[ind]]$right,drop=FALSE],1,max))
	col=node[[ind]]$right[col]
	row=node[[ind]]$left[row]

	#browser()
	hcl=funMerge(ind,row,col,hO,node,maxI,maxJ,cpp)
	
	
	dend=as.dendrogram(hcl)
	order=order.dendrogram(dend)
	hcl$order=order
	sum=0
	for( i in 1:(num-1))
	{
		sum=sum+dis[hc$order[i],hc$order[i+1]]
	}
	#browser()
	sum=0
	for( i in 1:(num-1))
	{
		sum=sum+dis[order[i],order[i+1]]
	}
	return(list(hcl=hcl,A=A,max=A[row,col]))

}

