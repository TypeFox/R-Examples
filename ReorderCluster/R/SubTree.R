SubTree <- function(ind,flag,node,hc,A,r,coef)
{
	if(flag[ind]==0)
	{
	return(list(hc=hc, node=node, leaf=node[[ind]][[1]][1],A=A,r=r))
	}
	for(i in 1:2)
	{
	if(hc$merge[ind,i]>0)
	{
		res=SubTree(hc$merge[ind,i],flag,node,hc,A,r,coef)
		hc=res$hc
		node=res$node
		A=res$A
		r=res$r
		if(res$leaf>0)
		{
			hc$merge[ind,i]=-res$leaf
			##browser()
			A[res$leaf,res$leaf]=length(node[[ind]][[i]])^coef
			r[res$leaf,res$leaf]=length(node[[ind]][[i]])
			node[[ind]][[i]]=res$leaf
			
		}
		else
		{
			node[[ind]][[i]]=c(node[[hc$merge[ind,i]]][[1]],node[[hc$merge[ind,i]]][[2]])
		}
	}
	}
	return(list(hc=hc, node=node,leaf=-1,A=A,r=r))
}
