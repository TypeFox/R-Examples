testBar<-function(hc)
{
	dd=dim(hc$merge)
	sss=list(left=NULL, right=NULL, LP=NULL, RP=NULL)
	sss=rep(list(sss),dd[1])	
	

	for(i in 1:dd[1])
	{
	for(j in 1:dd[2])
	{
		if(hc$merge[i,j]<0)
		{
			sss[[i]][[j]]=c(sss[[i]][[j]],-hc$merge[i,j])
		}
		else
		{
			sss[[i]][[j]]=c(sss[[i]][[j]],sss[[hc$merge[i,j]]][[1]],
			sss[[hc$merge[i,j]]][[2]])
		}
	}
	}
	sss
}