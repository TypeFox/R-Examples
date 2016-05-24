OrderingJoseph<-function(ind, hc, node, A, r, maxI, maxJ, class, coef)
{

	if(ind>=0)
	{
		res=OrderingJoseph(hc$merge[ind,1],hc, node, A, r, maxI, maxJ,class, coef)
		res=OrderingJoseph(hc$merge[ind,2],hc, node, res$A, res$r, res$maxI, res$maxJ, class, coef)
		##if(ind==dim(hc$merge)[1]) browser()

		A=res$A
		maxI=res$maxI
		maxJ=res$maxJ
		r=res$r
		for(i in 1:length(node[[ind]][[1]]))
		{
			if(ind==dim(hc$merge)[1])
			{
				#cat( "I = ", i,"\n")
				#cat("Elapsed time =",format(Sys.time()-ptm), "\n")
			}
			u=node[[ind]][[1]][i]
			if(length(node[[ind]][[1]])==1)
			{
			vremI=node[[ind]][[1]]
			}
			else
			{
			if(u%in%node[[hc$merge[ind,1]]][[1]])
			{
			vremI=node[[hc$merge[ind,1]]][[2]]
			}
			else
			{
			vremI=node[[hc$merge[ind,1]]][[1]]
			}
			}
			for(j in 1:length(node[[ind]][[2]]))
			{
				w=node[[ind]][[2]][j]

				if(length(node[[ind]][[2]])==1)
				{
				vremJ=node[[ind]][[2]]
				}
				else
				{
				if(w%in%node[[hc$merge[ind,2]]][[1]])
				{
				vremJ=node[[hc$merge[ind,2]]][[2]]
				}
				else
				{
				vremJ=node[[hc$merge[ind,2]]][[1]]
				}
				}
##algorithm code---------------
				max=-1
				for(i1 in 1:length(vremI))
				{
					for(j1 in 1:length(vremJ))
					{
					A[u,w]=A[u,vremI[i1]]+A[vremJ[j1],w]
					if(class[vremI[i1]]==class[vremJ[j1]])
					{
						A[u,w]=A[u,w]-r[u,vremI[i1]]^coef-r[w,vremJ[j1]]^coef+(r[u,vremI[i1]]+r[w,vremJ[j1]])^coef
					}	
						
					if(A[u,w]>max)
					{
						max=A[u,w]
						maxI[u,w]=vremI[i1]
						maxJ[u,w]=vremJ[j1]
					}
					if (all(class[vremJ]==class[vremJ[1]])) break
					}
				if (all(class[vremI]==class[vremI[1]])) break
				}
				A[u,w]=max
				A[w,u]=A[u,w]
				maxI[w,u]=maxJ[u,w]
				maxJ[w,u]=maxI[u,w]
				if((all(class[node[[ind]]$right]==class[node[[ind]]$right[1]]))&&(class[maxI[u,w]]==class[maxJ[u,w]]))
				{
					##r[u,w]=length(node[[ind]]$right)+r[u,maxI[u,w]]
					r[u,w]=r[maxJ[u,w],w]+r[u,maxI[u,w]]

				}
				else
				{
				r[u,w]=r[maxJ[u,w],w]
				}
				if((all(class[node[[ind]]$left]==class[node[[ind]]$left[1]]))&&(class[maxJ[w,u]]==class[maxI[w,u]]))
				{
					##r[w,u]=length(node[[ind]]$left)+r[w,maxI[w,u]]
					r[w,u]=r[maxJ[w,u],u]+r[w,maxI[w,u]]

				}
				else
				{
				r[w,u]=r[maxJ[w,u],u]
				}
			}
		}
	}
	return(list(A=A,maxI=maxI,maxJ=maxJ,r=r))
}

